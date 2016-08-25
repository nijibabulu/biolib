#! /usr/bin/env python

import docopt
import collections
from pylib import SeqFeatureIO

__doc__ = '''
Usage: gff3_to_bed.py [options] GFF3

This converts a gff3 file to a UCSC-style bed file, grouping
fetures by gene.

Options:
    --gene-name-attr ATTR   gene name attribute [default: Name]
    --gene-feature GENE     gene feature [default: mRNA]
    --exon-feature EXON     exon feature [default: exon]
    --cds-feature CDS       cds feature [default: CDS]
'''

args = docopt.docopt(__doc__)

gff3_strands = ['.','+','-']
class Gene(object):
    def __init__(self,name,exons,cds):
        self.cds = cds
        self.exons = exons
        self.name = name

        self.sizes = []
        self.starts = []

        exon_cmp = lambda a,b: cmp(a.location.start,b.location.start)
        self.exons = list(sorted(self.exons, cmp=exon_cmp))
        self.cds = list(sorted(self.cds, cmp=exon_cmp))
        self.start = int(self.exons[0].location.start)
        self.end = int(self.exons[0].location.start)
        self.chrom = self.exons[0].ref
        if len(self.cds):
            self.cds_start = int(self.cds[0].location.start)
            self.cds_end = int(self.cds[-1].location.end)
        else:
            self.cds_start = self.start
            self.cds_end = self.end
        self.strand = gff3_strands[self.exons[0].strand]
        for exon in self.exons:
            self.starts.append(int(exon.location.start)-self.start)
            self.sizes.append(len(exon))

    def __str__(self):
        return '\t'.join(
            str(x) for x in [self.chrom,self.start,self.end,self.name,
                             '0',self.strand,self.cds_start,self.cds_end,
                             '0',len(self.sizes),
                             ','.join(str(x) for x in self.sizes),
                             ','.join(str(x) for x in self.starts)])

    

class GeneFactory(object):
    def __init__(self):
        self.exons = []
        self.cds = []
    def add_exon(self,exon):
        if exon.type == args['--exon-feature']:
            self.exons.append(exon)
        if exon.type == args['--cds-feature']:
            self.cds.append(exon)
    def make_gene(self,name):
        return Gene(name,self.exons,self.cds)


gene_ids = {}
for sf in SeqFeatureIO.parse(open(args['GFF3']),'gff3'):
    if sf.type == args['--gene-feature']:
        gene_ids[sf.attributes['ID']] = sf.attributes[args['--gene-name-attr']]

facs = collections.defaultdict(GeneFactory)
for sf in SeqFeatureIO.parse(open(args['GFF3']),'gff3'):
    if sf.type == args['--exon-feature'] or sf.type ==args['--cds-feature']:
        if sf.attributes['Parent'] in gene_ids:
            facs[sf.attributes['Parent']].add_exon(sf)

for name,fac in facs.items():
    print fac.make_gene(gene_ids[name])

    
