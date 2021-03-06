"""Supplement the missing support for GFF3 from biopython.
Taken from FastaIO"""

import re
from Bio.SeqFeature import SeqFeature,FeatureLocation

_gff3_strand_to_numeric = { '+': 1, '-': -1, '.': 0}
_numeric_to_gff3_strand = { '1': '+', '-1': '-', '0': '.' }

#This is a generator function!
def GTFIterator(handle):
  """Generator function to iterate over Fasta records (as SeqRecord objects).

  handle - input file

  If this is not given, then the entire title line will be used
  as the description, and the first word as the id and name.

  Note that use of title2ids matches that of Bio.Fasta.SequenceParser
  but the defaults are slightly different.
  """
  line_no = 0
  #Skip any text before the first record (e.g. blank lines, comments)
  while True :
    line_no += 1
    line = handle.next()
    if line is None:
      return
    line = line.strip()
    if len(line) == 0 or line[0] == '#':
      continue
    try:
      ref,source,type,start,end,score,strand,frame,attributes = \
        line.split('\t')
    except:
      raise ValueError, 'Problem with line %d in %s.  Line was\n%s' %\
        (line_no,handle.name,line)

    attr_pairs = attributes.strip(';').split(';')
    attr_dict = {}
    for pair in attr_pairs:
        try:
            key,value = pair.strip().split(' ',1)
            value = value.strip('"')
            attr_dict[key] = value
        except:
            attr_dict['transcript_id']=pair
            

    #attr_dict = dict(map(lambda x: tuple([y.strip('"') 
                                          #for y in x.split(' ',1)]), attr_pairs))
    result = SeqFeature(location=FeatureLocation(int(start),int(end)),
      type=type,strand=_gff3_strand_to_numeric[strand],ref=ref,ref_db=source)
    result.name = result.id = attr_dict.get('transcript_id',None)
    result.attributes = attr_dict # not an official property of SeqFeature.
    result.frame = frame
    yield result

class GTFWriter:
    """Class to write Fasta format files."""
    def __init__(self, handle):
        """Create a Fasta writer.

        handle - Handle to an output file, e.g. as returned
                 by open(filename, "w")
        You can either use:

        myWriter = FastaWriter(open(filename,"w"))
        writer.write_file(myRecords)

        Or, follow the sequential file writer system, for example:

        myWriter = FastaWriter(open(filename,"w"))
        writer.write_header() # does nothing for Fasta files
        ...
        Multiple calls to writer.write_record() and/or writer.write_records()
        ...
        writer.write_footer() # does nothing for Fasta files
        writer.close()
        """
        self.handle = handle
        self._header_written = False
        self._feature_written = False

    def write_header(self):
      self._header_written = True

    def write_feature(self, feature):
      """Write a single GFF3 feature to the file."""
      if not self._header_written:
          self.write_header()
      self._feature_written = True
      source = getattr(feature,"ref_db","Unknown") 
      self.handle.write("\t".join([
        feature.ref,
        source,
        feature.type,
        str(int(feature.location.nofuzzy_start)+1), # gff3=1-base  SeqFeature=0-base
        str(int(feature.location.nofuzzy_end)+1),
        str(getattr(feature,"score",".") or '.'),
        _numeric_to_gff3_strand[str(feature.strand)],
        getattr(feature,"phase",".") or '.',
        ]))
      self.handle.write("\t")
      if hasattr(feature,"attributes"):
          for k,v in feature.attributes.items():
              self.handle.write('%s "%s"; ' % (k,v))
      self.handle.write("\n")

    def write_file(self, features): # should be in a superclass...
      written = 0
      for feature in features:
        self.write_feature(feature)
        written += 1
      return written
