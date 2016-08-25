"""Supplement the missing support for SeqFeatureIO from biopython.
Taken from SeqIO"""
import os
import GFF3IO
import GTFIO
from StringIO import StringIO
from Bio.SeqRecord import SeqRecord

#Convention for format names is "mainname-subtype" in lower case.
#Please use the same names as BioPerl where possible.
#
#Note that this simple system copes with defining
#multiple possible iterators for a given format/extension
#with the -subtype suffix
#
#Most alignment file formats will be handled via Bio.AlignIO

_FormatToIterator ={'gff3' : GFF3IO.GFF3Iterator,
                    'gtf'  : GTFIO.GTFIterator,
                    }

_FormatToWriter ={'gff3' : GFF3IO.GFF3Writer,
                  'gtf'  : GTFIO.GTFWriter,
                  }

readable_formats = _FormatToIterator.keys()
writeable_formats = _FormatToWriter.keys()

def write(features, handle, format) :
    """Write complete set of features to a file.

    features - A list (or iterator) of SeqFeature objects.
    handle    - File handle object to write to.
    format    - lower case string describing the file format to write.

    You should close the handle after calling this function.

    Returns the number of records written (as an integer).
    """

    #Try and give helpful error messages:
    if isinstance(handle, basestring) :
        raise TypeError('Need a file handle, not a string (i.e. not a filename)')
    if not isinstance(format, basestring) :
        raise TypeError('Need a string for the file format (lower case)')
    if not format :
        raise ValueError('Format required (lower case string)')
    if format != format.lower() :
        raise ValueError('Format string "%s" should be lower case' % format)
    if isinstance(features,SeqRecord):
        raise ValueError('Use a SeqRecord list/iterator, not just a single SeqRecord')

    #Map the file format to a writer class
    if format in _FormatToWriter :
        writer_class = _FormatToWriter[format]
        count = writer_class(handle).write_file(features)
    else:
        raise ValueError('Unknown format "%s"' % format)

    assert isinstance(count, int), 'Internal error - the underlying writer ' \
           + ' should have returned the record count, not %s' % repr(count)
    return count
    
def parse(handle, format):
    """Turns a sequence file into an iterator returning SeqRecords.

    handle   - handle to the file.
    format   - lower case string describing the file format.

    Typical usage, opening a file to read in, and looping over the record(s):

    """
    #Try and give helpful error messages:
    if isinstance(handle, basestring) :
        raise TypeError("Need a file handle, not a string (i.e. not a filename)")
    if not isinstance(format, basestring) :
        raise TypeError("Need a string for the file format (lower case)")
    if not format :
        raise ValueError("Format required (lower case string)")
    if format != format.lower() :
        raise ValueError("Format string '%s' should be lower case" % format)

    #Map the file format to a sequence iterator:    
    if format in _FormatToIterator :
        iterator_generator = _FormatToIterator[format]
        return iterator_generator(handle)
    else :
        raise ValueError("Unknown format '%s'" % format)

def to_dict(features, key_function=None) :
    """Turns a sequence iterator or list into a dictionary.

    features  - An iterator that returns SeqRecord objects,
                 or simply a list of SeqRecord objects.
    key_function - Optional function which when given a SeqRecord
                   returns a unique string for the dictionary key.

    e.g. key_function = lambda rec : rec.name
    or,  key_function = lambda rec : rec.description.split()[0]

    If key_function is ommitted then record.id is used, on the
    assumption that the records objects returned are SeqRecords
    with a unique id field.

    If there are duplicate keys, an error is raised.

    Example usage, defaulting to using the record.id as key:

    >>> from Bio import SeqIO
    >>> handle = open("GenBank/cor6_6.gb", "rU")
    >>> format = "genbank"
    >>> id_dict = SeqIO.to_dict(SeqIO.parse(handle, format))
    >>> print id_dict.keys()
    ['L31939.1', 'AJ237582.1', 'X62281.1', 'AF297471.1', 'X55053.1', 'M81224.1']
    >>> print id_dict["L31939.1"].description
    Brassica rapa (clone bif72) kin mRNA, complete cds.

    A more complex example, using the key_function argument in order to use
    a sequence checksum as the dictionary key:
    
    >>> from Bio import SeqIO
    >>> from Bio.SeqUtils.CheckSum import seguid
    >>> handle = open("GenBank/cor6_6.gb", "rU")
    >>> format = "genbank"
    >>> seguid_dict = SeqIO.to_dict(SeqIO.parse(handle, format), \
                                    key_function = lambda rec : seguid(rec.seq))
    >>> for key, record in seguid_dict.iteritems() :
    ...     print key, record.id
    SabZaA4V2eLE9/2Fm5FnyYy07J4 X55053.1
    l7gjJFE6W/S1jJn5+1ASrUKW/FA X62281.1
    /wQvmrl87QWcm9llO4/efg23Vgg AJ237582.1
    TtWsXo45S3ZclIBy4X/WJc39+CY M81224.1
    uVEYeAQSV5EDQOnFoeMmVea+Oow AF297471.1
    BUg6YxXSKWEcFFH0L08JzaLGhQs L31939.1

    """    
    if key_function is None :
        key_function = lambda sf : sf.name

    d = dict()
    for record in features :
        try:
            key = key_function(record)
        except:
            print record.attributes
            raise
        if key in d :
            d[key].append(record)
        else:
            d[key] = [record]
    return d


#if __name__ == "__main__":
#Run the doctests
#_test()
