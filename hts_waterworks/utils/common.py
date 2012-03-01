# coding=utf-8

""" Common resources and functions for working with pygr
"""

import string
import re
from pygr import worldbase, annotation, cnestedlist
import itertools
from itertools import ifilter, islice
import time

genome2resource = {'mm9':'Bio.Seq.Genome.MOUSE.mm9',
           'hg18':'Bio.Seq.Genome.HUMAN.hg18',
           'hg19':'Bio.Seq.Genome.HUMAN.hg19',
           'dm3':'Bio.Seq.Genome.DROME.dm3',
           'xenTro2':'Bio.Seq.Genome.XENTR.xenTro2'}

def getGenome(genome):
    if genome in genome2resource:
        genome = genome2resource[genome]
    return worldbase(genome, download=True)
def getFullGenomeName(genome):
    if genome in genome2resource:
        genome = genome2resource[genome]
    return genome



repeatSeq = '[atcg]'  # lowercase letters = repeat sequence
nSeq = '[nN]'      # don't allow N's in the sequence
degenTable = {'A':'A', 'C':'C', 'G':'G', 'T':'T',
              'R':'AG', 'Y':'CT', 'M':'CA', 'K':'TG', 'W':'TA', 'S':'CG',
              'B':'CTG', 'D':'ATG', 'H':'ATC', 'V':'ACG',
              'N':'ACGT'
              }

def consensus_to_regex(consensus):
    """Convert consensus sequence to regular expression"""
    return ''.join(['[%s]' % degenTable[l] for l in consensus])

iupac2matrix = dict(zip('ACGTRYMKWSBDHVN',
            [ [1,0,0,0], #A
                [0,1,0,0], #C
                [0,0,1,0], #G
            [0,0,0,1], #T
            [1,0,1,0], #R=AG
            [0,1,0,1], #Y=CT
            [1,1,0,0], #M=CA
            [0,0,1,1], #K=TG
            [1,0,0,1], #W=TA
            [0,1,1,0], #S=CG
            [0,1,1,1], #B=CTG
            [1,0,1,1], #D=ATG
            [1,1,0,1], #H=ATC
            [1,1,1,0], #V=ACG
            [1,1,1,1] #N=ACGT
            ]))



def parse_ucsc_range(rangeStr):
    """ Parse ucsc formatted range to its parts """
    chrom, start, stop = re.search(r'(\w+):([\d-]+)-([\d-]+)',
                                                    rangeStr).groups()
    start, stop = int(start), int(stop)
    return chrom, start, stop

def make_ucsc_range(chrom, start, stop):
    """ Create ucsc version of given coordinates """
    if start < 0 and stop < 0:
        # pygr-style coordinates on negative orientation
        start, stop = -stop, -start
    return '%s:%s-%s' % (chrom, start, stop)

def pygrSeqToBed(seq, name='.', score='.'):
    """Convert a pygr sequence to BED format"""
    return '\t'.join(map(str, [seq.id, seq.start, seq.stop, '.', '.',
                               '+' if seq.orientation==1 else '-']))

def bedCommentFilter(line):
    """Check if line is a comment or track line"""
    not_empty = len(line) > 0
    not_comment = line[0] not in ['#', '"']
    not_track = (len(line.split()) > 0 and line.split()[0] != 'track'
                                                and line.split()[0] != 'chr')
    return not_empty and not_comment and not_track


def readBedLines(bedlines, dataOnly=True):
    """Return data portion of a bed file"""
    filteredLines = itertools.ifilter(bedCommentFilter, bedlines)
    for line in filteredLines:
        features = line.strip().split('\t')
        if len(features) >= 6:
            chrom, start, stop, _, _, strand = features[:6]
        else:
            chrom, start, stop = features[:3]
            strand = '+'
        start, stop = int(start), int(stop)
        if dataOnly:
            yield chrom, start, stop, strand
        else:
            yield (chrom, start, stop) + tuple(features[3:5]) + (strand,) + \
                                                        tuple(features[6:])

def wrapBedToPygrSeqs(bedLines, genome):
    """Map bedLines to genome, as if they were coming from an AnnotDB"""
    for chrom, start, stop, strand in bedLines:
        slice = genome[chrom][start:stop]
        if strand in ['-', -1]:
            slice = -slice
        yield slice

def parseFastaLines(fastaLines):
    """Read lines from fasta file, one sequence at a time"""
    curID=None
    curSeq = ''
    for line in fastaLines:
        if line[0] == '>':
            if curID is not None:
                yield curID, curSeq
            curID = line[1:]
            curSeq = ''
        else:
            curSeq += line.strip()
    if curID is not None:
        yield curID, curSeq

_revCompTable = string.maketrans("ACGTNacgtn", "TGCANtgcan")
def reverseComplement(sequence):
    return sequence.translate(_revCompTable)[::-1]

def makeNormalSeq(seq, allowedSeq='acgtnACGTN'):
    "Replace non-DNA characters with N"
    return re.sub('[^%s]' % allowedSeq, 'N', seq)

def maskRepeats(seq):
    "Replace lower-case nucleotides in seq with N"
    return re.sub(repeatSeq, 'N', seq)

def peakIter(sequence):
    'peak at the first element, and chain it back so the iterable is untouched'
    seqIter = iter(sequence)
    firstElem = seqIter.next()
    return firstElem, itertools.chain([firstElem], seqIter)

class Bag(dict):
    """dict-like class that supports attribute access as well as getitem.
    >>> x = Bag()
    >>> x['foo'] = 'bar'
    >>> x.foo
    'bar'
    """
    def __init__(self, *args, **kw):
        dict.__init__(self, *args, **kw)
        for k in self.keys():
            self.__dict__[k] = self.__getitem__(k)
    def __setitem__(self, k, v):
        dict.__setitem__(self, k, v)
        self.__dict__[k] = v

def memoize(function):
    '''Decorator to memoize a function with immutable args'''
    cache = {}
    def decorated_function(*args):
        if args in cache:
            return cache[args]
        else:
            val = function(*args)
            cache[args] = val
            return val
    return decorated_function

def taketwo(sequence):
    '''Return two elements from sequence at a time
    >>>list(taketwo(range(10)))
    [(0, 1), (2, 3), (4, 5), (6, 7), (8, 9)]
    '''
    seqIter = iter(sequence)
    try:
        while True:
            part1 = seqIter.next()
            part2 = seqIter.next()
            yield (part1, part2)
    except StopIteration:
        pass

lastTime = None
def tic():
    global lastTime
    lastTime = time.time()

def toc():
    global lastTime
    if not lastTime:
        tic()
    return time.time() - lastTime

def breakpoint():
    """Drop-in debugger, preferring ipython pdb."""
    try:
        import ipdb as pdb
    except ImportError:
        import pdb
    pdb.set_trace()

def fastaToSequenceList(inputFile):
    """Convert fasta files to a list of sequences."""
    fasta_file = open(inputFile, "r")
    lines = fasta_file.readlines()
    sequences = []
    fasta_file.close()

    for line in lines:
        if line.startswith(">") != True:
            line = line.rstrip()
            sequences.append(line)
    return sequences

def parseFastq(fastqfile):
    "parse a fastq-formatted file, yielding a (header, sequence, quality) tuple"
    fastqiter = (l.strip('\n') for l in fastqfile)  # strip trailing newlines 
    fastqiter = ifilter(lambda l: l, fastqiter)  # skip blank lines
    while True:
        fqlines = list(islice(fastqiter, 4))
        if len(fqlines) == 4:
            header1,seq,header2,qual = fqlines
            if header1.startswith('@') and header2.startswith('+'):
                yield header1[1:], seq, qual
            else:
                raise ValueError("Invalid header lines: %s and %s" % (header1,
                                                                      header2))
        elif len(fqlines) == 0:
            raise StopIteration
        else:
            raise EOFError("Failed to parse four lines from fastq file!")

def parseGTF(gtfFile):
    """parse a GTF file and return the values in BED format, as:
chrom    source    feature    start    stop    score    strand    frame    {attributes}
    """
    for line in gtfFile:
        fields = line.rstrip('\n').split('\t')
        attrs = fields[8].strip().strip(';').split(';')
        attrs = dict(attr.strip().replace('"','').split(' ') for attr in attrs)
        yield tuple(fields[:8]) + (attrs,)


def flatten(l, ltypes=(list, tuple)):
    """Flatten an arbitrarily deep list of lists.
    from http://rightfootin.blogspot.com/2006/09/more-on-python-flatten.html
    """
    ltype = type(l)
    l = list(l)
    i = 0
    while i < len(l):
        while isinstance(l[i], ltypes):
            if not l[i]:
                l.pop(i)
                i -= 1
                break
            else:
                l[i:i + 1] = l[i]
        i += 1
    return ltype(l)

def bedToNLMSA(bedlines, genome, field_locations=dict(id=0, start=1, stop=2,
                                            name=3, score=4, orientation=-1)):
    "Build a pygr resource off of the BED file in_name"
    annotDB = annotation.AnnotationDB(None, genome, verbose=False,
                                      sliceAttrDict=field_locations)
    nlmsa = cnestedlist.NLMSA('tmp_bed', mode='memory', pairwiseMode=True,
                              bidirectional=False)
    index = 0
    skipped = 0
    for line in bedlines:
        if not line:
            continue
        fields = line.strip().split('\t')
        orientation = 1 if len(fields) < 6 or fields[5] == '+' else -1
        #print fields, orientation
        try:
            curAnnot = annotDB.new_annotation(index, fields + [orientation])
            nlmsa.addAnnotation(curAnnot)
            index += 1
        except KeyError as e:
            print ('Skipping row without matching chromosome: %s,' +\
                    'message: %s') % (row.id, e.message)
            skipped += 1
    #annotDB.close()
    nlmsa.build()
    return annotDB, nlmsa
