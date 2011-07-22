
import optparse
import re
import sys
from pygr import worldbase
from itertools import ifilter

from hts_waterworks.utils.common import getGenome, bedCommentFilter, maskRepeats

def main(argv=None):
    """ Output sequence of given genomic locations (ala BED or region) in FASTA format
    """
    if argv is None:
        argv = sys.argv[1:]

    parser = optparse.OptionParser(usage='usage: %prog [options] infile outfile\n' + main.__doc__)
    parser.add_option('-g', '--genome', dest='genome', type='string', default='mm9')
    opts,args = parser.parse_args(argv)
    
    if len(args) != 2:
        parser.print_help()
        print 'Please provide both an input and an output file!'
        sys.exit(-1)
    
    genome = getGenome(opts.genome)
    
    filein = open(args[0], 'r')
    fileout = open(args[1], 'w')  
    bedSeqs = getFastaFromBed(filein, genome)
    for fastaLine in bedSeqs:
        fileout.write(fastaLine)
    fileout.close()


def getFastaFromBed(fileLines, genome, maskRepeats=False):
    bedLines = ifilter(bedCommentFilter, fileLines)
    for index,line in enumerate(bedLines):
        fields = line.split('\t')
        if len(fields) >= 6 and re.search('[-|+]', fields[5]):
            chrom, start, stop, _, _, strand = fields[:6]
        else:
            chrom, start, stop = fields[:3]
            strand = '+'
        start = max(0,int(start))
        stop = max(int(stop), start)
        sequence = genome[chrom][start:stop]
        if strand == '-':
            sequence = -sequence
        seqStr = str(sequence)
        if maskRepeats:
            seqStr = maskRepeats(seqStr)
        fastaText = '>%s:%s-%s site_%s \n%s\n' % (chrom,start,stop,index,seqStr)
        yield fastaText

if __name__ == '__main__':
    main()

