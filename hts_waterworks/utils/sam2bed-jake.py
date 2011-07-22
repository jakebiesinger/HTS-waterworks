#!/usr/bin/env python
# encoding: utf-8
"""
sam2bed.py

Created by Aaron Quinlan on 2009-08-27.
Copyright (c) 2009 Aaron R. Quinlan. All rights reserved.

Modified by Jacob Biesinger 21 July 2011
"""

import sys
import optparse
import re


def main():
    """ Convert a SAM file to BED """

    parser = optparse.OptionParser(usage='%prog inputfile outputfile\n' + main.__doc__)
    opts, args = parser.parse_args()

    if len(args) != 2:
        parser.print_help()
        print 'Please specify an inputfile and an outputfile'
        sys.exit(-1)

    filein = open(args[0], 'r')
    fileout = open(args[1], 'w')
    for line in filein:
        if line[0] == '@':
            continue # skip comment characters
        fields = line.split('\t')
        #print fields
        samFlag = int(fields[1])
        if not (samFlag & 0x0004):  # Only create a BED entry if the read was aligned
            chrom, name, strand = fields[2], fields[0], getStrand(samFlag)
            start, end = str(int(fields[3])-1), str(int(fields[3]) + len(fields[9]) - 1)

            # Let's use the edit distance as the BED score.
            # Clearly alternatives exist.
            editPattern = re.compile('NM\:i\:(\d+)')
            editDistance = editPattern.findall(fields[12])
            editDistance = str(editDistance[0]) if editDistance else '0' # set to 0 if no field matches

            # Write out the BED entry
            fileout.write('\t'.join([chrom, start, end, name, editDistance, strand]) + '\n')

def getStrand(samFlag):
    return '-' if samFlag & 0x10 else '+'


if __name__ == "__main__":
    main()

