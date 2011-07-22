'''  bedOverlapRandomShuffle.py
        * Created by Jacob Biesinger
        * This version: Created October 2010
'''

import optparse
from collections import defaultdict
import sys
import random
import scipy
from bisect import insort

from hts_waterworks.utils.common import getGenome, readBedLines

def main():
    ''' Calculate the chance of observing a certain number of overlaps between two
        set of genomic regions.  Significance is estimated by randomly shuffling the
        positions (not changing the lengths or chromosomes) of one of the samples and reporting the number
        of overlapping sites in the shuffled sets.  '''
    usage = "%prog [options] bedFile1 bedFile2 \n" + main.__doc__
    parser = optparse.OptionParser(usage)
    parser.add_option('--genome', '-g', dest='genome', type='string', default=None,
                      help='The genome name the bed files come from, i.e, mm9 or hg19')
    parser.add_option('--num_shuffles', '-n', dest='num_shuffles', type='int', default=10000,
                      help='Number of times to shuffle bedFile1. default=%default')
    parser.add_option('--disjoint', '-d', action='store_true',
                      help='Make sure that there is no overlap in shuffled regions')
    parser.add_option('--quiet', '-q', action='store_true',
                      help='report only the overlap number (no messages)')
    parser.add_option('--unique_out', '-u', dest='unique_out', type='string', default=None,
                      help='print non-overlapping regions from bedfile1 to this file')
    parser.add_option('--report_col1', '-1', dest='report_col1', type='int', default=None,
                      help='bed column to use when reporting the overlap type. default:None')
    parser.add_option('--report_col2', '-2', dest='report_col2', type='int', default=None,
                      help='bed column to use when reporting the overlap type. default:None')
    parser.add_option('--file_report', '-f', dest='file_out', type='string', default=None,
                      help='where to file the overlap report')
    opts, args = parser.parse_args()
    if opts.genome is None:
        parser.print_help()
        print >>sys.stderr, 'You must specify a genome!'
        sys.exit(-1)
    if opts.num_shuffles < 0:
        parser.print_help()
        print >>sys.stderr, 'Must have a positive or 0 number of shuffles!'
        sys.exit(-1)
    genome = getGenome(opts.genome)
    chromSizes = dict((chrom, len(seq)) for chrom, seq in genome.iteritems())
    bedfile1 = open(args[0], 'r')
    bedfile2 = open(args[1], 'r')
    bedlines1 = sorted(readBedLines(bedfile1, dataOnly=False))
    bedlines2 = sorted(readBedLines(bedfile2, dataOnly=False))
    if opts.report_col1 is None:
        opts.report_col1 = 'none'
    if opts.report_col2 is None:
        opts.report_col2 = 'none'
    
    if not opts.quiet:
        print 'Original data:\t%s in %s\t%s in %s\t' % (args[0],len(bedlines1), args[1],len(bedlines2)),
    if opts.unique_out:
        originalOverlapCount, uniqueBed1 = getBedOverlap(bedlines1, bedlines2, alreadySorted=True, reportUnique=True, featureColumn1=opts.report_col1, featureColumn2=opts.report_col2)
        with open(opts.unique_out, 'w') as outfile:
            outfile.writelines('\n'.join('\t'.join(map(str, bedFields)) for bedFields in uniqueBed1))
    else:
        originalOverlapCount = getBedOverlap(bedlines1, bedlines2, alreadySorted=True, featureColumn1=opts.report_col1, featureColumn2=opts.report_col2)
    if not opts.quiet:
        print 'with %s overlaps or %s unique to bedfile1' % (originalOverlapCount, len(bedlines1) - originalOverlapCount)
    else:
        sys.stdout.write('\t' + str(originalOverlapCount))
    
    if opts.file_out:
        with open(opts.file_out, 'a') as outfile:
            print >> outfile, '\t'.join([args[0], str(len(bedlines1)), args[1],str(len(bedlines2)),
                                        'overlap: %s' % originalOverlapCount, 
                                        'unique to 1: %s' % (len(uniqueBed1) if opts.unique_out else len(bedlines1) - originalOverlapCount)])
    
    
    if opts.num_shuffles > 0:
        randOverlaps = [-1] * opts.num_shuffles  # preallocate
        
        print 'Generating %s random shuffles...' % opts.num_shuffles,
        for i in xrange(opts.num_shuffles):
            if i % 1000 == 0:
                print i,
                sys.stdout.flush()
            shuffledBeds1 = sorted(generateShuffledBed(bedlines1, chromSizes, disjoint=opts.disjoint))
            overlapCount = getBedOverlap(shuffledBeds1, bedlines2, alreadySorted=True)
            randOverlaps[i] = overlapCount
        print
        randomBetterCount = len(filter(lambda randVal: randVal >= originalOverlapCount, randOverlaps))
        randNumDistinctVals = len(set(randOverlaps))
        randHist, bins = scipy.histogram(randOverlaps, bins=min(randNumDistinctVals, 15))
        print 'Random overlap distribution is: \nbinCounts:\t%s\nbinEdges:  %s' % (randHist, bins)
        print 'Random shuffle:\t%s with at least as many overlaps, pvalue %s %s' % (randomBetterCount,
                                                '<' if randomBetterCount==0 else '=',
                                                max(1./opts.num_shuffles, float(randomBetterCount)/opts.num_shuffles))
        print 'Random mean:\t%s\tstdev:%s' % (scipy.mean(randOverlaps), scipy.std(randOverlaps))

def generateShuffledBed(bedLines, chromSizes, disjoint=False, numPerLine=1):
    '''Shuffle the location of the given regions without touching the length or chromosome.
        optionally make sure the samples are disjoint (don't overlap with each other).
        ensuring disjointness may take a while since each new sample has to check for
        overlap with previous samples.
        This function generates numPerLine random samples for each bedfile line 
        '''
    samplesSoFar = defaultdict(list)
    for line in bedLines:
        validSample = False
        chrom, start, stop = line[:3]
        chromLength = chromSizes[chrom]
        featureLength = abs(stop - start)
        for i in xrange(numPerLine):
            while not validSample:
                maxStart = chromLength - featureLength
                newStart = random.randrange(0, maxStart)  # bed files are 1-based
                newRegion = chrom, newStart, newStart + featureLength
                if disjoint:
                    if featureOverlapsSet(newRegion, samplesSoFar):
                        continue
                    else:
                        validSample = True
                        insort(samplesSoFar[chrom],newRegion)
                else:
                    validSample = True
            yield newRegion

def getBedOverlap(bed1, bed2, featureColumn1='none', featureColumn2='none', alreadySorted=False, reportUnique=False):
    '''Check the genomic overlap of bed1 with bed2. Input will be sorted by
    (chrom, start) unless sorted=True (meaning the data is already sorted).
    Optionally, include the index of a "feature" column, i.e., track which features
    overlap with which.  This latter option will return a dict(dict(overlap count)) with features as keys.
    reportUnique will report the lines of bed1 that don't overlap with bed2 as a second return value'''
    uniqueInBed1 = []
    if not alreadySorted:
        sortedBed1 = iter(sorted(bed1))
        sortedBed2 = iter(sorted(bed2))
    else:
        # still have to make a copy
        sortedBed1 = iter(bed1)
        sortedBed2 = iter(bed2)
    totalOverlap = defaultdict(lambda: defaultdict(int))  # dict(dict(int))
    # iterate over data, popping either region1 or region2
    try:  # make sure that they aren't empty
        region1 = sortedBed1.next()
        region2 = sortedBed2.next()
    except StopIteration:
        raise RuntimeError("There were no lines in one of the bed regions in getBedOverlap")
    try:  # go until one of the iterators is empty
        while True:
            reportVal1 = 'none' if featureColumn1 == 'none' else region1[featureColumn1]
            reportVal2 = 'none' if featureColumn2 == 'none' else region2[featureColumn2]
            overlap = checkOverlap(region1, region2)
            if overlap < 0:
                if reportUnique:
                    uniqueInBed1.append(region1)
                region1 = sortedBed1.next()  # advance region1
            elif overlap > 0:
                region2 = sortedBed2.next()  # advance region2
            else:
                totalOverlap[reportVal1][reportVal2] += 1
                region1 = sortedBed1.next()  # advance region1
    except StopIteration:
        if reportUnique:
            uniqueInBed1.extend(list(sortedBed1))  # any from bed1 that haven't been checked
    
    if featureColumn1 is 'none' and featureColumn2 is 'none':
        overlapReport = totalOverlap['none']['none']  # don't care about columns-- just get the number
    elif featureColumn2 is 'none':
        overlapReport = dict( (feature1, totalOverlap[feature1]['none']) for feature1 in totalOverlap)
    else:
        overlapReport = totalOverlap
    
    if reportUnique:
        return overlapReport, uniqueInBed1
    else:
        return overlapReport

def featureOverlapsSet(newFeature, featureSet):
    # can't do binary search since long features may overlap but have distant start/stops
    for region in featureSet[newFeature[0]]:
        direction = checkOverlap(newFeature, region)
        if direction == 0:
            return True
        elif direction == 1:  # this region is after all the reads in the featureSet
            return False
    return False

def checkOverlap(bedRegion1, bedRegion2):
    '''check if the given bed regions overlap, returning:
            0 if they overlap,
            -1 if region1 is before region2
            1 if region1 is after region2  '''
    chrom1, start1, stop1 = bedRegion1[0], int(bedRegion1[1]), int(bedRegion1[2])
    chrom2, start2, stop2 = bedRegion2[0], int(bedRegion2[1]), int(bedRegion2[2])
    if chrom1 < chrom2:
        return -1
    elif chrom1 > chrom2:
        return 1
    else:  # on same chrom
        if start1 < start2:
            if stop1 <= start2:  #should this be < ? no-- from FAQ: chromStart=0, chromEnd=100, spans bases numbered 0-99
                return -1
            else:
                return 0  # region1 overlaps or completely spans region2
        elif start1 < stop2:
            return 0  # region 1 starts after region2 start, but before region2 stop
        else:
            return 1  # region1 start is after region2 end

if __name__ == '__main__':
    main()