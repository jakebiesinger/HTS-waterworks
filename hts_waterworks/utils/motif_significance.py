#!/usr/bin/env python

import optparse
from pygr import worldbase
import motility
import sys
import cPickle as pickle

from hts_waterworks.utils.sampling import sample_resource, sample_genome
from hts_waterworks.utils.sequence_motif import (search_bed_file, Motif,
                        makePWMFromIUPAC, zscore_normal, pvalue_hypergeometric,
                        parseMotifsFromTransfac)
from hts_waterworks.utils.common import (getFullGenomeName, readBedLines,
                                     parseFastaLines)

def main(argv=None):
    """ Calculate significance of a motif in peaks with genomic background
    Can use restricted annotationDB, such as only promoter regions """

    parser = optparse.OptionParser("%prog [options] peaks.bed [outfile] \n"+main.__doc__)
    parser.add_option("--genome", '-g', dest="genome_resource", type="string",
                      help="""The pygr resource for the genome""")
    parser.add_option("--motif_file", '-m', dest="motif_file", type="string",
                      help="""The index file for all motifs, as a pickled dictionary, of pwm's or Motifs e.g.,
                      {"LRH_1":[[.25,.25,.1,.4],[.2,.2,.3,.3]]}""")
    parser.add_option("--consensus_file", '-c', dest="consensus_file", type="string",
                      help="""index file for consensus motifs (IUPAC format, one
                      per line in the file""")
    parser.add_option("--motif_key", '-k', dest="motif_key", type="string",
                      help="""The key for the current motif in motif_file, default=all""")
    parser.add_option('--zscore', '-z', dest='zscore', type='float', default=4.29,
                      help="""Calculate threshold score estimate from this Z-score. [default=%default]""")
    parser.add_option('--overlap_resource', dest='overlap_resource', type='string',
                      help="""Only count fg and bg that overlap with pygr resource""")
    parser.add_option('--bg_samples', dest='bg_samples', type='string',
                      help="""Pickled or Fasta file of background sequences to use instead of sampling the genome""")
    parser.add_option('--no_bg', dest='no_bg', action='store_true',
                      help="""skip sampling in the background""")
    parser.add_option('--report_region', type='string', help='Report the genomic regions of peaks with motif instances to this file')
    parser.add_option("--output_file", '-f', dest="output_file", type="string",
                      help="""Append the zscore information to the given file""")
    parser.add_option('--search_genome', action='store_true')
    if argv is None:
        argv = sys.argv[1:]
    opts, args = parser.parse_args(argv)
    if len(args) != 1:
        parser.print_help()
        print 'Specify the peaks bed file!'
        sys.exit(-1)
    if not opts.motif_file and not opts.consensus_file:
        parser.print_help()
        print 'Specify the motif file!'
        sys.exit(-1)

    updated_motifs = False
    print '# Loading resources...'
    opts.genome_resource = getFullGenomeName(opts.genome_resource)
    genome = worldbase(opts.genome_resource)
    if opts.overlap_resource:
        annotMap = worldbase(opts.overlap_resource)
        annotDB = worldbase(opts.overlap_resource + '_db')
    
    allMotifs = {}
    # load pickled dict of motifs
    if opts.motif_file:
        if opts.motif_file.endswith('.transfac'):
            allMotifs.update(parseMotifsFromTransfac(open(opts.motif_file, 'r').read()))
        else:
            allMotifs.update(pickle.load(open(opts.motif_file)))
    # create consensus dict of motifs
    if opts.consensus_file:
        with open(opts.consensus_file) as infile:
            for line in infile:
                name, consensus = line.strip().split('\t')
                allMotifs.update({name:makePWMFromIUPAC(consensus)})

    if opts.motif_key:
        allKeys = [opts.motif_key]
    else:
        allKeys = allMotifs.keys()
    
    # write a header
    if opts.output_file:
        outstr = '\t'.join(['peaks', 'motif', 'threshold_z', 'vs_bg_normal_Z',
                            'hypergeo_pvalue', 'fgMatches', 'fgSize',
                            'fgMatches/fgSize', 'bgMatches', 'bgSize'])
        open(opts.output_file, 'w').write(outstr)

    for motifKey in allKeys:
        print '# Loaded motif %s...' % motifKey
        pwm = allMotifs[motifKey]
        if isinstance(pwm, list):
            pwm = Motif(pwm)
            allMotifs[motifKey] = pwm
        if not pwm.bg_calculated():
            print '# Calculating motif background distribution...'
            pwm.calculate_background(genome)
            updated_motifs = True
        print 'motif %s: length=%s threshold=%s mean=%s sd=%s max_score=%s' % (motifKey, len(pwm), pwm.get_threshold(opts.zscore), pwm._mean, pwm._sd, pwm.max_score())
        
        if opts.search_genome and opts.report_region is not None:
            # search the genome with the motif
            print 'searching genome!'
            with open(opts.report_region, 'w') as outfile:
                for chrom in genome:
                    for match in pwm.find_in_region(genome[chrom]):
                        outstr = '{chrom}\t{start}\t{stop}\t{name}\t{score}\t{strand}\n'.format(chrom=chrom,
                                                                start=match[0], stop=match[1],
                                                                name=motifKey,score=pwm.calc_score(match[3]),
                                                                strand='+' if match[2] == 1 else '-')
                        outfile.write(outstr)
            continue
        
        allPeaks = open(args[0]).readlines()
        allPeaks = list(readBedLines(allPeaks))
        peakSizes = [stop - start for _, start, stop, _ in allPeaks]

        print '# Searching foreground sequence...'
        sys.stdout.flush()
        peakRegions = (genome[chrom][start:stop] for chrom, start, stop, _ in allPeaks)
        if opts.overlap_resource:
            # check to see if the bed line overlaps the resource
            overlappingRegions = [region for region in peakRegions \
                                        if len(annotMap[region]) > 0]
            # run a search in each of the overlapping regions
            motifInstancesInOverlap = [pwm.find_in_region(region, zscore=opts.zscore) \
                                        for region in overlappingRegions]
            fgSize = len(overlappingRegions)
            # count the number of peaks with at least one motif instance
            fgMatches = len(filter(lambda matches: len(matches) > 0, motifInstancesInOverlap))
        else:
            matchingPeaks = [region for region in peakRegions \
                                        if len(pwm.find_in_region(region, zscore=opts.zscore)) > 0]
            fgMatches = len(matchingPeaks)
            fgSize = len(allPeaks)
            
        if opts.report_region is not None:
            with open(opts.report_region, 'w') as outfile:
                outfile.writelines('%s\t%s\t%s\n' % (region.id, region.start, region.stop) for region in matchingPeaks)

        if opts.no_bg:
            outstr = '\t'.join([args[0], motifKey] + map(str, [opts.zscore, fgMatches, fgSize, 
                                                      float(fgMatches)/fgSize]))
            if opts.output_file:
                open(opts.output_file, 'a').write(outstr + '\n')
            else:
                print >>sys.stderr, outstr
        else:
            print '# Searching background sequence...'
            sys.stdout.flush()
            if opts.bg_samples:
                try:
                    bgSamples = pickle.load(open(opts.bg_samples))
                except:
                    try:
                        bgSamples = parseFastaLines(open(opts.bg_samples))
                    except:
                        raise RuntimeError("specified background samples file %s"
                                           "was niether a pickled file nor a fasta file!" %
                                           opts.bg_samples)
                
            elif opts.overlap_resource:
                bgSamples = sample_resource(annotDB, peakSizes, sampleSize=100000)
            else:
                bgSamples = sample_genome(genome, peakSizes, sampleSize=100000)
                #bgSamples = sample_genome(genome, peakSizes, sampleSize=100)
            bgSize = 0
            bgMatches = 0
            for region in bgSamples:
                bgSize += 1
                if len(pwm.find_in_region(region, zscore=opts.zscore)) > 0:
                    bgMatches += 1
    
            #calculate significance of foreground vs. background
            zscore = zscore_normal(fgMatches, fgSize, bgMatches, bgSize)
            pvalue = pvalue_hypergeometric(fgMatches, fgSize, bgMatches, bgSize)
            outstr = '\t'.join([args[0], motifKey] + map(str, ['thesh_z='+str(opts.zscore), zscore, pvalue, fgMatches, fgSize, float(fgMatches)/fgSize,bgMatches, bgSize]))
            if opts.output_file:
                open(opts.output_file, 'a').write(outstr + '\n')
            else:
                print >>sys.stderr, outstr
    if updated_motifs:
        print '# Saving motif info back to %s' % opts.motif_file
        pickle.dump(allMotifs, open(opts.motif_file, 'wb'))


if __name__ == '__main__':
    main()
