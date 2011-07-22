#!/usr/bin/env python

"""
    Calculate significance of the intersection between two sets of regions.
"""

import optparse
from pygr import worldbase
from pygr.metabase import WorldbaseNotFoundError
import motility
import sys
import cPickle as pickle
import os

import hts_waterworks.utils.sampling
from hts_waterworks.utils.common import getFullGenomeName
import hts_waterworks.utils.sequence_motif as sequence_motif
from build_annotations_from_csv import makeResourceFromBed

class Logger:
    """ simple logger that obeys quiet and verbose rules """
    def __init__(self, quietMode=False):
        self.quietMode = quietMode
    def log(self, outStr):
        if not self.quietMode:
            print outStr
    def error(self, outStr):
        print >>sys.stderr, outStr

# unescape table for python operators
escapedOperators = {
    '__lt__': '<',
    '__le__': '<=',
    '__eq__': '==',
    '__ne__': '!=',
    '__gt__': '>',
    '__ge__': '>=',
    '__sq__': '\'',
    '__dq__': '"' }            


def main():
    """ Calculate significance of the intersection between two sets of regions.
        Regions may be either BED files or pygr AnnotationDB's.
    """

    parser = optparse.OptionParser("%prog [options] resource1 resource2 \n"+main.__doc__)
    parser.add_option("--genome_resource", '-g', dest="genome_resource", type="string",
                      help="""The pygr resource for the genome""")
    parser.add_option('--filter_fxn1', dest='filter_fxn1', type='string', default='',
                      help="""Use the given function as a filter on what is considered a hit from resource1.
                      available variables are seq1,annot1, edge1.  e.g.,
                      --filter_fxn1="len(seq1) > 10" """)
    parser.add_option('--filter_fxn2', dest='filter_fxn2', type='string', default='',
                      help="""Use the given function as a filter on what is considered a hit from resource2.
                      available variables are seq2,annot2, edge2.  e.g.,
                      --filter_fxn2="float(annot2.FDR) < .25" """)
    parser.add_option("--format1", dest="format1", type="string", default='BED',
                      help="""Format of resource1. One of [bed, resource, file] corresponding
                      to a single BED file, a worldbase resource ID, or a list of IDs in a file. default:%default""")
    parser.add_option("--format2", dest="format2", type="string", default='BED',
                      help="""Format of resource2. See help for format1.""")
    parser.add_option("--name1", dest="name1", type="string", default='',
                      help="""Override the name for resource1.  Default=file or resource name""")
    parser.add_option("--name2", dest="name2", type="string", default='',
                      help="""Override the name for resource2.  Default=file or resource name""")
    parser.add_option('--overlap_resource', dest='overlap_resource', type='string', default='',
                      help="""Only count regions (both res1 and res2) that overlap with this worldbase ID""")
    
    #parser.add_option("--sample_size", '-s', dest="sample_size", type="int", default=10000,
    #                  help="""Total number of background samples to check for overlap""")
    parser.add_option("--output_file", '-f', dest="output_file", type="string",
                      help="""Write the significance calculation to the given file""")
    parser.add_option("--quiet", '-q', dest="quiet", action="store_true",
                      help="""Suppress progress reports from stdout""")
    
    opts, args = parser.parse_args()
    print opts, args
    log = Logger(opts.quiet)
    if len(args) != 2:
        parser.print_help()
        log.error('Need two genomic annotations! Please specify both resource1 and resource2 ')
        sys.exit(-1)
    print opts, args
    opts.genome_resource = getFullGenomeName(opts.genome_resource)
    log.log('# Loading genome resource %s' % opts.genome_resource)
    genome = worldbase(opts.genome_resource)
    if opts.overlap_resource:
        log.log('# Loading overlap resources %s and %s' % (opts.overlap_resource, opts.overlap_resource + '_db'))
        overlapMap = worldbase(opts.overlap_resource)
        overlapDB = worldbase(opts.overlap_resource + '_db')
    
    AllRes1Names, AllRes2Names = args
    for res1Name in AllRes1Names.split(','):
        
        if len(res1Name) == 0:
            continue
        opts.format1 = opts.format1.lower()
        if opts.format1 == 'bed':
        #if os.path.exists(res1Name):
            log.log('# Building resource1 from BED file %s' % res1Name)
            res1File = open(res1Name)
            res1Table, res1DB, res1Map = makeResourceFromBed(res1File, genome)
            res1File.close()
        elif opts.format1 == 'resource':
            log.log('# Loading resource1 %s from worldbase' % res1Name)
            res1Map = worldbase(res1Name)
        elif opts.format1 == 'file':
            res1_allVars = open(res1Name).readlines()
            log.log('# List for resource1 includes %s resources' % len(res1_allVars))
        else:
            parser.print_help()
            log.error('Unrecognized format specified for resource1: %s %s should be one of [bed, resource, file]' % (opts.format1, res1Name))
            sys.exit(-1)
        
        for res2Name in AllRes2Names.split(','):
            if len(res2Name) == 0:
                continue
            if opts.format2 == 'bed':
            #if os.path.exists(res2Name):
                log.log('# Building resource2 from BED file %s' % res2Name)
                res2File = open(res2Name)
                res2Table, res2DB, res2Map = makeResourceFromBed(res2File, genome)
                res2File.close()
            elif opts.format2 == 'resource':
                log.log('# Loading resource2 %s from worldbase' % res2Name)
                res2Map = worldbase(res2Name)
                try:
                    res2DB = worldbase(res2Name + '_db')
                except:
                    log.log('No DB found for resource2 at %s' % res2Name + '_db')
                    res2DB = None
            elif opts.format1 == 'file':
                log.error('several resource iteration not implemented yet...')
            else:
                parser.print_help()
                log.error('Unrecognized format specified for resource2: %s %s should be one of [bed, resource, file]' % (opts.format2, res2Name))
                sys.exit(-1)
        
            # Unescape if filter functions have been escaped
            for key, value in escapedOperators.items():
                if opts.filter_fxn1:
                    opts.filter_fxn1 = opts.filter_fxn1.replace( key, value )
                if opts.filter_fxn2:
                    opts.filter_fxn2 = opts.filter_fxn2.replace( key, value )
            
            res1Lengths = []
            res12Intersect = 0
            res2Count = 0
            
            #res1Size, res2Size, resIntersectSize = 0,0,0
            #res2SizeInBP = 0 
            log.log('# Calculating overlap between resources... Iterating over resource 1')
            sys.stdout.flush()
            for seq1, annot1, edge1  in res1Map.edges():
                if not opts.filter_fxn1 or eval(opts.filter_fxn1):      # no filter1 or passed it
                    if not opts.overlap_resource or len(list(get_overlap_edges_seq_msa(overlapMap, seq1))) > 0:  # no overlap req'd or seq1 overlaps
                        res1Lengths.append(len(annot1))
                        for seq2, annot2, edge2 in get_overlap_edges_seq_msa(seq1, res2Map):
                            if not opts.filter_fxn2 or eval(opts.filter_fxn2):  # no filter2 or passed it
                                if not opts.overlap_resource or len(list(get_overlap_edges_seq_msa(overlapMap, seq2))) > 0:  # no overlap req'd or seq2 overlaps
                                    #res12Intersect.append(len(annot2))  # only counting the bases that actually overlap
                                    res12Intersect += 1
            # only iterate over res2 if we don't have a db resource for it or there is some filtering necessary
            if not res2DB or opts.filter_fxn2 or opts.overlap_resource:
                log.log('# Iterating over resource 2')
                sys.stdout.flush()
                for seq2, annot2, edge2 in res2Map.edges():
                    #sys.stdout.flush()
                    #print '# iterating over res2 %s...' % res2Name, opts.overlap_resource,
                    if not opts.filter_fxn2 or eval(opts.filter_fxn2):  # no filter2 or passed it
                        if not opts.overlap_resource or len(list(get_overlap_edges_seq_msa(seq2, overlapMap))) > 0:
                            # instance of res2 found
                            #if res2Size % 1000 == 0:
                            #    print res2Size,
                            res2Count += 1
            else:
                res2Count = len(res2DB)
            log.log('# Calculating enrichment...')
            fgOverlap, fgSize = res12Intersect, sum(res1Lengths)
            bgOverlap, bgSize = res2Count, sum(len(chromSeq) for chromName, chromSeq in genome.iteritems() if '_' not in chromName)
            if fgSize == 0:
                log.error('ERROR: Empty resource1 or no hits passed filtering step!')
                log.error('fgOverlap, fgSize, bgOverlap, bgSize = %s %s %s %s' % (fgOverlap, fgSize, bgOverlap, bgSize))
            else:
                zscore = sequence_motif.zscore_hypergeometric(fgOverlap, fgSize, bgOverlap, bgSize)
                pvalue = sequence_motif.pvalue_hypergeometric(fgOverlap, fgSize, bgOverlap, bgSize)
                fold_enrichment = sequence_motif.fold_enrichment(fgOverlap, fgSize, bgOverlap, bgSize)
                if opts.name1:
                    curName1 = opts.name1
                else:
                    curName1 = res1Name
                if opts.name2:
                    curName2 = opts.name2
                else:
                    curName2 = res2Name
                outstr = '\t'.join(map(str, [curName1, curName2, zscore, pvalue, fold_enrichment, fgOverlap, fgSize, bgOverlap, bgSize]))
            
        
            
            
            #print '# Now sampling %s times...' % opts.sample_size
            #sys.stdout.flush()
            #bgMatches = 0
            #genomicSamples = sampling.sample_genome(genome, res1Lengths, sampleSize=opts.sample_size, excludeRepeat=False, excludeN=False)
            #for seq in genomicSamples:
            #    for seq2, annot2, edge2 in get_overlap_edges_seq_msa(seq, res2Map):
            #        if not opts.filter_fxn2 or eval(opts.filter_fxn2):  # no filter2 or passed it
            #            if not opts.overlap_resource or len(list(get_overlap_edges_seq_msa(seq, overlapMap))) > 0:
            #                # instance of res2 found
            #                #if res2Size % 1000 == 0:
            #                #    print res2Size,
            #                bgMatches += 1
            #zscore = sequence_motif.zscore_normal(resIntersectSize, res1Size, bgMatches, opts.sample_size)
            #pvalue = sequence_motif.pvalue_hypergeometric(resIntersectSize, res1Size, bgMatches, opts.sample_size)
            #outstr = '\t'.join(map(str, [res1Name, res2Name, zscore, pvalue, resIntersectSize, res1Size, bgMatches, opts.sample_size]))
            
            #print 'Iterating over resource 2'
            #for seq2, annot2, edge2 in res2Map.edges():
            #    #sys.stdout.flush()
            #    #print '# iterating over res2 %s...' % res2Name, opts.overlap_resource,
            #    if not opts.filter_fxn2 or eval(opts.filter_fxn2):  # no filter2 or passed it
            #        if not opts.overlap_resource or len(list(get_overlap_edges_seq_msa(seq2, overlapMap))) > 0:
            #            # instance of res2 found
            #            #if res2Size % 1000 == 0:
            #            #    print res2Size,
            #            res2Size += 1
            #            res2SizeInBP += len(seq2)
            #avgRes2Size = float(res2SizeInBP) / res2Size
            #genomeSize = sum(map(len, genome.itervalues()))
            #genomeTotalPartitions = float(genomeSize) / avgRes2Size
            #print '# Calculating enrichment significance...'
            #zscore = sequence_motif.zscore_normal(resIntersectSize, res1Size, res2Size, genomeTotalPartitions)
            #pvalue = sequence_motif.pvalue_hypergeometric(resIntersectSize, res1Size, res2Size, genomeTotalPartitions)
            #outstr = '\t'.join(map(str, [zscore, pvalue, resIntersectSize, res1Size, res2Size, genomeTotalPartitions]))
            
            print outstr
            if opts.output_file:
                open(opts.output_file, 'a').write(outstr + '\n')


def get_overlap_edges_seq_msa(seq, msa):
    """ safe method for getting annotation edges
    (just skips the edge if one of the chrom's isn't in the alignment)"""
    try:
        for seq, annot, edge in msa[seq].edges():
            yield seq, annot, edge
    except KeyError,e:
        pass

if __name__ == '__main__':
    main()
