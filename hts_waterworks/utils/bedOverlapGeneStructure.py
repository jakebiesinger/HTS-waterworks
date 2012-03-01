'''  bedOverlapRandomShuffle.py
        * Created by Jacob Biesinger
        * This version: Created October 2010
'''

import optparse
from collections import defaultdict
import sys
import os
from bisect import bisect
import itertools
import shlex
import subprocess
import tempfile
import scipy
from scipy.stats.kde import gaussian_kde
import matplotlib
from matplotlib import pyplot

from bedOverlapRandomShuffle import checkOverlap, generateShuffledBed

def main(argv=None):
    """ Report the overlap between a bed file and a set of refseq-formatted genes.
    The refseqGenes tab-separated file format may include a "bin" as the first
    column, and can have comments and extra columns, i.e.,
    [bin]   name    chrom   strand  txStart txEnd   cdsStart    cdsEnd  exonCount   exonStarts  exonEnds [...]
    
    For each region in bedFile, its placement on the genes is reported once, with
    priority in the following order:
        promoter, 5'UTR, 3'UTR, downstream, Exon, Intron, Intergenic (i.e., none of the above)
    """
    usage = "%prog [options] bedFile refseq_Genes [chrom.sizes] \n" + main.__doc__
    parser = optparse.OptionParser(usage)
    parser.add_option('--promoter_size', type='int', default=2000,
                      help='Size in bp of the upstream promoter region. default=%default')
    parser.add_option('--promoter_extend',type='int', default=0,
                      help='Extend the promoter regions past the Tx Start Site (into the 5utr or exon) by X bp. default=%default')
    parser.add_option('--downstream_size', type='int', default=500,
                      help='Size in bp of the downstream region. default=%default')
    parser.add_option('--downstream_extend', type='int', default=0,
                      help='Extend the downstream regions past the Tx Stop Site (into the 3utr or exon) by X bp. default=%default')
    parser.add_option('--random_shuffles', type='int', default=0,
                      help='Shuffle the location of the bed input X times, and report the mean overlaps. default=No shuffling')
    parser.add_option('--plot_dist_to_txStart', type='string', default=None,
                      help='Plot a histogram of the distances from the regions to the nearest txStart, possibly also plotting random shuffles.')
    parser.add_option('--pie_output', type='string', default=None,
                      help='Filename for a gene structure pie graph')
    parser.add_option('--bar_output', type='string', default=None,
                      help='Filename for a bar chart')
    parser.add_option('--genome', '-g', type='string', default=None,
                      help='Genome is specified to constrain sizes in random shuffling.  Not needed if there will be no shuffling.')
    parser.add_option('--output_file', '-f', type='string', default=None,
                      help='Save the overlap structure to the given file.')
    parser.add_option('--names_in_order', type='string', default="['promoter', 'utr5', 'exon', 'intron', 'utr3', 'down', 'intergenic']",
                      help='the order in which to count and report overlaps.')
    parser.add_option('--tss_hist', action='store_true', default=False,
                      help='plot the tss as a histogram instead of a density plot')
    if not argv:
        argv = sys.argv[1:]
    opts, args = parser.parse_args(argv)
    print argv, args
    if not (2 <= len(args) <= 3):  # two or three input args
        parser.error('Please specify both a bed file and a gene list as input!')
    if opts.random_shuffles > 0 and opts.genome is None:
        parser.error('If you want to do random shuffles, you need to also specify a genome!')
    old_size = matplotlib.rcParams['font.size']
    matplotlib.rcParams['font.size'] = 6
    bedfile = args[0]
    gene_base = args[1]
    namesInOrder = eval(opts.names_in_order)
    gene_files = []
    for name in namesInOrder:
        extension = '.' + name
        if name == 'promoter':
            extension += '%s_ext%s' % (opts.promoter_size, opts.promoter_extend)
        elif name == 'down':
            extension += '%s_ext%s' % (opts.downstream_size, opts.downstream_extend)
        if name != 'intergenic':
            gene_files.append(gene_base + extension)
    print gene_files
    #gene_files = [gene_base + extension for extension in ['.promoter%s_ext%s' % (opts.promoter_size, opts.promoter_extend),
    #                                                          '.utr5',
    #                                                          '.exon',
    #                                                          '.intron',
    #                                                          '.utr3',
    #                                                          '.down%s_ext%s' % (opts.downstream_size, opts.downstream_extend)]]
    
    for infile in gene_files:
        if not os.path.exists(infile):
            parser.error('Please run makeGeneStructure with your options first. Missing file %s' % infile)
    if not os.path.exists(gene_base + '.tss'):
        parser.error('Please run makeGeneStructure with your options first. Missing file %s' % gene_base + '.tss')
    if opts.random_shuffles > 0:
        if len(args) != 3:
            parser.error('Please specify the chrom sizes file for random shuffles')
        if not os.path.exists(args[2]):
            parser.error('Please specify a valid file for chrom sizes (no such file: %s)' % args[2])
        
    print '# Checking gene overlap...'
    regionCounts = getCategoriesForBed(bedfile, gene_files)
    counts = dict(zip(namesInOrder, regionCounts))
    if opts.plot_dist_to_txStart is not None:
        dists = list(getNearestTss(bedfile, gene_base+'.tss'))
        #print 'tss distances are:', dists
    
    print '%s total regions, distibuted in the genome as:\n%s' % (sum(regionCounts), counts)
    if opts.output_file is not None:
        with open(opts.output_file, 'a') as outfile:
            outfile.write('\t'.join([bedfile, str(sum(regionCounts)), str(counts)]) + '\n')
    
    expName = bedfile
    if opts.pie_output is not None:
        print '# Saving pie chart to %s' % opts.pie_output
        pyplot.figure(figsize=(8,8))
        pyplot.pie(regionCounts[::-1], labels=namesInOrder[::-1], colors=('r','y','g','c','chartreuse', 'b','m'))#, autopct='%.1f')
        #pyplot.pie(regionCounts[::-1], labels=namesInOrder[::-1])
        pyplot.title('%s %s regions overlap with Refseq Genes' % (expName, sum(regionCounts)))
        pyplot.savefig(opts.pie_output)
        pyplot.close()
    if opts.bar_output is not None:
        print '# Saving bar chart to %s' % opts.bar_output
        pyplot.plot((0,0),(0,len(regionCounts)))
        pyplot.title('%s %s regions overlap with Refseq Genes' % (expName, sum(regionCounts)))
        pyplot.bar(range(len(regionCounts)), regionCounts)
        pyplot.xticks(scipy.arange(len(regionCounts))+.4, namesInOrder, rotation=17)
        pyplot.savefig(opts.bar_output)
        pyplot.close()
    if opts.plot_dist_to_txStart and opts.random_shuffles <= 0:
        # only plot the TSS histogram
        if opts.tss_hist:
            pyplot.hist(dists, bins=30, normed=True)
        else:
            t_kde = gaussian_kde(dists)
            inds = scipy.linspace(-50000, 50000, 100)
            t_vals = t_kde.evaluate(inds)
            pyplot.plot(inds, t_vals, label=expName)
        pyplot.xlim([-50000, 50000])
        pyplot.title('Distance to TSS for %s regions from\n%s' % (sum(regionCounts), expName))
        pyplot.savefig(opts.plot_dist_to_txStart)
        pyplot.close()
    
    if opts.random_shuffles > 0:
        print '# Generating %s random shuffles...' % opts.random_shuffles
        with tempfile.NamedTemporaryFile(delete=False) as tmp:
            cmd = 'for i in `seq 1 %s`; do shuffleBed -i %s -g %s -chrom >> %s || break; done'
            cmd = cmd % (opts.random_shuffles, args[0], args[2], tmp.name)
            subprocess.check_call(cmd, shell=True)
            print '# Getting categories for shuffles...'
            randomCounts = getCategoriesForBed(tmp.name, gene_files)
            print randomCounts
            if opts.plot_dist_to_txStart is not None:
                print '# Getting nearest TSS for shuffles...'
                randomDists = list(getNearestTss(tmp.name, gene_base+'.tss'))
                #print 'random tss distances are:', randomDists
            
        averageRandCounts = dict(zip(namesInOrder, [float(count) / opts.random_shuffles for count in randomCounts]))
        print '%s total shuffles, average distribution in genes is:\n%s' % (opts.random_shuffles, averageRandCounts)
        if opts.output_file is not None:
            with open(opts.output_file, 'a') as outfile:
                outfile.write('\t'.join([bedfile + '_random', str(sum(averageRandCounts.values())), str(averageRandCounts)]) + '\n')
        randRegionCounts = [averageRandCounts[name] for name in namesInOrder]
        if opts.pie_output is not None:
            print 'saving random pie chart'
            pyplot.figure(figsize=(8,8))
            pyplot.axes()
            pyplot.pie(randRegionCounts[::-1], labels=namesInOrder[::-1])#, autopct='%.1f')#, autopct=lambda pct: str(round(pct, 1)))
            pyplot.title('Random regions overlap with Refseq Genes')
            pyplot.savefig(os.path.splitext(opts.pie_output)[0] + '_random.png')
            pyplot.close()
        if opts.bar_output is not None:
            print 'saving random bar chart'
            randBars = pyplot.bar(scipy.arange(len(randRegionCounts)), randRegionCounts, .4, color='y')
            realBars = pyplot.bar(scipy.arange(len(randRegionCounts))+.4, regionCounts, .4, color='r')
            pyplot.legend((randBars[0], realBars[0]), ('Random', expName))
            pyplot.title('Random regions overlap with Refseq Genes')
            pyplot.xticks(scipy.arange(len(randRegionCounts))+.2, namesInOrder, rotation=17)
            pyplot.savefig(os.path.splitext(opts.bar_output)[0] + '_random.png')
            pyplot.close()
        if opts.bar_output is not None:
            print 'saving enrichment chart'
            enrichment = scipy.array(regionCounts) / scipy.array(randRegionCounts)
            if len(enrichment<1) > 0:
                enrichment[enrichment < 1.] = -1./enrichment[enrichment<1]  #fractional enrichment becomes negative
            pyplot.barh(range(len(randRegionCounts)), enrichment)
            pyplot.plot((0,0),(0,len(randRegionCounts)))
            pyplot.plot((1,0),(1,len(randRegionCounts)))
            pyplot.plot((-1,0),(-1,len(randRegionCounts)))
            pyplot.title('Peak Enrichment in Refseq Genes vs. random')
            pyplot.yticks(scipy.arange(len(randRegionCounts))+.4, namesInOrder)
            pyplot.axes().xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=10, integer=True))  # more granular than default
            pyplot.xlabel('Fold Enrichment')
            pyplot.savefig(os.path.splitext(opts.bar_output)[0] + '_randomEnrichment.png')
            pyplot.close()
        if opts.plot_dist_to_txStart:
            # only plot the TSS histogram
            if len(dists) == 0 or len(randomDists) == 0:
                raise RuntimeError("dists is empty!", dists, len(randomDists))
            pyplot.hold(True)
            if opts.tss_hist:
                pyplot.hist(dists, bins=200, normed=True, color='black', alpha=.3)
            else:
                t_kde = gaussian_kde(dists)
                inds = scipy.linspace(-50000, 50000, 100)
                t_vals = t_kde.evaluate(inds)
                pyplot.plot(inds, t_vals, label=expName, color='black')
            pyplot.hold(True)
            if opts.tss_hist:
                pyplot.hist(randomDists, bins=200, normed=True, color='black', histtype='step')
            else:
                c_kde = gaussian_kde(randomDists)
                inds = scipy.linspace(-50000, 50000, 100)
                c_vals = c_kde.evaluate(inds)
                pyplot.plot(inds, c_vals, label='Random', color='blue')
            pyplot.xlim([-50000, 50000])
            pyplot.title('Distance to TSS for %s regions from\n%s' % (sum(regionCounts), expName))
            pyplot.xlabel('Distance to TSS in bp')
            pyplot.savefig(opts.plot_dist_to_txStart)
            pyplot.close()
    matplotlib.rcParams['font.size'] = old_size
            
def getCategoriesForBed(bedFile, geneFiles):
    '''Check for overlap of bedFile against geneFiles by progressively removing overlapping regions'''
    # get total line count in input
    total_count = sum(1 for i in open(bedFile))
    # build up a piped command, removing any overlaps from input to next job and counting the totals
    # counts are output to stderr, non-overlapping regions are output to next job
    cmd = shlex.split('intersectBed -v -a stdin -b stdout')
    countAndPass = '''import sys; count = 0;
outfile = open('%s', 'w')
for line in sys.stdin:
    count += 1
    print line,
    outfile.write(line)
print >> sys.stderr, count
outfile.close()
'''
    subsets = ' | intersectBed -v -a stdin -b %s | ' + ('python -c "%s" ' % countAndPass) 
    cmd = 'intersectBed -v -a %s -b %s ' % (bedFile, geneFiles[0])
    cmd += ' | python -c "%s" ' % (countAndPass % (bedFile + '.' + geneFiles[0]))
    for infile in geneFiles[1:]:
        cmd += subsets % (infile, bedFile + '.' + infile)
    cmd += ' > intergenic.sites'
    #print cmd
    p = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
    stderr_out = p.communicate()[1]
    remainingCounts = [line for line in stderr_out.split() if line]
    print remainingCounts
    remainingCounts = [int(line) for line in stderr_out.split() if line]
    origLeft = total_count
    overlapCounts = []
    for c in remainingCounts:
        overlap = origLeft - c
        overlapCounts.append(overlap)
        origLeft -= overlap
    
    return overlapCounts + [origLeft]

def getNearestTss(in_peaks, in_tss):
    "find the nearest TSS for each entry in bedname using the closestBed program"
    all_dists = []
    with tempfile.NamedTemporaryFile() as tmp_output:
        cmd = 'closestBed -a %s -b %s -t first -d > %s' % (in_peaks, in_tss, tmp_output.name)
        print cmd
        offset = len(open(in_peaks).readline().strip().split('\t'))
        p = subprocess.Popen(cmd, shell=True)
        p.communicate()
        with open(tmp_output.name) as infile:
            for line in infile:
                if not line:
                    continue
                fields = line.strip().split('\t')
                p_chrom, p_start, p_end = fields[0], int(fields[1]), int(fields[2])
                g_chrom, g_start, g_end = fields[offset], int(fields[offset+1]), int(fields[offset+2])
                if p_chrom == 'none' or g_chrom == 'none':
                    continue
                dist = min(g_start - p_start, g_start - p_end,
                           g_end - p_start, g_end - p_end, key=abs)  # key=abs => closest to 0
                all_dists.append(dist)
    return all_dists

if __name__ == '__main__':
    main()