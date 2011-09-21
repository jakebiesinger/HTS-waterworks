"""annotation.py
    Module to annotate feature locations with respect to genes and check for
    correlation with gene expression datasets

"""

#  Current Version: 0.1-13-gb6710ca
#  Last Modified: 2011-09-20 20:35

import shlex

from ruffus import (transform, suffix, files, follows, regex, split, collate,
                    add_inputs)
from ruffus.task import active_if
import pyper

from hts_waterworks.utils.ruffus_utils import sys_call, touch
from hts_waterworks.utils.common import parseGTF, flatten
import hts_waterworks.utils.makeGeneStructure as makeGeneStructure
import hts_waterworks.utils.bedOverlapGeneStructure as bedOverlapGeneStructure
from hts_waterworks.bootstrap import cfg, get_chrom_sizes
import hts_waterworks.call_peaks as call_peaks

@transform('*.gtfgenes', suffix('.gtfgenes'), '_genes')
def convert_gtf_genes_to_bed(in_gtf, out_bed):
    """convert the gene GTF files into "standard" UCSC format
    from multi-line GTF to:
    bin    name    chrom   strand  txStart txEnd   cdsStart    cdsEnd  exonCount   exonStarts  exonEnds [...]
    """
    all_genes = {}
    with open(in_gtf) as infile:
        for line in parseGTF(infile):
            chrom, src, feature, start, stop, score, strand, frame, attrs = line
            start, stop = int(start) - 1, int(stop) - 1  # shift to 0-based
            name = attrs['gene_id']
            if name not in all_genes:
                all_genes[name] = {}
                all_genes[name]['strand'] = '-'
                all_genes[name]['chrom'] = chrom
            if feature not in all_genes[name]:
                all_genes[name][feature] = []
            all_genes[name][feature].append((start, stop))
    
    with open(out_bed, 'w') as outfile:
        for gene_id, features in all_genes.iteritems():
            chrom, strand = features['chrom'], features['strand']
            all_posns = flatten([v for k, v in features.items() if k in
                                ['exon', 'CDS', 'start_codon', 'stop_codon']])
            txStart = min(all_posns)
            txEnd = max(all_posns)
            try:  # in xenTro2 JGI genes, 2987 don't have a start/stop_codon
                cdsStart = min(flatten([features['start_codon'],
                                               features['stop_codon']]))
                cdsEnd = max(flatten([features['start_codon'],
                                             features['stop_codon']]))
            except (KeyError, ValueError) as e:  # no start_codon or stop_codon
                # hack to indicate this is a non-coding gene
                # see makeGeneStructure for how this is set up
                cdsStart, cdsEnd = -1, -1
            exonCount = len(features['CDS'])
            exonStarts = ','.join(map(str, sorted(int(s[0])
                                            for s in features['CDS']))) + ','
            exonEnds = ','.join(map(str, sorted(int(s[1])
                                            for s in features['CDS']))) + ','
            # just make frames 0
            exonFrames = ','.join(['0'] * len(features['CDS'])) + ','
            outfile.write('\t'.join(map(str, [0, gene_id, chrom, strand,
                                    txStart, txEnd, cdsStart, cdsEnd,
                                    exonCount, exonStarts,
                                    exonEnds, exonFrames, gene_id,
                                    'none', 'none', ])) + '\n')


@active_if(cfg.getboolean('genes','download_refseq'))
@files(None, '%s.refseq_genes' % cfg.get('DEFAULT', 'genome'))
def get_refseq_genes(_, out_genes):
    """Download refseq genes from UCSC and reformat as BED"""
    url = 'http://hgdownload.cse.ucsc.edu/goldenPath/%s/database/refGene.txt.gz'
    url = url % cfg.get('DEFAULT', 'genome')
    sys_call('wget -N -P . %s' % url)
    sys_call('gunzip -f refGene.txt.gz')
    sys_call('mv refGene.txt %s' % out_genes)

@follows(get_refseq_genes, convert_gtf_genes_to_bed)
@transform('*_genes', suffix('_genes'), '_genes.all')
def refseq_genes_to_bed(in_genes, out_bed):
    """convert refseq genes file to BED format"""
    with open(in_genes) as infile:
        with open(out_bed, 'w') as outfile:
            for line in infile:
                fields = line.strip().split('\t')
                chrom, start, stop = fields[2], fields[4], fields[5]
                name, strand = fields[1], fields[3]
                outfile.write('\t'.join([chrom, start, stop, name,
                                         '0', strand]) + '\n')

@follows(get_refseq_genes, convert_gtf_genes_to_bed)
@split('*_genes', regex(r'(.*)_genes$'), r'\1_genes.*')
def refseq_genes_to_regions(in_genes, out_pattern):
    """make regions (promoter, downstream, 5UTR, etc) from refseq_genes"""
    args = shlex.split('''%s --promoter_size=%s --promoter_extend=%s
                          --downstream_size=%s --downstream_extend=%s
                          --with_gene_name''' % (
                            in_genes,
                            cfg.get('genes', 'promoter_size'),
                            cfg.get('genes', 'promoter_extend'),
                            cfg.get('genes', 'downstream_size'),
                            cfg.get('genes', 'downstream_extend')))
    makeGeneStructure.main(args)

@follows(get_refseq_genes, convert_gtf_genes_to_bed)
@collate(call_peaks.all_peak_caller_functions + ['*.custom.peaks'],
         regex(r'(.+)\.treat\.(.+)\.peaks'), 
         #add_inputs('*_genes.all'), r'\1.treat.\2.peaks.nearby.genes')
         add_inputs('*_genes.tss'), r'\1.treat.\2.peaks.nearby.genes')
def find_nearby_genes(in_files, out_genes):
    """report which genes are within a certain distance of a peak"""
    in_peaks, in_genes = in_files[0]
    tmp_output = in_peaks+'.nearby_genes_tmp'
    cmd = 'closestBed -a %s -b %s -t first -d > %s' % (in_peaks,
                                                       in_genes, tmp_output)
    sys_call(cmd)
    with open(tmp_output) as infile:
        with open(out_genes, 'w') as outfile:
            for line in infile:
                if not line:
                    continue
                fields = line.split('\t')
                dist = int(fields[12])
                if dist <= cfg.getint('genes', 'nearby_genes_max_dist'):
                    outfile.write(line)


#@jobs_limit(cfg.getint('DEFAULT', 'max_throttled_jobs'), 'throttled')
@follows(refseq_genes_to_regions, convert_gtf_genes_to_bed)
@split(call_peaks.all_peak_caller_functions + ['*.custom.peaks'] + ['*.merged.mapped_reads'],
         #regex(r'(.*).peaks'),
         regex(r'(.*)'),
         add_inputs('*_genes', get_chrom_sizes),
         [r'\1.genes.overlap',
          r'\1.genes.overlap.gene_structure_pie.png',
          r'\1.genes.overlap.gene_structure_pie_random.png',
          r'\1.genes.overlap.dist_to_txStart.png'])
def gene_overlap(in_files, out_files):
    """Check the overlap of peaks with a set of genes"""
    in_bed, in_gene_info, in_chrom_sizes = in_files
    out_gene_overlap = out_files[0]
    args = shlex.split('''%s %s %s --genome=%s --random_shuffles=0
                       --pie_output=%s.gene_structure_pie.png
                       --plot_dist_to_txStart=%s.dist_to_txStart.png
                       --bar_output=%s.fold_enrichment_vs_control.png
                       --output_file=%s
                       --promoter_size=%s --promoter_extend=%s
                       --downstream_size=%s --downstream_extend=%s
                       --names_in_order="%s"
                       ''' % (in_bed, in_gene_info, in_chrom_sizes,
                                cfg.get('DEFAULT', 'worldbase_genome'),
                                out_gene_overlap,out_gene_overlap,
                                out_gene_overlap, out_gene_overlap,
                                cfg.get('genes', 'promoter_size'),
                                cfg.get('genes', 'promoter_extend'),
                                cfg.get('genes', 'downstream_size'),
                                cfg.get('genes', 'downstream_extend'),
                                cfg.get('genes', 'gene_overlap_order')
                                ))
    #print args
    bedOverlapGeneStructure.main(args)







# expression data
@follows(find_nearby_genes)
@split('*.gene.expression', regex('(.*)'),
       add_inputs(refseq_genes_to_bed, find_nearby_genes),
       r'\1.with_peaks.*.ks_data')
def make_expression_ks(in_files, out_pattern):
    """Make a KS-test ready data file from gene expression and peak data"""
    in_expression, in_all_genes, in_peaks_list = (in_files[0], in_files[1],
                                                  in_files[2:])
    # get the gene list
    # default_val = cfg.getfloat('genes','ks_test_default_value')
    # gene_expr = dict((line.strip().split('\t')[3],default_val)
    #                               for line in open(in_all_genes))
    gene_expr = {}
    #print in_expression
    for line in open(in_expression):
        try:
            gene_id, expr_val = line.strip().split('\t')
        except:
            continue
        gene_expr[gene_id] = expr_val
    # sort the genes by expression value (low to high)
    gene_expr_sorted = sorted(gene_expr.items(), key=lambda x:float(x[1]))
    #print gene_expr_sorted[:20]
    #print gene_expr_sorted
    for in_peaks in in_peaks_list:
        #print in_peaks
        genes_with_peaks = {}
        for line in open(in_peaks):
            fields = line.strip().split('\t')
            peak_loc = '%s:%s-%s\t%s' % tuple(fields[:3] +
                                    [fields[4] if fields[4] != '.' else 'NA'])
            gene_id = fields[9]
            
            genes_with_peaks[gene_id] = peak_loc
        #print genes_with_peaks
        #print in_peaks, genes_with_peaks
        with open(in_expression + '.with_peaks.%s.ks_data' %
                                            in_peaks, 'w') as outfile:
            outfile.write('\t'.join(['gene_id', 'expression_val','has_peak',
                                     'peak_loc', 'peak_score']) + '\n')
            for gene_id, expr_val in gene_expr_sorted:
                has_peak = 1 if gene_id in genes_with_peaks else 0
                outfile.write('\t'.join(map(str, [gene_id, expr_val, has_peak,
                                genes_with_peaks[gene_id] if gene_id in
                                        genes_with_peaks else 'None\tNA'])) + '\n')
        # make data file with data reversed
        with open(in_expression + '.with_peaks.%s.reversed.ks_data' %
                                            in_peaks, 'w') as outfile:
            outfile.write('\t'.join(['gene_id', 'expression_val','has_peak',
                                     'peak_loc', 'peak_score']) + '\n')
            for gene_id, expr_val in reversed(gene_expr_sorted):
                has_peak = 1 if gene_id in genes_with_peaks else 0
                outfile.write('\t'.join(map(str, [gene_id, expr_val, has_peak,
                                genes_with_peaks[gene_id] if gene_id in
                                        genes_with_peaks else 'None\tNA'])) + '\n')

@transform(make_expression_ks, suffix('.ks_data'), '.expr_corr.png')
def draw_expression_correlation(in_data, out_png):
    """Correlation test to see if the correlation between expression values and
    peak quality (score column from peak file).
    """

    R_script = r"""
png('%(out_png)s')
d<-read.table(file="%(in_data)s", header=TRUE, sep="\t");
library(lattice);
r <- cor.test(d$expression_val, d$peak_score)
plot(d$expression_val, d$peak_score, xlab="expression value", ylab="peak score")
title(paste("R^2 = ", r$estimate, ", p-value = ", r$p.value));
dev.off()
""" % dict(in_data=in_data, out_png=out_png)
    #print R_script
    r = pyper.R()
    r(R_script)

@transform(make_expression_ks, suffix('.ks_data'), '.ks_plot.png')
def draw_expression_ks(in_data, out_png):
    """KS-test to see if the 1's and 0's from input is non-uniform"""

    R_script = r"""
d<-read.table(file="%(in_data)s", header=TRUE, sep="\t")
# This does the actual KS test
N<-dim(d)[1]
a<-d[1:N,3]
genes <- d[1:N, 1]
expr_vals <- d[1:N, 2]
q<-which(a==1)
t<-ks.test(q,"punif",0,N)
hits <- length(q)
# This is all for making the picture
png('%(out_png)s')
x<-array(1:N)
K<-sum(a)
G<-rep(-1*(K/(N-K)),times=N)
indexA <- a>0
G[indexA] <- 1  # G is 1 everywhere there is a hit, -K/N-K where a miss, sums to 0
rs<-cumsum(G)
maxRS = max(rs)
minRS = min(rs)
topY = 1.6*max(maxRS, 30)  # To give a bit of room for the hash marks

par(mar=c(4, 4, 6, 4) + 0.1)  # bit more room on the top and sides than usual
plot(NA,NA,xlim=c(0,N),ylim=c(minRS,topY),main=paste("%(plot_label)s", '\nHits: ', hits, ' pvalue: ', signif(t$p.value,2), '\n\n',sep=''),xlab='Gene Rank',ylab='Running Enrichment Score')
lines(x,rs)
if (maxRS > abs(minRS)){
	redX = which(rs == max(rs))
	redY = maxRS
} else {
	redX = which(rs == min(rs))
	redY = minRS
}
lines(c(0,N),c(0,0))
lines(c(redX,redX),c(0,redY),col='red')
# Hash marks
#rug(q,ticksize=.035,side=3)

library(beanplot)
#par(fig=c(0,1,.6,1), new=T)
beanwidth <- min(max(8, maxRS), 150)
beanplot(q, axes=T, horizontal=T, ylim=c(0,N), add=T, at=maxRS*1.3 + 20, boxwex=beanwidth)

#library(vioplot)
#par(fig=c(-.1,1.1,.6,1), new=T)
#vioplot(q, horizontal=TRUE)

# x-axis expression values
options(digits=3)  #stupid R.  This is supposed to make there be fewer digits used, but doesn't work in axis
axis(3, at=c(0, 1/4, 2/4, 3/4, 1) * N, labels=c(expr_vals[1], expr_vals[1*N/4],
                                                expr_vals[2*N/4], expr_vals[3*N/4],
                                                expr_vals[N]))

dev.off()
""" % dict(plot_label='KS-test for peaks within %sbp of genes\n%s'%(
                                    cfg.get('genes', 'nearby_genes_max_dist'),
                                    in_data), in_data=in_data, out_png=out_png)
    #print R_script
    r = pyper.R()
    r(R_script)



@active_if(False)
@split(call_peaks.all_peak_caller_functions, regex(r'(.*)'),
           [r'\1.go_genelist', r'\1.go_enrichedGO', r'\1.go_Rout'])
def gene_ontology(in_peaks, out_files):
    """Calculate the significance of the peaks near genes using BioConductor"""
    out_genes, out_go, out_raw = out_files
    cmd = """echo '
    peaks = read.table("%s", header=FALSE, sep="\t");
    peaks = data.frame(chr=as.factor(peaks[,1]), start=as.numeric(peaks[,2]),
                        end=as.numeric(peaks[,3]));
    peaks = RangedData(IRanges(start=peaks[,2], end=peaks[,3]), space=peaks[,1])
    source("http://bioconductor.org/biocLite.R");
    biocLite("ChIPpeakAnno");
    library(ChIPpeakAnno);
    mart<-useMart(biomart="ensembl",dataset="%s");
    tss = getAnnotation(mart, featureType="TSS");
    annopeaks = annotatePeakInBatch(peaks[, ], AnnotationData=tss);
    write.table(annopeaks, file="%s", sep="\t");
    ' | R --vanilla --slave > %s""" % (in_peaks, cfg.get('DEFAULT', 'R_mart'),
                                       out_genes, out_go, out_raw)
    print cmd
    touch(out_raw)
