'''  makeGeneStructure.py
        * Created by Jacob Biesinger
        * This version: Created October 2010
'''

import optparse
from collections import defaultdict
import sys
import os
from bisect import insort, bisect
import itertools


from hts_waterworks.utils.common import bedCommentFilter
from hts_waterworks.utils.bedOverlapRandomShuffle import checkOverlap, generateShuffledBed

def main(argv=None):
    """From an annotated refGene file, create a "gene structure."
    
    For each gene, an entry is placed in separate files representing these regions:
         promoter, 5'UTR, 3'UTR, downstream, Exon, Intron

    The refseqGenes tab-separated file format has the headers:
         bin   name    chrom   strand  txStart txEnd   cdsStart    cdsEnd  exonCount   exonStarts  exonEnds [...]
    and can be downloaded from:
    'http://hgdownload.cse.ucsc.edu/goldenPath/YOUR_GENOME/database/refGene.txt.gz'
    """
    usage = "%prog [options] refGene.txt\n" + main.__doc__
    parser = optparse.OptionParser(usage)
    parser.add_option('--promoter_size', type='int', default=2000,
                      help='Size in bp of the upstream promoter region. default=%default')
    parser.add_option('--promoter_extend',type='int', default=0,
                      help='Extend the promoter regions past the Tx Start Site (into the 5utr or exon) by X bp. default=%default')
    parser.add_option('--downstream_size', type='int', default=500,
                      help='Size in bp of the downstream region. default=%default')
    parser.add_option('--downstream_extend', type='int', default=0,
                      help='Extend the downstream regions past the Tx Stop Site (into the 3utr or exon) by X bp. default=%default')
    parser.add_option('--with_gene_name', action='store_true',
                      help='also report the gene name in the file')

    if not argv:
        argv = sys.argv[1:]
    global opts, args
    opts, args = parser.parse_args(argv)
    if len(args) != 1:
        parser.error('Please specify a gene list as input!')
    
    genefile = open(args[0], 'r')
    splitGeneStructure(genefile, args[0], opts.promoter_size, opts.promoter_extend,
                                    opts.downstream_size, opts.downstream_extend)
      
def splitGeneStructure(geneLines, outNameBase, promoterSize, promoterExtend, downstreamSize, downstreamExtend):
    '''parse genes from a refseqGene.txt file. Fields are tab-separated in the order:
    bin   name    chrom   strand  txStart txEnd   cdsStart    cdsEnd  exonCount   exonStarts  exonEnds [...] '''
    global opts, args
    geneLines = itertools.ifilter(bedCommentFilter, geneLines)  # filter comments and headers
    promoter_out = open(outNameBase + '.promoter%s_ext%s' % (promoterSize, promoterExtend), 'w')
    downstream_out = open(outNameBase + '.down%s_ext%s' % (downstreamSize, downstreamExtend), 'w')
    utr5_out = open(outNameBase + '.utr5', 'w')
    utr3_out = open(outNameBase + '.utr3', 'w')
    exon_out = open(outNameBase + '.exon', 'w')
    intron_out = open(outNameBase + '.intron', 'w')
    tss_out = open(outNameBase + '.tss', 'w')
    noncoding_out = open(outNameBase + '.noncoding', 'w')
    for line in geneLines:
        try:
            _ = int(line.strip().split('\t')[0])  # test if the first col is bin
            (name, chrom, strand, txStart, txEnd, cdsStart,
                cdsEnd, exons, name2, noncoding) = parse_gene_line(line, startCol=1)
        except:
            (name, chrom, strand, txStart, txEnd, cdsStart,
                cdsEnd, exons, name2, noncoding) = parse_gene_line(line, startCol=0)
        if noncoding:
            noncoding_out.write('\t'.join([chrom, str(txStart), str(txEnd)] + ([name, name2,strand] if opts.with_gene_name else [])) + '\n')
            continue
        if strand == '+':
            promoterLeft = max(0, txStart - promoterSize)
            promoterRight = txStart + promoterExtend 
            if promoterLeft == promoterRight:
                promoterRight += 1
            downstreamLeft = max(0, txEnd - downstreamExtend)
            downstreamRight = txEnd + downstreamSize
            if downstreamLeft == downstreamRight:
                downstreamRight += 1
            utr5Left = txStart
            utr5Right = cdsStart
            if utr5Left == utr5Right:
                utr5Right += 1
            utr3Left = cdsEnd
            utr3Right = txEnd
        else:
            promoterLeft = max(0, txEnd - promoterExtend)
            promoterRight = txEnd + promoterSize
            downstreamLeft = max(0, txStart - downstreamSize)
            downstreamRight = txStart + downstreamExtend
            utr5Left = cdsEnd
            utr5Right = txEnd
            utr3Left = txStart
            utr3Right = cdsStart
            
        if promoterLeft == promoterRight:
            promoterRight += 1
        if downstreamLeft == downstreamRight:
            downstreamRight += 1
        if utr5Left == utr5Right:
            utr5Right += 1
        if utr3Left == utr3Right:
            utr3Right += 1
        promoter_out.write('\t'.join([chrom, str(promoterLeft), str(promoterRight)] + ([name, name2,strand] if opts.with_gene_name else [])) + '\n')
        downstream_out.write('\t'.join([chrom, str(downstreamLeft), str(downstreamRight)] + ([name, name2,strand] if opts.with_gene_name else [])) + '\n')
        utr5_out.write('\t'.join([chrom, str(utr5Left), str(utr5Right)] + ([name, name2,strand] if opts.with_gene_name else [])) + '\n')
        utr3_out.write('\t'.join([chrom, str(utr3Left), str(utr3Right)] + ([name, name2,strand] if opts.with_gene_name else [])) + '\n')
        for exonLeft, exonRight in exons:
            if exonLeft == exonRight:
                exonRight += 1
            exon_out.write('\t'.join([chrom, str(exonLeft), str(exonRight)] + ([name, name2,strand] if opts.with_gene_name else [])) + '\n')
        for intronLeft, intronRight in [(exons[i][1], exons[i+1][0]) for i in range(len(exons)-1)]:
            if intronLeft == intronRight:
                intronRight += 1
            intron_out.write('\t'.join([chrom, str(intronLeft), str(intronRight)] + ([name, name2,strand] if opts.with_gene_name else []) ) + '\n')
        if strand == '+':
            tss_out.write('\t'.join([chrom, str(txStart), str(txStart + 1)] + ([name, name2,strand] if opts.with_gene_name else [])) + '\n')
        else:
            tss_out.write('\t'.join([chrom, str(txEnd), str(txEnd + 1)] + ([name, name2,strand] if opts.with_gene_name else [])) + '\n')

def parse_gene_line(geneline, startCol=1):
    """From a refseq-formatted gene file, parse a single gene's attributes"""
    fields = geneline.strip().split('\t')
    name, chrom, strand = fields[startCol: startCol + 3]
    txStart, txEnd, cdsStart, cdsEnd = map(int, fields[startCol+3:startCol+7])
    exons = zip(map(int, fields[startCol+8][:-1].split(',')), map(int, fields[startCol+9][:-1].split(',')))
    name2 = fields[startCol + 11]
    noncoding = name.startswith('NR_') or cdsStart < 0 or cdsEnd < 0
    return (name, chrom, strand, txStart, txEnd, cdsStart,
                cdsEnd, exons, name2, noncoding)

if __name__ == '__main__':
    main()
