

from hts_waterworks.utils.makeGeneStructure import parse_gene_line

def group_reads_by_gene(genefile, readsfile):
    """Group BED reads from readsfile by the gene they overlap with.
    
    readsfile must be sorted.
    """
    all_genes = map(parse_gene_line, open(genefile))
    all_genes.sort(key=lambda f: (f[1], int(f[3])))  # sort genes by chrom,start
    
    # get first gene
    all_genes = iter(all_genes)
    (g_name, g_chrom, g_strand, g_txStart, g_txEnd, g_cdsStart,
                g_cdsEnd, g_exons, g_name2, g_noncoding) = all_genes.next()
    
    # group reads by what gene they overlap with
    overlapping = []
    for line in open(readsfile):
        next_gene = False
        (r_chrom, r_start, r_end, r_name, r_score,
            r_strand) = line.strip().split('\t')
        r_start, r_end = int(r_start), int(r_end)
        # check for overlap...
        if r_chrom > g_chrom:
            next_gene = True
        elif r_chrom == g_chrom:
            if r_start >= g_end:
                next_gene = True
            else:
                overlaps = g_start <= r_start < g_end
                overlaps = overlaps or g_start < r_end <= g_end
                overlaps = overlaps or (r_start <= g_start and r_end >= g_end)
                if overlaps:
                    overlapping.append((r_chrom, r_start, r_stop, r_name,
                                        r_score, r_strand))
        if next_gene:
            yield g_name, overlapping
            overlapping = []
            (g_name, g_chrom, g_strand, g_txStart, g_txEnd, g_cdsStart,
                g_cdsEnd, g_exons, g_name2, g_noncoding) = all_genes.next()
    # reached the end of the reads, but there may still be genes...
    yield g_name, overlapping
    for (g_name, g_chrom, g_strand, g_txStart, g_txEnd, g_cdsStart,
                g_cdsEnd, g_exons, g_name2, g_noncoding) in all_genes:
        yield g_name, []


def group_adjacent_reads(read_groups, max_dist, use_count=False):
    """Agglomerate adjacent reads from different groups into a score matrix.
    
    Returns a tuple whose first value is the (chrom,start,stop) of the read
    inducing the grouping (left-most read) and the second value is a list
    containing the sum of scores (or count of reads if use_count is True)
    within max_dist of the inducing for each group.  Reads are induced from
    left-to-right and must be sorted by (chrom,start) within each group.
    
    For example, if there is a set of reads:
    group_1:
        chrX    500    501   read_1   8     +
        chrX    535    536   read_2   2     -
    group_2:
        chrX    480    481   read_3   3     +
    group_3:
        chrX    560    561  read_4    10    -
    then taking the list of the result would give:
    [
        (('chrX', 500, 501), [10, 3, 0]),
        (('chrX', 480, 481): [0, 0, 3])
    ]
    
    """
    inds = [0] * len(read_groups)
    # until inds[i] is at the end of the reads in all groups
    while any(inds[i] < len(read_groups[i] for i in range(len(read_groups)))):
        # find the group with the left-most read
        left_most = min([read_groups[inds[i]] for i in range(len(read_groups))
                         if inds[i] < len(read_groups[i])])
        # collect all the reads within max_dist of left_most for all groups
        totals = [0] * len(read_groups)
        for i in range(len(read_groups)):
            while (inds[i] < len(read_groups[i]) and  # group has more reads
                   read_groups[i][inds[i]][0] == left_most[0] and  # same chrom
                    read_groups[i][inds[i]][1] <= left_most[1] +
                                                    max_dist): # close enough
                if use_count:
                    totals[i] += 1
                else:
                    totals[i] += float(read_groups[i][inds[i]][4])  # add score
                inds[i] += 1
        yield (tuple(left_most), totals)
