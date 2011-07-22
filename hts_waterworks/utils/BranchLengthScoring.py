"""Branch Length Scoring.

Reference:
doi:10.1093/bioinformatics/btn605

Implemented by Kenny Daily 2009/2010
Institute of Genomics and Bioinformatics
University of California, Irvine

All software and related material contained 
herein can be downloaded freely for academic,
non-commercial, research use only.

For any other use, please contact Pierre Baldi at pfbaldi _at_ ics uci edu. 

Copyright 2009-2010

Contact: kdaily _at_ ics uci edu

"""

from collections import defaultdict

## REQUIRED PACKAGE, version 1.2
## http://www.daimi.au.dk/~mailund/newick.html
## http://www.daimi.au.dk/~mailund/newick/newick-1.2.tar.gz

import newick
import newick.tree

def parse_tree(treestring):
    """Wrapper on the newick.parse_tree to add parent edges and a leaf lookup.
    """

    my_tree = newick.parse_tree(treestring)
    newick.tree.add_parent_links(my_tree)

    return my_tree

def getMaxBBLS(tree):
    """Find max score by assuming P(motif) = 1 for all leaves.

    This is the sum of all the branches (max BLS score too)

    """
    
    leafNames = [leaf.id for leaf in tree.leaves]
    score_dict = dict([(name, 1.0) for name in leafNames])
    return BBLS(tree, score_dict)

class BaseVisitor(newick.tree.TreeVisitor):

    def get_sister(self, node):
        parent = node.parent
        left, right = self.get_children(parent)
        sister = None
        if left != node:
            sister = left
        else:
            sister = right          
        return sister

    def get_children(self, node):
        left, right = node.edges
        bl, sl, ll = left
        br, sr, lr = right
        return bl, br

class ProbVisitor(BaseVisitor):
    
    def __init__(self, scoredict):
        BaseVisitor.__init__(self)
        self.scoredict = scoredict
    
    def visit_leaf(self, leaf):
        leaf.prob = self.scoredict.get(leaf.id, 0.0)

    def post_visit_tree(self, tree):
        left, right = self.get_children(tree)
        tree.prob = 1 - ((1 - left.prob) * (1 - right.prob))

class EffectiveLengthVisitor(BaseVisitor):
    
    def __init__(self, tree):
        BaseVisitor.__init__(self)
        self.tree = tree
        self.tree.edge_length = 0.0
        
    def visit_leaf(self, l):
        l.effective_length = l.edge_length

    def pre_visit_edge(self, src, b, l, dst):
        dst.edge_length = l
                
    def post_visit_tree(self, t):
        if t != self.tree:
            my_length = t.edge_length
            left, right = self.get_children(t)    
            try:
                weighted_len = ((left.prob * left.effective_length) + (right.prob * right.effective_length)) / t.prob
                t.effective_length = my_length + weighted_len
            except ZeroDivisionError:
                t.effective_length = 0

        else:
            t.effective_length = 0

class TreeComplementVistor(BaseVisitor):
    
    def __init__(self, tree):
        BaseVisitor.__init__(self)
        self.tree = tree
        self.bbls = 0.0

    def post_visit_tree(self, t):
        left, right = self.get_children(t)
        self.bbls += (right.prob * left.prob * t.tree_complement_prob * (right.effective_length + left.effective_length))

    def get_bbls(self):
        return self.bbls
    
    def pre_visit_tree(self, t):
        if t == self.tree:
            t.tree_complement_prob = 1
        else:
            parent = t.parent
            sister = self.get_sister(t)
            t.tree_complement_prob = parent.tree_complement_prob * (1 - sister.prob)
        
def BBLS(tree, motifs):       
    """Convenience function to compute the BBLS.
    
    """

    ProbVis = ProbVisitor(motifs)
    ELVis = EffectiveLengthVisitor(tree)
    TCVis = TreeComplementVistor(tree)

    tree.dfs_traverse(ProbVis)
    tree.dfs_traverse(ELVis)
    tree.dfs_traverse(TCVis)

    return TCVis.get_bbls()

def BLS(tree, motifs):
    """Convenience function to compute the BLS.
    
    """

    new_motifs = {}
    for key, value in motifs.iteritems():
        if value > 0:
            new_motifs[key] = 1

    return BBLS(tree, new_motifs)

def _test():
    """A few test cases.

    """
    
    treestring = '(((aa:1,bb:3):10,cc:5):100,dd:7)'
    mmnewick = parse_tree(treestring)

    leaves = dict(aa=0.5, bb=0.5, cc=0.5, dd=0.5)
    bls = BLS(mmnewick, leaves)
    assert bls == 126.0

    bbls = BBLS(mmnewick, leaves)
    assert bbls == 56.375
    
    leaves = dict(aa=0.5)
    bls = BLS(mmnewick, leaves)
    assert bls == 0.0

    bbls = BBLS(mmnewick, leaves)
    assert bbls == 0.0

def main():
    pass

if __name__ == "__main__":
    _test()
    main()
