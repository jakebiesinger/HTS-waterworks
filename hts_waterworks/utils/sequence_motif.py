#!/usr/bin/env python

from pygr import worldbase
import motility
import cPickle as pickle
import random
import math
import scipy
import re
import scipy
import scipy.stats
import scipy.special


from sampling import calc_pwm_bg_dist, weighted_sample
from hts_waterworks.utils.common import makeNormalSeq, iupac2matrix

class Motif(motility.PWM):
    """ Extension of general PWM, adding pygr-style search and background calculation"""

    _sd = scipy.nan
    _mean = scipy.nan
    _zscore = 4.265
    _threshold = scipy.nan
    _attrs = {}

    # Pickle helper functions
    def __getstate__(self):
        result = self.__dict__.copy()
        del result['_m']
        return result
    def __setstate__(self, savedDict):
        self.__dict__ = savedDict
        motility.PWM.__init__(self, savedDict['matrix'])

    def calculate_background(self, genome, samples=None, **kwargs):
        """ Calculate the background distribution for this motif """
        if samples:
            mean, sd = calc_pwm_bg_dist(self, genome, samples=samples)
        elif len(kwargs) > 0:
            mean, sd = calc_pwm_bg_dist(self, genome, keepScores = False, **kwargs)
        else:
            mean, sd = calc_pwm_bg_dist(self, genome, keepScores = False)
        self._mean, self._sd = mean, sd
        return self

    def find_in_region(self, region, zscore=None):
        """ return motif instances at the precalculated threshold """
        if not self.bg_calculated():
            raise RuntimeError('Need to calculate the bacgkround distribution first!')
        threshold = self.get_threshold(zscore)
        seq = str(region) # get sequence
        seq = makeNormalSeq(seq)
        motifHits = self.find(seq, threshold) # get hits in sequence
        return motifHits

    def bg_calculated(self):
        "Check if the background distribution for this motif has been calculated"
        return not scipy.isnan(self._threshold) or (not scipy.isnan(self._mean) and not scipy.isnan(self._sd))
    def get_threshold(self, zscore=None):
        if not scipy.isnan(self._threshold):
            return self._threshold
        if not zscore:
            zscore = self._zscore
        thresh_for_z = self._mean + self._sd * zscore
        return min(self.max_score(), thresh_for_z)
    
    def generate_sites(self, to_generate):
        "Generate random sites from the Motif"
        columns = []
        for row in self.matrix:
            print row
            if min(row) < 0:
                row = [2**(el + math.log(.25, 2)) for el in row]
                row = [el / sum(row) for el in row]
            print row
            l = list(weighted_sample(zip(row, 'ACGT'), to_generate))
            random.shuffle(l)
            columns.append(l)
        sites = [''.join(letter) for letter in zip(*columns)]
        return sites


def makePWMFromIUPAC(iupac_motif):
    pwm = [iupac2matrix[l] for l in iupac_motif]
    pwm = Motif(pwm)
    pwm._threshold = len(iupac_motif)  #perfect match required
    return pwm
    


def search_bed_file(bedLines, pwm, genome, reportLocations=False):
    """ return lines of a bed file with instances of given motif """
    print 'searching'
    for line in bedLines:
        if line[0] == '#' or line.split()[0] == 'track':
            # skip header and comment lines
            continue
        chrom, start, stop = line.split()[:3]
        start, stop = int(start), int(stop)
        if start >= stop or start >= len(genome[chrom]):
            continue
        #print start, stop
        region = genome[chrom][start:stop]
        instances = pwm.find_in_region(region)
        motifFound = len(instances) > 0
        if motifFound:
            if reportLocations:
                yield line, instances
            else:
                yield line


""" Utilities for calculating the significance of motif enrichment
"""

def fold_enrichment(fgMatches, fgSize, bgMatches, bgSize):
    """ return fold-enrichment of fg vs. background """
    return (float(fgMatches) / fgSize) / (float(bgMatches) / bgSize)

def zscore_normal(fgMatches, fgSize, bgMatches, bgSize):
    """ Calculate zscore for an enriched foreground, assuming normal distribution
    """
    if bgSize == 0 or fgSize == 0:
        raise RuntimeError('bgSize and fgSize must not be 0 for calculation!')
    bgRate = float(bgMatches) / float(bgSize)
    denom = scipy.sqrt(bgRate * (1 - bgRate) * fgSize)
    if denom == 0.0:
        zscore = scipy.inf
    else:
        zscore = (fgMatches - bgRate * fgSize) / denom
    return zscore

def zscore2pval_normal(zscore):
    return scipy.stats.norm.sf(zscore, loc=0, scale=1)

def pvalue_binomial(fgMatches, fgSize, bgMatches, bgSize):
    """ Calculate binomial pvalue for enriched foreground
    """
    if bgSize == 0 or fgSize == 0:
        raise RuntimeError('bgSize and fgSize must not be 0 for calculation!')
    bgRate = float(bgMatches) / float(bgSize)
    pvalue = scipy.stats.binom.sf(fgMatches, fgSize, bgRate)
    return pvalue

def zscore_hypergeometric(fgMatches, fgSize, bgMatches, bgSize):
    """ Get the Zscore for foreground enrichment, using hypergeometric distribution """
    if bgSize == 0 or fgSize == 0:
        raise RuntimeError('bgSize and fgSize must not be 0 for calculation!')
    fgMatches, fgSize, bgMatches, bgSize = map(scipy.float_, [fgMatches, fgSize, bgMatches, bgSize])
    mean = fgSize * (float(bgMatches) / bgSize)
    var = (bgMatches * fgSize * (float(bgSize) - bgMatches) * (bgSize - fgSize)) / (float(bgSize) * bgSize * (bgSize-1))
    zscore = (fgMatches - mean) / scipy.sqrt(var)
    return zscore

def pvalue_hypergeometric(fgMatches, fgSize, bgMatches, bgSize, log=True):
    """ Get the p-value for foreground enrichment, using hypergeometric distribution.
        If log is True, the calculation is performed in log space (won't underflow very very small values)"""
    # convert from natural log to log10
    #survival_values = [log_hypergeom(x, fgSize, bgSize, bgMatches) / scipy.log(10) for x in xrange(fgMatches, min(fgSize, bgMatches))]
    #pvalue = sum(scipy.power(10, survival_values))
    fgMatches, fgSize, bgMatches, bgSize = map(scipy.float_, [fgMatches, fgSize, bgMatches, bgSize])
    
    survival_values = [hypergeom(x, fgSize, bgMatches, bgSize, log=log) for x in scipy.arange(fgMatches, min(fgSize+1, bgMatches+1))]
    pvalue = sum(survival_values if not log else scipy.exp(survival_values))
    return pvalue

def hypergeom(x, fgTotal, bgMatching, bgTotal, log=True):
    if log:
        return logchoose(fgTotal,x) + logchoose(bgTotal-fgTotal,bgMatching-x) - logchoose(bgTotal,bgMatching)
    else:
        return scipy.comb(fgTotal,x) * scipy.comb(bgTotal-fgTotal,bgMatching-x) / scipy.comb(bgTotal,bgMatching)

def logchoose(n,k):
   lgn1 = scipy.special.gammaln(n+1)
   lgk1 = scipy.special.gammaln(k+1)
   lgnk1 = scipy.special.gammaln(n-k+1)
   return lgn1 - (lgnk1 + lgk1)



def search_annotDB(annotDB, pwm, reportLocations=False):
    """ search an annotation database (i.e., promoter regions) for instances of a motif """
    for id, annot in annotDB.iteritems():
        motifHits = pwm.find_in_region(annot.sequence)
        if reportLocations:
            yield id, motifHits
        else:
            yield id

def alignAndCompareMotifs(motif1, motif2, reportAll=False, tryAllAlignments=True, reverseComp=True, quitThreshold=None, normalizeRows=True, fillValue=.25):
    """ Compare the PWM's for two motifs by calculating their correlation coefficient.
    By default, all possible alignments and orientations will be tried and the top coefficient will be reported.
    fillValue may be a number, or a 4-element array with nuc frequencies
    Returns (corrCoef, motif2_relative_posn, motif2_orientation) form best alignment, or the entire list if reportAll=True.
    """
    pwm1,pwm2 = motif1.matrix, motif2.matrix
    if normalizeRows:  # make sum in each row = 1
        pwm1, pwm2 = map(normalizePwmRows, [pwm1, pwm2])
    alignsToTry = xrange(-len(motif2) + 1, len(motif1)-1) if tryAllAlignments else [0]  # all possible shifts or no shifting
    results = []
    for curOffset in alignsToTry:
        curPwm1, curPwm2 = map(scipy.array, extendPWMs(pwm1, pwm2, curOffset, fillValue))
        # flatten arrays and take 1-dimensional correlation between them
        corrCoef = scipy.corrcoef(curPwm1.ravel(), curPwm2.ravel())[0,1] # top-right is correlation between matrices
        results.append([corrCoef, curOffset, 1])
        if quitThreshold is not None and corrCoef > quitThreshold:
            # return immediately if quit threshold has been passed
            break
        if reverseComp:
            curPwm2 = scipy.array(reverseComplement(curPwm2))
            corrCoef = scipy.corrcoef(curPwm1.ravel(), curPwm2.ravel())[0,1] # top-right is correlation between matrices
            results.append([corrCoef, curOffset, -1])
        if quitThreshold is not None and corrCoef > quitThreshold:
            # return immediately if quit threshold has been passed
            break
    if reportBest:
        results = scipy.array(results)
        best = results[results[:,0].argmax(), :] # choose the result (row) with the best corrCoef
        return best
    else:
        return results

def extendPWMs(pwm1, pwm2, offset2=0, fillValue=.25):
    """Extend both PWMs so that they are the same length by filling in values from fillValue.
    Optionally, a positive or negative offset for motif2 can be specified and both motifs will be filled.
    fillValue may be a number, or a 4-element array with nuc frequencies
    Returns (extendedPwm1, extendedPwm2) as 2D lists
    """
    # check for errors, convert pwms to list if necessary
    if type(fillValue) != list:
        fillValue = [fillValue] * 4 # extend to 4 nucleotides
    elif len(fillValue) != 4:
        raise RuntimeError('fillValue for extendPWMs must be a single number or a 4-element list!')
    if type(pwm1) == scipy.ndarray:
        pwm1 = pwm1.tolist()
    if type(pwm2) == scipy.ndarray:
        pwm2 = pwm2.tolist()
    
    if offset2 < 0:
        # prepend filler for pwm1
        pwm1 = [fillValue] * scipy.absolute(offset2) + pwm1
    elif offset2 > 0:
        # prepend filler for pwm2
        pwm2 = [fillValue] * scipy.absolute(offset2) + pwm2
    # extend the pwms as necessary on the right side
    pwm1 = pwm1 + [fillValue] * (len(pwm2) - len(pwm1))
    pwm2 = pwm2 + [fillValue] * (len(pwm1) - len(pwm2))
    
    return pwm1, pwm2

def get_LOD_matrix(oldPwm):
    # convert to scipy array and/or change to floats
    prior_min = 0.01
    bg_freqs =  [0.236212624584718,
                 0.263413621262458,
                 0.257475083056478,
                 0.242898671096346]
    if type(oldPwm) != scipy.ndarray:
        newPwm = scipy.array(oldPwm, dtype=float)
    else:
        newPwm = oldPwm.astype(float)
    # normalize each row
    for row in newPwm:
        # shift rows s.t. minimum = 0 then normalize by total in row.
            row[:] = [scipy.log((row[i] / row.sum() + prior_min) / bg_freqs[i]) for i in xrange(4) ]
    return newPwm if type(oldPwm) == scipy.ndarray else newPwm.tolist() # return in same format as incoming


def normalizePwmRows(oldPwm):
    """Return normalized PWM such that all rows sum to 1.
    Format (list vs. array) is preserved.
    """
    # convert to scipy array and/or change to floats
    if type(oldPwm) != scipy.ndarray:
        newPwm = scipy.array(oldPwm, dtype=float)
    else:
        newPwm = oldPwm.astype(float)
    # normalize each row
    for row in newPwm:
        # shift rows s.t. minimum = 0 then normalize by total in row.
        row[:] = row - row.min()  # Editing newPwm in-place via row[:]
        if not row.any():
            # all 0's in this row => uniform originally
            row[:] = [.25] * 4
        else:
            row[:] = row / row.sum()
    return newPwm if type(oldPwm) == scipy.ndarray else newPwm.tolist() # return in same format as incoming

def reverseComplement(matrix):
    'Return a new PWM (2D list) that is the reverse-complement of matrix'
    return [ row[::-1] for row in matrix[::-1] ]
    
def parseMotifsFromTransfac(transfacStr):
    'Parse a transfac string and return a dictionary of motifs'
    blockStart = re.compile(r'^AC\s+([\w\d.\- ]+)')
    blockEnd = re.compile('^//')
    blankLine = re.compile(r'(^\s*$)|(^XX)')
    matrixStart = re.compile('^P0')
    matrixLine = re.compile('^\d+')
    transfacLines = transfacStr.split('\n')
    motifID = None
    inMatrix = False
    motifDict = {}
    for line in transfacLines:
        #print line
        if blockEnd.search(line):  # reached end of an entry- add previous results
            if motifID:
                newMotif = Motif(curMotif['pwm'])
                newMotif._mean = curMotif['mean'] if 'mean' in curMotif else scipy.nan
                newMotif._sd = curMotif['stdev'] if 'stdev' in curMotif else scipy.nan
                newMotif._threshold = curMotif['threshold'] if 'threshold' in curMotif else scipy.nan
                newMotif._attrs = curMotif
                motifDict[motifID] = newMotif
            motifID = None
            curMotif = None
            inMatrix = False
        else:
            startMatch = blockStart.search(line)  # AC  ####   start of an entry
            if startMatch:
                motifID = startMatch.groups()[0].replace(' ', '_')
                curMotif = {'AC':[motifID]}
            elif not blankLine.search(line):  # not a blank line or just XX
                if matrixStart.search(line):
                    curMotif['pwm'] = []
                    inMatrix = True
                elif inMatrix and matrixLine.search(line):
                    curMotif['pwm'].append([float(val) for val in line.split()[1:] if re.search('[0-9.]+', val)])
                else:  # some other line match-- like DE or BF, etc
                    inMatrix = False
                    key, value = line[:2], line[2:].strip()
                    if key not in curMotif:
                        curMotif[key] = []
                    curMotif[key].append(value)
                    if 'mean:' in line:
                        curMotif['mean'] = float(line.split(':')[1].strip())
                    elif 'stdev:' in line:
                        curMotif['stdev'] = float(line.split(':')[1].strip())
                    elif 'threshold:' in line:
                        curMotif['threshold'] = float(line.split(':')[1].strip())
    if motifID:
        newMotif = Motif(curMotif['pwm'])
        newMotif._mean = curMotif['mean'] if 'mean' in curMotif else scipy.nan
        newMotif._sd = curMotif['stdev'] if 'stdev' in curMotif else scipy.nan
        newMotif._threshold = curMotif['threshold'] if 'threshold' in curMotif else scipy.nan
        newMotif._attrs = curMotif
        motifDict[motifID] = newMotif
    return motifDict 


def motifsToTransfac(motifDict):
    outstr = ''
    for key, motif in motifDict.items():
        outstr += '//\n'
        outstr += 'XX\n'.join(['AC  %s\n' % key, 'NA  %s\n' % key, 'DE  mean:%s\nDE  stdev:%s\nDE  threshold:%s\n' % (motif._mean, motif._sd, motif._threshold)])
        outstr += 'DE  Motif discovered by Meme\n'
        outstr += 'XX\nP0      A      C      G      T\n'
        outstr += '\n'.join(('%02d    ' % (rowIndex+1)) + '   '.join(str(weight) for weight in row) for rowIndex, row in enumerate(motif.matrix))
        outstr += '\nXX\n'
    return outstr



def pwmToConsensus(pwm):
    "convert a pwm to a consensus sequence"
    return ''.join(map(pwmRowToConsensus, pwm))
    
def pwmRowToConsensus(row):
    "convert one row of a pwm to a consensus letter"
    pA, pC, pG, pT = [float(e)/sum(row) for e in row]
    if pA >= .7:
        return 'A'
    elif pC >= .7:
        return 'C'
    elif pG >= .7:
        return 'G'
    elif pT >= .7:
        return 'T'
    elif pA + pT > .8:
        return 'W'
    elif pA + pC > .8:
        return 'M'
    elif pA + pG > .8:
        return 'R'
    elif pT + pC > .8:
        return 'Y'
    elif pT + pG > .8:
        return 'K'
    elif pC + pG > .8:
        return 'S'
    elif pA + pT > .8:
        return 'W'
    elif pA < .1:
        return 'B'
    elif pC < .1:
        return 'D'
    elif pG < .1:
        return 'H'
    elif pT < .1:
        return 'V'
    else:
        return 'N'


def parseMemeMotifs(meme_file, logOdds=True):
    if logOdds:
        block_start = 'Motif (\d+) position-specific scoring matrix'
    else:
        block_start = 'Motif (\d+) position-specific probability matrix'
    block_end = '-------------'
    threshold_str = 'bayes= (\d+\.*\d*)'
    inBlock = False
    inHeader = False
    skip=0
    allMatrices = {}
    curMatrix = []
    with open(meme_file) as meme_data:
        for line in meme_data:
            if skip:  #skip over current line
                skip -= 1
                continue
            if inBlock:
                end_match = re.search(block_end, line)
                if end_match:
                    inBlock = False
                    curMotif = Motif(curMatrix)
                    curMotif._threshold = thresholdVal
                    consensus = pwmToConsensus(curMatrix)
                    motif_name = motifID + '_' + consensus
                    allMatrices[motif_name] =  curMotif
                    #print 'added', curMotif, motif_name
                    curMatrix = []
                    motifID = 'None'
                    continue
                letterFreqs = [ float(column) for column in line.split() ]
                curMatrix.append(letterFreqs)
                continue
            start_match = re.search(block_start, line)
            if start_match:
                motifID = 'Motif_' + start_match.groups()[0]
                if logOdds:
                    motifID += '_LOD'
                else:
                    motifID += '_Freq'
                #print 'now reading', motifID
                skip = 1
                inHeader = True
                continue
            if inHeader:
                inHeader = False
                #if logOdds:
                #    thresholdVal = float(re.search(threshold_str, line).groups()[0])
                #else:
                #    thresholdVal = scipy.nan
                thresholdVal = scipy.nan
                inBlock = True
        # end parsing
    if len(allMatrices) == 0:
        raise RuntimeError("Didn't parse any motifs from %s", meme_file)
    return allMatrices

def parseDremeMotifs(dreme_str, logOdds=True):
    raise NotImplementedError("parseDremeMotifs is not complete!")
    if logOdds:
        blockStart = 'log-odds matrix:'
    else:
        blockStart = 'letter-probability matrix:'
    motif_start = re.compile(
r'''^MOTIF (?P<consensus>\w+).*$
.*alength= (?P<alength>\d+) w= (?P<w>\d+) nsites= (?P<nsites>\d+) E= (?P<e_val>.*)$
(^\d[\d. ]+$)+''' % blockStart, re.MULTILINE)
    inMotif = False
    inBlock = False
    inHeader = False
    allMatrices = {}
    curMatrix = []
    for line in dreme_str.split('\n'):
        if not inMotif:
            m = motif_start.search(line)

def parse_xms_motifs(filename):
    next = False
    write = False
    letterFreqs = []
    motifId = ''

    def startElement(name, attrs):
        if re.search("name", name): #if at the start of a new motif get the ID
            next=True
        if re.search("weight", name): #start writing to the matrices if in the block
            if not re.search("weightmatrix", name):
                write=True            
    def charElement(data):
        if write:
            letterFreqs.append(data)
            #print 'writing Motif Weight ' +data
        if next:
            motifId=data
            #print 'Now Parsing ', motifId
            next=False
    def endElement(name):
        if re.search("column", name):
            curMatrix.append(map(float,letterFreqs))
            #print repr(letterFreqs)
            letterFreqs=[]
        if re.search("weightmatrix", name):
            allMatrices[motifId]=curMatrix
        if re.search("weight", name):
            write=False
    
    p=parser.ParserCreate()
    p.StartElementHandler=startElement
    p.CharacterDataHandler=charElement
    p.EndElementHandler=endElement
    p.Parse(data)
    allMatrices = {}
    for key in allMatrices:
        motif=Motif(allMatrices[key])
        allMatrices[key]=motif
    return allMatrices

