#!/usr/bin/python

sizes = dict((f.split('\t')[0], int(f.split('\t')[1])) for f in open('hg18.chrom.sizes'))
binsize = 200
window_size = 50

from cookb_signalsmooth import smooth



from scipy import signal

def gauss_kern(size):
    """ Returns a normalized gauss kernel array for convolutions """
    size = int(size)
    x = sp.mgrid[-size:size+1]
    g = sp.exp(-(x**2/float(size) ))
    return g / g.sum()

def blur_image(im, n) :
    """ blurs the image by convolving with a gaussian kernel of typical
        size n.
    """
    g = gauss_kern(n)
    improc = signal.convolve(im, g, mode='same')
    return(improc)




import scipy as sp

import sys, os
from optparse import OptionParser
from MAT import AffyFileParser
def parse(fname, level, step, window):
    
    if not os.path.isfile(fname):
        print 'Error: could not open %s' % fname
        sys.exit('Stopping')

    oname = fname + '.wig'
    
    fp = open(oname, 'wt')
    fp.write( 'track type=wiggle_0\n' )
    
    fedges = open(fname + '.edges.bed', 'w')

    # mat file parsing parsing
    from MAT import AffyFileParser
    bar  =  AffyFileParser.CBARFileWriter()
    print fname
    bar.SetFileName(fname)
    if not bar.Read():
        raise Exception('Unable to properly read %s' % fname )

    seq  = AffyFileParser.CGDACSequenceResultItem()
    data = AffyFileParser.BarSequenceResultData()

    seq_range = range(bar.GetNumberSequences() - 1)
    for s in seq_range:
        bar.GetResults(s, seq)
        name  = seq.GetName()
        if name in ['XIST', 'chloroplast', 'NC_000964.1']:
            continue
        
        read_counts = sp.zeros(sizes[name] // binsize + 1)

        for i in range( 0, seq.GetNumberDataPoints()-1, step ):
            seq.GetData(i, 0, data)
            x = data.iValue
            seq.GetData(i, 1, data)
            y = data.fValue
            out ='%s\t%s\t%s\t%s\n' % ( name, x-1, x, y )
            fp.write( out )
            
            read_counts[x // binsize] += y
        
        d = sp.gradient(read_counts)
        d2 = sp.gradient(d)
        d = sp.array([0] + list(d[:-1]))
        d2 = sp.array([0,0] + list(d2[:-2]))
        
        #smoothed = smooth(read_counts, window_len=3)
        #read_counts[(read_counts > -.5) & (read_counts < .5)] = 0
        smoothed = blur_image(read_counts, window_size)
        #smooth_d = sp.gradient(smoothed)
        #smooth_d2 = sp.gradient(smooth_d)
        
        pos_to_neg = sp.diff(sp.sign(smoothed)) < 0
        #decreasing = (sp.diff(smoothed) < 0)
        
        #for x_d in ((smoothed > 1) & (smooth_d < 1)).nonzero()[0]:
        #for x_d in sp.nonzero(pos_to_neg & decreasing)[0]:
        for x_d in sp.nonzero(pos_to_neg)[0]:
            fedges.write('\t'.join(map(str, [name, (x_d+1) * binsize, (x_d + 2) * binsize])) + '\n')

    print '...saved bar data to: %s' % oname
                    

def usage():
    "Prints usage information"
    
    help = """
    
    USAGE: 
    
    bar2wig -t threshold -s step input_file

    PARAMETERS:
    
        - threshold is the minimal signal strength (default value 0.0)
        - step specifies the increment of the indices when sampling the  
          data and must be an integer >= 1 (default 1)
        
    Unspecified parameters will take the default values.

    EXAMPLE:
        
        bar2wig -t 1.5 -s 5 somedata.bar

    Will create the somedata.bar.wig file with threshold of 1.5
    while examining every 5th measurement. 

    SEE ALSO: wig2peak, wigalign
    """
   
    return help

if __name__ == '__main__':
    #sys.argv.extend( [ 'bar2text.py', '2' ] )
    try:
        import psyco
        psyco.full()
    except:
        pass

    parser = OptionParser()
    parser.set_usage( usage() )

    parser.add_option("-t", type="float", dest="level", default=0 )
    parser.add_option("-s", type="int", dest="step", default=1)
    parser.add_option("-w", type="int", dest="window", default=5)
    (options, args) = parser.parse_args()

    if len(args)!=1:
        parser.print_usage()
        sys.exit(-1)
    else:
        fname = args[-1]

    if options.step < 1:
        parser.print_usage()
        print 'Incorrect step option %s. It must be greater or equal to 1' % options.step
        sys.exit()

    print '\n...executing bar2wig with threshold=%s, step=%d' % (options.level, options.step) 

    parse(fname=fname, level=options.level, step=options.step, window=options.window)

    print '...done\n'
