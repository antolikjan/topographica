import numpy
import topo
import pickle
from param import normalize_path

def circular_dist(a,b,period):
    """
    The distance between a and b (scalars) in the periodic. 
    a,b have to be in (0,period)
    """
    return  numpy.minimum(abs(a-b), period - abs(a-b))

def rad_to_complex(vector):
    """
    Converts a vector/matrix of angles (0,2*pi) to vector/matrix of complex numbers (that will lie on the unit circle)
    """
    return numpy.cos(vector)+1j*numpy.sin(vector)
    

def angle_to_pi(array):
    """
    returns angles of complex numbers in array but in (0,2*pi) interval unlike numpy.angle the returns it in (-pi,pi)
    """
    return (numpy.angle(array)+ 4*numpy.pi) % (numpy.pi*2)
    

def circ_mean(matrix, weights=None, axis=None, low=0, high=numpy.pi*2,
              normalize=False):
    """
    Circular mean of matrix. Weighted if weights are not none.

    matrix     - matrix of data. Mean will be computed along axis axis.
    weights    - if not none, matrix of the same size as matrix
    low, high  - the min and max values that will be mapped onto the periodic
                 interval of (0, 2pi)
    axis       - axis along which to compute the circular mean. default = 0 (columns)
    normalize  - if True weights will be normalized along axis. If any weights
                 that are to be jointly normalized are all zero they will be
                 kept zero!

    return (angle, length) - where angle is the circular mean, and len is the
                           length of the resulting mean vector
    """

    # check whether matrix and weights are ndarrays
    if weights != None:
       assert matrix.shape == weights.shape 
    
    if axis == None:
        axis == 0

    # convert the periodic matrix to corresponding complex numbers
    m = rad_to_complex((matrix - low)/(high - low) * numpy.pi*2)

    # normalize weights
    if normalize:
        row_sums = numpy.sum(numpy.abs(weights), axis=axis)
        row_sums[numpy.where(row_sums == 0)] = 1.0
        if axis == 1:
            weights = weights / row_sums[:, numpy.newaxis]
        else:
            weights = numpy.transpose(weights) / row_sums[:, numpy.newaxis]
            weights = weights.T
            
    if weights == None:
        m = numpy.mean(m, axis=axis)
    else:
        z = numpy.multiply(m,weights)
        m = numpy.mean(z, axis=axis)

    return ((angle_to_pi(m) / (numpy.pi*2))*(high-low) + low, abs(m))


def compute_local_homogeneity_index(or_map,sigma):
    (xsize,ysize) = or_map.shape 
    
    lhi = numpy.zeros(or_map.shape) 
    
    for sx in xrange(0,xsize):
        for sy in xrange(0,ysize):
            lhi_current=[0,0]
            for tx in xrange(0,xsize):
                for ty in xrange(0,ysize):
                    lhi_current[0]+=numpy.exp(-((sx-tx)*(sx-tx)+(sy-ty)*(sy-ty))/(2*sigma*sigma))*numpy.cos(2*or_map[tx,ty])
                    lhi_current[1]+=numpy.exp(-((sx-tx)*(sx-tx)+(sy-ty)*(sy-ty))/(2*sigma*sigma))*numpy.sin(2*or_map[tx,ty])
            lhi[sx,sy]= numpy.sqrt(lhi_current[0]*lhi_current[0] + lhi_current[1]*lhi_current[1])/(2*numpy.pi*sigma*sigma)
            
    return lhi


def analyse_connectivity(sheet_name,proj_name,lhi):
    or_bins = numpy.linspace(0,numpy.pi/2,30)
    projection = topo.sim[sheet_name].projections()[proj_name]
    or_pref = topo.sim[projection.src.name].sheet_views['OrientationPreference'].view()[0]*numpy.pi
    or_pref_target = topo.sim[projection.dest.name].sheet_views['OrientationPreference'].view()[0]*numpy.pi
    
    domain_orientation_connection_strength = numpy.zeros((30,1))
    pinwheel_orientation_connection_strength = numpy.zeros((30,1))
    
    pinwheels = numpy.where(lhi.flatten()<0.1)[0]
    domains = numpy.where(lhi.flatten()>0.9)[0]
    
    for i in pinwheels:
        cf = projection.cfs.flatten()[i]
        this_or = or_pref_target.flatten()[i]
        ors = cf.input_sheet_slice.submatrix(or_pref).flatten()
        weights = numpy.multiply(cf.weights,cf.mask).flatten()
        
        for j,k in enumerate(numpy.digitize(circular_dist(ors,this_or,numpy.pi),or_bins)):
            print k    
            pinwheel_orientation_connection_strength[k-1] += weights[j]

    for i in domains:
        cf = projection.cfs.flatten()[i]
        this_or = or_pref_target.flatten()[i]
        ors = cf.input_sheet_slice.submatrix(or_pref).flatten()
        weights = numpy.multiply(cf.weights,cf.mask).flatten()
        
        for j,k in enumerate(numpy.digitize(circular_dist(ors,this_or,numpy.pi),or_bins)):
            domain_orientation_connection_strength[k-1] += weights[j]

        
    import pylab
    pylab.figure()
    pylab.subplot(2,1,1)
    pylab.hold('on')
    pylab.plot(numpy.linspace(0,numpy.pi/2,30,endpoint=False)+numpy.pi/2/60,pinwheel_orientation_connection_strength,'k',linewidth=2.0)
    pylab.plot(numpy.linspace(0,numpy.pi/2,30,endpoint=False)+numpy.pi/2/60,pinwheel_orientation_connection_strength,'ko')
    pylab.title(proj_name,fontsize=20)
    pylab.autoscale(tight=True)
    pylab.title('pinwheels')
    
    pylab.subplot(2,1,2)
    pylab.hold('on')
    pylab.plot(numpy.linspace(0,numpy.pi/2,30,endpoint=False)+numpy.pi/2/60,domain_orientation_connection_strength,'k',linewidth=2.0)
    pylab.plot(numpy.linspace(0,numpy.pi/2,30,endpoint=False)+numpy.pi/2/60,domain_orientation_connection_strength,'ko')
    pylab.title(proj_name,fontsize=20)
    pylab.autoscale(tight=True)
    pylab.title('domains')
    
    from param import normalize_path
    pylab.savefig(normalize_path('PPconnectivity: ' + proj_name + str(topo.sim.time()) + '.png'));

Analyse_connectivity():
    if True:
        e_lhi = compute_local_homogeneity_index(topo.sim['V1Complex'].sheet_views['OrientationPreference'].view()[0]*numpy.pi,2.0)    
        i_lhi = compute_local_homogeneity_index(topo.sim['V1ComplexInh'].sheet_views['OrientationPreference'].view()[0]*numpy.pi,2.0)    
        f = open('lhi2.0_E.pickle','wb')            
        pickle.dump(e_lhi,f)
        f.close()
        f = open('lhi2.0_I.pickle','wb')            
        pickle.dump(i_lhi,f)
        f.close()
    else:        
        f = open('lhi2.0_E.pickle','rb')            
        e_lhi = pickle.load(f)
        f = open('lhi2.0_I.pickle','rb')            
        i_lhi = pickle.load(f)

    
    analyse_connectivity('V1Complex','LongEE',e_lhi)
    analyse_connectivity('V1ComplexInh','LongEI',i_lhi)
