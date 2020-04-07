""" Wrapping C function Bragg using
    the Numpy declarations from Cython """

# cimport the Cython declarations for numpy
cimport numpy as np
import numpy as np
import datetime,time

np.import_array()
# cdefine the signature of our c function




# create the wrapper code, with numpy type annotations
def MonteCarloCCD(repeat, **kwargs):
    #Monte Carlo body CCD
    #Initialization
    size = kwargs['size']
    Bragg_peak_rel = kwargs['Bragg_peak_rel']
    Bloch_peak_rel = kwargs['Bloch_peak_rel']
    xc = kwargs['xc']
    vc = kwargs['vc']
    v_2nd_Bragg = kwargs['v_2nd_Bragg']
    g = kwargs['g']
    sigma_x = kwargs['sigma_x']
    sigma_v = kwargs['sigma_v']
    t = datetime.datetime.now()
    np.random.seed(int(time.mktime(t.timetuple())))
    
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] xy_t = np.zeros([size, 2], dtype=np.double)
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] z_t = np.zeros([size, 2], dtype=np.double)
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] v_t = np.zeros([size, 2], dtype=np.double)
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] vT = np.zeros([size, 2], dtype=np.double)
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] amp = np.zeros([size, 2], dtype=np.double)
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] phase = np.zeros([size, 2], dtype=np.double)
    
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] z0 = np.zeros([2], dtype=np.double)
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] v0 = np.ones([2], dtype=np.double)*(v_2nd_Bragg+g*(1.32-1.12))
    cdef int i = 0
    
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] random_var = np.zeros(size*2, dtype=np.double)
    
    random_var = np.random.normal(0, sigma_x, size*2)
    xy_t[:,0] = random_var[0:size]
    xy_t[:,1] = random_var[size:size*2]
    random_var = np.random.normal(0, sigma_v, size*2)
    vT[:,0] = random_var[:size]
    vT[:,1] = random_var[size:]
    
    for i in range(size):
        xy_t[i] += xc
        vT[i] += vc
    
    return xy_t


    #phase0 = np.average(phase[:,0],weights = amp[:,0])
    #phase1 = np.average(phase[:,1],weights = amp[:,1])
    #return (phase0+phase1)/16/n/(n+N)/2