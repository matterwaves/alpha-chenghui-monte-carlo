""" Wrapping C function Bragg using
    the Numpy declarations from Cython """

# cimport the Cython declarations for numpy
cimport numpy as np
import numpy as np

from ctypes import CDLL
from ctypes import RTLD_GLOBAL

gslcblas = CDLL('libgslcblas.dylib',mode=RTLD_GLOBAL)
gsl = CDLL('libgsl.dylib')


np.import_array()

# cdefine the signature of our c function
cdef extern from "Bragg.h":
     void Bragg (double * psi_real, double * psi_imag, double init_v, int picture0, int type0, int n0, double Omega10, double Omega20, double Delta0, double delta0, double pulsewidth0, double rampon0, double rampoff0, double omega_m0, double t0, double phic0, double stepsize, double * pulseshape, double * rampshape, int pulseshape_size, int rampshape_size, double * aoi, int aoi_size)

# create the wrapper code, with numpy type annotations
def CyBragg(args):

    
    cdef int picture, pulsetype, n, pulseshape_size, aoi_size
    cdef double Omega1, Omega2, Delta, delta, omega_m, t0, phic, init_v, pulsewidth, stepsize
    picture = args['picture']
    pulsetype = args['pulsetype']
    n = args['n']
    Omega1 = args['Omega1']
    Omega2 = args['Omega2']
    Delta = args['Delta']
    delta = args['delta']
    omega_m = args['omega_m']
    t0 = args['t0']
    phic = args['phic']
    init_v = args['init_v']
    stepsize = args['stepsize']
    pulsewidth = args['pulsewidth']
    HSize = n*3+11
    pulseshape_size = args['pulseshape_size']
    rampshape_size = args['rampshape_size']
    pulseshape0 = args['pulseshape']
    rampshape0 = args['rampshape']
    ramprate = args['ramprate']
    aoi0 = args['aoi']
    aoi_size = args['aoi_size']
    rampon = 0.01
    rampoff = 0.01
    if pulsetype > 1: 
        rampon = args['rampon']
        rampoff = args['rampoff']
        pulsewidth = args['Blochpulsewidth']+rampon+rampoff
        
    
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] psi_real = np.empty(HSize*aoi_size, dtype=np.double)
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] psi_imag = np.empty(HSize*aoi_size, dtype=np.double)
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] pulseshape = np.empty(pulseshape_size, dtype=np.double)
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] rampshape = np.empty(rampshape_size, dtype=np.double)
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] aoi = np.empty(aoi_size, dtype=np.double)

    for i in range(pulseshape_size):
        pulseshape[i] = pulseshape0[i]
	
    for i in range(rampshape_size):
        rampshape[i] = rampshape0[i]*ramprate
    
    for i in range(aoi_size):
        aoi[i] = aoi0[i]

    Bragg(<double*> np.PyArray_DATA(psi_real), <double*> np.PyArray_DATA(psi_imag), init_v, picture, pulsetype, n, Omega1, Omega2, Delta, delta, pulsewidth, rampon, rampoff, omega_m, t0, phic, stepsize, <double*> np.PyArray_DATA(pulseshape), <double*> np.PyArray_DATA(rampshape), pulseshape_size, rampshape_size, <double*> np.PyArray_DATA(aoi), aoi_size)

    return psi_real+1j*psi_imag