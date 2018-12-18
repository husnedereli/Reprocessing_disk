import numpy
cimport numpy as np
import ctypes

np.import_array()

cdef extern from "accretion_disk.h":
    int make_computation(int Nfilter, long int *computed_filter, double *time, double *flux, double *ratio, double *tau_time, int Ntime, int Ntau)

# create the wrapper code, with numpy type annotations
#def Rescal_filters_py( binningnumber, the string):
#    return rescal_filters(the string, binningnumber)



# create the wrapper code, with numpy type annotations
def AD_py( Nfilter, np.ndarray[long int, ndim=1, mode="c"] computed_filter not None, np.ndarray[double, ndim=1, mode="c"] time not None, np.ndarray[double, ndim=1, mode="c"] flux not None, np.ndarray[double, ndim=1, mode="c"] ratio not None, np.ndarray[double, ndim=1, mode="c"] tau_time not None, Ntime, Ntau):
    #cf = np.ascontiguousarray(computed_filter, dtype=ctypes.c_int)
    return make_computation(Nfilter, <long int*> np.PyArray_DATA(computed_filter), <double*> np.PyArray_DATA(time), <double*> np.PyArray_DATA(flux), <double*> np.PyArray_DATA(ratio), <double*> np.PyArray_DATA(tau_time), Ntime, Ntau)
