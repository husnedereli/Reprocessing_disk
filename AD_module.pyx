import numpy
cimport numpy as np
import ctypes

np.import_array()

cdef extern from "accretion_disk.h":
    int make_computation(int Nfilter, long int *computed_filter, double *time, double *flux, double *ratio, double *tau_time, int Ntime, int Ntau)

    int compute_LC(int filter_name, double *time_ILC, double *flux_ILC, int Ntime_ILC, double *t, double *flux, int Nt, double multiplicator)

    int compute_LC_no_integration(double *time_ILC, double *flux_ILC, int Ntime_ILC, double *t, double *flux, int Nt, double multiplicator)

    int disk_spectrum(double multiplicator, double *wavelength, double *spec, int len_spec)



def AD_py( Nfilter, np.ndarray[long int, ndim=1, mode="c"] computed_filter not None, np.ndarray[double, ndim=1, mode="c"] time not None, np.ndarray[double, ndim=1, mode="c"] flux not None, np.ndarray[double, ndim=1, mode="c"] ratio not None, np.ndarray[double, ndim=1, mode="c"] tau_time not None, Ntime, Ntau):
    return make_computation(Nfilter, <long int*> np.PyArray_DATA(computed_filter), <double*> np.PyArray_DATA(time), <double*> np.PyArray_DATA(flux), <double*> np.PyArray_DATA(ratio), <double*> np.PyArray_DATA(tau_time), Ntime, Ntau)


def LightCurve( filter_name, np.ndarray[double, ndim=1, mode="c"] time_ILC not None, np.ndarray[double, ndim=1, mode="c"] flux_ILC not None, Ntime_ILC, np.ndarray[double, ndim=1, mode="c"] t not None, np.ndarray[double, ndim=1, mode="c"] flux not None, Nt, multiplicator):
    return compute_LC( filter_name, <double*> np.PyArray_DATA(time_ILC), <double*> np.PyArray_DATA(flux_ILC), Ntime_ILC, <double*> np.PyArray_DATA(t), <double*> np.PyArray_DATA(flux), Nt, multiplicator)

def LightCurve_no_integration(np.ndarray[double, ndim=1, mode="c"] time_ILC not None, np.ndarray[double, ndim=1, mode="c"] flux_ILC not None, Ntime_ILC, np.ndarray[double, ndim=1, mode="c"] t not None, np.ndarray[double, ndim=1, mode="c"] flux not None, Nt, multiplicator):
    return compute_LC_no_integration( <double*> np.PyArray_DATA(time_ILC), <double*> np.PyArray_DATA(flux_ILC), Ntime_ILC, <double*> np.PyArray_DATA(t), <double*> np.PyArray_DATA(flux), Nt, multiplicator)


def Spectrum( multiplicator, np.ndarray[double, ndim=1, mode="c"] wavelength not None, np.ndarray[double, ndim=1, mode="c"] spec not None, len_spec):
    return disk_spectrum(multiplicator, <double*> np.PyArray_DATA(wavelength), <double*> np.PyArray_DATA(spec), len_spec)
