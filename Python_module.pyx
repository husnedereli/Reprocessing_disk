#""" Example of wrapping a C function that takes C double arrays as input using
#    the Numpy declarations from Cython """

# cimport the Cython declarations for numpy
import numpy
cimport numpy as np

# if you want to use the Numpy-C-API from Cython
# (not strictly necessary for this example, but good practice)
np.import_array()

# cdefine the signature of our c function
cdef extern from "CC.h":
    int CCooper(int Ngamma, int Ne, double Dt, double *Ee_bound, double *Eg_bound, double *Eg, double *fgamma, double *fgammanp1, double *fgamma_eq, double *fe, double *C)

    int evolv_e(int Ngamma, int Ne, double Dt, double *Ee, double *Ee_bound, double *Eg_bound, double *Eg, double *fgamma, double *fgammanp1, double *fgammaeq, double *fe, double *fenp1, double *C)

    int evolv_e_implicit(int Ngamma, int Ne, double Dt, double *Ee, double *Ee_bound, double *Eg_bound, double *Eg, double *fgamma, double *fgammanp1, double *fgamma_eq, double *fe, double *fenp1, double *C)

#    int evolv_e_diffusion(int Ngamma, int Ne, double Dt, double *Ee, double *Ee_bound, double *Gam, double *Gam_bound, double *Eg_bound, double *Eg, double *fgamma, double *fgammanp1, double *fgamma_eq,double *fe, double *fenp1, double *C)

    int evolv_e_CC_simple(int Ngamma, int Ne, double Dt, double *Ee, double *Ee_bound, double *Gam, double *Gam_bound, double *Eg_bound, double *Eg, double *fgamma, double *fgammanp1, double *fgamma_eq, double *fe, double *fenp1, double *fe_eq, double *C, double Etotinit)

    int evolv_e_CC_simple_iteration(int Ngamma, int Ne, double Dt, double *Ee, double *Ee_bound, double *Gam, double *Gam_bound, double *Eg_bound, double *Eg, double *fgamma,
    double *fgammanp1, double *fgamma_eq, double *fe, double *fenp1, double *fe_eq, double *C, double Etotinit)

# create the wrapper code, with numpy type annotations
def ChangCooper(Ngamma, Ne, Dt, np.ndarray[double, ndim=1, mode="c"] Ee_bound not None , np.ndarray[double, ndim=1, mode="c"] Eg_bound not None , np.ndarray[double, ndim=1, mode="c"] Eg not None, np.ndarray[double, ndim=1, mode="c"] fgamma not None, np.ndarray[double, ndim=1, mode="c"] fgammanp1 not None, np.ndarray[double, ndim=1, mode="c"] fgammaeq not None, np.ndarray[double, ndim=1, mode="c"] fe not None, np.ndarray[double, ndim=1, mode="c"] C not None):
    return CCooper(Ngamma, Ne, Dt, <double*> np.PyArray_DATA(Ee_bound), <double*> np.PyArray_DATA(Eg_bound), <double*> np.PyArray_DATA(Eg), <double*> np.PyArray_DATA(fgamma), <double*> np.PyArray_DATA(fgammanp1), <double*> np.PyArray_DATA(fgammaeq), <double*> np.PyArray_DATA(fe), <double*> np.PyArray_DATA(C))



def Evolve_e(Ngamma, Ne, Dt, np.ndarray[double, ndim=1, mode="c"] Ee not None , np.ndarray[double, ndim=1, mode="c"] Ee_bound not None , np.ndarray[double, ndim=1, mode="c"] Eg_bound not None , np.ndarray[double, ndim=1, mode="c"] Eg not None, np.ndarray[double, ndim=1, mode="c"] fgamma not None, np.ndarray[double, ndim=1, mode="c"] fgammanp1 not None, np.ndarray[double, ndim=1, mode="c"] fgamma_eq not None, np.ndarray[double, ndim=1, mode="c"] fe not None, np.ndarray[double, ndim=1, mode="c"] fenp1 not None, np.ndarray[double, ndim=1, mode="c"] C not None):
    return evolv_e(Ngamma, Ne, Dt, <double*> np.PyArray_DATA(Ee), <double*> np.PyArray_DATA(Ee_bound), <double*> np.PyArray_DATA(Eg_bound), <double*> np.PyArray_DATA(Eg), <double*> np.PyArray_DATA(fgamma), <double*> np.PyArray_DATA(fgammanp1), <double*> np.PyArray_DATA(fgamma_eq), <double*> np.PyArray_DATA(fe),<double*> np.PyArray_DATA(fenp1), <double*> np.PyArray_DATA(C))

def Evolve_e_implicit(Ngamma, Ne, Dt, np.ndarray[double, ndim=1, mode="c"] Ee not None , np.ndarray[double, ndim=1, mode="c"] Ee_bound not None , np.ndarray[double, ndim=1, mode="c"] Eg_bound not None , np.ndarray[double, ndim=1, mode="c"] Eg not None, np.ndarray[double, ndim=1, mode="c"] fgamma not None, np.ndarray[double, ndim=1, mode="c"] fgammanp1 not None, np.ndarray[double, ndim=1, mode="c"] fgamma_eq not None, np.ndarray[double, ndim=1, mode="c"] fe not None, np.ndarray[double, ndim=1, mode="c"] fenp1 not None, np.ndarray[double, ndim=1, mode="c"] C not None):
    return evolv_e_implicit(Ngamma, Ne, Dt, <double*> np.PyArray_DATA(Ee), <double*> np.PyArray_DATA(Ee_bound), <double*> np.PyArray_DATA(Eg_bound), <double*> np.PyArray_DATA(Eg), <double*> np.PyArray_DATA(fgamma), <double*> np.PyArray_DATA(fgammanp1), <double*> np.PyArray_DATA(fgamma_eq), <double*> np.PyArray_DATA(fe),<double*> np.PyArray_DATA(fenp1), <double*> np.PyArray_DATA(C))

#def Evolve_e_diffusion(Ngamma, Ne, Dt, np.ndarray[double, ndim=1, mode="c"] Ee not None , np.ndarray[double, ndim=1, mode="c"] Ee_bound not None , np.ndarray[double, ndim=1, mode="c"] Gam not None ,np.ndarray[double, ndim=1, mode="c"] Gam_bound not None , np.ndarray[double, ndim=1, mode="c"] Eg_bound not None , np.ndarray[double, ndim=1, mode="c"] Eg not None, np.ndarray[double, ndim=1, mode="c"] fgamma not None, np.ndarray[double, ndim=1, mode="c"] fgammanp1 not None, np.ndarray[double, ndim=1, mode="c"] fgamma_eq not None, np.ndarray[double, ndim=1, mode="c"] fe not None, np.ndarray[double, ndim=1, mode="c"] fenp1 not None, np.ndarray[double, ndim=1, mode="c"] C not None):
#    return evolv_e_diffusion(Ngamma, Ne, Dt, <double*> np.PyArray_DATA(Ee), <double*> np.PyArray_DATA(Ee_bound), <double*> np.PyArray_DATA(Gam), <double*> np.PyArray_DATA(Gam_bound), <double*> np.PyArray_DATA(Eg_bound), <double*> np.PyArray_DATA(Eg), <double*> np.PyArray_DATA(fgamma), <double*> np.PyArray_DATA(fgammanp1), <double*> np.PyArray_DATA(fgamma_eq), <double*> np.PyArray_DATA(fe),<double*> np.PyArray_DATA(fenp1), <double*> np.PyArray_DATA(C))


def Evolve_e_CC_simple(Ngamma, Ne, Dt,  np.ndarray[double, ndim=1, mode="c"] Ee not None ,
                                        np.ndarray[double, ndim=1, mode="c"] Ee_bound not None ,
                                        np.ndarray[double, ndim=1, mode="c"] Gam not None ,
                                        np.ndarray[double, ndim=1, mode="c"] Gam_bound not None ,
                                        np.ndarray[double, ndim=1, mode="c"] Eg_bound not None ,
                                        np.ndarray[double, ndim=1, mode="c"] Eg not None,
                                        np.ndarray[double, ndim=1, mode="c"] fgamma not None,
                                        np.ndarray[double, ndim=1, mode="c"] fgammanp1 not None,
                                        np.ndarray[double, ndim=1, mode="c"] fgamma_eq not None,
                                        np.ndarray[double, ndim=1, mode="c"] fe not None,
                                        np.ndarray[double, ndim=1, mode="c"] fenp1 not None,
                                        np.ndarray[double, ndim=1, mode="c"] fe_eq not None,
                                        np.ndarray[double, ndim=1, mode="c"] C not None,
                                        Etotinit):
    return evolv_e_CC_simple(Ngamma, Ne, Dt,    <double*> np.PyArray_DATA(Ee),
                                                <double*> np.PyArray_DATA(Ee_bound),
                                                <double*> np.PyArray_DATA(Gam),
                                                <double*> np.PyArray_DATA(Gam_bound),
                                                <double*> np.PyArray_DATA(Eg_bound),
                                                <double*> np.PyArray_DATA(Eg),
                                                <double*> np.PyArray_DATA(fgamma),
                                                <double*> np.PyArray_DATA(fgammanp1),
                                                <double*> np.PyArray_DATA(fgamma_eq),
                                                <double*> np.PyArray_DATA(fe),
                                                <double*> np.PyArray_DATA(fenp1),
                                                <double*> np.PyArray_DATA(fe_eq),
                                                <double*> np.PyArray_DATA(C),
                                                Etotinit)

def Evolve_e_CC_simple_iteration(Ngamma, Ne, Dt,  np.ndarray[double, ndim=1, mode="c"] Ee not None ,
                                        np.ndarray[double, ndim=1, mode="c"] Ee_bound not None ,
                                        np.ndarray[double, ndim=1, mode="c"] Gam not None ,
                                        np.ndarray[double, ndim=1, mode="c"] Gam_bound not None ,
                                        np.ndarray[double, ndim=1, mode="c"] Eg_bound not None ,
                                        np.ndarray[double, ndim=1, mode="c"] Eg not None,
                                        np.ndarray[double, ndim=1, mode="c"] fgamma not None,
                                        np.ndarray[double, ndim=1, mode="c"] fgammanp1 not None,
                                        np.ndarray[double, ndim=1, mode="c"] fgamma_eq not None,
                                        np.ndarray[double, ndim=1, mode="c"] fe not None,
                                        np.ndarray[double, ndim=1, mode="c"] fenp1 not None,
                                        np.ndarray[double, ndim=1, mode="c"] fe_eq not None,
                                        np.ndarray[double, ndim=1, mode="c"] C not None,
                                        Etotinit):
    return evolv_e_CC_simple_iteration(Ngamma, Ne, Dt,    <double*> np.PyArray_DATA(Ee),
                                                <double*> np.PyArray_DATA(Ee_bound),
                                                <double*> np.PyArray_DATA(Gam),
                                                <double*> np.PyArray_DATA(Gam_bound),
                                                <double*> np.PyArray_DATA(Eg_bound),
                                                <double*> np.PyArray_DATA(Eg),
                                                <double*> np.PyArray_DATA(fgamma),
                                                <double*> np.PyArray_DATA(fgammanp1),
                                                <double*> np.PyArray_DATA(fgamma_eq),
                                                <double*> np.PyArray_DATA(fe),
                                                <double*> np.PyArray_DATA(fenp1),
                                                <double*> np.PyArray_DATA(fe_eq),
                                                <double*> np.PyArray_DATA(C),
                                                Etotinit)















