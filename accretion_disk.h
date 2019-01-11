#ifndef ACCRETION_DISK_h
#define ACCRETION_DISK_h


//int rescal_filters(char *filter, int binningnumber);

int make_computation(int Nfilter, long int *computed_filter, double *time, double *flux, double *ratio, double *tau_time, int Ntime, int Ntau);

int compute_LC(int filter_name, double *time_ILC, double *flux_ILC, int Ntime_ILC, double *t, double *flux, int Nt, double multiplicator);

int compute_LC_no_integration(double *time_ILC, double *flux_ILC, int Ntime_ILC, double *t, double *flux, int Nt, double multiplicator);

int disk_spectrum(double multiplicator, double *wavelength, double *spec, int len_spec);

#endif
