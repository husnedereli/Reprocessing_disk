#ifndef ACCRETION_DISK_h
#define ACCRETION_DISK_h


//int rescal_filters(char *filter, int binningnumber);

int make_computation(int Nfilter, long int *computed_filter, double *time, double *flux, double *ratio, double *tau_time, int Ntime, int Ntau);

#endif
