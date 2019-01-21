#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "accretion_disk.h"

/**  Code specific declaration */
//#define Nr 66 			/** layers **/
//#define Ntheta 66		/** zones per layers **/

/**  Code specific declaration */
#define Nr 22 			/** layers **/
#define Ntheta 22		/** zones per layers **/

/** Definition of the physical constants **/
/**  Constant CGS */
#define c 2.99792458e10                 /** cm.s^-1         */
#define Msun 1.989e33                   /** g               */
#define h 6.6260755e-27                 /** erg.s           */
#define hc 1.9864474610385790e-16
#define kb 1.380657e-16                 /** erg.K^-1        */
#define sigmaSB 5.67051e-5              /** erg.cm^-2 K^-4 s^-1 */
#define Ggrav 6.67259e-8                /** cm^3.g^-1.s^-2  */
#define pc 3.08568e18                   /** cm              */
#define nm 1e-7                         /** cm              */
#define angstrom 1e-8                   /** cm              */

#define second_to_day 1.157407407e-5

/**  Other constant */
#define pi 3.14159265358979
#define degree_to_radian 0.0174532925

/**  Husne, 09/10/2018
  *  Here I define all the function to create and fill the disk
  */

///** define the Function of the Flux **/
int find_index(int Ntime, double t, double *time,  double tau, double *flux){
    /**  Husne, 13/12/2018
     *   Find the position of the t-tau in the light curve using Binary search algorithm
     */
    double t_tau = t-tau;
    int nb_loop = 0;
    /** Step1 */
    int L = time [0];
    int R = Ntime-1;
    /** Step2 */
    int m = 0;
    /** Check boundaries for element ouside of the array */
    if((t_tau)< time[0] || (t_tau)>time[R]){
        // printf("Undefined flux t-tau = %g\n",t_tau);
        return -1;
    }
    /** Check the boundaries */
    if(time[0] == (t_tau)){
        return 0;
    }
    if(time[R] == (t_tau)){
        return -1;
    }
    /** The element we are looking for is not at the boundary */
    while (L < R){
        /** Step3 */
        m = floor((L+R)/2);
        nb_loop +=1;
        if(time[m] == (t_tau)){
            return m;
        }
        if(time[m+1] == (t_tau)){
            return m+1;
        }
        /** Step4 */
        if (time[m] < (t_tau)){
            L = m ;
        }
        /** Step5 */
        if (time[m] > (t_tau)){
            R = m ;
        }
        if(R == L+1){
            return L;
        }
    }
    // printf("Undefined flux outside \n");
    return -1;
}




///** define the Function of the luminosity **/
double L_star(double L_bol, int Ntime, double t, double *time,  double tau, double *flux){///double omega, double t){
    //omega = 3*c/R_out

    /**  Look for the index */
    int index_flux = find_index(Ntime, t, time,  tau, flux);
    if(index_flux == -1){
        return -1.0;
    }

    /**  It is ok it is inside the time */

    double f = flux[index_flux]+(flux[index_flux+1]-flux[index_flux])*(t-tau-time[index_flux])/(time[index_flux+1]-time[index_flux]);;
    ///return 0.15*L_bol*f;
    return L_bol*f;
}

/** define the Function of the distance from the central variable source to disk elements **/
double r_star(double r, double h_star){
    return sqrt(pow(h_star,2.0)+pow(r,2.0));
}


/** define the Function of the temperature profile **/
double temp_profile(double t, double r, double rstar, double tau, double theta, double M, double M_rate, double r_in, double A, double h_star, double inc_angle, double L_bol, int Ntime, double *time, double *flux){

    /// Compute the time lag up to the radius. For speed purposed, it is now computed only one time in the main code.
    // double tau = sqrt(pow(h_star,2.0)+pow(r,2.0))+h_star*cos(inc_angle)-r*cos(theta*0.0174532925)*sin(inc_angle);
    // tau = tau/c;
    double Lstar = L_star(L_bol, Ntime, t, time,  tau, flux);
    if(Lstar < 0.0){
        return -1.0;
    }
    //double rstar = r_star(r, h_star);
    //printf("Contrib 1 = %g\t contrib 2 = %g\n ",((3.0*Ggrav*M*M_rate)/(8.0*pi*sigmaSB*pow(r,3.0)))*(1.0-sqrt(r_in/r)),  ((1.0-A)*(h_star*Lstar/(4.0*pi*sigmaSB*pow(rstar,3.0)))));
    //getchar();
    //printf()
    return pow(((3.0*Ggrav*M*M_rate)/(8.0*pi*sigmaSB*pow(r,3.0)))*(1.0-sqrt(r_in/r)) +((1.0-A)*(h_star*Lstar/(4.0*pi*sigmaSB*pow(rstar,3.0)))) ,0.25);
}

/** define the Function of the temperature profile when there is no illumination **/
double Temp_profile_steady(double r, double rstar, double M, double M_rate, double r_in){

    return pow(((3.0*Ggrav*M*M_rate)/(8.0*pi*sigmaSB*pow(r,3.0)))*(1.0-sqrt(r_in/r))  ,0.25);

}




/**  Husne, 09/10/2018
  *  Here I define all the function required to compute the spectra
  */

/** define the Planck Function **/
double Planck_Function(double lambda3, double lambda, double temperature){
    return ((2.0*hc)/lambda3)/(exp(hc/(lambda*kb*temperature))-1.0);
}

/** define the Function of predicted spectrum (SED) **/
double spectrum(double cos_inc_angle, double D2, double theta_in, double theta_out, double R_in, double R_out, double lambda3, double lambda, double temperature){
    return (cos_inc_angle/D2)*Planck_Function(lambda3, lambda, temperature)*(theta_out-theta_in)*0.5*(R_out*R_out - R_in*R_in);
}



/** define new type as regions and its elements to create disk **/
typedef struct region {
    double radius;
    double theta;
    double temp;
    double rstar;
    double tau;

    double temp_t;
    double temp_tptau;
} region;

/************************************* COMPUTE THE COLOR VARIABILITY ********************************************************/

int make_computation(int Nfilter, long int *computed_filter, double *time, double *flux, double *ratio, double *tau_time, int Ntime, int Ntau){


/*
    int ii;
    FILE *output;
    output = fopen("lc.txt","a");
    for(ii = 0; ii < Ntime; ii++){
        fprintf(output, "%g\t%g\n", time[ii], flux[ii]);
    }
    fclose(output);
*/

    /**  Husne, 9/10/2018
      *  Here I create and fill the disk. I compute the temperature and settle all the regions of the disk.
      */


    double *r;                      /** radius which is from the center of disk to the center of any region**/
    double *theta;                  /** azimuth angle which is from the origine to the r for any region**/
    r = (double *) calloc(Nr,sizeof(double));
    theta = (double *) calloc(Ntheta,sizeof(double));

    double M = 3.2e7*Msun;              /** M_sun, the black hole mass, converted to gr **/
    double Rg= (Ggrav*M)/(c*c);         /** gravitational radius **/
    double r_in= 6.0*Rg;                /** inner radius **/
    double r_out=10000*Rg;              /** outer radius **/


    double inc_angle = 45.0*0.0174532925;   /** inclination angle , converted to radian **/
    double cos_inc_angle = cos(inc_angle);  /** Cos of the inclination angle, avoid to recompute it all the time */
    double h_star = 10.0*Rg;                /** the vertical distance from the cetral variable source to disk **/
    double M_rate = 1.0*Msun/31557600.0;    /** M_sun yr^-1, the typical black hole accretion rate , converted to gr **/
                                            /** the numerical factor converts from year to second: we are working in cgs: cm gram second.*/
    double A = 0.5;                         /** the disk albedo **/
    double L_bol = 2.82e44;                 /** erg/s^-1, the bolometric luminosity **/

    /** Checking the values of the radii */
    // printf("Rg = %g\tR_int = %g\tR_out = %g \n", Rg, r_in, r_out);
    // getchar();

    /** create disks which contain the regions **/
    region *disk;
    disk = (region *) malloc(Nr*Ntheta*sizeof(region));

    /** the ratio of the outher and inner radius of each rings fixed **/
    double step = exp(log(r_out/r_in)/Nr);
    int i;
    for (i=0; i < Nr; i++){
        r[i] = r_in*pow(step,i);
    }

    for (i=0; i < Ntheta; i++){
        theta[i] = i*(360.0/Ntheta)*degree_to_radian;
    }

    /** fill the disks with elements (radius and theta) of regions **/
    int j;
    double tau;
    for (i=0; i < Nr; i++){
        for (j=0; j < Ntheta; j++){
            disk[i*Ntheta+j].radius = r[i];                 /** disk[0] region1, ... **/
            disk[i*Ntheta+j].theta = theta[j];
            disk[i*Ntheta+j].rstar = r_star(r[i], h_star);   /** disk[0] region1, ... **/

            /// Compute the time lag up to the radius.
            tau = sqrt(pow(h_star,2.0)+pow(r[i],2.0))+h_star*cos_inc_angle-r[i]*cos(theta[j])*sin(inc_angle);
            tau = (tau/c)*second_to_day;
            disk[i*Ntheta+j].tau = tau;
        }
    }




    //double L_star = 0.5*L_bol; /** the luminosity of central variable source **/
    //printf("Rg = %g\t r_star = %g\t M_rate = %g\t")








    /**  Husne,  9/10/2018
      *
      */
    /** for the computation of luminosity so it is for temperature **/
    ///double omega = 10.0*c/r_out;
    /** for the computation of the radiation from the disk. **/
    double D = 75.01*1e6*pc;                   /** Mpc distance from observer to the source, converted to cm **/
    double D2 = D*D;
    double R_in;
    double R_out;
    double theta_in;
    double theta_out;




    /** ************************************************
     * ************************************************
     * ************************************************
     * ************************************************
     * ************************************************
     */




    /**  Husne,  11/10/2018
     *  Convolotion with the filter bandpass.
     Read a txt file for U bandpass.
     */
    //double filtername[6] = {0, 1, 2, 3, 4, 5}; //* filter names: 0=UVW2, 1=UVM2, 2=UVW1, 3=U, 4=B, 5=V */
    //int Nfilter = 6;
    double **wavelength;
    double **wavelength3;
    wavelength = (double **) malloc(Nfilter*sizeof(double*)); //* create an array */
    wavelength3 = (double **) malloc(Nfilter*sizeof(double*)); //* create an array */
    double **transmission;
    transmission = (double **) malloc(Nfilter*sizeof(double*));
    double c1_filtername, c2_filtername;
    int numberofloop_filtername;
    int *numberofloop;
    numberofloop = (int*) calloc(Nfilter,sizeof(int)); //* create an array */

    for (j=0; j < Nfilter; j++){
        //printf("Begining of loop \t j = %d\tNfilter = %d\n", j, Nfilter);
        //getchar();
        FILE *input_filtername;
        numberofloop_filtername = 0.0;
        //*it is important, when the filter number given as a "0" make computation*/
        if (computed_filter[j] == 0){
            switch(j) {
                case 0 : //*it is UVW2 filter then*/
                    input_filtername=fopen("Filter/UVW2_binned5.txt","r");//* open a text file for reading */
                    /**  Here %lf means type double */
                    /// step 1 caunt the number of loop
                    while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                        if(c2_filtername > 0.01){
                            numberofloop_filtername = numberofloop_filtername + 1;                 //* caunt the number of loop */
                        }
                        /// Previous line is equivalent to      numberofloop_U = numberofloop_U + 1;                 //* caunt the number of loop */
                    }
                    numberofloop[j] = numberofloop_filtername;
                    fclose(input_filtername);
                    break;
                case 1 : //*it is UVM2 filter then*/
                    input_filtername=fopen("Filter/UVM2_binned5.txt","r");      //* open a text file for reading */
                    /**  Here %lf means type double */
                    /// step 1
                    while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                        if(c2_filtername > 0.01){
                            numberofloop_filtername = numberofloop_filtername + 1;                 //* caunt the number of loop */
                        }
                        /// Previous line is equivalent to      numberofloop_U = numberofloop_U + 1;                 //* caunt the number of loop */
                    }
                    numberofloop[j] = numberofloop_filtername;
                    fclose(input_filtername);
                    break;
                case 2 : //*it is UVW1 filter then*/
                    input_filtername=fopen("Filter/UVW1_binned5.txt","r");      //* open a text file for reading */
                    /**  Here %lf means type double */
                    /// step 1
                    while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                        if(c2_filtername > 0.01){
                            numberofloop_filtername = numberofloop_filtername + 1;                 //* caunt the number of loop */
                        }
                        /// Previous line is equivalent to      numberofloop_U = numberofloop_U + 1;                 //* caunt the number of loop */
                    }
                    numberofloop[j] = numberofloop_filtername;
                    fclose(input_filtername);
                    break;
                case 3 : //*it is U filter then*/
                    input_filtername=fopen("Filter/U_binned5.txt","r");      //* open a text file for reading */
                    /**  Here %lf means type double */
                    /// step 1
                    while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                        if(c2_filtername > 0.01){
                            numberofloop_filtername = numberofloop_filtername + 1;                 //* caunt the number of loop */
                        }
                        /// Previous line is equivalent to      numberofloop_U = numberofloop_U + 1;                 //* caunt the number of loop */
                    }
                    numberofloop[j] = numberofloop_filtername;
                    fclose(input_filtername);
                    break;
                case 4 : //*it is B filter then*/
                    input_filtername=fopen("Filter/B_binned5.txt","r");      //* open a text file for reading */
                    /**  Here %lf means type double */
                    /// step 1
                    while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                        if(c2_filtername > 0.01){
                            numberofloop_filtername = numberofloop_filtername + 1;                 //* caunt the number of loop */
                        }
                        /// Previous line is equivalent to      numberofloop_U = numberofloop_U + 1;                 //* caunt the number of loop */
                    }
                    numberofloop[j] = numberofloop_filtername;
                    fclose(input_filtername);
                    break;
                case 5 : //*it is V filter then*/
                    input_filtername=fopen("Filter/V_binned5.txt","r");      //* open a text file for reading */
                    /**  Here %lf means type double */
                    /// step 1
                    while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                        if(c2_filtername > 0.01){
                            numberofloop_filtername = numberofloop_filtername + 1;                 //* caunt the number of loop */
                        }
                        /// Previous line is equivalent to      numberofloop_U = numberofloop_U + 1;                 //* caunt the number of loop */
                    }
                    numberofloop[j] = numberofloop_filtername;
                    fclose(input_filtername);
                    break;
            }

        }

        /// step 2 to create arrays
        wavelength[j] = (double *) calloc(numberofloop[j],sizeof(double)); //* create an array */
        wavelength3[j] = (double *) calloc(numberofloop[j],sizeof(double)); //* create an array */
        transmission[j] = (double *) calloc(numberofloop[j],sizeof(double));

        /// step 3 fill the arrays
        //*it is important, when the filter number given as a "0" make computation*/
        if (computed_filter[j] == 0){
            switch(j) {
                case 0 : //*it is UVW2 filter then*/
                    input_filtername=fopen("Filter/UVW2_binned5.txt","r");//* open a text file for reading */
                    /**  Here %lf means type double */
                    /// step 1
                    i = 0;
                    while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                        /// i = i + 1 ;
                        if(c2_filtername > 0.01){
                            wavelength[j][i]= c1_filtername*angstrom;
                            wavelength3[j][i]= wavelength[j][i]*wavelength[j][i]*wavelength[j][i];
                            transmission[j][i]= c2_filtername;
                            i += 1;
                        }
                    }
                    fclose(input_filtername);
                    break;
                case 1 : //*it is UVM2 filter then*/
                    input_filtername=fopen("Filter/UVM2_binned5.txt","r");      //* open a text file for reading */
                    /**  Here %lf means type double */
                    /// step 1
                    i = 0;
                    while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                        /// i = i + 1 ;
                        if(c2_filtername > 0.01){
                            wavelength[j][i]= c1_filtername*angstrom;
                            wavelength3[j][i]= wavelength[j][i]*wavelength[j][i]*wavelength[j][i];
                            transmission[j][i]= c2_filtername;
                            i += 1;
                        }
                    }
                    fclose(input_filtername);
                    break;
                case 2 : //*it is UVW1 filter then*/
                    input_filtername=fopen("Filter/UVW1_binned5.txt","r");      //* open a text file for reading */
                    /**  Here %lf means type double */
                    /// step 1
                    i = 0;
                    while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                        /// i = i + 1 ;
                        if(c2_filtername > 0.01){
                            wavelength[j][i]= c1_filtername*angstrom;
                            wavelength3[j][i]= wavelength[j][i]*wavelength[j][i]*wavelength[j][i];
                            transmission[j][i]= c2_filtername;
                            i += 1;
                        }
                    }
                    fclose(input_filtername);
                    break;
                case 3 : //*it is U filter then*/
                    input_filtername=fopen("Filter/U_binned5.txt","r");      //* open a text file for reading */
                    /**  Here %lf means type double */
                    /// step 1
                    i = 0;
                    while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                        /// i = i + 1 ;
                        if(c2_filtername > 0.01){
                            wavelength[j][i]= c1_filtername*angstrom;
                            wavelength3[j][i]= wavelength[j][i]*wavelength[j][i]*wavelength[j][i];
                            transmission[j][i]= c2_filtername;
                            i += 1;
                        }
                    }
                    fclose(input_filtername);
                    break;
                case 4 : //*it is B filter then*/
                    input_filtername=fopen("Filter/B_binned5.txt","r");      //* open a text file for reading */
                    /**  Here %lf means type double */
                    /// step 1
                    i = 0;
                    while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                        /// i = i + 1 ;
                        if(c2_filtername > 0.01){
                            wavelength[j][i]= c1_filtername*angstrom;
                            wavelength3[j][i]= wavelength[j][i]*wavelength[j][i]*wavelength[j][i];
                            transmission[j][i]= c2_filtername;
                            i += 1;
                        }
                    }
                    fclose(input_filtername);
                    break;
                case 5 : //*it is V filter then*/
                    input_filtername=fopen("Filter/V_binned5.txt","r");      //* open a text file for reading */
                    /**  Here %lf means type double */
                    /// step 1
                    i = 0;
                    while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                        /// i = i + 1 ;
                        if(c2_filtername > 0.01){
                            wavelength[j][i]= c1_filtername*angstrom;
                            wavelength3[j][i]= wavelength[j][i]*wavelength[j][i]*wavelength[j][i];
                            transmission[j][i]= c2_filtername;
                            i += 1;
                        }
                    }
                    fclose(input_filtername);
                    break;
            }
        }
    }


    /*for (j=0;j<Nfilter;j++){
        printf("j = %d\t  %d\n",j, numberofloop[j]);
    }*/






    /** ************************************************
     * ************************************************
     * ************************************************
     * ************************************************
     * ************************************************
     */





    /**  Husne, 20/10/2018
     * I define the time and tau_time as arrays
     */
    // 1) SET the tau so give the value of tau
    //const int Ntau = 7;
    //double tau_time[7] = {3.0, 6.0, 10.0, 20.0, 40.0, 100.0, 200.0};

    //printf("%g\t\n",tau_time);  //* print the arrays */


    int Nt = 2200;           /** time[Ntime-1]-time[0]=1048.1724000000004 **/
    double *t;
    t = (double *) calloc(Nt,sizeof(double));
    for (i=0; i<Nt; i++){
        t[i] = ((double) i)/2.0;
        //printf("t=%g\n", t[i]);
    }

    // a) make the sum of all the S(tau, t)*1/Nt
    // i) compute all S for all t
    /**  Husne,  19/10/2018 *  compute the color variability and plot them */
    double avarage_SBU = 0.0;
    double S_BU;            /**color_variation of BU */



    /**  Husne,  18/10/2018 *  Now compute the integral for U band. */

    double deltaLambda_U;

    /**  Husne,  19/10/2018 *  Now compute the integral for B band.*/

    double deltaLambda_B;



    double flux_t_U;
    double flux_t_B;
    double flux_tptau_U;
    double flux_tptau_B;

    double Temperature_t;
    double Temperature_tptau;
    double Integral;

    double f_U_im1 = 0.0;
    double f_B_im1 = 0.0;
    double f_U_i = 0.0;
    double f_B_i = 0.0;


    int m;
    int n;
    int l;
    int k;
    double stepR = exp(log(r_out/r_in)/Nr);
    double sqrt_stepR = sqrt(stepR);
    double stepT = 360.0*degree_to_radian/Ntheta;

    int can_do_computation;
    int nb_computation;
    double stau;
    /** Loop for the Number of filter, because I need to compute the average S for all tau by using two bands */
    for (m=0;m<Nfilter;m++){
        for (n=m+1;n<Nfilter;n++){
            //*it is important, when the filter number given as a "0" make computation for all casses between two filters */
            if(computed_filter[m] == 0 && computed_filter[n] == 0){
                printf("\n");
                /** Loop for the tau, because I need to compute the average S for all tau */
                for (l=0; l < Ntau; l++){
                    S_BU = 0.0;
                    nb_computation = 0;
                    /** Loop for the time, because I need to compute an average with respect to time */
                    /// to put back k < Nt
                    for (k=0; k < Nt; k++){
                        if(k % 100 == 0){
                            printf("m = %d\tn=%d\tl=%d\ttau=%g\tk = %d\n", m, n, l, tau_time[l], k);
                        }
                        /**  I need to compute S for all t
                         *  I need to compute f(t+tau), f(t) for U band
                         */

                        flux_t_U = 0.0;
                        flux_t_B = 0.0;
                        flux_tptau_U = 0.0;
                        flux_tptau_B = 0.0;

                        /**  I need to check that for all elements of the disk I never obtained a negative temperature, which means that the computation can be done. */
                        can_do_computation = 0;
                        for (j=0; j < Nr*Ntheta; j++){
                            /// temperature at time t in U
                            Temperature_t = temp_profile(t[k], disk[j].radius, disk[j].rstar, disk[j].tau, disk[j].theta, M, M_rate, r_in, A, h_star, inc_angle, L_bol, Ntime, time, flux);
                            disk[j].temp_t = Temperature_t;
                            /// temperature at time t + tau in U
                            Temperature_tptau = temp_profile(t[k] + tau_time[l], disk[j].radius, disk[j].rstar, disk[j].tau, disk[j].theta, M, M_rate, r_in, A, h_star, inc_angle, L_bol, Ntime, time, flux);
                            disk[j].temp_tptau = Temperature_tptau;
                            if(Temperature_t < 0.0 || Temperature_tptau < 0.0){
                                can_do_computation = 1;
                            }
                        }

                        /** Loop for the radius and theta, because I need to compute the temparature and spectrum of disk */
                        /// f is the summ of contribution from all the disk elements.

                        if(can_do_computation == 0){
                            //nb_computation += 1;
                            for (j=0; j < Nr*Ntheta; j++){

                                R_in = disk[j].radius/sqrt_stepR;            /** from the center to the first layer of any region **/
                                R_out = disk[j].radius*sqrt_stepR;           /** from the center to the last layer of any region **/
                                theta_in = disk[j].theta - 0.5*stepT;         /** from the origine to the first layer of any region on the bottom**/
                                theta_out = disk[j].theta + 0.5*stepT;        /** from the origine to the last layer of any region on the top**/

                                /**  Now I compute the integral for the U-band */
                                /// temperature at time t in U
                                // stored in the disk structure for every time
                                //Temperature_t = temp_profile(t[k], disk[j].radius,disk[j].rstar, disk[j].tau, disk[j].theta, M, M_rate, r_in, A, h_star, inc_angle, L_bol, Ntime, time, flux);
                                Temperature_t = disk[j].temp_t;
                                // printf("temperature 1 = %g\t", Temperature_t);
                                /// Initialization of the sum to compute the integral over the filter
                                Integral = 0.0;

                                /** Loop for the band, because I need to compute the integral over bandpass */
                                ///numberofloop_U = numberofloop[m]
                                for(i = 1; i < numberofloop[m] ; i++){

                                    deltaLambda_U = (wavelength[m][i] - wavelength[m][i-1]);
                                    f_U_i = spectrum(cos_inc_angle, D2, theta_in, theta_out, R_in, R_out, wavelength3[m][i], wavelength[m][i], Temperature_t)*transmission[m][i];
                                    f_U_im1 = spectrum(cos_inc_angle, D2, theta_in, theta_out, R_in, R_out, wavelength3[m][i-1], wavelength[m][i-1], Temperature_t)*transmission[m][i-1];

                                    Integral += (f_U_i+f_U_im1)*deltaLambda_U/2.0;
                                }

                                flux_t_U = flux_t_U + Integral;
                                if(flux_t_U != flux_t_U){
                                    printf("Flux_t_u = NAN \t j = %d\t Integral = %g\n", j, Integral);
                                    getchar();
                                }

                                /// temperature at time t + tau in U
                                //Temperature_tptau = temp_profile(t[k] + tau_time[l], disk[j].radius, disk[j].rstar, disk[j].tau, disk[j].theta, M, M_rate, r_in, A, h_star, inc_angle, L_bol, Ntime, time, flux);
                                Temperature_tptau = disk[j].temp_tptau;

                                // printf("temperature 2 = %g\t", Temperature_tptau);
                                /// Initialization of the sum to compute the integral over the filter
                                Integral = 0.0;
                                /** Loop for the one band, because I need to compute the integral over bandpass */
                                ///numberofloop_U = numberofloop[m]
                                for(i = 1; i < numberofloop[m] ; i++){

                                    deltaLambda_U = (wavelength[m][i] - wavelength[m][i-1]);
                                    f_U_i    = spectrum(cos_inc_angle, D2, theta_in, theta_out, R_in, R_out, wavelength3[m][i], wavelength[m][i], Temperature_tptau)*transmission[m][i];
                                    f_U_im1  = spectrum(cos_inc_angle, D2, theta_in, theta_out, R_in, R_out, wavelength3[m][i-1], wavelength[m][i-1], Temperature_tptau)*transmission[m][i-1];

                                    Integral += (f_U_i+f_U_im1)*deltaLambda_U/2.0;
                                }
                                flux_tptau_U = flux_tptau_U + Integral;
                                if(flux_tptau_U != flux_tptau_U){
                                    printf("Flux_t_uptau = NAN \t j = %d\t Integral = %g\n", j, Integral);
                                    getchar();
                                }




                                /**  Now I compute the integral for the B-band */
                                /// temperature at time t in B, no need to recompute it
                                //Temperature_t = temp_profile(t[k], disk[j].radius, disk[j].rstar, disk[j].tau, disk[j].theta, M, M_rate, r_in, A, h_star, inc_angle, L_bol, Ntime, time, flux);
                                // printf("temperature 3 = %g\t", Temperature_t);
                                ///if(Temperature_t!=Temperature_t){
                                ///printf("%.13g\t\n", Temperature_t);
                                ///getchar();
                                ///}

                                /// Initialization of the sum to compute the integral over the filter
                                Integral = 0.0;
                                /** Loop for the another band, because I need to compute the integral over bandpass */
                                ///numberofloop_B = numberofloop[n]
                                for(i = 1; i < numberofloop[n] ; i++){

                                    deltaLambda_B = (wavelength[n][i] - wavelength[n][i-1]);
                                    f_B_i    = spectrum(cos_inc_angle, D2, theta_in, theta_out, R_in, R_out, wavelength3[n][i], wavelength[n][i], Temperature_t)*transmission[n][i];
                                    f_B_im1  = spectrum(cos_inc_angle, D2, theta_in, theta_out, R_in, R_out, wavelength3[n][i-1], wavelength[n][i-1], Temperature_t)*transmission[n][i-1];

                                    Integral += (f_B_i+f_B_im1)*deltaLambda_B/2.0;
                                }

                                flux_t_B = flux_t_B + Integral;
                                if(flux_t_B != flux_t_B){
                                    printf("Flux_t_b = NAN \t j = %d\t Integral = %g\n", j, Integral);
                                    getchar();
                                }

                                /// temperature at time t + tau in U
                                // No need to recompute
                                //Temperature_tptau = temp_profile(t[k] + tau_time[l], disk[j].radius, disk[j].rstar, disk[j].tau, disk[j].theta, M, M_rate, r_in, A, h_star, inc_angle, L_bol, Ntime, time, flux);
                                // printf("temperature 4 = %g\n", Temperature_tptau);
                                /// Initialization of the sum to compute the integral over the filter
                                Integral = 0.0;
                                /** Loop for the another band, because I need to compute the integral over bandpass */
                                ///numberofloop_B = numberofloop[n]
                                for(i = 1; i < numberofloop[n] ; i++){

                                    deltaLambda_B = (wavelength[n][i] - wavelength[n][i-1]);
                                    f_B_i    = spectrum(cos_inc_angle, D2, theta_in, theta_out, R_in, R_out, wavelength3[n][i], wavelength[n][i], Temperature_tptau)*transmission[n][i];
                                    f_B_im1  = spectrum(cos_inc_angle, D2, theta_in, theta_out, R_in, R_out, wavelength3[n][i-1], wavelength[n][i-1], Temperature_tptau)*transmission[n][i-1];

                                    Integral += (f_B_i+f_B_im1)*deltaLambda_B/2.0;
                                }
                                flux_tptau_B = flux_tptau_B + Integral;
                                if(flux_tptau_B != flux_tptau_B){
                                    printf("Flux_tptau_b = NAN \t j = %d\t Integral = %g\n", j, Integral);
                                    getchar();
                                }
                            }
                            stau =  (flux_tptau_B-flux_t_B) / (flux_tptau_U-flux_t_U);
                            if(stau > 0.0){
                                S_BU = S_BU + stau;
                                nb_computation +=1;
                            }
                            //printf("%.13g\t\n", flux_t_B);
                            //printf("%.13g\t\n", flux_t_U);
                        }
                    }
                    //avarage_SBU = S_BU/Nt;
                    avarage_SBU = S_BU/ ((double) nb_computation);
                    ratio[l] = avarage_SBU;
                    printf("%.13g\t\n", avarage_SBU);
                }
            }
        }
    }

    for (j=0; j < Nfilter; j++){
        free(wavelength[j]);
        free(wavelength3[j]);
        free(transmission[j]);
    }

    free(r);
    free(theta);
    free(disk);
    free(wavelength);
    free(wavelength3);
    free(transmission);


    return 0;
}



/******************************************** COMPUTE LIGHT CURVE ****************************************************************/


/**  What are the parameters ?
  *     1- filter_name : name of the filter
  *     2- time_ILC : time of the illuminating light-curve
  *     3- flux_ILC : flux of the illuminating light-curve
  *     4- Ntime_ILC : length of array time_ILC and flux_ILC
  *     5- t : the time at which the LC should be computed
  *     6- flux :  the vector to be filled by the light-curve
  *     7- Nt :  the length of those arrays
  */
int compute_LC(int filter_name, double *time_ILC, double *flux_ILC, int Ntime_ILC, double *t, double *flux, int Nt, double multiplicator){


    /**  Husne, 9/10/2018
      *  Here I create and fill the disk. I compute the temperature and settle all the regions of the disk.
      */


    double *r;                        /** radius which is from the center of disk to the center of any region**/
    double *theta;                    /** azimuth angle which is from the origine to the r for any region**/
    r = (double *) calloc(Nr,sizeof(double));
    theta = (double *) calloc(Ntheta,sizeof(double));

    double M = 3.2e7*Msun;              /** M_sun, the black hole mass, converted to gr **/
    double Rg= (Ggrav*M)/(c*c);         /** gravitational radius **/
    double r_in= 6.0*Rg;                /** inner radius **/
    double r_out=10000*Rg;              /** outer radius **/


    double inc_angle = 45.0*degree_to_radian;   /** inclination angle , converted to radian **/
    double cos_inc_angle = cos(inc_angle);  /** Cos of the inclination angle, avoid to recompute it all the time */
    double h_star = 10.0*Rg;                /** the vertical distance from the cetral variable source to disk **/
    double M_rate = 1.0*Msun/31557600.0;    /** M_sun yr^-1, the typical black hole accretion rate , converted to g.s^-1 **/
    /** the numerical factor converts from year to second: we are working in cgs: cm gram second.*/
    double A = 0.5;                         /** the disk albedo **/
    double L_bol = multiplicator*2.82e44;                 /** erg/s^-1, the bolometric luminosity **/

    /** Checking the values of the radii */
    // printf("Rg = %g\tR_int = %g\tR_out = %g \n", Rg, r_in, r_out);
    // getchar();

    /** create disks which contain the regions **/
    region *disk;
    disk = (region *) malloc(Nr*Ntheta*sizeof(region));

    /** the ratio of the outher and inner radius of each rings fixed **/
    double step = exp(log(r_out/r_in)/Nr);
    int i;
    for (i=0; i < Nr; i++){
        r[i] = r_in*pow(step,(double) i);
    }

    for (i=0; i < Ntheta; i++){
        theta[i] = ((double) i)*(360.0/Ntheta)*degree_to_radian;
    }

    /** fill the disks with elements (radius and theta) of regions **/
    int j;
    double tau;
    for (i=0; i < Nr; i++){
        for (j=0; j < Ntheta; j++){
            disk[i*Ntheta+j].radius = r[i];                 /** disk[0] region1, ... **/
            disk[i*Ntheta+j].theta = theta[j];
            disk[i*Ntheta+j].rstar = r_star(r[i], h_star);   /** disk[0] region1, ... **/

            /// Compute the time lag up to the radius.
            tau = sqrt(pow(h_star,2.0)+pow(r[i],2.0))+h_star*cos_inc_angle-r[i]*cos(theta[j])*sin(inc_angle);
            tau = (tau/c)*second_to_day;
            disk[i*Ntheta+j].tau = tau;
        }
    }


    /**  Husne,  9/10/2018
     *
     */
    /** for the computation of luminosity so it is for temperature **/
    ///double omega = 10.0*c/r_out;
    /** for the computation of the radiation from the disk. **/
    double D = 75.01*1e6*pc;                   /** Mpc distance from observer to the source, converted to cm **/
    double D2 = D*D;
    double R_in;
    double R_out;
    double theta_in;
    double theta_out;


    /** ************************************************
     * ************************************************
     * ************************************************
     * ************************************************
     * ************************************************
     */




    /**  Husne,  11/10/2018
     *  Convolotion with the filter bandpass.
     Read a txt file for U bandpass.
     */
    //double filtername[6] = {0, 1, 2, 3, 4, 5}; //* filter names: 0=UVW2, 1=UVM2, 2=UVW1, 3=U, 4=B, 5=V */
    //int Nfilter = 6;
    double *wavelength;
    double *wavelength3;
    double *transmission;

    double c1_filtername, c2_filtername;
    int numberofloop = 0;

    FILE *input_filtername;
    //*it is important, when the filter number given as a "0" make computation*/

    switch(filter_name) {
        case 0 : //*it is UVW2 filter then*/
            input_filtername=fopen("Filter/UVW2_binned5.txt","r");//* open a text file for reading */
            /**  Here %lf means type double */
            /// step 1 caunt the number of loop
            while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                if(c2_filtername > 0.01){
                    numberofloop = numberofloop + 1;                 //* caunt the number of loop */
                }
                /// Previous line is equivalent to      numberofloop_U = numberofloop_U + 1;                 //* caunt the number of loop */
            }
            fclose(input_filtername);
            break;
        case 1 : //*it is UVM2 filter then*/
            input_filtername=fopen("Filter/UVM2_binned5.txt","r");      //* open a text file for reading */
            /**  Here %lf means type double */
            /// step 1
            while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                if(c2_filtername > 0.01){
                    numberofloop = numberofloop + 1;                 //* caunt the number of loop */
                }
                /// Previous line is equivalent to      numberofloop_U = numberofloop_U + 1;                 //* caunt the number of loop */
            }
            fclose(input_filtername);
            break;
        case 2 : //*it is UVW1 filter then*/
            input_filtername=fopen("Filter/UVW1_binned5.txt","r");      //* open a text file for reading */
            /**  Here %lf means type double */
            /// step 1
            while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                if(c2_filtername > 0.01){
                    numberofloop = numberofloop + 1;                 //* caunt the number of loop */
                }
                /// Previous line is equivalent to      numberofloop_U = numberofloop_U + 1;                 //* caunt the number of loop */
            }
            fclose(input_filtername);
            break;
        case 3 : //*it is U filter then*/
            input_filtername=fopen("Filter/U_binned5.txt","r");      //* open a text file for reading */
            /**  Here %lf means type double */
            /// step 1
            while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                if(c2_filtername > 0.01){
                    numberofloop = numberofloop + 1;                 //* caunt the number of loop */
                }
                /// Previous line is equivalent to      numberofloop_U = numberofloop_U + 1;                 //* caunt the number of loop */
            }
            fclose(input_filtername);
            break;
        case 4 : //*it is B filter then*/
            input_filtername=fopen("Filter/B_binned5.txt","r");      //* open a text file for reading */
            /**  Here %lf means type double */
            /// step 1
            while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                if(c2_filtername > 0.01){
                    numberofloop = numberofloop + 1;                 //* caunt the number of loop */
                }
                /// Previous line is equivalent to      numberofloop_U = numberofloop_U + 1;                 //* caunt the number of loop */
            }
            fclose(input_filtername);
            break;
        case 5 : //*it is V filter then*/
            input_filtername=fopen("Filter/V_binned5.txt","r");      //* open a text file for reading */
            /**  Here %lf means type double */
            /// step 1
            while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                if(c2_filtername > 0.01){
                    numberofloop = numberofloop + 1;                 //* caunt the number of loop */
                }
                /// Previous line is equivalent to      numberofloop_U = numberofloop_U + 1;                 //* caunt the number of loop */
            }
            fclose(input_filtername);
            break;
    }
    /**  DAMIEN : WE STOP LOOKING AT THE PB HERE */

    /// step 2 to create arrays
    wavelength = (double *) calloc(numberofloop,sizeof(double)); //* create an array */
    wavelength3 = (double *) calloc(numberofloop,sizeof(double)); //* create an array */
    transmission = (double *) calloc(numberofloop,sizeof(double));

    /// step 3 fill the arrays
    //*it is important, when the filter number given as a "0" make computation*/
    switch(filter_name) {
        case 0 : //*it is UVW2 filter then*/
            input_filtername=fopen("Filter/UVW2_binned5.txt","r");//* open a text file for reading */
            /**  Here %lf means type double */
            /// step 1
            i = 0;
            while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                /// i = i + 1 ;
                if(c2_filtername > 0.01){
                    wavelength[i]= c1_filtername*angstrom;
                    wavelength3[i]= wavelength[i]*wavelength[i]*wavelength[i];
                    transmission[i]= c2_filtername;
                    i += 1;
                }
            }
            fclose(input_filtername);
            break;
        case 1 : //*it is UVM2 filter then*/
            input_filtername=fopen("Filter/UVM2_binned5.txt","r");      //* open a text file for reading */
            /**  Here %lf means type double */
            /// step 1
            i = 0;
            while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                /// i = i + 1 ;
                if(c2_filtername > 0.01){
                    wavelength[i]= c1_filtername*angstrom;
                    wavelength3[i]= wavelength[i]*wavelength[i]*wavelength[i];
                    transmission[i]= c2_filtername;
                    i += 1;
                }
            }
            fclose(input_filtername);
            break;
        case 2 : //*it is UVW1 filter then*/
            input_filtername=fopen("Filter/UVW1_binned5.txt","r");      //* open a text file for reading */
            /**  Here %lf means type double */
            /// step 1
            i = 0;
            while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                /// i = i + 1 ;
                if(c2_filtername > 0.01){
                    wavelength[i]= c1_filtername*angstrom;
                    wavelength3[i]= wavelength[i]*wavelength[i]*wavelength[i];
                    transmission[i]= c2_filtername;
                    i += 1;
                }
            }
            fclose(input_filtername);
            break;
        case 3 : //*it is U filter then*/
            input_filtername=fopen("Filter/U_binned5.txt","r");      //* open a text file for reading */
            /**  Here %lf means type double */
            /// step 1
            i = 0;
            while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                /// i = i + 1 ;
                if(c2_filtername > 0.01){
                    wavelength[i]= c1_filtername*angstrom;
                    wavelength3[i]= wavelength[i]*wavelength[i]*wavelength[i];
                    transmission[i]= c2_filtername;
                    i += 1;
                }
            }
            fclose(input_filtername);
            break;
        case 4 : //*it is B filter then*/
            input_filtername=fopen("Filter/B_binned5.txt","r");      //* open a text file for reading */
            /**  Here %lf means type double */
            /// step 1
            i = 0;
            while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                /// i = i + 1 ;
                if(c2_filtername > 0.01){
                    wavelength[i]= c1_filtername*angstrom;
                    wavelength3[i]= wavelength[i]*wavelength[i]*wavelength[i];
                    transmission[i]= c2_filtername;
                    i += 1;
                }
            }
            fclose(input_filtername);
            break;
        case 5 : //*it is V filter then*/
            input_filtername=fopen("Filter/V_binned5.txt","r");      //* open a text file for reading */
            /**  Here %lf means type double */
            /// step 1
            i = 0;
            while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                /// i = i + 1 ;
                if(c2_filtername > 0.01){
                    wavelength[i]= c1_filtername*angstrom;
                    wavelength3[i]= wavelength[i]*wavelength[i]*wavelength[i];
                    transmission[i]= c2_filtername;
                    i += 1;
                }
            }
            fclose(input_filtername);
            break;
    }


    /*for (j=0;j<Nfilter;j++){
     printf("j = %d\t  %d\n",j, numberofloop[j]);
     }*/






    /** ************************************************
     * ************************************************
     * ************************************************
     * ************************************************
     * ************************************************
     */





    /**  Husne, 20/10/2018
     * I define the time and tau_time as arrays
     */
    // 1) SET the tau so give the value of tau
    //const int Ntau = 7;
    //double tau_time[7] = {3.0, 6.0, 10.0, 20.0, 40.0, 100.0, 200.0};

    //printf("%g\t\n",tau_time);  //* print the arrays */




    /**  Husne, 15/11/2018
      *  Compute all the other light curve in different filters (UVM2, UVW1, U, B, V)
      *  by using the illuminating light curves of UVW2
      */

    /**  Husne,  18/10/2018 *  Now compute the integral for U band. */
    double deltaLambda_U;
    double *flux_t_U;
    double Temperature_t;
    double Integral;
    double f_U_im1 = 0.0;
    double f_U_i = 0.0;
    ///int Ntime = 0;




    FILE *output;
    int m;
    int k;
    double stepR = exp(log(r_out/r_in)/Nr);
    double sqrt_stepR = sqrt(stepR);
    double stepT = 360.0*degree_to_radian/Ntheta;

    int can_do_computation;
    ///int nb_computation;
    /** Loop for the Number of filter, because I need to compute the flux for each band */

    /** Loop for the time */
    for (k=0; k < Nt; k++){
        if(k % 100 == 0){
            printf("m = %d\t\tk = %d\n", m, k);
        }

        /**
         *  I need to compute f(t) in the required band
         */

        /**  I need to check that for all elements of the disk I never obtained a negative temperature, which means that the computation can be done. */
        can_do_computation = 0;
        for (j=0; j < Nr*Ntheta; j++){
            /// temperature at time t in U
            Temperature_t = temp_profile(   t[k], disk[j].radius, disk[j].rstar, disk[j].tau, disk[j].theta, M, M_rate, r_in, A, h_star, inc_angle, L_bol, Ntime_ILC, time_ILC, flux_ILC);
            disk[j].temp_t = Temperature_t;
            if(Temperature_t < 0.0){
                can_do_computation = 1;
            }
        }

        /** Loop for the radius and theta, because I need to compute the temparature and spectrum of disk */
        /// f is the summ of contribution from all the disk elements.

        if(can_do_computation == 0){
            //nb_computation += 1;
            for (j=0; j < Nr*Ntheta; j++){

                R_in = disk[j].radius/sqrt_stepR;            /** from the center to the first layer of any region **/
                R_out = disk[j].radius*sqrt_stepR;           /** from the center to the last layer of any region **/
                theta_in = disk[j].theta - 0.5*stepT;         /** from the origine to the first layer of any region on the bottom**/
                theta_out = disk[j].theta + 0.5*stepT;        /** from the origine to the last layer of any region on the top**/

                /**  Now I compute the integral for the U-band */
                /// temperature at time t in U
                Temperature_t = disk[j].temp_t;
                /// Initialization of the sum to compute the integral over the filter
                Integral = 0.0;

                /** Loop for the band, because I need to compute the integral over bandpass */
                ///numberofloop_U = numberofloop[m]
                for(i = 1; i < numberofloop ; i++){

                    deltaLambda_U = (wavelength[i] - wavelength[i-1]);
                    f_U_i = spectrum(cos_inc_angle, D2, theta_in, theta_out, R_in, R_out, wavelength3[i], wavelength[i], Temperature_t)*transmission[i];
                    f_U_im1 = spectrum(cos_inc_angle, D2, theta_in, theta_out, R_in, R_out, wavelength3[i-1], wavelength[i-1], Temperature_t)*transmission[i-1];

                    Integral += (f_U_i+f_U_im1)*deltaLambda_U/2.0;
                }
                flux[k] = flux[k] + Integral;


                ///nb_computation +=1;
                ///printf("%.13g\t\n", flux_t_U[k]);
                /*if(flux_t_U[k] != flux_t_U[k]){
                    printf("Flux_t_u = NAN \t j = %d\t Integral = %g\n", j, Integral);
                    getchar();
                }*/
            }

        }


        /*printf("time = %g\t\tflux = %.13g\t\n", t[k], flux_t_U[k]);
        output = fopen("lc_U_disk.txt","a");
        ///if(flux_t_U[k] > 0){
        fprintf(output, "%g\t%g\n", t[k], flux_t_U[k]);
        ///} */

    }




//    fclose(output);


    free(r);
    free(theta);
    free(disk);
    free(wavelength);
    free(wavelength3);
    free(transmission);

    return 0;
}


/******************************************** COMPUTE The LIGHT CURVE Without INTEGRATION OVER BAND PASSES **************************/



/**  What are the parameters ?
 *     1- time_ILC : time of the illuminating light-curve
 *     2- flux_ILC : flux of the illuminating light-curve
 *     3- Ntime_ILC : length of array time_ILC and flux_ILC
 *     4- t : the time at which the LC should be computed
 *     5- flux :  the vector to be filled by the light-curve
 *     6- Nt :  the length of those arrays
 */
int compute_LC_no_integration(double *time_ILC, double *flux_ILC, int Ntime_ILC, double *t, double *flux, int Nt, double multiplicator){


    /**  Husne, 9/10/2018
     *  Here I create and fill the disk. I compute the temperature and settle all the regions of the disk.
     */


    double *r;                        /** radius which is from the center of disk to the center of any region**/
    double *theta;                    /** azimuth angle which is from the origine to the r for any region**/
    r = (double *) calloc(Nr,sizeof(double));
    theta = (double *) calloc(Ntheta,sizeof(double));

    double M = 3.2e7*Msun;              /** M_sun, the black hole mass, converted to gr **/
    double Rg= (Ggrav*M)/(c*c);         /** gravitational radius **/
    double r_in= 6.0*Rg;                /** inner radius **/
    double r_out=10000*Rg;              /** outer radius **/


    double inc_angle = 45.0*degree_to_radian;   /** inclination angle , converted to radian **/
    double cos_inc_angle = cos(inc_angle);  /** Cos of the inclination angle, avoid to recompute it all the time */
    double h_star = 10.0*Rg;                /** the vertical distance from the cetral variable source to disk **/
    double M_rate = 1.0*Msun/31557600.0;    /** M_sun yr^-1, the typical black hole accretion rate , converted to g.s^-1 **/
    /** the numerical factor converts from year to second: we are working in cgs: cm gram second.*/
    double A = 0.5;                         /** the disk albedo **/
    double L_bol = multiplicator*2.82e44;                 /** erg/s^-1, the bolometric luminosity **/

    /** Checking the values of the radii */
    // printf("Rg = %g\tR_int = %g\tR_out = %g \n", Rg, r_in, r_out);
    // getchar();

    /** create disks which contain the regions **/
    region *disk;
    disk = (region *) malloc(Nr*Ntheta*sizeof(region));

    /** the ratio of the outher and inner radius of each rings fixed **/
    double step = exp(log(r_out/r_in)/Nr);
    int i;
    for (i=0; i < Nr; i++){
        r[i] = r_in*pow(step,(double) i);
    }

    for (i=0; i < Ntheta; i++){
        theta[i] = ((double) i)*(360.0/Ntheta)*degree_to_radian;
    }

    /** fill the disks with elements (radius and theta) of regions **/
    int j;
    double tau;
    for (i=0; i < Nr; i++){
        for (j=0; j < Ntheta; j++){
            disk[i*Ntheta+j].radius = r[i];                 /** disk[0] region1, ... **/
            disk[i*Ntheta+j].theta = theta[j];
            disk[i*Ntheta+j].rstar = r_star(r[i], h_star);   /** disk[0] region1, ... **/

            /// Compute the time lag up to the radius.
            tau = sqrt(pow(h_star,2.0)+pow(r[i],2.0))+h_star*cos_inc_angle-r[i]*cos(theta[j])*sin(inc_angle);
            tau = (tau/c)*second_to_day;
            disk[i*Ntheta+j].tau = tau;
        }
    }


    /**  Husne,  9/10/2018
     *
     */
    /** for the computation of luminosity so it is for temperature **/
    ///double omega = 10.0*c/r_out;
    /** for the computation of the radiation from the disk. **/
    double D = 75.01*1e6*pc;                   /** Mpc distance from observer to the source, converted to cm **/
    double D2 = D*D;
    double R_in;
    double R_out;
    double theta_in;
    double theta_out;


    /** ************************************************
     * ************************************************
     * ************************************************
     * ************************************************
     * ************************************************
     */




    /**  Husne,  11/01/2019
     *  Central value of UVW2 = 212 nm = 212*1e-7 cm
     *  Central value of UVM2 = 231 nm = 231*1e-7 cm
     *  Central value of UVW1 = 291 nm = 291*1e-7 cm
     *  Central value of U = 344 nm = 344*1e-7 cm
     *  Central value of B = 450 nm = 450*1e-7 cm
     *  Central value of V = 543 nm = 543*1e-7 cm
     */

    //double wavelength_UVW2_3 = 212*nm*212*nm*212*nm;
    //double wavelength_UVW2 = 212*nm;
    double wavelength_UVW2_3 = 231*nm*231*nm*231*nm;
    double wavelength_UVW2 = 231*nm;

    /** ************************************************
     * ************************************************
     * ************************************************
     * ************************************************
     * ************************************************
     */





    /**  Husne, 20/10/2018
     * I define the time and tau_time as arrays
     */
    // 1) SET the tau so give the value of tau
    //const int Ntau = 7;
    //double tau_time[7] = {3.0, 6.0, 10.0, 20.0, 40.0, 100.0, 200.0};

    //printf("%g\t\n",tau_time);  //* print the arrays */




    /**  Husne, 15/11/2018
     *  Compute all the other light curve in different filters (UVM2, UVW1, U, B, V)
     *  by using the illuminating light curves of UVW2
     */

    /**  Husne,  18/10/2018 *  Now compute the integral for U band. */
    double deltaLambda_U;
    double *flux_t_U;
    double Temperature_t;
    double Integral;
    double f_U_im1 = 0.0;
    double f_U_i = 0.0;
    ///int Ntime = 0;




    FILE *output;
    int m;
    int k;
    double stepR = exp(log(r_out/r_in)/Nr);
    double sqrt_stepR = sqrt(stepR);
    double stepT = 360.0*degree_to_radian/Ntheta;

    int can_do_computation;
    ///int nb_computation;
    /** Loop for the Number of filter, because I need to compute the flux for each band */

    /** Loop for the time */
    for (k=0; k < Nt; k++){
        if(k % 100 == 0){
            printf("m = %d\t\tk = %d\n", m, k);
        }

        /**
         *  I need to compute f(t) in the required band
         */

        /**  I need to check that for all elements of the disk I never obtained a negative temperature, which means that the computation can be done. */
        can_do_computation = 0;
        for (j=0; j < Nr*Ntheta; j++){
            /// temperature at time t in U
            Temperature_t = temp_profile(   t[k], disk[j].radius, disk[j].rstar, disk[j].tau, disk[j].theta, M, M_rate, r_in, A, h_star, inc_angle, L_bol, Ntime_ILC, time_ILC, flux_ILC);
            disk[j].temp_t = Temperature_t;
            if(Temperature_t < 0.0){
                can_do_computation = 1;
            }
        }

        /** Loop for the radius and theta, because I need to compute the temparature and spectrum of disk */
        /// f is the summ of contribution from all the disk elements.

        if(can_do_computation == 0){
            //nb_computation += 1;
            for (j=0; j < Nr*Ntheta; j++){

                R_in = disk[j].radius/sqrt_stepR;            /** from the center to the first layer of any region **/
                R_out = disk[j].radius*sqrt_stepR;           /** from the center to the last layer of any region **/
                theta_in = disk[j].theta - 0.5*stepT;         /** from the origine to the first layer of any region on the bottom**/
                theta_out = disk[j].theta + 0.5*stepT;        /** from the origine to the last layer of any region on the top**/

                /**  Now I compute the integral for the U-band */
                /// temperature at time t in U
                Temperature_t = disk[j].temp_t;
                /// Initialization of the sum to compute the integral over the filter

                /** Loop for the band, because I need to compute the integral over bandpass */
                ///numberofloop_U = numberofloop[m]
                ////for(i = 1; i < numberofloop ; i++){

                f_U_i = spectrum(cos_inc_angle, D2, theta_in, theta_out, R_in, R_out, wavelength_UVW2_3, wavelength_UVW2, Temperature_t);
                    ///f_U_im1 = spectrum(cos_inc_angle, D2, theta_in, theta_out, R_in, R_out, wavelength_UVW2_3, wavelength_UVW2, Temperature_t);
                flux[k] = flux[k] + f_U_i;
                ///nb_computation +=1;
                ///printf("%.13g\t\n", flux_t_U[k]);
                /*if(flux_t_U[k] != flux_t_U[k]){
                 printf("Flux_t_u = NAN \t j = %d\t Integral = %g\n", j, Integral);
                 getchar();
                 }*/
            }
        }
        /*printf("time = %g\t\tflux = %.13g\t\n", t[k], flux_t_U[k]);
         output = fopen("lc_U_disk.txt","a");
         ///if(flux_t_U[k] > 0){
         fprintf(output, "%g\t%g\n", t[k], flux_t_U[k]);
         ///} */
    }
    
    free(r);
    free(theta);
    free(disk);
    return 0;
}





/******************************************** COMPUTE THE SPECTRUM **************************************************************/

int disk_spectrum(double multiplicator, double *wavelength, double *spec, int len_spec){



    double *r;                        /** radius which is from the center of disk to the center of any region**/
    double *theta;                    /** azimuth angle which is from the origine to the r for any region**/
    r = (double *) calloc(Nr,sizeof(double));
    theta = (double *) calloc(Ntheta,sizeof(double));

    double M = 3.2e7*Msun;              /** M_sun, the black hole mass, converted to gr **/
    double Rg= (Ggrav*M)/(c*c);         /** gravitational radius **/
    double r_in= 6.0*Rg;                /** inner radius **/
    double r_out=10000*Rg;              /** outer radius **/


    double inc_angle = 45.0*degree_to_radian;   /** inclination angle , converted to radian **/
    double cos_inc_angle = cos(inc_angle);  /** Cos of the inclination angle, avoid to recompute it all the time */
    double h_star = 10.0*Rg;                /** the vertical distance from the cetral variable source to disk **/
    double M_rate = 1.0*Msun/31557600.0;    /** M_sun yr^-1, the typical black hole accretion rate , converted to g.s^-1 **/
    /** the numerical factor converts from year to second: we are working in cgs: cm gram second.*/
    double A = 0.5;                         /** the disk albedo **/
    double L_bol = multiplicator*0.0;                 /** erg/s^-1, the bolometric luminosity **/

    /** Checking the values of the radii */
    // printf("Rg = %g\tR_int = %g\tR_out = %g \n", Rg, r_in, r_out);
    // getchar();

    /** create disks which contain the regions **/
    region *disk;
    disk = (region *) malloc(Nr*Ntheta*sizeof(region));

    /** the ratio of the outher and inner radius of each rings fixed **/
    double step = exp(log(r_out/r_in)/Nr);
    int i;
    for (i=0; i < Nr; i++){
        r[i] = r_in*pow(step,(double) i);
    }

    for (i=0; i < Ntheta; i++){
        theta[i] = ((double) i)*(360.0/Ntheta)*degree_to_radian;
    }

    /** fill the disks with elements (radius and theta) of regions **/
    int j;
    double tau;
    for (i=0; i < Nr; i++){
        for (j=0; j < Ntheta; j++){
            disk[i*Ntheta+j].radius = r[i];                 /** disk[0] region1, ... **/
            disk[i*Ntheta+j].theta = theta[j];
            disk[i*Ntheta+j].rstar = r_star(r[i], h_star);   /** disk[0] region1, ... **/

            /// Compute the time lag up to the radius.
            tau = sqrt(pow(h_star,2.0)+pow(r[i],2.0))+h_star*cos_inc_angle-r[i]*cos(theta[j])*sin(inc_angle);
            tau = (tau/c)*second_to_day;
            disk[i*Ntheta+j].tau = tau;
        }
    }


    /**  Husne,  9/10/2018
     *
     */
    /** for the computation of luminosity so it is for temperature **/
    ///double omega = 10.0*c/r_out;
    /** for the computation of the radiation from the disk. **/
    double D = 75.01*1e6*pc;                   /** Mpc distance from observer to the source, converted to cm **/
    double D2 = D*D;
    double R_in;
    double R_out;
    double theta_in;
    double theta_out;


    /**  Husne, 15/11/2018
      *  Compute all the other light curve in different filters (UVM2, UVW1, U, B, V)
      *  by using the illuminating light curves of UVW2
      */

    /**  Husne,  18/10/2018 *  Now compute the integral for U band. */
    double Temperature_t;
    double Integral;
    double f_U_i = 0.0;
    ///int Ntime = 0;




    int k;
    double stepR = exp(log(r_out/r_in)/Nr);
    double sqrt_stepR = sqrt(stepR);
    double stepT = 360.0*degree_to_radian/Ntheta;

    ///int nb_computation;
    /** Loop for the Number of filter, because I need to compute the flux for each band */

    for (j=0; j < Nr*Ntheta; j++){
        Temperature_t = Temp_profile_steady( disk[j].radius, disk[j].rstar,  M, M_rate, r_in);
        disk[j].temp = Temperature_t;
        /*if(j % 10 == 0){
            printf("j = %d\t Temp = %g\n", j , Temperature_t);
        }*/
    }

//    printf("Len_spec = %d\n", len_spec);
    for(k = 0 ; k < len_spec ; k ++){
        Integral = 0.0;
        for (j=0; j < Nr*Ntheta; j++){

     //       Temperature_t = Temp_profile_steady( disk[j].radius, disk[j].rstar,  M, M_rate, r_in);

            R_in = disk[j].radius/sqrt_stepR;            /** from the center to the first layer of any region **/
            R_out = disk[j].radius*sqrt_stepR;           /** from the center to the last layer of any region **/
            theta_in = disk[j].theta - 0.5*stepT;         /** from the origine to the first layer of any region on the bottom**/
            theta_out = disk[j].theta + 0.5*stepT;        /** from the origine to the last layer of any region on the top**/
            f_U_i = spectrum(cos_inc_angle, D2, theta_in, theta_out, R_in, R_out, wavelength[k]*wavelength[k]* wavelength[k], wavelength[k], disk[j].temp);
            Integral += f_U_i;
 /*           if(k == 10){
                if(j % 10 == 0){
                    printf("j = %d\t integral = %g\n", j , f_U_i);
                }
            }*/
        }
        spec[k] = Integral;

    }


    free(r);
    free(theta);
    free(disk);

    return 0;
}




//
//int main(){
//
//
//    /** Define argument for the filters */
//    int Nfilter = 6;
//    int computed_filter[6] = {1, 1, 1, 1, 0, 0};
//
//
//
//    make_computation(Nfilter, computed_filter);
//
//    return 0;
//
//}
//
//


