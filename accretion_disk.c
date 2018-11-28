#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "accretion_disk.h"

/**  Code specific declaration */
#define Nr 66 			/** layers **/
#define Ntheta 66		/** zones per layers **/

/** Definition of the physical constants **/
/**  Constant CGS */
#define c 2.99792458e10                 /** cm.s^-1         */
#define Msun 1.989e33                   /** g               */
#define h 6.6260755e-27                 /** erg.s           */
#define hc 1.9864474610385790e-16
#define kb 1.380657e-16                 /** erg.K^-1        */
#define sigmaSB 5.67051e-5
#define Ggrav 6.67259e-8                /** cm^3.g^-1.s^-2  */
#define pc 3.08568e18                   /** cm              */
#define nm 1e-7                         /** cm              */
#define angstrom 1e-8                   /** cm              */


/**  Other constant */
#define pi 3.14159265358979


/**  Husne, 09/10/2018
  *  Here I define all the function to create and fill the disk
  */



///** define the Function of the luminosity **/
double L_star(double L_bol, double time, double flux){///double omega, double t){
    //omega = 3*c/R_out
    return 0.8*L_bol*time*flux;///(1.0+sin(omega*t)); ///** 0.15*L_bol **/
}

/** define the Function of the distance from the central variable source to disk elements **/
double r_star(double r, double h_star){
    return sqrt(pow(h_star,2.0)+pow(r,2.0));
}


/** define the Function of the temperature profile **/
double temp_profile(double time, double r, double rstar, double tau, double theta, double M, double M_rate, double r_in, double A, double h_star, double inc_angle, double L_bol, double flux){

    /// Compute the time lag up to the radius. For speed purposed, it is now computed only one time in the main code.
    // double tau = sqrt(pow(h_star,2.0)+pow(r,2.0))+h_star*cos(inc_angle)-r*cos(theta*0.0174532925)*sin(inc_angle);
    // tau = tau/c;
    double Lstar = L_star(L_bol, flux, time - tau);
    //double rstar = r_star(r, h_star);
    //printf("Contrib 1 = %g\t contrib 2 = %g\n ",((3.0*Ggrav*M*M_rate)/(8.0*pi*sigmaSB*pow(r,3.0)))*(1.0-sqrt(r_in/r)),  ((1.0-A)*(h_star*Lstar/(4.0*pi*sigmaSB*pow(rstar,3.0)))));
    //getchar();
    //printf()
    return pow(((3.0*Ggrav*M*M_rate)/(8.0*pi*sigmaSB*pow(r,3.0)))*(1.0-sqrt(r_in/r)) +((1.0-A)*(h_star*Lstar/(4.0*pi*sigmaSB*pow(rstar,3.0)))) ,0.25);
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
} region;


int make_computation(int Nfilter, long int *computed_filter, double *time, double *flux){


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
    double L_bol = 2.82e46; /** the bolometric luminosity **/

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
        theta[i] = i*(360.0/Ntheta);
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
            tau = sqrt(pow(h_star,2.0)+pow(r[i],2.0))+h_star*cos_inc_angle-r[i]*cos(theta[j]*0.0174532925)*sin(inc_angle);
            tau = tau/c;
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
                    break;
            }
            fclose(input_filtername);
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
                    break;
            }
            fclose(input_filtername);
        }
    }


    for (j=0;j<Nfilter;j++){
        printf("j = %d\t  %d\n",j, numberofloop[j]);
    }






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
    const int Ntau = 7;
    double tau_time[7] = {3.0, 6.0, 10.0, 20.0, 40.0, 100.0, 200.0};

    //printf("%g\t\n",tau_time);  //* print the arrays */


    int Ntime = 0;            /** 5000days*86400 = seconds **/
    ///double *time;
    ///time = (double *) calloc(Ntime,sizeof(double));
    ///for (i=0; i<Ntime; i++){
   ///     time[i] = i/2.0;
   ///     //printf("time=%g\n", t[i]);
   /// }

    // a) make the sum of all the S(tau, t)*1/Ntime
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
    double stepT = 360.0/Ntheta;
    /** Loop for the Number of filter, because I need to compute the average S for all tau by using two bands */
    for (m=0;m<Nfilter;m++){
        for (n=m+1;n<Nfilter;n++){
            //*it is important, when the filter number given as a "0" make computation for all casses between two filters */
            if(computed_filter[m] == 0 && computed_filter[n] == 0){
                printf("\n");
                /** Loop for the tau, because I need to compute the average S for all tau */
                for (l=0; l < Ntau; l++){
                    S_BU = 0.0;
                    /** Loop for the time, because I need to compute an average with respect to time */
                    for (k=0; k < Ntime; k++){
                        if(k % 10 == 0){
                            printf("k = %d\n", k);
                        }
                        /**  I need to compute S for all t
                         *  I need to compute f(t+tau), f(t) for U band
                         */

                        flux_t_U = 0.0;
                        flux_t_B = 0.0;
                        flux_tptau_U = 0.0;
                        flux_tptau_B = 0.0;

                        /** Loop for the radius and theta, because I need to compute the temparature and spectrum of disk */
                        /// f is the summ of contribution from all the disk elements.

                        for (j=0; j < Nr*Ntheta; j++){

                            R_in = disk[j].radius/sqrt_stepR;            /** from the center to the first layer of any region **/
                            R_out = disk[j].radius*sqrt_stepR;           /** from the center to the last layer of any region **/
                            theta_in = disk[j].theta - 0.5*stepT;         /** from the origine to the first layer of any region on the bottom**/
                            theta_out = disk[j].theta + 0.5*stepT;        /** from the origine to the last layer of any region on the top**/

                            /**  Now I compute the integral for the U-band */
                            /// temperature at time t in U
                            Temperature_t = temp_profile(time[k], disk[j].radius,disk[j].rstar, disk[j].tau, disk[j].theta, M, M_rate, r_in, A, h_star, inc_angle, L_bol, flux[k]);
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


                            /// temperature at time t + tau in U
                            Temperature_tptau = temp_profile(time[k] + tau_time[l], disk[j].radius, disk[j].rstar, disk[j].tau, disk[j].theta, M, M_rate, r_in, A, h_star, inc_angle, L_bol, flux[k]);
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





                            /**  Now I compute the integral for the B-band */
                            /// temperature at time t in B
                            Temperature_t = temp_profile(time[k], disk[j].radius, disk[j].rstar, disk[j].tau, disk[j].theta, M, M_rate, r_in, A, h_star, inc_angle, L_bol, flux[k]);
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


                            /// temperature at time t + tau in U
                            Temperature_tptau = temp_profile(time[k] + tau_time[l], disk[j].radius, disk[j].rstar, disk[j].tau, disk[j].theta, M, M_rate, r_in, A, h_star, inc_angle, L_bol, flux[k]);
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
                        }

                        S_BU = S_BU + (flux_tptau_B-flux_t_B) / (flux_tptau_U-flux_t_U);
                        printf("%.13g\t\n", flux_t_B);
                        printf("%.13g\t\n", flux_t_U);
                        

                    }
                    avarage_SBU = S_BU/Ntime;
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


