#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/**  Code specific declaration */
#define Nr 66 			/** layers **/
#define Ntheta 66		/** zones per layers **/

/** Definition of the physical constants **/
/**  Constant CGS */
#define c 2.99792458e10                 /** cm.s^-1         */
#define Msun 1.989e33                   /** g               */
#define h 6.6260755e-27                 /** erg.s           */
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

/** define the Function of the lag **/
double lag_tao(double r, double theta, double inc_angle, double h_star){
    return sqrt(pow(h_star,2.0)+pow(r,2.0))+h_star*cos(inc_angle)-r*cos(theta)*sin(inc_angle);
}

///** define the Function of the luminosity **/
//double L_star(double t, double r, double theta, double inc_angle, double h_star){
//    //return t-lag_tao(r, theta, inc_angle, h_star)/c;
//    return 0;
//}

double L_star(double L_bol, double omega, double t){
    //omega = 3*c/R_out
    return 0.15*L_bol*(1+sin(omega*t));
}



/** define the Function of the distance from the central variable source to disk elements **/
double r_star(double r, double h_star){
    return sqrt(pow(h_star,2.0)+pow(r,2.0));
}

/** define the Function of the temperature profile **/
double temp_profile(double t, double r, double theta, double M, double M_rate, double r_in, double A, double h_star, double inc_angle, double L_bol, double omega){

    //double Lstar = L_star(t, r, theta, inc_angle, h_star);
    /// Compute the time lag up to the radius.
    double tau = sqrt(pow(h_star,2.0)+pow(r,2.0))+h_star*cos(inc_angle)-r*cos(theta*0.0174532925)*sin(inc_angle);
    tau = tau/c;
    double Lstar = L_star(L_bol, omega, t - tau);
    double rstar = r_star(r, h_star);
    //printf("Contrib 1 = %g\t contrib 2 = %g\n ",((3.0*Ggrav*M*M_rate)/(8.0*pi*sigmaSB*pow(r,3.0)))*(1.0-sqrt(r_in/r)),  ((1.0-A)*(h_star*Lstar/(4.0*pi*sigmaSB*pow(rstar,3.0)))));
    //getchar();
    //printf()
    return pow(((3.0*Ggrav*M*M_rate)/(8.0*pi*sigmaSB*pow(r,3.0)))*(1.0-sqrt(r_in/r)) +((1.0-A)*(h_star*Lstar/(4.0*pi*sigmaSB*pow(rstar,3.0)))) ,0.25);
}




/**  Husne, 09/10/2018
  *  Here I define all the function required to compute the spectra
  */

/** define the Planck Function **/
double Planck_Function(double lambda, double temperature){
    return (2.0*h*c)/pow(lambda,3.0)/(exp(h*c/(lambda*kb*temperature))-1.0);
}

/** define the Function of predicted spectrum (SED) **/
double spectrum(double inc_angle, double D, double theta_in, double theta_out, double R_in, double R_out, double lambda, double temperature){
    return (cos(inc_angle)/pow(D,2.0))*Planck_Function(lambda, temperature)*(theta_out-theta_in)*(pow(R_out,2.0)/2.0-pow(R_in,2.0)/2.0);
}



/** define new type as regions and its elements to create disk **/
typedef struct region {
    double radius;
    double theta;
    double temp;
    double temp_tplustau;
} region;



int main(){


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
    for (i=0; i < Nr; i++){
        for (j=0; j < Ntheta; j++){
            disk[i*Ntheta+j].radius = r[i];  /** disk[0] region1, ... **/
            disk[i*Ntheta+j].theta = theta[j];
        }
    }



    double inc_angle = 45.0*0.0174532925;   /** inclination angle , converted to radian **/
    double h_star = 10.0*Rg;                /** the vertical distance from the cetral variable source to disk **/
    double M_rate = 1.0*Msun/31557600.0;    /** M_sun yr^-1, the typical black hole accretion rate , converted to gr **/
                                            /** the numerical factor converts from year to second: we are working in cgs: cm gram second.*/
    double A = 0.5;                         /** the disk albedo **/
    double L_bol = 2.82e46; /** the bolometric luminosity **/
    //double L_star = 0.5*L_bol; /** the luminosity of central variable source **/

    //printf("Rg = %g\t r_star = %g\t M_rate = %g\t")

   
    /**  Husne, 20/10/2018
     * I define the time and tau_time as arrays
     */
    
    //double t;
    int Ntime = 5000;            /** 5000days*86400 = seconds **/
    double *t;
    t = (double *) calloc(Ntime,sizeof(double));
    for (i=0; i<Ntime; i++){
        t[i] = i/2.0;
        //printf("t=%g\n", t[i]);
    }

    //getchar();
    
    double Ntau = 94.8;            /** days*86400 = seconds deconvolotion timescale **/
    double *tau_time;
    t = (double *) calloc(Ntau,sizeof(double));
    for (i=0; i<Ntau; i++){
        tau_time[i] = i/2.0;
        //printf("tau_time=%g\n", t[i]);
    }
    
    
    
    
    

    double lag;
    double Lstar;
    double rstar;
    double temperature;
    double omega = 10.0*c/r_out;
    /** call the functions **/
    //FILE *test;
    //test = fopen("temperature.txt","a");
    int k;
    for (i=0; i < Nr; i++){
        for (j=0; j < Ntheta; j++){
            //t = 10.0;
            for (k=0; k < Ntime; k++){
                temperature = temp_profile (t[k], r[i], theta[j], M, M_rate, r_in, A, h_star, inc_angle, L_bol, omega);
                //printf("Temperature[%d]: %g\n",i, temperature);
                //fprintf(test, "%g\t%g\t%g\n", r[i], theta[j], temperature);
                /** fill the disks with elements (temp) of regions **/
                //disk[i*Ntheta+j].temp = temperature;
                disk[i*Ntheta+j].temp = temperature/k;
            }
        }
        //fprintf(test, "\n");
    }
    //fclose(test);

    /**  Husne, 22/10/2018
     * Compute the temperature for t+tau
     */
    double *Time;
    for (i=0; i<Ntime; i++){
        for (i=0; i<Ntau; i++){
            Time[i] = t[i]+tau_time[j];
        }
    }


    double temperature_tplustau;
    /** call the functions **/
    //FILE *test;
    //test = fopen("temperature.txt","a");
    for (i=0; i < Nr; i++){
        for (j=0; j < Ntheta; j++){
            //t = 10.0;
            for (k=0; k < Ntime; k++){
                temperature_tplustau = temp_profile (Time[k], r[i], theta[j], M, M_rate, r_in, A, h_star, inc_angle, L_bol, omega);
                //printf("Temperature[%d]: %g\n",i, temperature);
                //fprintf(test, "%g\t%g\t%g\n", r[i], theta[j], temperature);
                /** fill the disks with elements (temp) of regions **/
                //disk[i*Ntheta+j].temp_tplustau = temperature_tplustau;
                disk[i*Ntheta+j].temp_tplustau = temperature_tplustau/k;
            }
        }
        //fprintf(test, "\n");
    }
    //fclose(test);
    
    
    
    

    /**  Husne,  9/10/2018
      *  Now I compute the radiation from such disk.
      */


    double D = 75.01*1e6*pc;                   /** Mpc distance from observer to the source, converted to cm **/
    double R_in;
    double R_out;
    double theta_in;
    double theta_out;
    double SED;
    double lambda;
    /** call the functions **/
    for (j=0; j < Nr*Ntheta; j++){


        R_in = disk[j].radius/sqrt(step);        /** from the center to the first layer of any region **/
        R_out = disk[j].radius*sqrt(step);       /** from the center to the last layer of any region **/
        theta_in = disk[j].theta-(step/2.0);    /** from the origine to the first layer of any region on the bottom**/
        theta_out = disk[j].theta+(step/2.0);   /** from the origine to the last layer of any region on the top**/
        SED = spectrum(inc_angle, D, theta_in, theta_out, R_in, R_out, lambda, disk[j].temp);
        //printf("SED[%d]: %g\n",j, SED);
    }






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

    FILE *input_U;
    double c1_U, c2_U;
    int numberofloop_U = 0;

    input_U=fopen("bess_U.txt","r");      //* open a text file for reading */

    /**  Here %lf means type double */
    /// step 1
    while(fscanf(input_U,"%lf%lf", &c1_U, &c2_U) !=EOF ){
        numberofloop_U++;                 //* caunt the number of loop */
        /// Previous line is equivalent to      numberofloop_U = numberofloop_U + 1;                 //* caunt the number of loop */
    }
    fclose(input_U);

    /// step 2
    double *wavelength_U;                 //* create an array */
    wavelength_U = (double *) calloc(numberofloop_U,sizeof(double));
    double *transmission_U;
    transmission_U = (double *) calloc(numberofloop_U,sizeof(double));


    /// step 3
    input_U=fopen("bess_U.txt","r"); //* open a text file for reading */
    i = 0;
    while(fscanf(input_U,"%lf%lf", &c1_U, &c2_U) !=EOF ){

        wavelength_U[i] = c1_U;             //* fill the array */
        transmission_U[i] = c2_U;
        i += 1 ;
        /// i = i + 1 ;
    }
    fclose(input_U);

    for(i = 0; i < numberofloop_U ; i++){
         //printf("%g\t%g\n",wavelength_U[i], transmission_U[i]);  //* print the arrays */
    }


    /**  Husne,  18/10/2018
     *  Now compute the integral for U band.
     */
    double compute_integral_U = 0.0;
    double compute_integral_U_tplustau = 0.0;
    double deltaLambda_U;
    double summ_region_with_i_U;
    double summ_region_with_im1_U;
    double summ_region_with_i_U_tplustau;
    double summ_region_with_im1_U_tplustau;
    for(i = 1; i < numberofloop_U ; i++){

        deltaLambda_U = (wavelength_U[i]*angstrom-wavelength_U[i-1]*angstrom);
        summ_region_with_i_U = 0.0;
        summ_region_with_im1_U = 0.0;

        for (j=0; j < Nr*Ntheta; j++){
            R_in = disk[j].radius/sqrt(step);        /** from the center to the first layer of any region **/
            R_out = disk[j].radius*sqrt(step);       /** from the center to the last layer of any region **/
            theta_in = disk[j].theta-(step/2.0);    /** from the origine to the first layer of any region on the bottom**/
            theta_out = disk[j].theta+(step/2.0);   /** from the origine to the last layer of any region on the top**/

            summ_region_with_i_U += spectrum(inc_angle, D, theta_in, theta_out, R_in, R_out, wavelength_U[i]*angstrom, disk[j].temp);
            summ_region_with_im1_U += spectrum(inc_angle, D, theta_in, theta_out, R_in, R_out, wavelength_U[i-1]*angstrom, disk[j].temp);
            summ_region_with_i_U_tplustau += spectrum(inc_angle, D, theta_in, theta_out, R_in, R_out, wavelength_U[i]*angstrom, disk[j].temp_tplustau);
            summ_region_with_im1_U_tplustau += spectrum(inc_angle, D, theta_in, theta_out, R_in, R_out, wavelength_U[i-1]*angstrom, disk[j].temp_tplustau);

        }
        compute_integral_U = compute_integral_U + deltaLambda_U*0.5*(transmission_U[i-1]*summ_region_with_im1_U + transmission_U[i]*summ_region_with_i_U ) ;
        compute_integral_U_tplustau = compute_integral_U_tplustau + deltaLambda_U*0.5*(transmission_U[i-1]*summ_region_with_im1_U_tplustau + transmission_U[i]*summ_region_with_i_U_tplustau) ;

    }
    //printf("%g\t\n",compute_integral_U);  //* print the arrays */



    /**  Husne,  19/10/2018
     *  Convolotion with the filter bandpass.
     Read a txt file for B bandpass.
     */

    FILE *input_B;
    double c1_B, c2_B;
    int numberofloop_B = 0;

    input_B=fopen("bess_B.txt","r");      //* open a text file for reading */

    /**  Here %lf means type double */
    /// step 1
    while(fscanf(input_B,"%lf%lf", &c1_B, &c2_B) !=EOF ){
        numberofloop_B++;                 //* caunt the number of loop */
        /// Previous line is equivalent to      numberofloop_B = numberofloop_B + 1;                 //* caunt the number of loop */
    }
    fclose(input_B);

    /// step 2
    double *wavelength_B;                 //* create an array */
    wavelength_B = (double *) calloc(numberofloop_B,sizeof(double));
    double *transmission_B;
    transmission_B = (double *) calloc(numberofloop_B,sizeof(double));


    /// step 3
    input_B=fopen("bess_B.txt","r"); //* open a text file for reading */
    i = 0;
    while(fscanf(input_B,"%lf%lf", &c1_B, &c2_B) !=EOF ){

        wavelength_B[i] = c1_B;             //* fill the array */
        transmission_B[i] = c2_B;
        i += 1 ;
        /// i = i + 1 ;
    }
    fclose(input_B);

    for(i = 0; i < numberofloop_B ; i++){
        //printf("%g\t%g\n",wavelength_B[i], transmission_B[i]);  //* print the arrays */
    }


    /**  Husne,  19/10/2018
     *  Now compute the integral for B band.
     */
    double compute_integral_B = 0.0;
    double compute_integral_B_tplustau = 0.0;
    double deltaLambda_B;
    double summ_region_with_i_B;
    double summ_region_with_im1_B;
    double summ_region_with_i_B_tplustau;
    double summ_region_with_im1_B_tplustau;
    for(i = 1; i < numberofloop_B ; i++){

        deltaLambda_B = (wavelength_B[i]*angstrom-wavelength_B[i-1]*angstrom);
        summ_region_with_i_B = 0.0;
        summ_region_with_im1_B = 0.0;

        for (j=0; j < Nr*Ntheta; j++){
            R_in = disk[j].radius/sqrt(step);        /** from the center to the first layer of any region **/
            R_out = disk[j].radius*sqrt(step);       /** from the center to the last layer of any region **/
            theta_in = disk[j].theta-(step/2.0);    /** from the origine to the first layer of any region on the bottom**/
            theta_out = disk[j].theta+(step/2.0);   /** from the origine to the last layer of any region on the top**/

            summ_region_with_i_B += spectrum(inc_angle, D, theta_in, theta_out, R_in, R_out, wavelength_B[i]*angstrom, disk[j].temp);
            summ_region_with_im1_B += spectrum(inc_angle, D, theta_in, theta_out, R_in, R_out, wavelength_B[i-1]*angstrom, disk[j].temp);
            summ_region_with_i_B_tplustau += spectrum(inc_angle, D, theta_in, theta_out, R_in, R_out, wavelength_B[i]*angstrom, disk[j].temp_tplustau);
            summ_region_with_im1_B_tplustau += spectrum(inc_angle, D, theta_in, theta_out, R_in, R_out, wavelength_B[i-1]*angstrom, disk[j].temp_tplustau);


        }
        compute_integral_B = compute_integral_B + deltaLambda_B*0.5*(transmission_B[i-1]*summ_region_with_im1_B + transmission_B[i]*summ_region_with_i_B);
        compute_integral_B_tplustau = compute_integral_B_tplustau + deltaLambda_B*0.5*(transmission_B[i-1]*summ_region_with_im1_B_tplustau + transmission_B[i]*summ_region_with_i_B_tplustau);

    }
    printf("%g\t\n",compute_integral_B);  //* print the arrays */



    /** ************************************************
     * ************************************************
     * ************************************************
     * ************************************************
     * ************************************************
     */


    /**  Husne,  19/10/2018
     *  compute the color variability and plot them
     */
    double *color_variation_BU;
    color_variation_BU = (double *) calloc(numberofloop_B,sizeof(double));
    for(i = 0; i < numberofloop_B ; i++){
        //for(j = 0; j < numberofloop_U ; j++){
        color_variation_BU[i] = (compute_integral_B_tplustau-compute_integral_B)/(compute_integral_U_tplustau-compute_integral_U);
        printf("%g\t\n",color_variation_BU[i]);  //* print the arrays */
       // }
    }


    free(r);
    free(theta);
    free(disk);
    free(wavelength_U);
    free(transmission_U);
    free(wavelength_B);
    free(transmission_B);


    return 0;
}







