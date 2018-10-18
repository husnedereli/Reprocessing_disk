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

/**  Other constant */
#define pi 3.14159265358979


/**  Husne, 09/10/2018
  *  Here I define all the function to create and fill the disk
  */

/** define the Function of the lag **/
double lag_tao(double r, double theta, double inc_angle, double h_star){
    return sqrt(pow(h_star,2.0)+pow(r,2.0))+h_star*cos(inc_angle)-r*cos(theta)*sin(inc_angle);
}

/** define the Function of the luminosity **/
double L_star(double t, double r, double theta, double inc_angle, double h_star){
    //return t-lag_tao(r, theta, inc_angle, h_star)/c;
    return 0;
}


/** define the Function of the distance from the central variable source to disk elements **/
double r_star(double r, double h_star){
    return sqrt(pow(h_star,2.0)+pow(r,2.0));
}

/** define the Function of the temperature profile **/
double temp_profile(double t, double r, double theta, double M, double M_rate, double r_in, double A, double h_star, double inc_angle){

    double Lstar = L_star(t, r, theta, inc_angle, h_star);
    double rstar = r_star(r, h_star);
    return pow(((3.0*Ggrav*M*M_rate)/(8.0*pi*sigmaSB*pow(r,3.0)))*(1.0-sqrt(r_in/r)) +((1.0-A)*(h_star*Lstar/(4.0*pi*sigmaSB*pow(rstar,3.0)))),0.25);
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


/**  Husne, 18/10/2018
 *  Here I define all the function required for the convolotion with the filter bandpass.
 */

/** define the response Function **/
/*double response(double R, double lambda_max, double lambda_min){
    return (pow(R(lambda_max),2)/2)-(pow(R(lambda_min),2)/2);
}*/

/** define the convolation Function **/
/*double convolation(double inc_angle, double D, double theta_in, double theta_out, double R_in, double R_out, double lambda, double temperature, double R, double lambda_max, double lambda_min){
    double spectra = spectrum(inc_angle, D, theta_in, theta_out, R_in, R_out, lambda, temperature);
    return spectra*response(R, lambda_max, lambda_min);
}*/




/** define new type as regions and its elements to create disk **/
typedef struct region {
    double radius;
    double theta;
    double temp;
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
    double M_rate = 1.0*Msun;               /** M_sun yr^-1, the typical black hole accretion rate , converted to gr **/
    double A = 0.5;                         /** the disk albedo **/
    //double L_bol = 2.82e44; /** the bolometric luminosity **/
    //double L_star = 0.5*L_bol; /** the luminosity of central variable source **/


    double t;
    double lag;
    double Lstar;
    double rstar;
    double temperature;
    /** call the functions **/
    for (i=0; i < Nr; i++){
        for (j=0; j < Ntheta; j++){
            t = 10.0;
            temperature = temp_profile (t, r[i], theta[j], M, M_rate, r_in, A, h_star, inc_angle);
            // printf("Temperature[%d]: %g\n",i, temperature);
            /** fill the disks with elements (temp) of regions **/
            disk[i*Ntheta+j].temp = temperature;
        }
    }

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
      *  Read a txt file for U bandpass.
      */

    FILE *input;
    double c1, c2;
    int numberofloop = 0;

    input=fopen("Bessel_U-1.txt","r");      //* open a text file for reading */

    /**  Here %lf means type double */
    /// step 1
    while(fscanf(input,"%lf%lf", &c1, &c2) !=EOF ){
        numberofloop++;                 //* caunt the number of loop */
        /// Previous line is equivalent to      numberofloop = numberofloop + 1;                 //* caunt the number of loop */
    }
    fclose(input);

    /// step 2
    double *wavelength;                 //* create an array */
    wavelength = (double *) calloc(numberofloop,sizeof(double));
    double *transmission;
    transmission = (double *) calloc(numberofloop,sizeof(double));


    /// step 3
    input=fopen("Bessel_U-1.txt","r"); //* open a text file for reading */
    i = 0;
    while(fscanf(input,"%lf%lf", &c1, &c2) !=EOF ){

        wavelength[i] = c1;             //* fill the array */
        transmission[i] = c2;
        i += 1 ;
        /// i = i + 1 ;
    }
    fclose(input);

    for(i = 0; i < numberofloop ; i++){
        // printf("%g\t%g\n",wavelength[i], transmission[i]);  //* print the arrays */
    }


    /**  Husne,  18/10/2018
     *  Now compute the integral.
     */
    double compute_integral = 0.0;
    //double N = 100;
    double deltaLambda;
    double summ_region_with_i;
    double summ_region_with_im1;
    for(i = 1; i < numberofloop ; i++){

        deltaLambda = (wavelength[i]-wavelength[i-1]);
        summ_region_with_i = 0.0;
        summ_region_with_im1 = 0.0;

        for (j=0; j < Nr*Ntheta; j++){
            R_in = disk[j].radius/sqrt(step);        /** from the center to the first layer of any region **/
            R_out = disk[j].radius*sqrt(step);       /** from the center to the last layer of any region **/
            theta_in = disk[j].theta-(step/2.0);    /** from the origine to the first layer of any region on the bottom**/
            theta_out = disk[j].theta+(step/2.0);   /** from the origine to the last layer of any region on the top**/

            summ_region_with_i += spectrum(inc_angle, D, theta_in, theta_out, R_in, R_out, wavelength[i], disk[j].temp)
            summ_region_with_im1 += spectrum(inc_angle, D, theta_in, theta_out, R_in, R_out, wavelength[i-1], disk[j].temp)

        }
        compute_integral = compute_integral + deltaLambda*0.5*(transmission[i-1]*summ_region_with_im1 + transmission[i]*summ_region_with_i ) ;

    }
    printf("%g\t\n",compute_integral);  //* print the arrays */


    free(wavelength);
    free(transmission);

    free(r);
    free(theta);
    free(disk);

    return 0;
}







