#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/**  Code specific declaration */
#define Nr 66 			/** layers **/
#define Ntheta 66		/** zones per layers **/

/** Definition of the physical constants **/
/**  Constant CGS */
#define c 2.99792458e10                 /** cm.s^-1 */
#define Msun 1.989e33                   /** gr       */
#define h 6.6260755e-27                 /** erg.s   */
#define kb 1.380657e-16                 /** erg.K^-1  */
#define sigmaSB 5.67051e-5
#define Ggrav 6.67259e-8                /** cm.s^-2 */
#define pc 3.08568e18                   /** cm      */

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
    return ((3.0*Ggrav*M*M_rate)/(8.0*pi*sigmaSB*pow(r,3.0)))*(1.0-sqrt(r_in/r))+((1.0-A)*(h_star*L_star(t, r, theta, inc_angle, h_star)/(4.0*pi*sigmaSB*pow(r_star(r, h_star),3.0))));
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
} region;



int main(){


    /**  Husne, 9/10/2018
      *  Here I create and fill the disk. I compute the temperature and settle all the regions of the disk.
      */


    double *r;                      /** radius which is from the center of disk to the center of any region**/
    double *theta;                  /** azimuth angle which is from the origine to the r for any region**/
    r = (double *) calloc(Nr,sizeof(double));
    theta = (double *) calloc(Ntheta,sizeof(double));

    double M = 3.2e7*Msun;               /** M_sun, the black hole mass, converted to gr **/
    double Rg= (Ggrav*M)/(c*c);         /** gravitational radius **/
    double r_in= 6.0*Rg;                /** inner radius **/
    double r_out=10000*Rg;              /** outer radius **/

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
            disk[i*Ntheta+j].theta=theta[j];
        }
    }



    double inc_angle = 45.0*0.0174532925; /** inclination angle , converted to radian **/
    double h_star = 10.0*Rg; /** the vertical distance from the cetral variable source to disk **/
    double M_rate = 1.0*Msun; /** M_sun yr^-1, the typical black hole accretion rate , converted to gr **/
    double A = 0.5; /** the disk albedo **/
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
            printf("Temperature[%d]: %g\n",i, temperature);
            /** fill the disks with elements (temp) of regions **/
            disk[i*Ntheta+j].temp = temperature;
        }
    }

    /**  Husne,  9/10/2018
      *  Now I compute the radiation from such disk.
      */

    double D;
    double R_in;
    double R_out;
    double theta_in;
    double theta_out;
    double SED;
    double lambda;
    /** call the functions **/
    for (j=0; j < Nr*Ntheta; j++){
        D = 75.01*1e6*pc;                   /** Mpc distance from observer to the source, converted to cm **/
        lambda = h*c/(kb*disk[j].temp);
        R_in = disk[j].radius/sqrt(step);        /** from the center to the first layer of any region **/
        R_out = disk[j].radius*sqrt(step);       /** from the center to the last layer of any region **/
        theta_in = disk[j].theta-(step/2.0);    /** from the origine to the first layer of any region on the bottom**/
        theta_out = disk[j].theta+(step/2.0);   /** from the origine to the last layer of any region on the top**/
        SED = spectrum(inc_angle, D, theta_in, theta_out, R_in, R_out, lambda, disk[j].temp);
        //printf("SED[%d]: %g\n",j, SED);
    }
    
    

    
    
    
    free(r);
    free(theta);
    free(disk);

    return 0;
}







