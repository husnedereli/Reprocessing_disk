#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define Nr 66 /** layers **/
#define Ntheta 66 /** zones per layers **/

/** Definition of the physical constants **/
#define pc 3.08567758e13 /** km **/
#define G 4.302e-3 /** pc M_sun^-1 (km/s)^-1 **/
#define c 3e5 /** km/s **/
#define sigma1 5.670373e-8 /** W m^-2 K^4 the Stefan–Boltzmann constant  **/
#define pi 3.14





/** Function of the lag **/
double lag_tao (double r, double theta, double inc_angle, double h_star){
    return sqrt(pow(h_star,2.0)+pow(r,2.0))+h_star*cos(inc_angle)-r*cos(theta)*sin(inc_angle);
}

/** Function of the distance from the central variable source to disk elements **/
double r_star(double r, double h_star){
    return sqrt(pow(h_star,2.0)+pow(r,2.0));
}


/** Function of the temperature profile **/
double temp_profile (double t, double r, double theta, double M, double M_rate, double r_in, double A, double h_star, double L_star, double inc_angle, double lag_tao, double r_star){
    return ((3.0*G*M*M_rate)/(8.0*pi*sigma1*pow(r,3.0)))*(1.0-sqrt(r_in/r))+(1.0-A)*((h_star*L_star*(t-lag_tao))/(4.0*pi*sigma1*pow(r_star,3.0)));
}

/** define new type as region **/
typedef struct region {
    double r;
    double theta;
    double temp;
} region;



int main()
{
    double *r; 				/** radius **/
    double *theta; 			/** azimuth angle **/
    r = (double *) calloc(Nr,sizeof(double));
    theta = (double *) calloc(Ntheta,sizeof(double));

    double M = 3.2e7; /** M_sun, the black hole mass **/
    double Rg= (G*pc*M)/(c*c); /** gravitational radius **/
    double r_in= 6.0 * Rg; /** inner radius **/
    double r_out=3200.0*Rg; /** outer radius **/

    /** create disks **/
    region *disk;
    disk = (region *) malloc(Nr*Ntheta*sizeof(region));

    /** the ratio of the outher and inner radius of each rings fixed **/
    double step = exp(log(r_out/r_in)/Nr);
    int i;
    for (i=0; i < Nr; i++){
        r[i] = r_in*pow(step,i);
        //create_disk.r1 = r[i];
    }

    for (i=0; i < Ntheta; i++){
        theta[i] = i*(360.0/Ntheta);
        //create_disk.theta1 = theta[i];
    }
    
    /** fill the disks with regions **/
    int j;
    for (i=0; i < Nr; i++){
        r[i] = r_in*pow(step,i);
        for (j=0; j < Ntheta; j++){
            theta[j] = j*(360.0/Ntheta);
            disk[i j]= r[i]*theta[j];
        }
    }
    
    double inc_angle = 45.0; /** inclination angle **/
    double h_star = 10.0*Rg; /** the vertical distance from the cetral variable source to disk **/
    double M_rate = 100.0; /** ??? the accretion rate **/
    double A = 0.5; /** the disk albedo **/
    double L_bol = 2.82e44; /** the bolometric luminosity **/
    double L_star = 0.5*L_bol; /** the luminosity of central variable source **/


    
    
    for (i=0; i < Nr; i++){
        double t = 10.0;
        double lag = lag_tao (r[i], theta[i], inc_angle, h_star);
        double rstar = r_star(r[i], h_star);
        double temperature = temp_profile (t, r[i], theta[i], M, M_rate, r_in, A, h_star, L_star, inc_angle, lag, rstar);
        //printf ("Output: %f", temperature);
        //create_disk.temp1 = temperature;
    }


    
    free(r);
    free(theta);
    free(disk);

    return 0;
}







