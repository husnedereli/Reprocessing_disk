#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define Nr 66 /** layers **/
#define Ntheta 66 /** zones per layers **/

/** Definition of the physical constants **/
#define pc 3.08567758e13 /** km **/
#define G 4.302e-3 /** pc M_sun^-1 (km/s)^-1 **/
#define c 3e5 /** km/s **/
#define sigma = 5.670373e-8 /** W m^-2 K^4 the Stefan–Boltzmann constant  **/
#define pi 3.14


/** Function of the lag **/
double inc_angle = 45; /** inclination angle **/
int lag_tao (double r, double theta, int inc_angle){
    sqrt(h_star**2+r**2)+h_star*cos(inc_angle)-r*cos(theta)*sin(inc_angle)
}

/** Function of the distance from the cetral variable source to disk elements **/
double r_star (double r){
    sqrt(h_star**2+r**2)
}


/** Function of the temperature profile **/
double temp_profile (int t, double r, double theta ){
    return ((3.0*G*M*M_rate)/(8.0*pi*sigma*pow(r,3.0)))*(1.0-sqrt(r_in/r))+(1.0-A)*((h_star*L_star(t-lag_tao(int inc_angle)))/(4.0*pi*sigma*pow(r_star,3.0)));
}
    

    





int main()
{
    double *r; 				/** radius **/
    double *theta; 			/** azimuth angle **/
    r = (double *) calloc(Nr,sizeof(double));
    theta = (double *) calloc(Ntheta,sizeof(double));
    
    
    double M = 3.2e7; /** M_sun, the black hole mass **/
    double Rg= (G*pc*M)/(c*c); /** gravitational radius **/
    double r_in= 6 * Rg; /** inner radius **/
    double r_out=3200*Rg; /** outer radius **/
    
    /** the ratio of the outher and inner radius of each rings fixed **/
    double step = exp(log(r_out/r_in)/Nr);
    int i;
    for (i=0; i < Nr; i++){
        r[i] = r_in*pow(step,i);
    }
    
    for (i=0; i < Ntheta; i++){
        theta[i] = i*(360/Ntheta);
    }

    
    free(r);
    free(theta);
    /** the temperature profile **/
    double M_rate = 100; /** ??? the accretion rate **/
    double A = 0.5; /** the disk albedo **/
    double h_star = 10*Rg; /** the vertical distance from the cetral variable source to disk **/
    double L_bol = 2.82e44; /** the bolometric luminosity **/
    double L_star = 0.5*L_bol; /** the luminosity of central variable source **/
    
    
    
    
    
    
    
    return 0;
}







