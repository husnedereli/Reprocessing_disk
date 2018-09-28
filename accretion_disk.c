#include <stdio.h>
#include <stdlib.h>


#define Nr 66
#define Ntheta 66

/** Definition of the physical constants **/
#define pc 3.08567758e13; /** km **/
#define G 4.302e-3; /** pc M_sun^-1 (km/s)^-1 **/
#define c 3e5 /** km/s **/


int main()
{
    double *r; /** radius **/
    double *theta; /** azimuth angle **/
    r = (double *) calloc(Nr,sizeof(double));
    theta = (double *) calloc(Ntheta,sizeof(double));
    
    
    double M = 3.2e7; /** M_sun, the black hole mass **/
    double Rg= (G*pc*M)/(c*c);
    double r_min= 6 * Rg;
    double r_max=1000;
    
    
    double step = exp(log(r_max/r_min)/Nr);
    int i;
    for (i=0; i < Nr; i++){
        r[i] = r_min*pw(step,i);
    }
    
    int i;
    for (i=0; i < Ntheta; i++){
        theta[i] = i*(360/Ntheta);
    }

    
    free(r);
    free(theta);
    return 0;
}