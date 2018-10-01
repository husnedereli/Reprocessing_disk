#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define Nr 66 /** layers **/
#define Ntheta 66 /** zones per layers **/

/** Definition of the physical constants **/
#define pc 3.08567758e13 /** km **/
#define G 4.302e-3 /** pc M_sun^-1 (km/s)^-1 **/
#define c 3e5 /** km/s **/


int main()
{
    double *r; /** radius **/
    double *theta; /** azimuth angle **/
    r = (double *) calloc(Nr,sizeof(double));
    theta = (double *) calloc(Ntheta,sizeof(double));
    
    
    double M = 3.2e7; /** M_sun, the black hole mass **/
    double Rg= (G*pc*M)/(c*c); /** gravitational radius **/
    double r_in= 6 * Rg; /** inner radius **/
    double r_out=1000; /** outer radius **/
    
    
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
    return 0;
}