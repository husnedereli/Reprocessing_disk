#include <stdio.h>
#include <stdlib.h>


#define Nr 66
#define Ntheta 66

/** Definition of the physical constants **/
#define pc 3.08567758e13; /** km **/
#define G 4.302e-3; /** pc M_sun^-1 (km/s)^-1 **/

int main()
{
    double *r; /** radius **/
    double *theta; /** azimuth angle **/
    r = (double *) calloc(Nr,sizeof(double));
    theta = (double *) calloc(Ntheta,sizeof(double));
    
    double G = 4.302;
    double M = 20;
    double c = 10;
    double Rg= (G*M)/(c*c);
    double r_min= 6*Rg;
    double r_max=1000;
    
    
    free(r);
    free(theta);
    return 0;
}