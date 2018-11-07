#include <stdio.h>
#include <stdlib.h>
#include <math.h>



int main(){
/**  Husne,  11/10/2018
 *Read a txt file for U bandpass.
 */
    int i;
    FILE *input_U;
    double c1_U, c2_U;
    int numberofloop_U = 0;

    input_U=fopen("U.txt","r");      //* open a text file for reading */

    /**  Here %lf means type double */
    /// step 1
    while(fscanf(input_U,"%lf%lf", &c1_U, &c2_U) !=EOF ){
        if(c2_U > 0.01){
            numberofloop_U++;                 //* caunt the number of loop */
        /// Previous line is equivalent to      numberofloop_U = numberofloop_U + 1;                 //* caunt the number of loop */
        }
    }
    fclose(input_U);

    /// step 2
    double *wavelength_U;                 //* create an array */
    wavelength_U = (double *) calloc(numberofloop_U,sizeof(double));
    double *transmission_U;
    transmission_U = (double *) calloc(numberofloop_U,sizeof(double));


    /// step 3
    input_U=fopen("U.txt","r"); //* open a text file for reading */
    i = 0;
    while(fscanf(input_U,"%lf%lf", &c1_U, &c2_U) !=EOF ){
        if(c2_U > 0.01){
            wavelength_U[i] = c1_U;             //* fill the array */
            transmission_U[i] = c2_U;
            i += 1 ;
        /// i = i + 1 ;
        }
    }
    fclose(input_U);

    for(i = 0; i < numberofloop_U ; i++){
        //printf("%g\t%g\n",wavelength_U[i], transmission_U[i]);  //* print the arrays */
    }


    int avarage_wavelength_U = 0.0;
    int j = 0.0;
    while (5+5*j< numberofloop_U) {
        int sum = 0.0;
        for(i = 0; i < 5; i++){
            sum = sum + wavelength_U[i+5*j];
        }
        avarage_wavelength_U = sum/5;
        printf("Average is: %d",avarage_wavelength_U);
    }
    
    


    return 0;
}
