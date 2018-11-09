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

    input_U=fopen("Filter/B.txt","r");      //* open a text file for reading */

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
    input_U=fopen("Filter/B.txt","r"); //* open a text file for reading */
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

    //for(i = 0; i < numberofloop_U ; i++){
    //    printf("i = %d\t %g\t%g\n",i, wavelength_U[i], transmission_U[i]);  //* print the arrays */
    //}


    double *average_wavelength_U;
    average_wavelength_U = (double *) calloc(numberofloop_U,sizeof(double));
    double *average_transmission_U;
    average_transmission_U = (double *) calloc(numberofloop_U,sizeof(double));

    int binningnumber = 5;
    double sumW;
    double sumT;
    int j = 0;

    FILE *output_U;
    output_U=fopen("Filter/B_binned5.txt", "a");

    while (binningnumber+binningnumber*j< numberofloop_U) {
        sumW = 0.0;
        sumT = 0.0;
        for(i = 0; i < binningnumber; i++){
            sumW = sumW + wavelength_U[i+binningnumber*j];
            sumT = sumT + transmission_U[i+binningnumber*j];
        //    printf("binningnumber = %d\t sumW = %g\t sumT = %g\n", binningnumber, sumW, sumT);
        }
        average_wavelength_U[j] = sumW/binningnumber;
        average_transmission_U[j] = sumT/binningnumber;
        //printf("T = %g\t W = %g\n", average_transmission_U[j], average_wavelength_U[j]);

        fprintf (output_U, "%g\t%g\n", average_wavelength_U[j], average_transmission_U[j]);
        j=j+1;
    }
    fclose(output_U);

    
    free(wavelength_U);
    free(average_wavelength_U);
    free(transmission_U);
    free(average_transmission_U);

    return 0;
}
