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
double L_star(double L_bol, double omega, double t){
    //omega = 3*c/R_out
    return 0.15*L_bol*(1.0+sin(omega*t));
}

/** define the Function of the distance from the central variable source to disk elements **/
double r_star(double r, double h_star){
    return sqrt(pow(h_star,2.0)+pow(r,2.0));
}


/** define the Function of the temperature profile **/
double temp_profile(double t, double r, double theta, double M, double M_rate, double r_in, double A, double h_star, double inc_angle, double L_bol, double omega){
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
} region;


int make_computation(int Nfilter, int *computed_filter){


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








    /**  Husne,  9/10/2018
      *
      */
    /** for the computation of luminosity so it is for temperature **/
    double omega = 10.0*c/r_out;
    /** for the computation of the radiation from the disk. **/
    double D = 75.01*1e6*pc;                   /** Mpc distance from observer to the source, converted to cm **/
    double R_in;
    double R_out;
    double theta_in;
    double theta_out;




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
    //double filtername[6] = {0, 1, 2, 3, 4, 5}; //* filter names: 0=UVW2, 1=UVM2, 2=UVW1, 3=U, 4=B, 5=V */
    //int Nfilter = 6;
    double **wavelength;
    wavelength = (double **) malloc(Nfilter*sizeof(double*)); //* create an array */
    double **transmission;
    transmission = (double **) malloc(Nfilter*sizeof(double*));
    double c1_filtername, c2_filtername;
    int numberofloop_filtername;
    int *numberofloop;
    numberofloop = (int*) calloc(Nfilter,sizeof(int)); //* create an array */

    for (j=0; j < Nfilter; j++){
        //printf("Begining of loop \t j = %d\tNfilter = %d\n", j, Nfilter);
        //getchar();
        FILE *input_filtername;
        numberofloop_filtername = 0.0;
        //*it is important, when the filter number given as a "0" make computation*/
        if (computed_filter[j] == 0){
            switch(j) {
                case 0 : //*it is UVW2 filter then*/
                    input_filtername=fopen("UVW2.txt","r");//* open a text file for reading */
                    /**  Here %lf means type double */
                    /// step 1 caunt the number of loop
                    while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                        numberofloop_filtername = numberofloop_filtername + 1;                 //* caunt the number of loop */
                        /// Previous line is equivalent to      numberofloop_U = numberofloop_U + 1;                 //* caunt the number of loop */
                    }
                    numberofloop[j] = numberofloop_filtername;
                    break;
                case 1 : //*it is UVM2 filter then*/
                    input_filtername=fopen("UVM2.txt","r");      //* open a text file for reading */
                    /**  Here %lf means type double */
                    /// step 1
                    while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                        numberofloop_filtername++;                 //* caunt the number of loop */
                        /// Previous line is equivalent to      numberofloop_U = numberofloop_U + 1;                 //* caunt the number of loop */
                    }
                    numberofloop[j] = numberofloop_filtername;
                    break;
                case 2 : //*it is UVW1 filter then*/
                    input_filtername=fopen("UVW1.txt","r");      //* open a text file for reading */
                    /**  Here %lf means type double */
                    /// step 1
                    while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                        numberofloop_filtername++;                 //* caunt the number of loop */
                        /// Previous line is equivalent to      numberofloop_U = numberofloop_U + 1;                 //* caunt the number of loop */
                    }
                    numberofloop[j] = numberofloop_filtername;
                    break;
                case 3 : //*it is U filter then*/
                    input_filtername=fopen("U.txt","r");      //* open a text file for reading */
                    /**  Here %lf means type double */
                    /// step 1
                    while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                        numberofloop_filtername++;                 //* caunt the number of loop */
                        /// Previous line is equivalent to      numberofloop_U = numberofloop_U + 1;                 //* caunt the number of loop */
                    }
                    numberofloop[j] = numberofloop_filtername;
                    break;
                case 4 : //*it is B filter then*/
                    input_filtername=fopen("B.txt","r");      //* open a text file for reading */
                    /**  Here %lf means type double */
                    /// step 1
                    while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                        numberofloop_filtername++;                 //* caunt the number of loop */
                        /// Previous line is equivalent to      numberofloop_U = numberofloop_U + 1;                 //* caunt the number of loop */
                    }
                    numberofloop[j] = numberofloop_filtername;
                    break;
                case 5 : //*it is V filter then*/
                    input_filtername=fopen("V.txt","r");      //* open a text file for reading */
                    /**  Here %lf means type double */
                    /// step 1
                    while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                        numberofloop_filtername++;                 //* caunt the number of loop */
                        /// Previous line is equivalent to      numberofloop_U = numberofloop_U + 1;                 //* caunt the number of loop */
                    }
                    numberofloop[j] = numberofloop_filtername;
                    break;
            }
            fclose(input_filtername);
        }

        /// step 2 to create arrays
        wavelength[j] = (double *) calloc(numberofloop[j],sizeof(double)); //* create an array */
        transmission[j] = (double *) calloc(numberofloop[j],sizeof(double));

        /// step 3 fill the arrays
        //*it is important, when the filter number given as a "0" make computation*/
        if (computed_filter[j] == 0){
            switch(j) {
                case 0 : //*it is UVW2 filter then*/
                    input_filtername=fopen("UVW2.txt","r");//* open a text file for reading */
                    /**  Here %lf means type double */
                    /// step 1
                    i = 0;
                    while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                        /// i = i + 1 ;
                        wavelength[j][i]= c1_filtername;
                        transmission[j][i]= c2_filtername;
                        i += 1;
                    }
                    break;
                case 1 : //*it is UVM2 filter then*/
                    input_filtername=fopen("UVM2.txt","r");      //* open a text file for reading */
                    /**  Here %lf means type double */
                    /// step 1
                    i = 0;
                    while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                        /// i = i + 1 ;
                        wavelength[j][i]= c1_filtername;
                        transmission[j][i]= c2_filtername;
                        i += 1;
                    }
                    break;
                case 2 : //*it is UVW1 filter then*/
                    input_filtername=fopen("UVW1.txt","r");      //* open a text file for reading */
                    /**  Here %lf means type double */
                    /// step 1
                    i = 0;
                    while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                        /// i = i + 1 ;
                        wavelength[j][i]= c1_filtername;
                        transmission[j][i]= c2_filtername;
                        i += 1;
                    }
                    break;
                case 3 : //*it is U filter then*/
                    input_filtername=fopen("U.txt","r");      //* open a text file for reading */
                    /**  Here %lf means type double */
                    /// step 1
                    i = 0;
                    while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                        /// i = i + 1 ;
                        wavelength[j][i]= c1_filtername;
                        transmission[j][i]= c2_filtername;
                        i += 1;
                    }
                    break;
                case 4 : //*it is B filter then*/
                    input_filtername=fopen("B.txt","r");      //* open a text file for reading */
                    /**  Here %lf means type double */
                    /// step 1
                    i = 0;
                    while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                        /// i = i + 1 ;
                        wavelength[j][i]= c1_filtername;
                        transmission[j][i]= c2_filtername;
                        i += 1;
                    }
                    break;
                case 5 : //*it is V filter then*/
                    input_filtername=fopen("V.txt","r");      //* open a text file for reading */
                    /**  Here %lf means type double */
                    /// step 1
                    i = 0;
                    while(fscanf(input_filtername,"%lf%lf", &c1_filtername, &c2_filtername) !=EOF ){
                        /// i = i + 1 ;
                        wavelength[j][i]= c1_filtername;
                        transmission[j][i]= c2_filtername;
                        i += 1;
                    }
                    break;
            }
            fclose(input_filtername);
        }
    }


    //for (j=0;j<Nfilter;j++){
    //    for (i=0;i<numberofloop[j];i++){
    //        printf("j = %d\t i = %d\t %g\t%g\n",j, i, wavelength[j][i], transmission[j][i]);
   //     }
   //     getchar();
   //     printf("\n");
   // }







    /** ************************************************
     * ************************************************
     * ************************************************
     * ************************************************
     * ************************************************
     */





    /**  Husne, 20/10/2018
     * I define the time and tau_time as arrays
     */
    // 1) SET the tau so give the value of tau
    const int Ntau = 7;
    double tau_time[7] = {3.0, 6.0, 10.0, 20.0, 40.0, 100.0, 200.0};

    //printf("%g\t\n",tau_time);  //* print the arrays */


    int Ntime = 10;            /** 5000days*86400 = seconds **/
    double *time;
    time = (double *) calloc(Ntime,sizeof(double));
    for (i=0; i<Ntime; i++){
        time[i] = i/2.0;
        //printf("time=%g\n", t[i]);
    }

    // a) make the sum of all the S(tau, t)*1/Ntime
    // i) compute all S for all t
    /**  Husne,  19/10/2018 *  compute the color variability and plot them */
    double avarage_SBU = 0.0;
    double S_BU;            /**color_variation of BU */



    /**  Husne,  18/10/2018 *  Now compute the integral for U band. */

    double deltaLambda_U;

    /**  Husne,  19/10/2018 *  Now compute the integral for B band.*/

    double deltaLambda_B;



    double flux_t_U;
    double flux_t_B;
    double flux_tptau_U;
    double flux_tptau_B;

    double Temperature_t;
    double Temperature_tptau;
    double Integral;

    double f_U_im1 = 0.0;
    double f_B_im1 = 0.0;
    double f_U_i = 0.0;
    double f_B_i = 0.0;





    /** Loop for the Number of filter, because I need to compute the average S for all tau by using two bands */
    int m;
    int n;
    for (m=0;m<Nfilter;m++){
        for (n=m+1;n<Nfilter;n++){
            //*it is important, when the filter number given as a "0" make computation for all casses between two filters */
            if(computed_filter[m] == 0 && computed_filter[n] == 0){
                /** Loop for the tau, because I need to compute the average S for all tau */
                int l;
                int k;
                for (l=0; l < Ntau; l++){
                    S_BU = 0.0;
                    printf("\n");
                    /** Loop for the time, because I need to compute an average with respect to time */
                    for (k=0; k < Ntime; k++){

                        /**  I need to compute S for all t
                         *  I need to compute f(t+tau), f(t) for U band
                         */

                        flux_t_U = 0.0;
                        flux_t_B = 0.0;
                        flux_tptau_U = 0.0;
                        flux_tptau_B = 0.0;

                        /** Loop for the redius and theta, because I need to compute the temparature and spectrum of disk */
                        /// f is the summ of constribution from all the disk elements.
                        for (j=0; j < Nr*Ntheta; j++){
                            R_in = disk[j].radius/sqrt(step);            /** from the center to the first layer of any region **/
                            R_out = disk[j].radius*sqrt(step);           /** from the center to the last layer of any region **/
                            theta_in = disk[j].theta-(step/2.0);         /** from the origine to the first layer of any region on the bottom**/
                            theta_out = disk[j].theta+(step/2.0);        /** from the origine to the last layer of any region on the top**/

                            /**  Now I compute the integral for the U-band */
                            /// temperature at time t in U
                            Temperature_t = temp_profile(time[k], disk[j].radius, disk[j].theta, M, M_rate, r_in, A, h_star, inc_angle, L_bol, omega);
                            /// Initialization of the sum to compute the integral over the filter
                            Integral = 0.0;

                            /** Loop for the band, because I need to compute the integral over bandpass */
                            ///numberofloop_U = numberofloop[m]
                            for(i = 1; i < numberofloop[m] ; i++){

                                deltaLambda_U = (wavelength[m][i] - wavelength[m][i-1])*angstrom;
                                f_U_i = spectrum(inc_angle, D, theta_in, theta_out, R_in, R_out, wavelength[m][i]*angstrom, Temperature_t)*transmission[m][i];
                                f_U_im1 = spectrum(inc_angle, D, theta_in, theta_out, R_in, R_out, wavelength[m][i-1]*angstrom, Temperature_t)*transmission[m][i-1];

                                Integral += (f_U_i+f_U_im1)*deltaLambda_U/2.0;
                            }

                            flux_t_U = flux_t_U + Integral;


                            /// temperature at time t + tau in U
                            Temperature_tptau = temp_profile(time[k] + tau_time[l], disk[j].radius, disk[j].theta, M, M_rate, r_in, A, h_star, inc_angle, L_bol, omega);
                            /// Initialization of the sum to compute the integral over the filter
                            Integral = 0.0;
                            /** Loop for the one band, because I need to compute the integral over bandpass */
                            ///numberofloop_U = numberofloop[m]
                            for(i = 1; i < numberofloop[m] ; i++){

                                deltaLambda_U = (wavelength[m][i] - wavelength[m][i-1])*angstrom;
                                f_U_i    = spectrum(inc_angle, D, theta_in, theta_out, R_in, R_out, wavelength[m][i]*angstrom, Temperature_tptau)*transmission[m][i];
                                f_U_im1  = spectrum(inc_angle, D, theta_in, theta_out, R_in, R_out, wavelength[m][i-1]*angstrom, Temperature_tptau)*transmission[m][i-1];

                                Integral += (f_U_i+f_U_im1)*deltaLambda_U/2.0;
                            }
                            flux_tptau_U = flux_tptau_U + Integral;





                            /**  Now I compute the integral for the B-band */
                            /// temperature at time t in B
                            Temperature_t = temp_profile(time[k], disk[j].radius, disk[j].theta, M, M_rate, r_in, A, h_star, inc_angle, L_bol, omega);
                            /// Initialization of the sum to compute the integral over the filter
                            Integral = 0.0;
                            /** Loop for the another band, because I need to compute the integral over bandpass */
                            ///numberofloop_B = numberofloop[n]
                            for(i = 1; i < numberofloop[n] ; i++){

                                deltaLambda_B = (wavelength[n][i] - wavelength[n][i-1])*angstrom;
                                f_B_i    = spectrum(inc_angle, D, theta_in, theta_out, R_in, R_out, wavelength[n][i]*angstrom, Temperature_t)*transmission[n][i];
                                f_B_im1  = spectrum(inc_angle, D, theta_in, theta_out, R_in, R_out, wavelength[n][i-1]*angstrom, Temperature_t)*transmission[n][i-1];

                                Integral += (f_B_i+f_B_im1)*deltaLambda_B/2.0;
                            }

                            flux_t_B = flux_t_B + Integral;


                            /// temperature at time t + tau in U
                            Temperature_tptau = temp_profile(time[k] + tau_time[l], disk[j].radius, disk[j].theta, M, M_rate, r_in, A, h_star, inc_angle, L_bol, omega);
                            /// Initialization of the sum to compute the integral over the filter
                            Integral = 0.0;
                            /** Loop for the another band, because I need to compute the integral over bandpass */
                            ///numberofloop_B = numberofloop[n]
                            for(i = 1; i < numberofloop[n] ; i++){

                                deltaLambda_B = (wavelength[n][i] - wavelength[n][i-1])*angstrom;
                                f_B_i    = spectrum(inc_angle, D, theta_in, theta_out, R_in, R_out, wavelength[n][i]*angstrom, Temperature_tptau)*transmission[n][i];
                                f_B_im1  = spectrum(inc_angle, D, theta_in, theta_out, R_in, R_out, wavelength[n][i-1]*angstrom, Temperature_tptau)*transmission[n][i-1];

                                Integral += (f_B_i+f_B_im1)*deltaLambda_B/2.0;
                            }
                            flux_tptau_B = flux_tptau_B + Integral;
                        }

                        S_BU = S_BU + (flux_tptau_B-flux_t_B) / (flux_tptau_U-flux_t_U);

                    }
                    avarage_SBU = S_BU/Ntime;
                    printf("%.13g\t\n", avarage_SBU);
                }
            }
        }
    }

    free(r);
    free(theta);
    free(disk);
    free(wavelength);
    free(transmission);


    return 0;
}





int main(){


    /** Define argument for the filters */
    int Nfilter = 6;
    //int i;
    //int j;
    int computed_filter[6] = {1, 0, 1, 0, 0, 1};
    //for(i = 0; i < Nfilter ; i++){
        //for(j = i+1; j < Nfilter ; j++){
            //if(computed_filter[i] == 0 && computed_filter[j] == 0){

            //}
        //}
   // }


    make_computation(Nfilter, computed_filter);

    return 0;

}




