/*--------------------------------------------------------------------------
! -- PHYS704
!-- Program : Advection (wave) equation by FTCS scheme
!--------------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <algorithm>

/***********************************************************************
 * Globals
 ***********************************************************************/ 
const double PI = atan(1.0)*4.0;

/* timestep*/
const double c = 1.;
const double dx = 5;
const double dt = .2;
const double cdt = c*dt;

const int m = 81;
const int nmax = 100; // (m-1)*dx/cdt;
const int nout=nmax/20;
const int mx = 4;           // number of grids to skip in output 

class Wave {
    public:
        double f0[m], f1[m], gc[m];
        Wave();
        double gF(int);         //general wave function -> F = c*f
        double gF_prime(int, int); 
        double gF_prime2(int, int);
        double f_prime(int);         // derivative of wave functin f'
        double v(int);
        double v_prime(int);       //derivative of c         
};

Wave::Wave() {
    f0 = {0.0};         // wave from previous timestep
    f1 = {0.0};         // calculated wave for current timestep
    std::fill_n(gc, m, c); // fill generalized c array
}

// Generalized wave function 
double Wave::gF(int i) {
   return f0[i] * v(i); 
}

double Wave::f_prime(int i) {
    return (f0[i] - f0[i-1])/dx;
}

double Wave::v(int i) {
    return 25.*(1-f0[i]);
}

double Wave::gF_prime(int i, int frac) {
    switch (frac) {
        case 1:
            return (v(i+1)*f0[i+1] - v(i)*f0[i]) / (f0[i+1] - f0[i]);
            break;
        case -1:
            return (v(i-1)*f0[i-1] - v(i)*f0[i]) / (f0[i-1] - f0[i]);
            break;
        case 0:
            return c;
            break;
        default:
            printf("you did something stupid\n");
            return 0.0;
    }
}

double Wave::gF_prime2(int i, int frac) {
    switch (frac) {
        case 1:
            return 25.*(1 - 2*(f0[i-1] + f0[i])/2);
            break;
        case -1:
            return 25.*(1 - 2*(f0[i+1] + f0[i])/2);
            break;
        case 0:
            return c;
            break;
        default:
            printf("you did something stupid\n");
            return 0.0;
    }
}

/* Initialisation of wave */
int initialize_sin_wave(Wave *wave) {

    for(int i=5;i<=55;i++){
        wave->f0[i]=sin(2*PI*(i-5)/50.0);
    }

    return 0;
}

int initialize_square_wave(Wave *wave) {

    int i = 5;

    // make square wave
    for(i=20;i<=39;i++){
        wave->f0[i]=1;
    }

    wave->f0[40] = .5;

    /*
    double p = 1;

    // make small slope in the front
    while (p > 0) {
        wave->f0[i] = p;
        p -= .1;
        i++;
    }
    */

    return 0;
}

/* FTCS scheme */
int ftcs_increment(Wave *wave) {

    for(int i=1;i<m-1;i++) {
        wave->f1[i] = wave->f0[i] - cdt/(2.0*dx) * (wave->f0[i+1] - wave->f0[i-1]);   
    }

    return 0;
}

/* Lax scheme */
int lax_increment(Wave *wave) {

    for(int i=1;i<m-1;i++) {
        wave->f1[i] = .5*(wave->f0[i+1] + wave->f0[i-1]) - cdt/(2.0*dx) * (wave->f0[i+1] - wave->f0[i-1]);   
    }

    return 0;
}

int lax_wendroff_increment(Wave *wave) {

    for(int i=1;i<m-1;i++) {
        wave->f1[i] = wave->f0[i] 
                      - cdt/(2.0*dx) * (wave->f0[i+1] - wave->f0[i-1]) 
                      + 2.0 * (cdt/(2.0*dx))*(cdt/(2.0*dx)) * (wave->f0[i+1] - 2.0*wave->f0[i] + wave->f0[i-1]);   
    }

    return 0;
}

/* Lax-Wendroff scheme */
int general_lax_wendroff_increment(Wave *wave) {

    for(int i=1;i<m-1;i++) {
        wave->f1[i] = wave->f0[i] 
                      - dt/(2.0*dx) * ( wave->gF(i+1) - wave->gF(i-1) ) 
                      + (dt)*(dt) / (2.0*dx) * (wave->gF_prime(i, 0) * ( wave->gF(i+1) - wave->gF(i) )/dx
                                                           - wave->gF_prime(i, 0) * ( wave->gF(i) - wave->gF(i-1))/dx);
    }

    return 0;
}

int move_traffic(Wave *wave) {

    for(int i=1;i<m-1;i++) {
        wave->f1[i] = wave->f0[i] 
                      - dt/(2.0*dx) * ( wave->gF(i+1) - wave->gF(i-1) ) 
                      + (dt)*(dt) / (2.0*dx) * (wave->gF_prime2(i, 1) *  ( wave->gF(i+1) - wave->gF(i))/dx
                                                           - wave->gF_prime2(i, -1) * ( wave->gF(i) - wave->gF(i-1))/dx);
    }

    wave->f1[m-1] = wave->f1[m-2];
    wave->f1[0] = wave->f1[m-1];

    return 0;
}


/* Output to file ___ */
int write_to_file(char label[], double f[], double t) {
    FILE *file;
    file=fopen(label,"a");
    double x = 0.0;

    for(int i=0; i<m; i+=mx){
        fprintf(file,"%16.8e  %16.8e  %16.8e\n", t, x, f[i]);
        x += dx*mx;
    }

    fprintf(file,"\n");
    fclose(file);

    return 0;
}

int main(void) {
    char label[20];
    class Wave wave;

    // parameters of the system 
    printf("\nnumber of timestep   = %d \n", nmax);
    printf("interval of snapshot = %e \n", nout*dt);
    printf("h = %16.8e \n c*tau = %16.8e \n\n", dx, cdt);

    // initialization 
    double t = 0.0;
    int nn = 0;                  // output counter
    //initialize_sin_wave(&wave);
    initialize_square_wave(&wave);

    // Output the initial values 
    //  (snapshot)
    sprintf(label,"advec_%05d.dat",nn);
    write_to_file(label, wave.f0, t);
    //  (t_x diagram)
    sprintf(label,"t_x_f.dat");
    write_to_file(label, wave.f0, t);
 
    /* loop for timestep */
    for(int n=1; n<=nmax; n++) { 

        t += dt;

        //ftcs_increment(&wave);
        //lax_increment(&wave);
        //lax_wendroff_increment(&wave);
        //general_lax_wendroff_increment(&wave);
        move_traffic(&wave);

        // Output information if at resolution
        if(n%nout == 0) {  
            nn += 1;
            
            // OUTPUT (snapshot)
            sprintf(label,"advec_%05d.dat",nn);
            write_to_file(label, wave.f1, t);

            // OUTPUT (t_x diagram)
            sprintf(label,"t_x_f.dat");
            write_to_file(label, wave.f1, t);
        }

        // update array f0 with new values
        for(int i=0; i<m; i++) { 
            wave.f0[i]=wave.f1[i];
        }
    }

    return 0;
}
