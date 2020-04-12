#include <stdio.h>
#include <string>
#include <iostream>
#include <math.h>

/***********************************************************************
 * Project 01: Orbit of satellite
 ***********************************************************************/ 
/* 
 * To change between Euler, RK4, and adaptive RK: change in function 'increment'
 */

/***********************************************************************
 * Globals
************************************************************************/

/* gravitational constant g=9.8 [kg/s^2] */
const double GM = 1.0;

/* pi = 3.141592...*/
const double pi = atan(1.0)*4.0;

/*  maximum time-step, time-step (dt [s]) */

const int maxstep = 40;  
const double dt = 0.1;

/* acceptable error */
const double err = 0.001;

/* Safty Step */
const double S1 = 0.9;
const double S2 = 4.;
  

/***********************************************************************
 * functions for parameters
 ***********************************************************************/ 
double get_acceleration(double x, double y, double *ax, double *ay) {
    double r = sqrt(x*x + y*y);
    *ax = -(GM/(r*r*r)) * x;
    *ay = -(GM/(r*r*r)) * y;
    return 0;
}

double get_velocity(double initial, double acceleration, double tau) {
    return initial + acceleration * tau;
}

double get_kinetic_energy(double vx, double vy){
    return (vx*vx+vy*vy)/2.0;
}

double get_potential_energy(double x, double y){
    double r = sqrt(x*x + y*y);
    return -GM/r;
}

double rk_integration(double temp[], double x, double y, double vx, double vy, double tau) {
    double F1[4], F2[4], F3[4], F4[4];

    get_acceleration(x, y, &F1[2], &F1[3]);
    F1[1] = vy;
    F1[0] = vx;

    get_acceleration(x + F1[0]*.5*tau, y + F1[1]*.5*tau, &F2[2], &F2[3]);
    F2[1] = get_velocity(vy, F1[3], .5*tau);
    F2[0] = get_velocity(vx, F1[2], .5*tau);

    get_acceleration(x + F2[0]*.5*tau, y + F2[1]*.5*tau, &F3[2], &F3[3]);
    F3[1] = get_velocity(vy, F2[3], .5*tau);
    F3[0] = get_velocity(vx, F2[2], .5*tau);

    get_acceleration(x + F3[0]*tau, y + F3[1]*tau, &F4[2], &F4[3]);
    F4[1] = get_velocity(vy, F3[3], tau);
    F4[0] = get_velocity(vx, F3[2], tau);
    
    temp[3] = vy + 1/6. * (F1[3] + 2 * F2[3] + 2 * F3[3] + F4[3]) * tau;
    temp[2] = vx + 1/6. * (F1[2] + 2 * F2[2] + 2 * F3[2] + F4[2]) * tau;
    temp[1] = y  + 1/6. * (F1[1] + 2 * F2[1] + 2 * F3[1] + F4[1]) * tau;
    temp[0] = x  + 1/6. * (F1[0] + 2 * F2[0] + 2 * F3[0] + F4[0]) * tau;

    return 0;
}

/***********************************************************************
 * Main Loop
 ***********************************************************************/ 
int main(void) {
    double ek,ep,et;
    double delta_c, relerr;
    double t;
    int i;
  
    /* define arrays */ 
    double x[maxstep+1],y[maxstep+1],vx[maxstep+1],vy[maxstep+1];
  
    FILE *fp;
  
    /* open an output file */
    fp=fopen("orbit.dat","w");
  
    /* initial position (0.0,1.0m) & velocity */
    x[0]=1.; y[0]=0.0; vx[0]= 0.0; vy[0]=0.25; 
  
    ek = get_kinetic_energy(vx[0], vy[0]);
    ep = get_potential_energy(x[0], y[0]);
    et = ek + ep;
  
    /* output the initial values */
    i=0;  t=0.0;
    fprintf(fp,"%16.8e %16.8e %16.8e %16.8e %16.8e  %16.8e  %16.8e  %16.8e\n",t,x[0],y[0],vx[0],vy[0],ek,ep,et);

    // time step that is updated each iteration
    double tau = dt;
    double tau_old = dt;

    for(i=0;i<=maxstep;i++){

        /* 
         * temp = { x, y, vx, vy }
         * for full step and half steps
         *
         * */
        double temp_full[4]  = {0};
        double temp_half[4]  = {0};

        //full step
        rk_integration(temp_full, x[i], y[i], vx[i], vy[i], tau);

        //Two half steps
        rk_integration(temp_half, x[i], y[i], vx[i], vy[i], .5*tau);
        rk_integration(temp_half, temp_half[0], temp_half[1], temp_half[2], temp_half[3], .5*tau);

        //  errors
        relerr = sqrt( pow(temp_full[0], 2)  + 
                       pow(temp_full[1], 2)  +
                       pow(temp_full[2], 2) + 
                       pow(temp_full[3], 2));
       
        delta_c = sqrt( pow(temp_full[0] - temp_half[0], 2) +
                        pow(temp_full[1] - temp_half[1], 2) +
                        pow(temp_full[2] - temp_half[2], 2) + 
                        pow(temp_full[3] - temp_half[3], 2)) / relerr;

        // increment time and calculate the time step for next iteration
        t += tau;
        tau_old = tau;
        tau = S1 * dt * pow( fabs( err / ( delta_c + pow(10,-10) ) ), .2);
        
        // Check for drastic changes
        if (tau > (S2 * tau_old) ) {
            tau = S2 * tau_old;
        }
        if (tau < ( tau_old / S2) ) {
            tau = tau_old / S2;
        }

        // update position
        // update velocity
        if (delta_c < err) {
            x[i+1]  = temp_full[0];
            y[i+1]  = temp_full[1];
            vx[i+1] = temp_full[2];
            vy[i+1] = temp_full[3];
        }
        else {
            x[i+1]  = temp_half[0];
            y[i+1]  = temp_half[1];
            vx[i+1] = temp_half[2];
            vy[i+1] = temp_half[3];
        }

        /*  energy calculations */
        ek = get_kinetic_energy(vx[i+1], vy[i+1]);
        ep = get_potential_energy(x[i+1], y[i+1]);
        et = ek + ep;
    
        /* ouput the updated values */
        fprintf(fp,"%16.8e %16.8e %16.8e %16.8e %16.8e  %16.8e  %16.8e  %16.8e\n",t,x[i+1],y[i+1],vx[i+1],vy[i+1],ek,ep,et);
    }
  
    fclose(fp);
    
    return 0;
}
