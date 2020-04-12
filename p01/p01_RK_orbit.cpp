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
double GM = 1.0;

/* pi = 3.141592...*/
double pi = atan(1.0)*4.0;

/*  maximum time-step, time-step (dt [s]) */
int maxstep = 200000;  
double dt = 0.1;
  

/***********************************************************************
 * functions for parameters
 ***********************************************************************/ 
double get_acceleration(double x, double y, double *ax, double *ay) {
    double r = sqrt(x*x + y*y);
    *ax = -(GM/(r*r*r)) * x;
    *ay = -(GM/(r*r*r)) * y;
    return 0;
}

double get_velocity(double initial, double acceleration, double tau = dt) {
    return initial + acceleration*tau;
}

double get_kinetic_energy(double vx, double vy){
    return (vx*vx+vy*vy)/2.0;
}

double get_potential_energy(double x, double y){
    double r = sqrt(x*x + y*y);
    return -GM/r;
}

double rk_integration(double x[], double y[], double vx[], double vy[], int i) {
    double F1[4], F2[4], F3[4], F4[4];

    get_acceleration(x[i], y[i], &F1[2], &F1[3]);
    F1[1] = vy[i];
    F1[0] = vx[i];

    get_acceleration(x[i] + F1[0]*.5*dt, y[i] + F1[1]*.5*dt, &F2[2], &F2[3]);
    F2[1] = get_velocity(vy[i], F1[3], .5*dt);
    F2[0] = get_velocity(vx[i], F1[2], .5*dt);

    get_acceleration(x[i] + F2[0]*.5*dt, y[i] + F2[1]*.5*dt, &F3[2], &F3[3]);
    F3[1] = get_velocity(vy[i], F2[3], .5*dt);
    F3[0] = get_velocity(vx[i], F2[2], .5*dt);

    get_acceleration(x[i] + F3[0]*dt, y[i] + F3[1]*dt, &F4[2], &F4[3]);
    F4[1] = get_velocity(vy[i], F3[3], dt);
    F4[0] = get_velocity(vx[i], F3[2], dt);
    
    vy[i+1] = vy[i] + 1/6. * (F1[3] + 2 * F2[3] + 2 * F3[3] + F4[3]) * dt;
    vx[i+1] = vx[i] + 1/6. * (F1[2] + 2 * F2[2] + 2 * F3[2] + F4[2]) * dt;

    y[i+1] = y[i] + 1/6. * (F1[1] + 2 * F2[1] + 2 * F3[1] + F4[1]) * dt;
    x[i+1] = x[i] + 1/6. * (F1[0] + 2 * F2[0] + 2 * F3[0] + F4[0]) * dt;

    return 0;
}

/***********************************************************************
 * Main Loop
 ***********************************************************************/ 
int main(void) {
    double ax,ay,ek,ep,et;
    double t;
    int i;
  
    /* define arrays */ 
    double x[maxstep+1],y[maxstep+1],vx[maxstep+1],vy[maxstep+1];
  
    FILE *fp;
  
    /* open an output file */
    fp=fopen("orbit_1_4.dat","w");
  
    /* initial position (0.0,1.0m) & velocity */
    x[0]=1; y[0]=0.0; vx[0]= 0.0; vy[0]=1.4;
  
    ek = get_kinetic_energy(vx[i], vy[i]);
    ep = get_potential_energy(x[i], y[i]);
    et = ek + ep;
  
    /* output the initial values */
    i=0;  t=0.0;
    fprintf(fp,"%16.8e %16.8e %16.8e %16.8e %16.8e  %16.8e  %16.8e  %16.8e\n",t,x[0],y[0],vx[0],vy[0],ek,ep,et);

    for(i=0;i<=maxstep;i++){

        rk_integration(x, y, vx, vy, i);
 
        //x[i+1] = x[i] + F1[1] *dt
        /*  energy calculations */
        ek = get_kinetic_energy(vx[i+1], vy[i+1]);
        ep = get_potential_energy(x[i+1], y[i+1]);
        et = ek + ep;
    
        t = t + dt;
    
        /* ouput the updated values */
        fprintf(fp,"%16.8e %16.8e %16.8e %16.8e %16.8e  %16.8e  %16.8e  %16.8e\n",t,x[i+1],y[i+1],vx[i+1],vy[i+1],ek,ep,et);
    }
  
    fclose(fp);
    
    return 0;
}
