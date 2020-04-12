#include <stdio.h>
#include <fstream>
#include <iostream>         // std::cout
#include <math.h>
#include <cmath>        // std::abs
#include <algorithm>    // std::random_shuffle
#include <vector>           // std::vector

const double PI = atan(1.0)*4;      // pi -> 3.1415926 etc
const int particle_num = 10000;       // number of particles in the simulation 
const int steps = 100;              // number of steps the program runs for
const double th = 0.0001;           // thresh hold to consider a collision
const int n_out = 20;

const int max_bin = 1000;

// bin size of count, range of count (vmin,vmax)
const double bin = 0.1;
const double vmax = 8.0;
const double vmin = -8.0;
const double nbin = round( (vmax-vmin)/bin );

// Particle objects
// contains information on velocity in each direction
// return the energy of the particle
class Particle {
    public:
        double vx;  // velocity in x
        double vy;  // velocity in y
        double vz;  // velocity in z
        double m;   // mass
        int tracer;
        double kinetic_energy();
        double x_momentum();
        double y_momentum();
        double z_momentum();
        Particle();
};

// initial values for a particle 
// vx,vy,vz are assigned a random value (-1,1) in the initial distribution function
// mass is set to 1
Particle::Particle() {
    vx = 0;
    vy = 0;
    vz = 0;
    m = 1;
    tracer = 0;
}

double Particle::kinetic_energy() {
    return .5 * m * (vx*vx + vy*vy + vz*vz);
}

double Particle::x_momentum() {
    return m * vx;
}

double Particle::y_momentum() {
    return m * vy;
}

double Particle::z_momentum() {
    return m * vz;
}

// Random number generator; Uniform dist. in [0,1)
double randm( int *seed ) {
    // Input
    //   seed    Integer seed (DO NOT USE A SEED OF ZERO)
    // Output
    //	 rand    Random number uniformly distributed in [0,1)

    const double a = 16807.0;
    const double m = 2147483647.0;
    double temp = a * (*seed);
    *seed = (int)(fmod(temp,m));
    double randm = *seed/m;

    return randm;
}

// generate an array of tuples [(A0,B0), (A1,B1), ...]
int zip(int A[], Particle B[], std::vector<std::pair<int,Particle>> &zipped){
    
    for (int i=0; i<particle_num; i++){
        zipped.push_back(std::make_pair(A[i], B[i]));
    }

    return 0;
}

// Take the values from the sorted tuple array and assign them to the initial arrays 
int unzip(std::vector<std::pair<int, Particle>> &zipped, int A[], Particle B[]){

    for(int i=0; i<particle_num; i++){
        A[i] = zipped[i].first;
        B[i] = zipped[i].second;
    }

    return 0;
}

// Conditional to sort by the first element of the pair 
bool pairCompare(const std::pair<int, Particle>& firstElem, const std::pair<int, Particle>& secondElem) {
      return firstElem.first < secondElem.first;
}

// Return True if an element is in a given array
// !!!!!!!!!! NOT USED RIGHT NOW !!!!!!!!
bool in_array(int x, int array[]) {

    // number of elements in array
    int n = sizeof(array)/sizeof(int);

    for(int i=0; i<n; i++){
        if (x == array[i]) {
            return true;
        }
    }

    return false; 
}

// randomly shuffle the indices in a given Particle array
int shuffle(Particle particles[]) {
    

    int rand_index[particle_num]={0};
    std::vector<std::pair<int, Particle>> zipped;
    
    // fill r with random integers between 0 and the number of particles.
    for (int i=0; i < particle_num; i++) {
        rand_index[i] = (rand() % particle_num);
    }

    zip(rand_index, particles, zipped);
    std::sort(zipped.begin(), zipped.end(), pairCompare);
    unzip(zipped, rand_index, particles);

    return 0;
}

// Generate a randomized distribution of velocities
// generates random value between vmin and vmax for a velocity in each direction 
int initial_distribution(Particle particles[], int *seed) {
    /*
    for (int i=0; i<particle_num; i++) {

        particles[i].vx = 2*vmax*(randm(seed) - .5);
        particles[i].vy = 2*vmax*(randm(seed) - .5);
        particles[i].vz = 2*vmax*(randm(seed) - .5);
        particles[i].tracer = i;
    }
    */

    for(int i=0;i<particle_num/2;i++){
        particles[i].vx = sqrt(-2*log(1 - randm(seed))) * cos(2*PI*randm(seed)) + 2;
        particles[i].vy = sqrt(-2*log(1 - randm(seed))) * cos(2*PI*randm(seed));
        particles[i].vz = sqrt(-2*log(1 - randm(seed))) * cos(2*PI*randm(seed));
    }

    for(int i=particle_num/2;i<particle_num;i++){
        particles[i].vx = sqrt(-2*log(1 - randm(seed))) * cos(2*PI*randm(seed)) - 2;
        particles[i].vy = sqrt(-2*log(1 - randm(seed))) * cos(2*PI*randm(seed));
        particles[i].vz = sqrt(-2*log(1 - randm(seed))) * cos(2*PI*randm(seed));
    }

    /*
    for(int i=0;i<particle_num/3;i++){
        particles[i].vx = sqrt(-2*log(1 - randm(seed))) * cos(2*PI*randm(seed)) + 2;
        particles[i].vy = sqrt(-2*log(1 - randm(seed))) * cos(2*PI*randm(seed));
        particles[i].vz = sqrt(-2*log(1 - randm(seed))) * cos(2*PI*randm(seed)) - 2;
    }

    for(int i=particle_num/3;i<particle_num*2./3;i++){
        particles[i].vx = sqrt(-2*log(1 - randm(seed))) * cos(2*PI*randm(seed)) - 2;
        particles[i].vy = sqrt(-2*log(1 - randm(seed))) * cos(2*PI*randm(seed));
        particles[i].vz = sqrt(-2*log(1 - randm(seed))) * cos(2*PI*randm(seed)) - 2;
    }

    for(int i=particle_num*2./3;i<particle_num;i++){
        particles[i].vx = sqrt(-2*log(1 - randm(seed))) * cos(2*PI*randm(seed));
        particles[i].vy = sqrt(-2*log(1 - randm(seed))) * cos(2*PI*randm(seed));
        particles[i].vz = sqrt(-2*log(1 - randm(seed))) * cos(2*PI*randm(seed)) + 2;
    }
    */

    return 0;
}

// Checks to see if two particles will collide with each other 
bool does_collide(double max_relv, Particle A, Particle B) {

    double v = 0;

    // magnitude of relative velocity  
    v = std::sqrt((A.vx - B.vx) * (A.vx - B.vy) + 
                  (A.vy - B.vy) * (A.vx - B.vy) + 
                  (A.vz - B.vz) * (A.vx - B.vy)); 

    // The pair is accepted as collision partner if relative speed is greater than the set threshold times a
    // random number [0,1)
    return v > max_relv * (rand() % 100)/100.;
}

double rel_v(Particle *A, Particle *B) {

    double vrx = A->vx - B->vx;
    double vry = A->vy - B->vy;
    double vrz = A->vz - B->vz;

    return sqrt(vrx*vrx + vry*vry + vrz*vrz);
}

// calculate exchange of energy
int exchange_energy(Particle *A, Particle *B) {

    double cm[3] = { (A->vx + B->vx)*.5, (A->vy + B->vy)*.5, (A->vz + B->vz)*.5};
    double vr2[3] = {0};

    double phi = 2.*PI*(rand() % 100)/100.;
    double costheta = (rand() % 100)/100.;
    double sintheta = sqrt(1 - costheta*costheta);

    double vr = rel_v(A, B);

    vr2[0] = vr*sintheta*cos(phi);
    vr2[1] = vr*sintheta*sin(phi);
    vr2[2] = vr*costheta;


    A->vx = cm[0] + .5*vr2[0];
    A->vy = cm[1] + .5*vr2[1];
    A->vz = cm[2] + .5*vr2[2];
    
    B->vx = cm[0] - .5*vr2[0];
    B->vy = cm[1] - .5*vr2[1];
    B->vz = cm[2] - .5*vr2[2];

    return 1;
}



// check the collision of particles in the randomly shuffled array 
int check_collisions(double max_relv, Particle particles[]) {

    for (int i=1; i < particle_num; i+=2){

        // if the particles collide, exchange momentum/energy
        // else move to next pair
        if ( does_collide(max_relv, particles[i-1], particles[i]) ){
            exchange_energy(&particles[i-1], &particles[i]);
        }
    }

    return 1;
}

double get_max_relativle_velocity(Particle particles[]){
    
    double max_v = 0;
    double vr = 0;

    for (int i=1; i < particle_num; i+=2) {

        vr = rel_v(&particles[i], &particles[i-1]);

        if ( vr > max_v) {
            max_v = vr;
        }
    }

    return max_v;
}


// return the momentum of the particle 
double momentum(double p[]) {
    return sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
}

// return Kinetic Energy of the particle 
double energy(Particle particle) {

    return .5 * ( particle.vx * particle.vx 
                + particle.vy * particle.vy
                + particle.vz * particle.vz );
}

// Record information about the system every n_out steps 
void diag(Particle particles[]) {
    int ip;
    static int i_time = 0;
    double counter[max_bin] = {0};
    double total_E = 0, total_P = 0, Z = 0;
    double p[3] = {0};
    //  char fname[15];
    char *fname = new char[32];
    FILE *fp;
    //std::fstream fs;

    std::cout << i_time << std::endl;

    // Create phase file
    sprintf(fname, "phase_%05d.dat", i_time);
    fp = fopen(fname, "w");

    //fs.open (("phase_%05d.dat", i_time), std::fstream::in);

    for (int i=0 ; i<particle_num ; i++) {

        fprintf(fp, "%5d  %16.8e  %16.8e  %16.8e\n", i, particles[i].vx, particles[i].vy, particles[i].vz);

        //fs << "%5d  %16.8e  %16.8e  %16.8e\n", i, particles[i].vx, particles[i].vy, particles[i].vz << endl;

        // Counting for px histogram
        ip = round(particles[i].vx/bin) + nbin/2;

        //if ((ip < 40) && (i_time == 0)) printf("%6d  %16.8e", ip, counter[ip]);

        if ((ip>=0) && (ip<=nbin)) {
            //if ((ip < 40) && (i_time == 0)) printf("  in  ");
            counter[ip] += 1.0;
        }
        //if ((ip < 40) && (i_time == 0)) printf("%6d  %16.8e\n", ip, counter[ip]);

        double e = energy(particles[i]);
        total_E += e;
        Z += exp(-e);
        p[0] += particles[i].vx;
        p[1] += particles[i].vy;
        p[2] += particles[i].vz;
    }

    //fs.close();

    fclose(fp);

    // Create count file
    sprintf(fname, "count_%05d.dat", i_time);
    fp = fopen(fname, "w");

    for (int i=0 ; i<=nbin ; i++) {
        //if (i_time == 0) printf("%6d  %16.8e\n", i, counter[i]);
        fprintf(fp, "%16.8e  %16.8e\n", vmin+bin*i, counter[i]);
    }

    fclose(fp);


    // Conservation of Energy and Momentum diag

    double A = -log(Z);
    double S = total_E - A;
    total_P = momentum(p);
    sprintf(fname, "conservation.dat", i_time);
    fp = fopen(fname, "a");
    fprintf(fp, "%5d  %16.8e  %16.8e  %16.8e  %16.8e\n", i_time, total_E, total_P, A, S);  
    fclose(fp);

    i_time++;
    delete[] fname;
}

int main(void) {
    Particle particles[particle_num];
    int seed = 3;

    // generate <particle_num> number of particles with velocities in x,y,z
    initial_distribution(particles, &seed);

    /*
    for (int j=0; j<particle_num; j++){
        std::cout << particles[j].num << ", ";
    }
    std::cout << std::endl;
    */

    diag(particles);

    // iterate through number of steps in the simulation
    // shuffled the array randomly, and compare each particle to its neighbor
    // if there is a collision - calculate transfer of energy
    for(int i=0; i < steps; i++) {
        double max_relv = 0;

        shuffle(particles);

        max_relv = get_max_relativle_velocity(particles);
        
        // print sorted random array
        /*
        std::cout << "Particles: ";
        for (int j=0; j<particle_num; j++){
            std::cout << particles[j].num << ", ";
        }
        std::cout << std::endl;
        */

        // record information
        if ( i % (steps / n_out) == 0 ){
            diag(particles);
        }

        // check collision for all neighboring particles
        // if collision, calculate exchange of energy 
        check_collisions(max_relv, particles);
    }

    return 1;
}
