#include <iostream>
#include "constants.h"
#include <armadillo>
#include <cmath>
#include <math.h>
#include <fstream>

using namespace std;
using namespace arma;

void RK4(vec *r, vec *v, double distance, double dt);
vec f(vec r, double distance);



/*
 1. Constants
 2. Implement RK4 in main
 3. Solve earth-sun
 4. Make planet class
 5. Solve earth-sun
 6. Make Solver class
 7. Solve earth-sun
 8. Device tests for the system
 9. Manually implement Mars
 10. Generally implement Mars
 11. Generally implement solar system
 12. If time: Open GL
*/

int main()
{
    double Tfinal = 40; //years
    double delta_t = 0.0001; //years
    int N = Tfinal/delta_t;

    // Initializing earth-sun
    vec X = zeros<vec>(N); vec Vx = zeros<vec>(N);
    vec Y = zeros<vec>(N); vec Vy = zeros<vec>(N);
    X(0) = 1;   Vx(0) = 0; // Initial conditions for x-position and velocity
    Y(0) = 0;   Vy(0) = 10; // Initial conditions for y-position and vecolity


    for (int i=0;i<N-1;i++){
        vec r = zeros<vec>(2); vec v = zeros<vec>(2);
        r << X(i) << Y(i); v << Vx(i) << Vy(i);
        double distance = sqrt(r(0)*r(0) + r(1)*r(1));
        RK4(&r,&v,distance, delta_t);
        X(i+1) = r(0); Y(i+1) = r(1); Vx(i+1) = v(0); Vy(i+1) = v(1);
    }



    // Writing to file. This is working ok.
    ofstream myfile;
    myfile.open("earth_sun.txt");
    for (int i=0;i<N;i++){
        myfile << X(i) << " " << Y(i) << " " << i*delta_t << endl;
    }
    myfile.close();

}

void RK4(vec *r, vec *v, double distance, double dt){
    vec k1 = zeros<vec>(2); vec k2 = zeros<vec>(2); vec k3 = zeros<vec>(2); vec k4 = zeros<vec>(2);

    k1 = dt * f(*r, distance);
    k2 = dt * f(*r + k1/2, distance);
    k3 = dt * f(*r + k2/2, distance);
    k4 = dt * f(*r + k3, distance);
    *v = *v + (k1 + 2*k2 + 2*k3 + k4) / 6;
    *r = *v * dt;
}

vec f(vec r, double distance){
    return -4*pi*pi / (distance*distance*distance) * r;
}


/*
vec f(vec r, vec Rsun, double m){ // r is the position [x,y]. m is not the moving planets mass, but the interracting object.
    double R = Dist(DIST(Rsun,r));
    return 4*3.14*3.14 / (R*R*R) * DIST(r,Rsun);
}

vec inline DIST(vec r1, vec r2){ return r2 - r1; } // Returning distance between two planets
vec inline Dist(vec r){ return sqrt(r(0)*r(0) + r(1)*r(1)); }
*/
