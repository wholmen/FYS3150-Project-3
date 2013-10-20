#include <iostream>
#include "constants.h"
#include <armadillo>
#include <cmath>
#include <math.h>
#include <fstream>

using namespace std;
using namespace arma;

void RK4(vec *r, vec *v, double distance, double dt);
vec f(vec r, double dist);
double Distance(vec r1, vec r2);

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
    double Tfinal = 10; //years
    double delta_t = 0.1; //years
    int N = Tfinal/delta_t;

    // Initializing earth-sun
    vec X = zeros<vec>(N); vec Vx = zeros<vec>(N);
    vec Y = zeros<vec>(N); vec Vy = zeros<vec>(N);
    X(0) = 1;   Vx(0) = 0; // Initial conditions for x-position and velocity
    Y(0) = 0;   Vy(0) = 2*pi; // Initial conditions for y-position and vecolity


    for (int i=0;i<N-1;i++){
        vec r = zeros<vec>(2); vec v = zeros<vec>(2); vec r2 = zeros<vec>(2);
        r << X(i) << Y(i); v << Vx(i) << Vy(i);
        double distance = Distance(r,r2);
        RK4(&r,&v,distance,delta_t);
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
    k1 = dt * f(*r,           distance);
    k2 = dt * f(*r + dt*k1/2, distance);
    k3 = dt * f(*r + dt*k2/2, distance);
    k4 = dt * f(*r + dt*k3,   distance);
    *v = *v + (k1 + 2*k2 + 2*k3 + k4) / 6;
    *r = *r + *v * dt;
}

vec f(vec r, double d){
    return -4*pi*pi / (d*d*d) * r;
}

double Distance(vec r1, vec r2){
    vec R = r1-r2;
    return sqrt(R(0)*R(0) + R(1)*R(1));
}

