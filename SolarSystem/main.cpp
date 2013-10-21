#include <iostream>
#include <armadillo>
#include <cmath>
#include <math.h>
#include <fstream>
#include "solver.h"
#include "constants.h"
#include "planet.h"

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
    vec r = zeros<vec>(2); vec v = zeros<vec>(2); vec r2 = zeros<vec>(2);

/*
    for (int i=0;i<N-1;i++){
        vec r = zeros<vec>(2); vec v = zeros<vec>(2);
        r << X(i) << Y(i); v << Vx(i) << Vy(i);
        Solver earth_sun(&r,&v);
        earth_sun.RK4(&r,&v,distance,delta_t);
        X(i+1) = r(0); Y(i+1) = r(1); Vx(i+1) = v(0); Vy(i+1) = v(1);
    }
*/
    Planet earth(Mearth, N);
    earth.initial_condition(1,0,0,2*pi);
    earth.position(0).print();

    for (int i=0;i<N-1;i++){
        r = earth.position(i); v = earth.velocity(i);
        double distance = Distance(r,r2);
        Solver earth_sun(&r,&v);
        earth_sun.RK4(&r,&v,distance,delta_t);
        earth.new_position(r(0),r(1),i); earth.new_velocity(v(0),v(1),i);
    }

/*
    // Writing to file. This is working ok.
    ofstream myfile;
    myfile.open("earth_sun.txt");
    for (int i=0;i<N-2;i++){
        r = earth.position(i);
        myfile << r(0) << " " << r(1) << " " << i*delta_t << endl;
    }
    myfile.close();
*/
}


double Distance(vec r1, vec r2){
    vec R = r1-r2;
    return sqrt(R(0)*R(0) + R(1)*R(1));
}

