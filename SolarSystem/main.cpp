#include <iostream>
#include "constants.h"
#include <armadillo>
#include <cmath>

using namespace std;
using namespace arma;

void RK4(double *x, double *y);
double Force(double x, double y);

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
    double delta_t = 0.1; //years
    int N = Tfinal/delta_t;
    vec X = zeros<vec>(N,1); vec Vx = zeros<vec>(N,1);
    vec Y = zeros<vec>(N,1); vec Vy = zeros<vec>(N,1);
    X(0) = Rearth;   Vx(0) = 0; // Initial conditions for x-position and velocity
    Y(0) = 0;        Vy(0) = sqrt(Gconst * Msun / (X(0)**2+Y(0)**2)); // Initial conditions for y-position and vecolity


    RK4(&X(0),&Y(0));


}

void RK4(double *x, double *y){
    *x = *x + 2;
    *y = *y *3;



}

double Force(double x, double y){
    double F = -Gconst * Mearth * Msun / (x*x + y*y);
    return F;
}
