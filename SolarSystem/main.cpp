#include <iostream>
#include "constants.h"
#include <armadillo>
#include <cmath>
#include <math.h>

using namespace std;
using namespace arma;

void RK4(double *x, double *y, double *vx, double *vy, double T, double dt, int N);
double Force(double x, double y,int o);

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

    // Initializing earth-sun
    vec X = zeros<vec>(N); vec Vx = zeros<vec>(N);
    vec Y = zeros<vec>(N); vec Vy = zeros<vec>(N);
    X(0) = Rearth;   Vx(0) = 0; // Initial conditions for x-position and velocity
    Y(0) = 0;        Vy(0) = sqrt(Gconst * Msun) / (X(0)*X(0) + Y(0)*Y(0)); // Initial conditions for y-position and vecolity

    RK4(&X(0),&Y(0),&Vx(0),&Vy(0),Tfinal,delta_t,N);


}

void RK4(double *x, double *y, double *vx, double *vy, double T, double dt, int N){


    // Euler-Chromer
    double ax,ay,Fx,Fy;
    for (int i=0; i<N; i++){
        Fx = Force(x[i],y[i],0);
        Fy = Force(x[i],y[i],1);
        ax = Fx/Mearth;
        ay = Fy/Mearth;
        vx[i+1] = vx[i] + ax*dt;
        vy[i+1] = vy[i] + ay*dt;
        x[i+1] = x[i] + vx[i+1]*dt;
        y[i+1] = y[i] + vy[i+1]*dt;
    }
}

double Force(double x, double y, int o){
    double F;
    if (o==0){ //x-direction
        F = -Gconst * Mearth * Msun / (x*x*x + y*y*y) * x;
    }else if (o==1){ //y-direction
        F = -Gconst * Mearth * Msun / (x*x*x + y*y*y) * y;
    }
    return F;
}
