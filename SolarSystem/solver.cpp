#include "solver.h"

Solver:: Solver(vec *r, vec *v, vec r2, double dt){
    *r = *r; *v = *v; r2 = r2; dt = dt;
    RK4();
}

void Solver::RK4(){
    vec k1 = zeros<vec>(2); vec k2 = zeros<vec>(2); vec k3 = zeros<vec>(2); vec k4 = zeros<vec>(2);
    k1 = dt * f(*r,           distance(*r,r2));
    k2 = dt * f(*r + dt*k1/2, distance(*r,r2));
    k3 = dt * f(*r + dt*k2/2, distance(*r,r2));
    k4 = dt * f(*r + dt*k3,   distance(*r,r2));
    *v = *v + (k1 + 2*k2 + 2*k3 + k4) / 6;
    *r = *r + *v * dt;
}

vec Solver::f(vec r, double d){
    return -4*pi*pi / (d*d*d) * r;
}

double Solver::distance(vec r1, vec r2){
    vec R = r1 - r2;
    return sqrt(R(0)*R(0) + R(1)*R(1));
}
