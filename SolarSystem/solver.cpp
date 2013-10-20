#include "solver.h"

Solver:: Solver(vec *r, vec *v){
}

void Solver::RK4(vec *r, vec *v, double distance, double dt){
    vec k1 = zeros<vec>(2); vec k2 = zeros<vec>(2); vec k3 = zeros<vec>(2); vec k4 = zeros<vec>(2);
    k1 = dt * f(*r,           distance);
    k2 = dt * f(*r + dt*k1/2, distance);
    k3 = dt * f(*r + dt*k2/2, distance);
    k4 = dt * f(*r + dt*k3,   distance);
    *v = *v + (k1 + 2*k2 + 2*k3 + k4) / 6;
    *r = *r + *v * dt;
}

vec Solver::f(vec r, double d){
    return -4*pi*pi / (d*d*d) * r;
}

