#include "solver.h"

Solver:: Solver(double GCONST, int P){
    Gconst = GCONST; p = P;
}

void Solver::RK4(vec *r, vec *v, mat pos, vec mas, double dt){
    vec k1 = zeros<vec>(2); vec k2 = zeros<vec>(2); vec k3 = zeros<vec>(2); vec k4 = zeros<vec>(2);

    k1 = dt * f(*r          ,pos,mas);
    k2 = dt * f(*r + dt*k1/2,pos,mas);
    k3 = dt * f(*r + dt*k2/2,pos,mas);
    k4 = dt * f(*r + dt*k3  ,pos,mas);
    *v = *v + (k1 + 2*k2 + 2*k3 + k4) / 6;
    *r = *r + *v * dt;
}

vec Solver::f(vec r, mat pos, vec mas){
    vec a = zeros<vec>(2);
    int size = 10;
    for (int pl=0;pl<size;pl++){
        if (pl!=p){
            vec r2 = zeros<vec>(2); r2 = pos.col(pl);
            a = a - mas(pl) / (distance(r,r2)*distance(r,r2)*distance(r,r2)) * vector(r,r2);
        }
    }
    return a*Gconst;
}

double Solver::distance(vec r1, vec r2){
    vec R = r1 - r2;
    return sqrt(R(0)*R(0) + R(1)*R(1));
}

vec Solver::vector(vec r1, vec r2){
    return r2 - r1;
}
