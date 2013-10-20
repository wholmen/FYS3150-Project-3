#ifndef SOLVER_H
#define SOLVER_H
#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;


class Solver
{

    double mass; double T;
    vec *r, *v;
    double pi = 3.14;

public:
    Solver(vec *r, vec *v);
    void RK4(vec *r, vec *v, double distance, double dt);
    vec f(vec r, double dist);
};

#endif // SOLVER_H
