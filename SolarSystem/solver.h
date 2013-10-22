#ifndef SOLVER_H
#define SOLVER_H
#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;


class Solver
{
    vec *r, *v, r2;
    double pi = 3.14; double dt;

public:
    Solver(vec *r, vec *v, vec r2, double dt);
    void RK4();
    vec f(vec r, double dist);
    double distance(vec r1, vec r2);
};

#endif // SOLVER_H
