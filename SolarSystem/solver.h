#ifndef SOLVER_H
#define SOLVER_H
#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;


class Solver
{
    double Gconst; int p;

public:
    Solver(double GCONST, int P);
    void RK4(vec *r, vec *v, mat pos, vec mas, double dt);
    vec f(vec r, mat pos, vec mas);
    double distance(vec r1, vec r2);
    vec vector(vec r1, vec r2);
};

#endif // SOLVER_H
