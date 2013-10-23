#ifndef PLANET_H
#define PLANET_H
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

class Planet
{
    double mass; int N;
    vec r, v, X, Y, Vx, Vy;

public:
    Planet(double m, int n);
    void initialize(int N);
    void initial_condition(double x, double y, double vx, double vy);
    vec position(int i);
    vec velocity(int i);
    void new_position(vec r, int i);
    void new_velocity(vec v, int i);
    double MASS();
};

#endif // PLANET_H
