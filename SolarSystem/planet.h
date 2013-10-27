#ifndef PLANET_H
#define PLANET_H
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

class Planet
{
    double mass; int N;
    vec X, Y, Vx, Vy, Ep, Ek, Etot;

public:
    Planet(double m, int n);
    void initialize(int N);
    void initial_condition(double x, double y, double vx, double vy);
    vec position(int i);
    vec velocity(int i);
    void new_position(vec r, int i);
    void new_velocity(vec v, int i);
    double MASS();
    void energy(double EP,int i);
    void Energytot();
    double Etotal(int i);
};

#endif // PLANET_H
