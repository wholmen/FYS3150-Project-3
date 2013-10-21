#ifndef PLANET_H
#define PLANET_H
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

class Planet
{
    double mass; int N;
    vec X = zeros<vec>(N); vec Y = zeros<vec>(N);
    vec Vx = zeros<vec>(N); vec Vy = zeros<vec>(N);
    vec r = zeros<vec>(2); vec v = zeros<vec>(2);
public:
    Planet(double mass, int N);
    void initial_condition(double x, double y, double vx, double vy);
    vec position(int i);
    vec velocity(int i);
    void new_position(double x, double y, int i);
    void new_velocity(double vx, double vy, int i);
};

#endif // PLANET_H
