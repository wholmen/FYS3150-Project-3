#include "planet.h"

Planet::Planet(double m, int n)
{
    N = n; mass = m;
    Planet::initialize(N);
}

void Planet::initialize(int N){
    X = zeros<vec>(N); Y = zeros<vec>(N);
    Vx = zeros<vec>(N); Vy = zeros<vec>(N);
    Ek = zeros<vec>(N); Ep = zeros<vec>(N); Etot = zeros<vec>(N);
}

void Planet::initial_condition(double x, double y, double vx, double vy)
{
    X(0) = x; Y(0) = y; Vx(0) = vx; Vy(0) = vy;
}

void Planet::new_position(vec r, int i)
{
    X(i+1) = r(0); Y(i+1) = r(1);
}

void Planet::new_velocity(vec v, int i)
{
    Vx(i+1) = v(0); Vy(i+1) = v(1);
}

vec Planet::position(int i)
{
    vec r = zeros<vec>(2);
    r(0) = X(i); r(1) = Y(i);
    return r;
}

vec Planet::velocity(int i)
{
    vec v = zeros<vec>(2);
    v(0) = Vx(i); v(1) = Vy(i);
    return v;
}

double Planet::MASS(){
    return mass;
}

void Planet::energy(double EP,int i){
    Ep(i) = EP;
    Ek(i) = (Vx(i)*Vx(i) + Vy(i)*Vy(i)) * 0.5 * mass;
}

void Planet::Energytot(){
    Etot = Ek + Ep;
}

double Planet::Etotal(int i){
    return Etot(i);
}


