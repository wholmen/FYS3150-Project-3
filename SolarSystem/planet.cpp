#include "planet.h"

Planet::Planet(double mass, int N)
{
    Planet::N = N; Planet::mass = mass;
}

void Planet::initial_condition(double x, double y, double vx, double vy)
{
    Planet::X(0) = x; Planet::Y(0) = y; Planet::Vx(0) = vx; Planet::Vy(0) = vy;
}

void Planet::new_position(double x, double y, int i)
{
    Planet::X(i+1) = x; Planet::Y(i+1) = y;
}

void Planet::new_velocity(double vx, double vy, int i)
{
    Planet::Vx(i+1) = vx; Planet::Vy(i+1) = vy;
}

vec Planet::position(int i)
{
    Planet::r(0) = Planet::X(i); Planet::r(1) = Planet::Y(i);
    return Planet::r;
}

vec Planet::velocity(int i)
{
    Planet::v(0) = Planet::Vx(i); Planet::v(1) = Planet::Vy(i);
    return Planet::v;
}
