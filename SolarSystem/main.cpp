#include <iostream>
#include <armadillo>
#include <cmath>
#include <math.h>
#include <fstream>
#include "solver.h"
#include "constants.h"
#include "planet.h"

using namespace std;
using namespace arma;

void RK4(vec *r, vec *v, double distance, double dt);
vec f(vec r, double dist);
double Distance(vec r1, vec r2);

int main()
{
    double Tfinal = 10; //years
    double dt = 0.01; //years
    int N = Tfinal/dt;

    // Initializing earth-sun
    vec r_earth = zeros<vec>(2); vec v_earth = zeros<vec>(2); vec r_sun = zeros<vec>(2); vec v_sun;

    Planet earth(Mearth, N); earth.initial_condition(1,0,0,2*pi);
    Planet sun(Msun, N); sun.initial_condition(0,0,0,0);


    Planet **planets = new Planet*[2];
    int p = 0;
    //                      Mass                                        x y vx vy
    planets[p] = new Planet(Mearth, N); planets[p] -> initial_condition(1,0,0,2*pi); p++;
    planets[p] = new Planet(Msun, N);   planets[p] -> initial_condition(0,0,0,2*pi); p++;
    vec pos[p]; vec vel[p];
    vec force_r[p]; double force_m[p];

    for (int i=0; i<N-1; i++){

        // Extracting all positions and velocities from class Planet.
        for (int j=0;j<p;j++){
            pos[j] = planets[j] -> position(j); vel[j] = planets[j] -> velocity(j);
        }

        for (int j=0;j<p;j++){
            //Creating acceleration function for planet: p
            for (int e=o;e<p;e++){
                if (e != j){
                    force_r[e] = pos[e];
                    force_m[e] = planets[e] -> MASS();
                }}
            Solver eval(&pos[j],&vel[j],force_e,force_m,dt); // must create a function creator inside Solver based on force_e and force_m. These together can make accelleration.
            eval.RK4();

            // Must insert new pos[j] into planet-class!
        }



    for (int i=0;i<N-1;i++){

        r_earth = earth.position(i); v_earth = earth.velocity(i);
        r_sun = sun.position(i); v_sun = sun.velocity(i);

        Solver earth_sun(&r_earth,&v_earth,r_sun,dt);// Solver sun_earth(&r_sun,&v_sun);
        earth_sun.RK4();
        //sun_earth.RK4(&r_sun,&v_sun,Distance(r_earth,r_sun), dt);

        earth.new_position(r_earth,i); earth.new_velocity(v_earth,i);
        sun.new_position(r_sun,i); sun.new_velocity(v_sun, i);
    }


    // Writing to file. This is working ok.
    ofstream myfile; ofstream myfile2;
    myfile.open("earth_sun.txt");
    myfile2.open("sun_earth.txt");
    for (int i=0;i<N;i++){
        myfile << earth.position(i)(0) << " " << earth.position(i)(1) << " " << i*dt << endl;
        myfile2 << sun.position(i)(0) << " " << sun.position(i)(1) << " " << i*dt << endl;
    }
    myfile.close();

}


double Distance(vec r1, vec r2){
    vec R = r1-r2;
    return sqrt(R(0)*R(0) + R(1)*R(1));
}

