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

int main()
{
    double Tfinal = 1000; //years
    double dt = 0.01; //years
    int N = Tfinal/dt; //iteration points

    // Feeded planet-data into the planet class.
    Planet **planets = new Planet*[10]; int p = 0;
    //                      Mass                                               x        y      vx        vy
    planets[p] = new Planet(Msun,     N); planets[p] -> initial_condition(0        ,0      ,0       ,auy*13.1*Mjupiter/Msun); p++;
    planets[p] = new Planet(Mmercury, N); planets[p] -> initial_condition(-Rmercury,0      ,0       ,-auy*47.9); p++;
    planets[p] = new Planet(Mvenus,   N); planets[p] -> initial_condition(0        ,Rvenus ,auy*35.0,0        ); p++;
    planets[p] = new Planet(Mearth,   N); planets[p] -> initial_condition(-Rearth  ,0      ,0       ,auy*29.8 ); p++;
    planets[p] = new Planet(Mmars,    N); planets[p] -> initial_condition(0        ,Rmars  ,auy*24.1,0        ); p++;
    planets[p] = new Planet(Mjupiter, N); planets[p] -> initial_condition(-Rjupiter,0      ,0       ,-auy*13.1); p++;
    planets[p] = new Planet(Msaturn,  N); planets[p] -> initial_condition(Rsaturn  ,0      ,0       ,-9.7*auy ); p++;
    planets[p] = new Planet(Muranus,  N); planets[p] -> initial_condition(Ruranus  ,0      ,0       ,6.8*auy  ); p++;
    planets[p] = new Planet(Mneptun,  N); planets[p] -> initial_condition(0        ,Rneptun,-5.4*auy,0        ); p++;
    planets[p] = new Planet(Mpluto,   N); planets[p] -> initial_condition(-Rpluto  ,0      ,0       ,auy*4.7  ); p++;

    // Creating one solver instance for each planet.
    Solver **solvers = new Solver*[p];
    for (int pl=0;pl<p;pl++){
        solvers[pl] = new Solver(Gconst,pl);
    }

    // Starting the loop over each timestep i.
    for (int i=0; i<N-1; i++){

        // Retrieving planet-data from respective classes and inserting the info into matrixes.
        mat pos = zeros<mat>(2,p); mat vel=zeros<mat>(2,p); vec mas = zeros<vec>(p);
        for (int j=0;j<p;j++){
            pos.col(j) += planets[j]->position(i);      //(2xP) - matrix with each column representing coordinates for respective planet.
            vel.col(j) += planets[j]->velocity(i);      //(2xP) - matrix with each column representing velocity for respective planet.
            mas(j) = planets[j] -> MASS();              //(1xP) - vector where each element is representing mass for respective planet.
        }

        for (int pl=0;pl<p;pl++){                       // Looping over each planet.
            vec r = zeros<vec>(2); r = pos.col(pl);     // Coordinates for planet pl.
            vec v = zeros<vec>(2); v = vel.col(pl);     // Velocity for planet pl.
            solvers[pl]->RK4(&r,&v,pos,mas,dt);         // RK4 calculates new position and velocity for planet pl.

            planets[pl]->energy(solvers[pl]->pot(),i);
            planets[pl]->new_position(r,i); planets[pl]->new_velocity(v,i); // Feed new position/velocity into respective planet class.
        } // end of planet loop

    } // end of time loop

    for (int pl=0;pl<p;pl++){planets[pl]->Energytot();} // Creating Etot inside each instance of planets.
    vec L = zeros<vec>(N);
    for (int i=0;i<N;i++){
        double r = sqrt(planets[3]->position(i)(0)*planets[3]->position(i)(0) + planets[3]->position(i)(1)*planets[3]->position(i)(1));
        double v = sqrt(planets[3]->velocity(i)(0)*planets[3]->velocity(i)(0) + planets[3]->velocity(i)(1)*planets[3]->velocity(i)(1));
        L(i) = r*v;
    }


    // -----------------------------------------------------------------------
    // ------------------------ Writing to files -----------------------------
    // -----------------------------------------------------------------------

    ofstream myfile0; ofstream myfile1; ofstream myfile2; ofstream myfile3; ofstream myfile4;
    ofstream myfile5; ofstream myfile6; ofstream myfile7; ofstream myfile8; ofstream myfile9;
    myfile0.open("sun.txt"); myfile1.open("mercury.txt"); myfile2.open("venus.txt"); myfile3.open("earth.txt");
    myfile4.open("mars.txt"); myfile5.open("jupiter.txt"); myfile6.open("saturn.txt"); myfile7.open("uranus.txt");
    myfile8.open("neptun.txt"); myfile9.open("pluto.txt");

    for (int i=0;i<N;i++){
        myfile0  << planets[0] ->position(i)(0) << " " << planets[0] ->position(i)(1) << " " << planets[0]->Etotal(i) << " " << i*dt << endl;
        myfile1  << planets[1] ->position(i)(0) << " " << planets[1] ->position(i)(1) << " " << planets[1]->Etotal(i) << " " << i*dt << endl;
        myfile2  << planets[2] ->position(i)(0) << " " << planets[2] ->position(i)(1) << " " << planets[2]->Etotal(i) << " " << i*dt << endl;
        myfile3  << planets[3] ->position(i)(0) << " " << planets[3] ->position(i)(1) << " " << planets[3]->Etotal(i) << " " << i*dt << " " << L(i) << endl;
        myfile4  << planets[4] ->position(i)(0) << " " << planets[4] ->position(i)(1) << " " << planets[4]->Etotal(i) << " " << i*dt << endl;
        myfile5  << planets[5] ->position(i)(0) << " " << planets[5] ->position(i)(1) << " " << planets[5]->Etotal(i) << " " << i*dt << endl;
        myfile6  << planets[6] ->position(i)(0) << " " << planets[6] ->position(i)(1) << " " << planets[6]->Etotal(i) << " " << i*dt << endl;
        myfile7  << planets[7] ->position(i)(0) << " " << planets[7] ->position(i)(1) << " " << planets[7]->Etotal(i) << " " << i*dt << endl;
        myfile8  << planets[8] ->position(i)(0) << " " << planets[8] ->position(i)(1) << " " << planets[8]->Etotal(i) << " " << i*dt << endl;
        myfile9  << planets[9] ->position(i)(0) << " " << planets[9] ->position(i)(1) << " " << planets[9]->Etotal(i) << " " << i*dt << endl;
    }
    myfile0.close(); myfile1.close(); myfile2.close(); myfile3.close(); myfile4.close();
    myfile5.close(); myfile6.close(); myfile7.close(); myfile8.close(); myfile9.close();

    return 0;
}
