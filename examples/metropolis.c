// File name: main.c
// Description: The main file for running the MC algorithm. 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "../lib/lib.h"

int main(void){
  
  srand(time(0)); // generate seed (this has to be done once per program)
  
  int num_parts = 100; // the number of particles per type
  int num_types = 2; // the number of types

  // here, we're modelling a binary fluid
  int total = 2*num_parts; // the number of particles in the system

  // the dimensions of the problem
  double temp = 1.00;
  double xdim = 10.0;
  double ydim = xdim;
  double zdim = xdim;  

  // parameters
  double cutoff = 0.5*xdim;
  double mvd = 0.0005;
  double tol = 5e-3;

  // Global parameters
  int p; // random number for particles
  double energy, en_change, rnum, weight;

  // Monte Carlo statistics
  int accepted=0;

  // File pointers
  FILE *fp; // position files
  FILE *ep; // energy file
  FILE *info; // information file

  // move details
  int num_moves = 1000000;
  int every = 1;
  int after = 1000;

  // ------------------------------------------------------------------------
  // MONTE CARLO SIMULATION
  // ------------------------------------------------------------------------
  
  // initialize the type matrix.
  struct type_matrix type_data = init_types(num_types);
  // initialize the energy data vector
  struct energy_data ehis = energy_his(num_moves, every, after);  
  
  // initialize the energy storage
  //  double en_data[(num_moves-after)/every] = {0.0};  
  
  set_sigma(type_data.sigma, 1,1, 1.0);
  set_sigma(type_data.sigma, 0,0, 1.0);
  set_sigma(type_data.sigma, 1,0, 1.0);

  set_epsilon(type_data.epsilon, 1, 1, 1.0);
  set_epsilon(type_data.epsilon, 0, 0, 1.0);
  set_epsilon(type_data.epsilon, 1, 0, 1.0);

  // the system: all particle information is contained here.
  struct particle sys[total];
  // the test state that is used in energy evaluations.
  struct particle test[total]; 

  // add particles to the system
  for (int i=0; i<num_parts; i++){
    add_particle(sys, i, 0, xdim, ydim, zdim);
  }
  for (int j=num_parts; j<total; j++){
    add_particle(sys, j, 1, xdim, ydim, zdim);
  }

  // pinfo(sys, total);
  gnupos(fp, sys, total, "initial");

  double test_en;
  // calculate the initial energy
  energy = system_energy(sys,
			 type_data,
			 total,
			 cutoff,
			 xdim, ydim, zdim);

  printf("Initial energy: %f\n", energy);
  
  // assign the system particles to the test system
  for (int k=0; k<=num_moves; k++){
    // pick a random particle
    p = prand(total);
    
    // assign the system in use to the test system
    copy_system(test, sys, total);

    // move only the test system!
    //    en_change = particle_move(test,
    en_change = particle_move(test,
			      type_data,
			      p,
			      mvd,
			      total,
			      cutoff,
			      xdim, ydim, zdim);

    rnum = urand();    
    weight = exp(-en_change/(kb*temp));    
    if (weight >= rnum){ 
      copy_system(sys, test, total);
      energy += en_change;
      accepted +=1;
    }

    // automatically adds energy data to history
    add_en(&ehis, k, energy);
    
  }

  // SANITY CHECK
  sanity_check(sys, energy, tol, type_data, total, cutoff, xdim, ydim, zdim);
  
  // print final positions of particles to visualise
  gnupos(fp, sys, total, "final");
  
  logfile(info,
	  sys,
	  accepted,
	  num_moves,
	  total,
	  num_types,
	  xdim,
	  ydim,
	  zdim,
	  temp,
	  mvd);

  dump_en(&ehis, ep, "energy.txt");

  // free the type matrices.
  free_types(type_data);
  free_energy(ehis);
  return 0;
}
  

