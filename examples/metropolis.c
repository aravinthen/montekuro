// File name: main.c
// Description: The main file for running the MC algorithm. 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "../lib/lib.h"

int main(void){
  srand(time(0)); // generate seed (this has to be done once per program)
  
  int num_parts = 5; // the number of particles per type
  int num_types = 2; // the number of types

  // here, we're modelling a binary fluid
  int total = 2*num_parts; // the number of particles in the system
  int num_moves = 100;

  // the dimensions of the problem
  double temp = 1.0;
  double xdim = 1.0;
  double ydim = 1.0;
  double zdim = 1.0;  

  // parameters
  double cutoff=1.5;
  double mvd = 0.1;

  // Global parameters
  int p;
  double energy, en_change, rnum, weight;
  
  // initialize the type matrix.

  struct type_matrix type_data = init_types(num_types);  
  set_sigma(type_data.sigma, 1,1, 1.0);
  set_sigma(type_data.sigma, 0,0, 1.0);
  set_sigma(type_data.sigma, 1,0, 1.0);

  set_epsilon(type_data.epsilon, 1,1, 1.0);
  set_epsilon(type_data.epsilon, 0,0, 0.5);
  set_epsilon(type_data.epsilon, 1,0, 0.2);

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


  double test_en;
  // calculate the initial energy
  energy = system_energy(sys,
			 type_data,
			 total,
			 cutoff,
			 xdim, ydim, zdim);
  
  // assign the system particles to the test system
  for (int k=0; k<num_moves; k++){
    // pick a random particle
    p = prand(total);
    
    // assign the system in use to the test system
    copy_system(test, sys, total);

    // move only the test system!
    //    en_change = particle_move(test,
    particle_move(test,
		  type_data,
		  p,
		  mvd,
		  total,
		  cutoff,
		  xdim, ydim, zdim);
    
    rnum = urand();
    test_en = system_energy(test,
			    type_data,
			    total,
			    cutoff,
			    xdim, ydim, zdim);
    
    weight = exp(-(test_en - energy)/(kb*temp));

    if (weight > rnum){ 
      copy_system(sys, test, total);
      energy = test_en;
    }
    printf("%f\n", energy);
  };
  
  printf("Success!\n");
  
  // free the type matrices.
  free_types(type_data);
}
