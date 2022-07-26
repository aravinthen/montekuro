// Program Name: system.h
// Author: Aravinthen Rajkumar
// Description: main code for monte carlo simulation of a Lennard-Jones fluid.

#include "lib.h"

struct type_matrix init_types(int num_types){
	// Inititalises a type matrix. There are actually two matrices
	// within this structure, but both must be square matrices of 
	// dimension "num_types".
	// The information on how types interact with each other is
	// given by the indices of each matrix.
	struct type_matrix types;
	types.num_types = num_types;

	// initialise an array of pointers to double pointers with [num_types] elements
	types.epsilon = malloc(sizeof(double * )*num_types);
	types.sigma = malloc(sizeof(double * )*num_types);

	// populate the array of pointers to double pointers with arrays of double pointers
	for (int i=0;i<num_types;i++){
	  types.epsilon[i] = malloc(sizeof(double)*num_types);
	  types.sigma[i] = malloc(sizeof(double)*num_types);
	}

	for (int i=0;i<num_types; i++)
	  for (int j=0; j<num_types; j++)
            types.epsilon[i][j] = 0.0;

	return types;
};

// The functions below act on a type matrix object through which data
// about the Lennard-Jones parameters will be obtained.
// Within these functions, the type matrix itself is always referred 
// to as "types". It may be called something else in a program.

void set_sigma( double * sigmatrix[] , int i, int j, double val){
   sigmatrix[i][j] = val;
   sigmatrix[j][i] = val;
}

void set_epsilon( double * epsmatrix[] , int i, int j, double val){
  epsmatrix[i][j] = val;
  epsmatrix[j][i] = val;
}

void free_types(struct type_matrix types){
  free(types.epsilon);
  free(types.sigma);
}

// ---------------------------------------------------------------------

void add_particle(struct particle system[],
		  int label, int type, 
		  double xrange, double yrange, double zrange){
  
  // creates a random configuration of particles
  system[label].id = label;
  system[label].species = type;
  system[label].x = xrange*(double)rand()/RAND_MAX;
  system[label].y = yrange*(double)rand()/RAND_MAX;
  system[label].z = zrange*(double)rand()/RAND_MAX;
}

void copy_system(struct particle copy[],
		 struct particle original[],
		 int total){

  // copies the particles in a system
  
  for (int i=0; i<total; i++){
    copy[i].id = i;
    copy[i].species = original[i].species;
    copy[i].x = original[i].x;
    copy[i].y = original[i].y;
    copy[i].z = original[i].z;
  }
}

