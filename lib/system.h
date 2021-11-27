// Program Name: main.c
// Author: Aravinthen Rajkumar
// Description: main program for monte carlo simulation of a Lennard-Jones fluid.
// I'm assuming that the mass of the atoms in the system are unity in Lennard-Jones units. 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define pi 3.14159

// The functions below act on a type matrix object through which data
// about the Lennard-Jones parameters will be obtained.
// Within these functions, the type matrix itself is always referred 
// to as "types". It may be called something else in a program.

struct type_matrix {
	// epsilon and sigma are pointers to a region of memory for doubles.
	// they will be allocated memory using malloc()
	int num_types;
	double * epsilon;
	double * sigma;
};

struct type_matrix init_types(int num_types){
	// Inititalises a type matrix. There are actually two matrices
	// within this structure, but both must be square matrices of 
	// dimension "num_types".
	// The information on how types interact with each other is
	// given by the indices of each matrix.
	struct type_matrix types;
	types.num_types = num_types;
	types.epsilon = malloc(sizeof(double)*num_types*num_types);
	types.sigma = malloc(sizeof(double)*num_types*num_types);
	
	return types;
}

void define_sigma(struct type_matrix * types, 
									int type1, int type2, double sigma){
}

void free_types(struct type_matrix types){
	free(types.epsilon);
	free(types.sigma);
}

// ---------------------------------------------------------------------

struct particle {
  // this is the basic particle structure
  int id;   // the "id" of the particle
	int species; 
  double x;  // x coordinate 
  double y;  // y coordinate
  double z;  // z coordinate
};

void add_particle(struct particle system[],
									int label, int type, 
									double xrange, double yrange, double zrange){

  // creates a random configuration of particles
  // random number generator
	system[label].id = label;
	system[label].species = type;
	system[label].x = xrange*(double)rand()/RAND_MAX;
	system[label].y = yrange*(double)rand()/RAND_MAX;
	system[label].z = zrange*(double)rand()/RAND_MAX;
}

