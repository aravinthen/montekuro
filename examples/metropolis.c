// File name: main.c
// Description: The main file for running the MC algorithm. 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "../lib/system.h"

int main(void){
	int num_parts = 100; // the number of particles
	int num_types = 3; // the number of types 
	
	// the dimensions of the problem
	double xdim = 1.0;
	double ydim = 1.0;
	double zdim = 1.0;

	// initialize the type matrix.
	struct type_matrix type_data = init_types(num_types);
	define_sigma(&type_data, 1, 1, 1.0);

	// the system: all particle information is contained here. 
	struct particle system[num_parts];	

	// add particles to the system
	for (int i=0; i++; i<num_parts){
		add_particle(system, i, 1, xdim, ydim, zdim);
	}

	// free the type matrices.
	free_types(type_data);
}
