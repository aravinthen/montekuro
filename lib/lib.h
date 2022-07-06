// Program name: lib.h
// Author: Aravinthen Rajkumar
// Description: Header file for library functions

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define pi 3.14159

// --------------------------------------------------------------------------
struct type_matrix {
	// epsilon and sigma are pointers to a region of memory for doubles.
	// they will be allocated memory using malloc()
	int num_types;
	double ** epsilon;
	double ** sigma;
};

struct particle {
  // this is the basic particle structure
  int id;   // the "id" of the particle
  int species; 
  double x;  // x coordinate 
  double y;  // y coordinate
  double z;  // z coordinate
};

// --------------------------------------------------------------------------
struct type_matrix init_types(int num_types);

double particle_move(struct particle system[],
		     struct type_matrix types,
		     int label,
		     double magnitude,
		     int nump,
		     double cutoff,
		     double xrange, double yrange, double zrange);

double system_energy(struct particle system[],
		     struct type_matrix types,
		     int num_parts,
		     double cutoff,
		     double xrange, double yrange, double zrange);

void pinfo(struct particle system[],
	   int num);

void set_sigma(double * sigmatrix[] , int i, int j, double val);
void set_epsilon(double * epsmatrix[] , int i, int j, double val);
void free_types(struct type_matrix types);

void add_particle(struct particle system[],
		  int label, int type, 
		  double xrange, double yrange, double zrange);

void copy_system(struct particle system1[],
		 struct particle system2[],
		 int total);
