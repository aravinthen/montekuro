// Program name: lib.h
// Author: Aravinthen Rajkumar
// Description: Header file for library functions

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define pi 3.14159
#define kb 1

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

struct energy_data {
  int num_moves;
  int every;
  int after;
  int point;
  double * history;
};

// --------------------------------------------------------------------------
struct type_matrix init_types(int num_types);
struct energy_data energy_his(int num_moves, int every, int after);

double particle_move(struct particle system[],
		     struct type_matrix types,
		     int label,
		     double magnitude,
		     int nump,
		     double cutoff,
		     double xrange, double yrange, double zrange);

double urand();
int prand(int total);

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
void free_energy(struct energy_data energy);

void add_particle(struct particle system[],
		  int label, int type, 
		  double xrange, double yrange, double zrange);

void copy_system(struct particle system1[],
		 struct particle system2[],
		 int total);

void gnupos(FILE*fp, struct particle system[], int num, char fname[]);

void logfile(FILE*info,
	     struct particle system[],
	     int accepted,
	     int moves,
	     int total,
	     int ntypes,
	     double xlen,
	     double ylen,
	     double zlen,
	     double T,
	     double step);

int type_count(struct particle system[], int total, int type);

void add_en(struct energy_data * edata, int k, double energy);

void dump_en(struct energy_data * edata,
	     FILE*en,
	     char fname[]);

int sanity_check(struct particle system[],
		 double final_en,
		 double tol,
		 struct type_matrix type_data,
		 int total,
		 double cutoff,
		 double xdim, double ydim, double zdim);

