// Program Name: main.c
// Author: Aravinthen Rajkumar
// Description: main program for monte carlo simulation of a Lennard-Jones fluid.
// I'm assuming that the mass of the atoms in the system are unity in Lennard-Jones units. 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

struct particle {
  // this is the basic particle structure
  int id;   // the "id" of the particle
  double x;  // x coordinate 
  double y;  // y coordinate
  double z;  // z coordinate
};


// Main algorithm functions ----------------------------------------------------------------------
void initialise(struct particle system[], int num, double xrange, double yrange, double zrange){
  // creates a random configuration of particles
  // random number generator
  for(int i=0;i<num;i++){
    system[i].id = i;
    system[i].x = xrange*(double)rand()/RAND_MAX;
    system[i].y = yrange*(double)rand()/RAND_MAX;
    system[i].z = zrange*(double)rand()/RAND_MAX;
  }
}


long double total_energy(struct particle system[],double sigma, double epsilon,
			 int num, double xrange, double yrange, double zrange,  double cutoff){
  // calculates the potential energy of a single particle.
  // this is done by summing over all the other atoms in the system within a cutoff.
  // The potential used here is the Lennard-Jones cut-off potential.
  // The cutoff demands the use of periodic boundary conditions.

  double dx, dy, dz;
  double distance;
  double energy = 0.0;

  for (int i=0; i<num; i++){
    for (int j=i+1; j<num; j++){
      if (i==j) continue; // ignores self-interaction
      else {
	// minimum image convention for periodic systems
	dx = system[i].x - system[j].x;
	dy = system[i].y - system[j].y;
	dz = system[i].z - system[j].z;

	// correct x for periodicity
	if (dx > 0.5*xrange) dx= dx-xrange;
	if (dx <= -0.5*xrange) dx=dx+xrange;
	// correct y for periodicity
	if (dy > 0.5*yrange) dy= dy-yrange;
	if (dy <= -0.5*yrange) dy=dy+yrange;
	// correct z for periodicity
	if (dz > 0.5*zrange) dz=dz-zrange;
	if (dz <= -0.5*zrange) dz=dz+zrange;

	distance = sqrt(dx*dx + dy*dy + dz*dz);
      
	if (distance <= cutoff){
	  energy += 4*epsilon*(pow(sigma/distance, 12) - pow(sigma/distance, 6));
	}
      }    
    }
  }
  return energy;
  
}

long double particle_energy(struct particle system[],struct particle p, double sigma, double epsilon,
		      int num, double xrange, double yrange, double zrange,  double cutoff){
  // calculates the potential energy of a single particle.
  // this is done by summing over all the other atoms in the system within a cutoff.
  // The potential used here is the Lennard-Jones cut-off potential.
  // The cutoff demands the use of periodic boundary conditions.

  int id = p.id;
  double dx, dy, dz;
  double distance;
  double energy = 0.0;

  for (int i=0; i<num; i++){
    if (i==id) continue; // ignores the particle itself
    
    else {
      // minimum image convention for periodic systems
      dx = p.x - system[i].x;
      dy = p.y - system[i].y;
      dz = p.z - system[i].z;

      // correct x for periodicity
      if (dx > 0.5*xrange) dx= dx-xrange;
      if (dx <= -0.5*xrange) dx=dx+xrange;
      // correct y for periodicity
      if (dy > 0.5*yrange) dy= dy-yrange;
      if (dy <= -0.5*yrange) dy=dy+yrange;
      // correct z for periodicity
      if (dz > 0.5*zrange) dz=dz-zrange;
      if (dz <= -0.5*zrange) dz=dz+zrange;

      distance = sqrt(dx*dx + dy*dy + dz*dz);
      
      if (distance <= cutoff){
	energy += 4*epsilon*(pow(sigma/distance, 12) - pow(sigma/distance, 6));
      }
    }    
  }

  return energy;
}

long double trial_energy(struct particle system[], int index, double trialx, double trialy, double trialz,
		    double sigma, double epsilon, double cutoff,
		    int num, double xrange, double yrange, double zrange){
  // calculates the potential energy of a trial position.
  // an index must be supplied so that the position ignores the particle that 

  double dx, dy, dz;
  double distance;
  double energy = 0.0;

  for (int i=0; i<num; i++){
    if (i==index) continue; // ignores the particle itself
    
    else {
      // minimum image convention for periodic systems
      dx = trialx - system[i].x;
      dy = trialy - system[i].y;
      dz = trialz - system[i].z;

      // correct x for periodicity
      if (dx > 0.5*xrange) dx= dx-xrange;
      if (dx <= -0.5*xrange) dx=dx+xrange;
      // correct y for periodicity
      if (dy > 0.5*yrange) dy= dy-yrange;
      if (dy <= -0.5*yrange) dy=dy+yrange;
      // correct z for periodicity
      if (dz > 0.5*zrange) dz=dz-zrange;
      if (dz <= -0.5*zrange) dz=dz+zrange;

      distance = sqrt(dx*dx + dy*dy + dz*dz);
      
      if (distance <= cutoff){
	energy += 4*epsilon*(pow(sigma/distance, 12) - pow(sigma/distance, 6));
      }
    }    
  }

  return energy;
}

int move(struct particle system[], double delta, double sigma, double epsilon,
	  int num, double temp, double xrange, double yrange, double zrange,  double cutoff){
  // displaces a particle by a set amount within the system.
  // app
  // 
  int rand_id;
  long double o_energy, n_energy;
  long double accept; // the acceptance probability
  long double beta = 1/(temp*1.38064e-23);
  long double prob;
  double trialx, trialy, trialz;

  // pick a particle
  rand_id = rand() % (num + 1);
  // calculate old energy
  o_energy = particle_energy(system, system[rand_id], sigma, epsilon,
			     num, xrange, yrange, zrange, cutoff);

  // create a trial position
  trialx = fmod(system[rand_id].x + delta*rand()/RAND_MAX, xrange);
  trialy = fmod(system[rand_id].y + delta*rand()/RAND_MAX, yrange);
  trialz = fmod(system[rand_id].z + delta*rand()/RAND_MAX, zrange);

  // calculate energy of trial position
  n_energy = trial_energy(system, rand_id, trialx, trialy, trialz,
			  sigma, epsilon, cutoff,
			  num, xrange, yrange, zrange);
  
  // Pick an acceptance number to sample the probability distribution
  accept = (double)rand()/RAND_MAX;

  // if accepted, change rand_id's positions to the trial position.
  // otherwise, do nothing.  

  prob = exp(-beta*(n_energy - o_energy));

  if (accept < prob){
    system[rand_id].x = trialx;
    system[rand_id].y = trialy;
    system[rand_id].z = trialz;
    return 1; // move made.
  } else {
    return 0; // move not made.
  }  
  
}

// -----------------------------------------------------------------------------------------------

void print_info(FILE *f, const char * filename,
		struct particle system[],
		int num, double sigma, double epsilon,
		double xrange, double yrange, double zrange, double cutoff){
  // produces a file in which the system information is inputted
  f = fopen(filename,"w");
  fprintf(f, "ID\t\tX\t\tY\t\tZ\t\tEnergy\n");  
  long double energy;
  long double total;
  for(int i=0;i<num;i++){
    energy = particle_energy(system, system[i], sigma, epsilon, num, xrange, yrange, zrange, cutoff);    
    fprintf(f, "%d\t\t%f\t%f\t%f\t%Lf\n", system[i].id, system[i].x, system[i].y, system[i].z, energy);
  }

  total = total_energy(system, sigma, epsilon, num, xrange, yrange, zrange, cutoff);
  fprintf(f, "\t\t\t\t\t\t\t\t%Lf\n", total);
  fclose(f);
}

int main(void){
  FILE *fp;
  // lennard jones parameters
  double sigma = 1.0;
  double epsilon = 1.0;
  double cutoff = 2.5;

  // simulation box
  double xrange = 5.0;
  double yrange = 5.0;
  double zrange = 5.0;
  double volume = xrange*yrange*zrange;

  long double initial_e;
  long double final_e;
  // specific parameters for monte carlo simulation  
  int num = 500; // num. particles
  double delta = 0.2; // translation distance
  double temp = 1.0;
  int ncycle = 100000; // number of mc cycles
  long int scycle = 0; // number of successful cycles.
  long double energy;
  
  // printf("BASIC METROPOLIS MONTE CARLO ALGORITHM:\n");

  // set random number seed. This only needs to be done once.
  srand((unsigned int) time(0));
  
  struct particle sample[num];
  initialise(sample, num, xrange, yrange, zrange);
  
  initial_e = total_energy(sample, sigma, epsilon, num, xrange, yrange, zrange, cutoff);  
  print_info(fp, "p_initial.txt", sample, num, sigma, epsilon, xrange, yrange, zrange, cutoff);
  
  for (int cycle=0;cycle<=ncycle;cycle++){
    scycle += move(sample, delta, sigma, epsilon, num, temp, xrange, yrange, zrange, cutoff);
  }
  
  final_e = total_energy(sample, sigma, epsilon, num, xrange, yrange, zrange, cutoff);
  print_info(fp, "p_output.txt", sample, num, sigma, epsilon, xrange, yrange, zrange, cutoff);

  fp = fopen("sim_details.txt", "w");
  fprintf(fp, "Simulation of Lennard-Jones liquid.\n");
  fprintf(fp, "Sigma: \t\t\t\t\t%f\n", sigma);
  fprintf(fp, "Epsilon: \t\t\t\t%f\n", epsilon);
  fprintf(fp, "Lennard-Jones Cutoff: \t\t\t%f\n", cutoff);
  fprintf(fp, "Volume: \t\t\t\t%f\n", volume);
  fprintf(fp, "Number of atoms: \t\t\t%d\n", num);
  fprintf(fp, "Temperature: \t\t\t\t%f\n", temp);
  fprintf(fp, "Number of successful moves: \t\t%ld/%d\n", scycle, ncycle);
  fprintf(fp, "Average energy (initial): \t\t%Lf\n", initial_e/num);
  fprintf(fp, "Average energy (final): \t\t%Lf\n", final_e/num);
  fclose(fp);
}
