// Program Name: calcs.c
// Author: Aravinthen Rajkumar
// Description: the calculations required to make use of MC simulations

#include "lib.h"
#include "math.h"

int prand(int total){
  // produced a random number from 0 to 1
  return (rand() % (total)); ;
}

double urand(){
  // produced a random number from 0 to 1
  return (double)rand() / (double)((unsigned) RAND_MAX + 1) ;
}

double system_energy(struct particle system[],
		     struct type_matrix types,
		     int num_parts,
		     double cutoff,
		     double xrange, double yrange, double zrange){  
  // calculate the total system energy.
  // [P] THIS OUGHT TO BE PARALLELISED.

  double cx, cy, cz; // current particle
  double px, py, pz; // next particle
  double dx, dy, dz; // differences between particles
  double r; // distance for LJ formula. sticks with physics convention
  double sig, eps;

  double total_en = 0.0;
      
  for (int i=0;i<num_parts-1;i++){
    // obtain the position of the current particle once.
    cx = system[i].x;
    cy = system[i].y;
    cz = system[i].z;
    for (int j=i+1; j<num_parts;j++){
      // get type data
      sig = types.sigma[system[i].species][system[j].species];
      eps = types.epsilon[system[i].species][system[j].species];
      // Correct the compared position to account for periodicity
      // general rule: (px + 2*xrange) mod xrange
      px = fmod(system[j].x + 2*xrange, xrange);
      py = fmod(system[j].y + 2*yrange, yrange);
      pz = fmod(system[j].z + 2*zrange, zrange);

      dx = px-cx;
      dy = py-cy;
      dz = pz-cz;

      r = sqrt(dx*dx + dy*dy + dz*dz); 

      if (r < cutoff){	
	total_en += 4*eps*(pow((sig/r), 12) - pow((sig/r), 6));
      }      
    }
  }
  
  return total_en;
  
}
