// Program Name: moves.c
// Author: Aravinthen Rajkumar
// Description: various moves used in a Monte-Carlo simulation.

#include "math.h"
#include "lib.h"

double particle_move(struct particle system[],
		     struct type_matrix types,
		     int label,
		     double magnitude,
		     int nump,
		     double cutoff,
		     double xrange, double yrange, double zrange){
  
  // Moves a particle with |system|, then calculates the associated change
  // in energy of the move.
  
  double newx, newy, newz;
  double norm;
  double old_energy=0.0;
  double energy=0.0;
  double cx, cy, cz; // current particle position
  double nx, ny, nz; // new particle position
  double px, py, pz; // particles being sumed over
  double dx, dy, dz; // change in position
  double r; // normed distance
  double sig, eps;
  
  int nspec = system[label].species; // species of the labelled particle  
  cx = system[label].x;
  cy = system[label].y;
  cz = system[label].z;
  
  for (int i=0;i<nump;i++){
    if (i!=label){
      sig = types.sigma[system[i].species][nspec];
      eps = types.epsilon[system[i].species][nspec];

      // correct new positon for periodicity
      px = fmod(system[i].x + 2*xrange, xrange);
      py = fmod(system[i].y + 2*yrange, yrange);
      pz = fmod(system[i].z + 2*zrange, zrange);

      dx = px-cx;
      dy = py-cy;
      dz = pz-cz;
      r = sqrt(dx*dx + dy*dy + dz*dz);
      
      if (r < cutoff){	
	old_energy += 4*eps*(pow((sig/r), 12) - pow((sig/r), 6));
      }      
    }
  }

  // generate position vector
  newx = (double)rand()/RAND_MAX;
  newy = (double)rand()/RAND_MAX;
  newz = (double)rand()/RAND_MAX;

  // normalise position vector
  norm = sqrt(newx*newx + newy*newy + newz*newz);
  newx = newx/norm;
  newy = newy/norm;
  newz = newz/norm;

  // add position vector to particle position and apply periodic boundaries
  nx = fmod((system[label].x + newx), xrange);
  ny = fmod((system[label].y + newy), yrange);
  nz = fmod((system[label].z + newz), zrange);

  // assign new value to system
  system[label].x = nx;
  system[label].y = ny;
  system[label].z = nz;
  
  for (int j=0;j<nump;j++){
    if (j!=label){
      sig = types.sigma[system[j].species][nspec];
      eps = types.epsilon[system[j].species][nspec];

      // correct new positon for periodicity
      px = fmod(system[j].x + 2*xrange, xrange);
      py = fmod(system[j].y + 2*yrange, yrange);
      pz = fmod(system[j].z + 2*zrange, zrange);

      dx = px-nx;
      dy = py-ny;
      dz = pz-nz;
      r = sqrt(dx*dx + dy*dy + dz*dz);
      
      if (r < cutoff){	
	energy += 4*eps*(pow((sig/r), 12) - pow((sig/r), 6));
      }      

    }
  }
  return energy-old_energy;
}


