// Program Name: info.h
// Author: Aravinthen Rajkumar
// This program is used to show various information about a system.

#include <stdio.h>
#include "lib.h"

void pinfo(struct particle system[], int num){
  // prints all the data of the particles within the input system 
  printf("\nID\tType\tx\t\ty\t\tz\n");
  for (int i=0;i<num;i++){
    printf("%d\t%d\t%f\t%f\t%f\n", 
	   system[i].id,
	   system[i].species,
	   system[i].x,
	   system[i].y,
	   system[i].z);
  }
}

// ------------------------------------------------------------------
// File print methods
// ------------------------------------------------------------------ 
void gnupos(FILE*fp, struct particle system[], int num, char fname[]){
  // prints all the data of the particles within the input system
  fp = fopen(fname, "w");
  for (int i=0;i<num;i++){
    fprintf(fp,
	    "%f\t%f\t%f\n", 
	    system[i].x,
	    system[i].y,
	    system[i].z);
  }
  fclose(fp);
}

int type_count(struct particle system[], int total, int type){
  int count=0;
  for (int i=0;i<total;i++){
    if (system[i].species == type) count++;
  }
  return count;
}

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
	     double step){

  int typenum;
  // prints all the data of the particles within the input system
  info = fopen("STATS", "w");
  fprintf(info, "MC Simulation with %d particles of %d different types.\n", total, ntypes);
  fprintf(info, "Volume: %f*%f*%f = %f\n", xlen, ylen, zlen, xlen*ylen*zlen);
  fprintf(info, "Temperature: %f\n", T);
  fprintf(info, "Step length: %f\n", step);

  // Count the particles:
  for (int i=0;i<ntypes;i++){
    typenum = type_count(system, total, i);
    fprintf(info, "Number fraction of Type %d: %f\n", i, (float)typenum/total);
  }
  
  fprintf(info, "Acceptance rate: (%d/%d) %f\n", accepted, moves, (float)accepted/moves);
  fclose(info);
}
