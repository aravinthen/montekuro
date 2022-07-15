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
