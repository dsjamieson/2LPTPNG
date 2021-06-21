#include <math.h>
#include "allvars.h"
#include "proto.h"

void checkchoose(void)
{

  if(WhichSpectrum != 0 && WhichTransfer !=0 ) {
      fprintf(stdout,"\n ERROR: You are running with two imput files: power and tranfer function \n Please select only one of them\n"); exit(2);
  } 
  
  if (Fnl != 0.) {
    if(WhichSpectrum != 0) {
		fprintf(stdout,"\n Fnl != 0. requires the transfer function as input\n switch WhichSpectrum to zero in the input parameter file\n"); 
		exit(2);
    }
  }

}

