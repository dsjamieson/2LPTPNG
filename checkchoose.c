#include <math.h>
#include "allvars.h"
#include "proto.h"

void checkchoose(void)
{

  // first check not both power and transfer are used

  if(WhichSpectrum != 0 && WhichTransfer !=0 ) {
      fprintf(stdout,"\n ERROR: You are running with two imput files: power and tranfer function \n Please select only one of both\n"); exit(2);}; 

  // second check if running ns!=1

  if(PrimordialIndex != 1.0) {

      fprintf(stdout,"\n You are running a tilted spectrum index\n");

#ifdef LOCAL_FNL 
  if(WhichSpectrum != 0) { fprintf(stdout,"\n Local ns with tilted power spectrum requires the transfer function as input\n switch WhichSpectrum to zero in the input parameter file\n"); exit(2);};
#endif 

#ifdef ORTOG_FNL 
      fprintf(stdout,"\n ERROR: For non-loccal non-Gaussian models\n The power spectrum index should be set to 1.0 \n Tilted power spectrum may be implemented in the next version\n"); exit(2); 
#endif 

#ifdef EQUIL_FNL 
      fprintf(stdout,"\n ERROR: For non-loccal non-Gaussian models\n The power spectrum index should be set to 1.0 \n Tilted power spectrum may be implemented in the next version\n"); exit(2); 
#endif 

  }


#ifdef LOCAL_FNL 
  if(WhichSpectrum != 0) { fprintf(stdout,"\n Local model requires the transfer function as input\n switch WhichSpectrum to zero in the input parameter file\n"); exit(2);};
#endif 

#ifdef EQUIL_FNL 
  if(WhichSpectrum != 0) { fprintf(stdout,"\n NonLocal Models require the input to be the Transfer Function,\n switch WhichSpectrum to zero in the input parameter file\n"); exit(2);};
#endif 
#ifdef ORTOG_FNL 
  if(WhichSpectrum != 0) { fprintf(stdout,"\n NonLocal Models require the input to be the Transfer Function,\n switch WhichSpectrum to zero in the input parameter file\n"); exit(2);};
#endif 

}

