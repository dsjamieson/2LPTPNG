#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <gsl/gsl_rng.h>
#include "allvars.h"
#include "proto.h"

#define ASSERT_ALLOC(cond) {                                                                                  \
  if(!cond)                                                                                                   \
    {                                                                                                         \
      printf("failed to allocate %g Mbyte on Task %d\n", bytes / (1024.0 * 1024.0), ThisTask);                \
      printf("bailing out.\n");                                                                               \
      FatalError(1);                                                                                          \
    }                                                                                                         \
}

int frequency_of_primes (int n) {
  int i,j;
  int freq=n-1;
  for (i=2; i<=n; ++i) for (j=sqrt(i);j>1;--j) if (i%j==0) {--freq; break;}
  return freq;
}

void print_timed_done(int n) {
  clock_t tot_time = clock() - start_time;
  int tot_hours = (int) floor(((double) tot_time) / 60. / 60. / CLOCKS_PER_SEC);
  int tot_mins = (int) floor(((double) tot_time) / 60. / CLOCKS_PER_SEC) - 60 * tot_hours;
  double tot_secs = (((double) tot_time) / CLOCKS_PER_SEC) - 60. * (((double) tot_mins) * 60. * ((double) tot_hours));
  clock_t diff_time = clock() - previous_time;
  int diff_hours = (int) floor(((double) diff_time) / 60. / 60. / CLOCKS_PER_SEC);
  int diff_mins = (int) floor(((double) diff_time) / 60. / CLOCKS_PER_SEC) - 60 *diff_hours;
  double diff_secs = (((double) diff_time) / CLOCKS_PER_SEC) - 60. * (((double) diff_mins) * 60. * ((double) diff_hours));
  for (int i = 0; i < n; i++)
   printf(" ");
  printf("Done [%02d:%02d:%05.2f, %02d:%02d:%05.2f]\n", diff_hours, diff_mins, diff_secs, tot_hours, tot_mins, tot_secs);
  previous_time = clock();
  return;
}

int main(int argc, char **argv) {

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);

  start_time = clock();
  previous_time = start_time;

  if(argc < 2) {
    if(ThisTask == 0) {
	  fprintf(stdout, "\nParameters are missing.\n");
	  fprintf(stdout, "Call with <ParameterFile>\n\n");
	}
    MPI_Finalize();
    exit(0);
  }

  read_parameterfile(argv[1]);
  checkchoose();  
  set_units();
  initialize_transferfunction(); 
  initialize_powerspectrum(); 
  initialize_ffts(); 
  read_glass(GlassFile);

  if (ThisTask == 0)
    print_setup();

  displacement_fields();

  if (ThisTask == 0) {printf("Writing initial conditions snapshot..."); fflush(stdout);};
  write_particle_data();
  if (ThisTask == 0 ) print_timed_done(10);
  if(NumPart)
    free(P);
  free_ffts();
  MPI_Barrier(MPI_COMM_WORLD);
  print_spec();
  MPI_Finalize();		/* clean up & finalize MPI */
  exit(0);
}

void print_setup(void) {
  char pstr[79];

  for (int i = 0; i < 79; i++)
    pstr[i] = '*';
  printf("%s\n", pstr);

#ifdef ONLY_GAUSSIAN
  char exec[] = "  2LPT  ";
#endif
#ifdef LOCAL_FNL  
  char exec[] = "2LPTNGLC";
#endif
#ifdef EQUIL_FNL 
  char exec[] = "2LPTNGEQ";
#endif
#ifdef ORTOG_FNL 
  char exec[] = "2LPTNGOR";
#endif
//wrc
#ifdef ORTOG_LSS_FNL 
  char exec[] = "2LPTNGORLSS";
#endif


  pstr[0] = '*';
  for (int i = 1; i < 79/2 - 3; i++)
    pstr[i] = ' ';
  sprintf(&pstr[79/2 - 4], "%s", exec);
  for (int i = 79/2 + 4; i < 78; i++)
    pstr[i] = ' ';
  pstr[78] = '*';
  printf("%s\n", pstr);
  for (int i = 0; i < 79; i++)
    pstr[i] = '*';
  printf("%s\n", pstr);

  printf("\n");

  printf("      Box = %.2e   Nmesh = %d   Nsample = %d\n", Box, Nmesh, Nsample);
  printf("                       Nglass = %d    GlassTileFac = %d\n\n", Nglass, GlassTileFac);
  printf("        Omega = %.4f         OmegaLambda = %.4f    OmegaBaryon = %.2e\n", Omega,  OmegaLambda, OmegaBaryon);
  printf("       sigma8 = %.4f     PrimoridalIndex = %.4f       Redshift = %.2e\n", Sigma8, PrimordialIndex, Redshift);
  printf("  HubbleParam = %.4f  OmegaDM_2ndSpecies = %.2e          fNL = %+.2e\n\n", HubbleParam, OmegaDM_2ndSpecies, Fnl); 
  printf("   FixedAmplitude = %d    PhaseFlip = % d   SphereMode = %d    Seed = %d\n", FixedAmplitude, PhaseFlip, SphereMode, Seed);   

  for (int i = 0; i < 79; i++)
    pstr[i] = '*';
  printf("\n%s\n\n", pstr);
  return;
}

void displacement_fields(void) {
  MPI_Request request;
  MPI_Status status;
  gsl_rng *random_generator;
  int i, j, k, ii, jj, kk, axes;
  int n;
  int sendTask, recvTask;
#ifdef ONLY_GAUSSIAN
  double p_of_k, delta;
#else
  double t_of_k, phig, Beta, twb;
  fftw_complex *(cpot);  /* For computing nongaussian fnl ic */
  fftw_real *(pot);
#ifndef LOCAL_FNL  
// *************************** DSJ *******************************
  double exp_1_over_3, exp_1_over_6, kmag_1_over_3, kmag_2_over_3;
// *************************** DSJ *******************************
  fftw_complex *(cpartpot); /* For non-local fluctuations */
  fftw_real *(partpot);
  fftw_complex *(cp1p2p3sym);
  fftw_real *(p1p2p3sym);
  fftw_complex *(cp1p2p3sca);
  fftw_real *(p1p2p3sca);
  fftw_complex *(cp1p2p3nab);
  fftw_real *(p1p2p3nab);
  fftw_complex *(cp1p2p3tre);
  fftw_real *(p1p2p3tre);
// ****  wrc ****
#ifdef ORTOG_LSS_FNL
  // Will contain phi_G/k (hence invs)
  fftw_complex *(cp1p2p3inv);
  fftw_real *(p1p2p3inv);
  // names correspond to kernels in Appendix A of Coulton et al.
  fftw_complex *(cp1p2p3_K12D);
  fftw_real *(p1p2p3_K12D);
  fftw_complex *(cp1p2p3_K12E);
  fftw_real *(p1p2p3_K12E);
  fftw_complex *(cp1p2p3_K12F);
  fftw_real *(p1p2p3_K12F);
  fftw_complex *(cp1p2p3_K12G);
  fftw_real *(p1p2p3_K12G);
  double orth_p = 27/(-21+743/(7*(20*pow(PI,2)-193)));
  double orth_t = (2+20/9*orth_p)/(6+10/3*orth_p);
#endif
// ****  wrc ****
#endif
#endif


  double fac, vel_prefac, vel_prefac2;
  double kvec[3], kmag, kmag2;
  double phase, ampl, hubble_a;
  double u, v, w;
  double f1, f2, f3, f4, f5, f6, f7, f8;
  double dis, dis2, maxdisp, max_disp_glob;
  unsigned int *seedtable;
// ******* FAVN *****
  double phase_shift; 
// ******* FAVN *****

  unsigned int bytes, nmesh3;
  int coord;
  fftw_complex *(cdisp[3]), *(cdisp2[3]) ; /* ZA and 2nd order displacements */
  fftw_real *(disp[3]), *(disp2[3]) ;

  fftw_complex *(cdigrad[6]);
  fftw_real *(digrad[6]);



#ifdef CORRECT_CIC
  double fx, fy, fz, ff, smth;
#endif

  hubble_a = Hubble * sqrt(Omega / pow(InitTime, 3) + (1 - Omega - OmegaLambda) / pow(InitTime, 2) + OmegaLambda);
  vel_prefac = InitTime * hubble_a * F_Omega(InitTime) / sqrt(InitTime);
  vel_prefac2 = InitTime * hubble_a * F2_Omega(InitTime) / sqrt(InitTime);

// ******************************************** FAVN **********************************************
  phase_shift = 0.0;
  if (PhaseFlip==1)
  	phase_shift = PI;
// ******************************************** FAVN **********************************************

  fac = pow(2 * PI / Box, 1.5);

  maxdisp = 0;

  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

  gsl_rng_set(random_generator, Seed);

  if(!(seedtable = malloc(Nmesh * Nmesh * sizeof(unsigned int))))
    FatalError(4);

  for(i = 0; i < Nmesh / 2; i++)
    {
      for(j = 0; j < i; j++)
	seedtable[i * Nmesh + j] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i + 1; j++)
	seedtable[j * Nmesh + i] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i; j++)
	seedtable[(Nmesh - 1 - i) * Nmesh + j] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i + 1; j++)
	seedtable[(Nmesh - 1 - j) * Nmesh + i] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i; j++)
	seedtable[i * Nmesh + (Nmesh - 1 - j)] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i + 1; j++)
	seedtable[j * Nmesh + (Nmesh - 1 - i)] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i; j++)
	seedtable[(Nmesh - 1 - i) * Nmesh + (Nmesh - 1 - j)] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i + 1; j++)
	seedtable[(Nmesh - 1 - j) * Nmesh + (Nmesh - 1 - i)] = 0x7fffffff * gsl_rng_uniform(random_generator);
    }


#ifdef ONLY_GAUSSIAN

  if (ThisTask == 0) {printf("Setting up gradient of Gaussian potential..."); fflush(stdout);};

  for(axes=0,bytes=0; axes < 3; axes++)
    {
      cdisp[axes] = (fftw_complex *) malloc(bytes += sizeof(fftw_real) * TotalSizePlusAdditional);
      disp[axes] = (fftw_real *) cdisp[axes];
    }

  ASSERT_ALLOC(cdisp[0] && cdisp[1] && cdisp[2]);


#if defined(MULTICOMPONENTGLASSFILE) && defined(DIFFERENT_TRANSFER_FUNC)
  for(Type = MinType; Type <= MaxType; Type++)
#endif
    {

      /* first, clean the array */
      for(i = 0; i < Local_nx; i++)
	for(j = 0; j < Nmesh; j++)
	  for(k = 0; k <= Nmesh / 2; k++)
	    for(axes = 0; axes < 3; axes++)
	      {
		cdisp[axes][(i * Nmesh + j) * (Nmesh / 2 + 1) + k].re = 0;
		cdisp[axes][(i * Nmesh + j) * (Nmesh / 2 + 1) + k].im = 0;
	      }

      for(i = 0; i < Nmesh; i++)
	{
	  ii = Nmesh - i;
	  if(ii == Nmesh)
	    ii = 0;
	  if((i >= Local_x_start && i < (Local_x_start + Local_nx)) ||
	     (ii >= Local_x_start && ii < (Local_x_start + Local_nx)))
	    {
	      for(j = 0; j < Nmesh; j++)
		{
		  gsl_rng_set(random_generator, seedtable[i * Nmesh + j]);
		  
		  for(k = 0; k < Nmesh / 2; k++)
		    {
		      phase = gsl_rng_uniform(random_generator) * 2 * PI;
//************ FAVN ***************
              phase += phase_shift;
//************ FAVN ***************
		      do
			ampl = gsl_rng_uniform(random_generator);
		      while(ampl == 0);
		      
		      if(i == Nmesh / 2 || j == Nmesh / 2 || k == Nmesh / 2)
			continue;
		      if(i == 0 && j == 0 && k == 0)
			continue;
		      
		      if(i < Nmesh / 2)
			kvec[0] = i * 2 * PI / Box;
		      else
			kvec[0] = -(Nmesh - i) * 2 * PI / Box;
		      
		      if(j < Nmesh / 2)
			kvec[1] = j * 2 * PI / Box;
		      else
			kvec[1] = -(Nmesh - j) * 2 * PI / Box;
		      
		      if(k < Nmesh / 2)
			kvec[2] = k * 2 * PI / Box;
		      else
			kvec[2] = -(Nmesh - k) * 2 * PI / Box;
		      
		      kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
		      kmag = sqrt(kmag2);
		      
		      if(SphereMode == 1)
			{
			  if(kmag * Box / (2 * PI) > Nsample / 2)	/* select a sphere in k-space */
			    continue;
			}
		      else
			{
			  if(fabs(kvec[0]) * Box / (2 * PI) > Nsample / 2)
			    continue;
			  if(fabs(kvec[1]) * Box / (2 * PI) > Nsample / 2)
			    continue;
			  if(fabs(kvec[2]) * Box / (2 * PI) > Nsample / 2)
			    continue;
			}
		      
		      p_of_k = PowerSpec(kmag);	      
// ************ FAVN/DSJ ************
			  if (!FixedAmplitude)
		        p_of_k *= -log(ampl);
// ************ FAVN/DSJ ************
		      
		      delta = fac * sqrt(p_of_k) / Dplus;	/* scale back to starting redshift */
		      
		      if(k > 0)
			{
			  if(i >= Local_x_start && i < (Local_x_start + Local_nx))
			    for(axes = 0; axes < 3; axes++)
			      {
				cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].re =
				  -kvec[axes] / kmag2 * delta * sin(phase);
				cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].im =
				  kvec[axes] / kmag2 * delta * cos(phase);
			      }
			}
		      else	/* k=0 plane needs special treatment */
			{
			  if(i == 0)
			    {
			      if(j >= Nmesh / 2)
				continue;
			      else
				{
				  if(i >= Local_x_start && i < (Local_x_start + Local_nx))
				    {
				      jj = Nmesh - j;	/* note: j!=0 surely holds at this point */
				      
				      for(axes = 0; axes < 3; axes++)
					{
					  cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].re =
					    -kvec[axes] / kmag2 * delta * sin(phase);
					  cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].im =
					    kvec[axes] / kmag2 * delta * cos(phase);
					  
					  cdisp[axes][((i - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k].re =
					    -kvec[axes] / kmag2 * delta * sin(phase);
					  cdisp[axes][((i - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k].im =
					    -kvec[axes] / kmag2 * delta * cos(phase);
					}
				    }
				}
			    }
			  else	/* here comes i!=0 : conjugate can be on other processor! */
			    {
			      if(i >= Nmesh / 2)
				continue;
			      else
				{
				  ii = Nmesh - i;
				  if(ii == Nmesh)
				    ii = 0;
				  jj = Nmesh - j;
				  if(jj == Nmesh)
				    jj = 0;
				  
				  if(i >= Local_x_start && i < (Local_x_start + Local_nx))
				    for(axes = 0; axes < 3; axes++)
				      {
					cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].re =
					  -kvec[axes] / kmag2 * delta * sin(phase);
					cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].im =
					  kvec[axes] / kmag2 * delta * cos(phase);
				      }
				  
				  if(ii >= Local_x_start && ii < (Local_x_start + Local_nx))
				    for(axes = 0; axes < 3; axes++)
				      {
					cdisp[axes][((ii - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) +
					      k].re = -kvec[axes] / kmag2 * delta * sin(phase);
					cdisp[axes][((ii - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) +
					      k].im = -kvec[axes] / kmag2 * delta * cos(phase);
				      }
				}
			    } 
			}
		    }
		}
	    }
	}

  if (ThisTask == 0 ) print_timed_done(4);
      
#else    /* non gaussian initial potential fnl type  */

      if (ThisTask == 0) {printf("Setting Gaussian potential..."); fflush(stdout);};

      bytes=0; /*initialize*/
      cpot = (fftw_complex *) malloc(bytes += sizeof(fftw_real) * TotalSizePlusAdditional);
      pot = (fftw_real *) cpot;

      ASSERT_ALLOC(cpot);


      /* first, clean the cpot array */
      for(i = 0; i < Local_nx; i++)
        for(j = 0; j < Nmesh; j++)
          for(k = 0; k <= Nmesh / 2; k++)
              {
                cpot[(i * Nmesh + j) * (Nmesh / 2 + 1) + k].re = 0;
                cpot[(i * Nmesh + j) * (Nmesh / 2 + 1) + k].im = 0;
              }

      /* Ho in units of UnitLength_in_cm and c=1, i.e., internal units so far  */
      /* Beta = 3/2 H(z)^2 a^3 Om(a) / D0 = 3/2 Ho^2 Om0 / D0 at redshift z = 0.0 */ 
      Beta = 1.5 * Omega / (2998. * 2998. / UnitLength_in_cm / UnitLength_in_cm * 3.085678e24 * 3.085678e24 ) / D0 ;        

      for(i = 0; i < Nmesh; i++)
        {
          ii = Nmesh - i;
          if(ii == Nmesh)
            ii = 0;
          if((i >= Local_x_start && i < (Local_x_start + Local_nx)) ||
             (ii >= Local_x_start && ii < (Local_x_start + Local_nx)))
            {
              for(j = 0; j < Nmesh; j++)
                {
                  gsl_rng_set(random_generator, seedtable[i * Nmesh + j]);

                  for(k = 0; k < Nmesh / 2; k++)
                    {
                      phase = gsl_rng_uniform(random_generator) * 2 * PI;
// ***************** FAVN *****************
                      phase += phase_shift;
// ***************** FAVN *****************
                      do
                        ampl = gsl_rng_uniform(random_generator);

                      while(ampl == 0);

                      if(i == Nmesh / 2 || j == Nmesh / 2 || k == Nmesh / 2)
                        continue;
                      if(i == 0 && j == 0 && k == 0)
                        continue;

                      if(i < Nmesh / 2)
                        kvec[0] = i * 2 * PI / Box;
                      else
                        kvec[0] = -(Nmesh - i) * 2 * PI / Box;

                      if(j < Nmesh / 2)
                        kvec[1] = j * 2 * PI / Box;
                      else
                        kvec[1] = -(Nmesh - j) * 2 * PI / Box;

                      if(k < Nmesh / 2)
                        kvec[2] = k * 2 * PI / Box;
                      else
                        kvec[2] = -(Nmesh - k) * 2 * PI / Box;

                      kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
                      kmag = sqrt(kmag2);

                      if(SphereMode == 1)
                        {
                          if(kmag * Box / (2 * PI) > Nsample / 2)       /* select a sphere in k-space */
                            continue;
                        }
                      else
                        {
                          if(fabs(kvec[0]) * Box / (2 * PI) > Nsample / 2)
                            continue;
                          if(fabs(kvec[1]) * Box / (2 * PI) > Nsample / 2)
                            continue;
                          if(fabs(kvec[2]) * Box / (2 * PI) > Nsample / 2)
                            continue;
                        }

                      phig = Anorm * exp( PrimordialIndex * log(kmag) );   /* initial normalized power */
// ************** FAVN/DSJ ***************
					  if (!FixedAmplitude)
					  	phig *= -log(ampl);
// ***************** FAVN/DSJ *************
                      
                      phig = sqrt(phig) * fac * Beta / kmag2;    /* amplitude of the initial gaussian potential */
               
                      if(k > 0)
                        {
                          if(i >= Local_x_start && i < (Local_x_start + Local_nx))
                               {

                                coord = ((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k;

                                cpot[coord].re = phig * cos(phase);
                                cpot[coord].im = phig * sin(phase);

                               }
                        }
                      else      /* k=0 plane needs special treatment */
                        {
                          if(i == 0)
                            {
                              if(j >= Nmesh / 2)
                                continue;
                              else
                                {
                                  if(i >= Local_x_start && i < (Local_x_start + Local_nx))
                                    {
                                      jj = Nmesh - j;   /* note: j!=0 surely holds at this point */

                                          coord = ((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k;

                                          cpot[coord].re = phig * cos(phase);
                                          cpot[coord].im = phig * sin(phase);


                                          coord = ((i - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k; 
                                          cpot[coord].re =  phig * cos(phase);
                                          cpot[coord].im = -phig * sin(phase);

                                    }
                                }
                            }
                          else  /* here comes i!=0 : conjugate can be on other processor! */
                            {
                              if(i >= Nmesh / 2)
                                continue;
                              else
                                {
                                  ii = Nmesh - i;
                                  if(ii == Nmesh)
                                    ii = 0;
                                  jj = Nmesh - j;
                                  if(jj == Nmesh)
                                    jj = 0;

                                  if(i >= Local_x_start && i < (Local_x_start + Local_nx))
                                      {

                                        coord = ((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k;
          
                                        cpot[coord].re = phig * cos(phase);
                                        cpot[coord].im = phig * sin(phase);

                                      }
                                  if(ii >= Local_x_start && ii < (Local_x_start + Local_nx))
                                      {
                                        coord = ((ii - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k;

                                        cpot[coord].re = phig * cos(phase);
                                        cpot[coord].im = -phig * sin(phase);

                                      }
                                }
                            }
                        }
                    }
                }
            }
        }

     if (ThisTask == 0 ) print_timed_done(19);

     /*** For non-local models it is important to keep all factors of SQRT(-1) as done below ***/
     /*** Notice also that there is a minus to convert from Bardeen to gravitational potential ***/

#ifdef LOCAL_FNL  

      if (ThisTask == 0) {printf("Computing local non-Gaussian potential..."); fflush(stdout);};

      /******* LOCAL PRIMORDIAL POTENTIAL ************/
      rfftwnd_mpi(Inverse_plan, 1, pot, Workspace, FFTW_NORMAL_ORDER);
      fflush(stdout);

      /* square the potential in configuration space */

      for(i = 0; i < Local_nx; i++)
        for(j = 0; j < Nmesh; j++)
          for(k = 0; k < Nmesh; k++)
            {
             coord = (i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k; 
   
             pot[coord] = pot[coord] + Fnl * pot[coord]*pot[coord];
   
            }

      MPI_Barrier(MPI_COMM_WORLD);
      rfftwnd_mpi(Forward_plan, 1, pot, Workspace, FFTW_NORMAL_ORDER);

       /* remove the N^3 I got by forwardfurier and put zero to zero mode */

      nmesh3 = ((unsigned int) Nmesh) * ((unsigned int) Nmesh ) * ((unsigned int) Nmesh);    
      for(i = 0; i < Local_nx; i++)
        for(j = 0; j < Nmesh; j++)
          for(k = 0; k <= Nmesh / 2 ; k++)
            {
              coord = (i * Nmesh + j) * (Nmesh / 2 + 1) + k;

// ****************************** DSJ *************************
		      if(SphereMode == 1)
			{
              ii = i + Local_x_start;

			  if(ii < Nmesh / 2)
				kvec[0] = ii * 2 * PI / Box;
			  else
				kvec[0] = -(Nmesh - ii) * 2 * PI / Box;

			  if(j < Nmesh / 2)
				kvec[1] = j * 2 * PI / Box;
			  else
				kvec[1] = -(Nmesh - j) * 2 * PI / Box;

			  if(k < Nmesh / 2)
				kvec[2] = k * 2 * PI / Box;
			  else
				kvec[2] = -(Nmesh - k) * 2 * PI / Box;

			  kmag = sqrt(kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2]);

			  if(kmag * Box / (2 * PI) > Nsample / 2) {	/* select a sphere in k-space */
				cpot[coord].re = 0.;
				cpot[coord].im = 0.;
				continue;
			  }
			}
// ****************************** DSJ *************************

              cpot[coord].re /= (double) nmesh3; 
              cpot[coord].im /= (double) nmesh3; 

            }
        
       if(ThisTask == 0) {
              cpot[0].re = 0.;
              cpot[0].im = 0.; 
           }


      if (ThisTask == 0 ) print_timed_done(7);

#else

#ifdef EQUIL_FNL
      if (ThisTask == 0) {printf("Computing equilateral non-Gaussian potential..."); fflush(stdout);};
#endif
#ifdef ORTOG_FNL
      if (ThisTask == 0) {printf("Computing orthogonal non-Gaussian potential... "); fflush(stdout);};
#endif

// ****  wrc ****
#ifdef ORTOG_LSS_FNL
      if (ThisTask == 0) {printf("Computing orthogonal-LSS non-Gaussian potential... "); fflush(stdout);};
#endif

// ****  wrc ****
      /**** NON-LOCAL PRIMORDIAL POTENTIAL **************/ 

          /* allocate partpotential */
	
	  cpartpot = (fftw_complex *) malloc(bytes = sizeof(fftw_real) * TotalSizePlusAdditional);
	  partpot = (fftw_real *) cpartpot;
	  ASSERT_ALLOC(cpartpot);

	  cp1p2p3sym = (fftw_complex *) malloc(bytes = sizeof(fftw_real) * TotalSizePlusAdditional);
	  p1p2p3sym = (fftw_real *) cp1p2p3sym;
	  ASSERT_ALLOC(cp1p2p3sym);

	  cp1p2p3sca = (fftw_complex *) malloc(bytes = sizeof(fftw_real) * TotalSizePlusAdditional);
	  p1p2p3sca = (fftw_real *) cp1p2p3sca;
	  ASSERT_ALLOC(cp1p2p3sca);

	  cp1p2p3nab = (fftw_complex *) malloc(bytes = sizeof(fftw_real) * TotalSizePlusAdditional);
	  p1p2p3nab = (fftw_real *) cp1p2p3nab;
	  ASSERT_ALLOC(cp1p2p3nab);

     	  cp1p2p3tre = (fftw_complex *) malloc(bytes = sizeof(fftw_real) * TotalSizePlusAdditional);
	  p1p2p3tre = (fftw_real *) cp1p2p3tre;
	  ASSERT_ALLOC(cp1p2p3tre);


// ****  wrc ****
#ifdef ORTOG_LSS_FNL
    cp1p2p3inv = (fftw_complex *) malloc(bytes = sizeof(fftw_real) * TotalSizePlusAdditional);
    p1p2p3inv  = (fftw_real *) cp1p2p3inv ;
    ASSERT_ALLOC(cp1p2p3inv);
    cp1p2p3_K12D = (fftw_complex *) malloc(bytes = sizeof(fftw_real) * TotalSizePlusAdditional);
    p1p2p3_K12D  = (fftw_real *) cp1p2p3_K12D ;
    ASSERT_ALLOC(cp1p2p3_K12D);
    cp1p2p3_K12E = (fftw_complex *) malloc(bytes = sizeof(fftw_real) * TotalSizePlusAdditional);
    p1p2p3_K12E  = (fftw_real *) cp1p2p3_K12E ;
    ASSERT_ALLOC(cp1p2p3_K12E);
    cp1p2p3_K12F = (fftw_complex *) malloc(bytes = sizeof(fftw_real) * TotalSizePlusAdditional);
    p1p2p3_K12F  = (fftw_real *) cp1p2p3_K12F ;
    ASSERT_ALLOC(cp1p2p3_K12F);
    cp1p2p3_K12G = (fftw_complex *) malloc(bytes = sizeof(fftw_real) * TotalSizePlusAdditional);
    p1p2p3_K12G  = (fftw_real *) cp1p2p3_K12G ;
    ASSERT_ALLOC(cp1p2p3_K12G);
#endif
// ****  wrc ****
// ******************* DSJ ***********************
	  exp_1_over_3 = (4. - PrimordialIndex) / 3.; // n_s modified exponent for generalized laplacian/inverse laplacian, exponent for k^2
	  exp_1_over_6 = (4. - PrimordialIndex) / 6.; // n_s modified exponent for generalized conjugate gradient magnitude and its inverse, exponent for k^2 
// ******************* DSJ ***********************

      /* first, clean the array */
      for(i = 0; i < Local_nx; i++)
	for(j = 0; j < Nmesh; j++)
	  for(k = 0; k <= Nmesh / 2; k++)
	      {
		coord = (i * Nmesh + j) * (Nmesh / 2 + 1) + k;
		cp1p2p3sym[coord].re = 0;
		cp1p2p3sym[coord].im = 0;
		cp1p2p3sca[coord].re = 0;
		cp1p2p3sca[coord].im = 0;
                cp1p2p3nab[coord].re = 0;
                cp1p2p3nab[coord].im = 0;
                cp1p2p3tre[coord].re = 0;
                cp1p2p3tre[coord].im = 0;

// ****  wrc ****
#ifdef ORTOG_LSS_FNL
                cp1p2p3inv[coord].re = 0;
                cp1p2p3inv[coord].im = 0;
#endif
// ****  wrc ****
		cpartpot[coord].re = 0;
		cpartpot[coord].im = 0;
	      }


      /* multiply by k */
          
      for(ii = 0; ii < Local_nx; ii++)
        for(j = 0; j < Nmesh; j++)
          for(k = 0; k <= Nmesh / 2 ; k++)
            {

              coord = (ii * Nmesh + j) * (Nmesh / 2 + 1) + k;
              i = ii + Local_x_start;
   
                   /* already are zero */ 
                   /*  if(i == 0 && j == 0 && k == 0); continue */

                      if(i < Nmesh / 2)
                        kvec[0] = i * 2 * PI / Box;
                      else
                        kvec[0] = -(Nmesh - i) * 2 * PI / Box;

                      if(j < Nmesh / 2)
                        kvec[1] = j * 2 * PI / Box;
                      else
                        kvec[1] = -(Nmesh - j) * 2 * PI / Box;

                      if(k < Nmesh / 2)
                        kvec[2] = k * 2 * PI / Box;
                      else
                        kvec[2] = -(Nmesh - k) * 2 * PI / Box;

                      kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
// ************************************ DSJ ********************************
					  kmag_2_over_3 = pow(kmag2, exp_1_over_3);                      
					  kmag_1_over_3 = pow(kmag2, exp_1_over_6);                      

                      cpartpot[coord].re = kmag_1_over_3 * cpot[coord].re;
                      cpartpot[coord].im = kmag_1_over_3 * cpot[coord].im;
                      
                      cp1p2p3nab[coord].re = kmag_2_over_3 * cpot[coord].re;
                      cp1p2p3nab[coord].im = kmag_2_over_3 * cpot[coord].im; 

// ****  wrc ****
#ifdef ORTOG_LSS_FNL
                      cp1p2p3inv[coord].re =  cpot[coord].re/kmag_1_over_3;
                      cp1p2p3inv[coord].im =  cpot[coord].im/kmag_1_over_3; 
#endif
// ****  wrc ****                  
// ************************************ DSJ ********************************
            }

       MPI_Barrier(MPI_COMM_WORLD);

       /* Fourier back to real */
      rfftwnd_mpi(Inverse_plan, 1, pot, Workspace, FFTW_NORMAL_ORDER);
      rfftwnd_mpi(Inverse_plan, 1, partpot, Workspace, FFTW_NORMAL_ORDER);
      rfftwnd_mpi(Inverse_plan, 1, p1p2p3nab, Workspace, FFTW_NORMAL_ORDER);

// ****  wrc ****
#ifdef ORTOG_LSS_FNL
      rfftwnd_mpi(Inverse_plan, 1, p1p2p3inv, Workspace, FFTW_NORMAL_ORDER);
#endif
// ****  wrc ****
      MPI_Barrier(MPI_COMM_WORLD);

      /* multiplying terms in real space  */


      for(i = 0; i < Local_nx; i++)
        for(j = 0; j < Nmesh; j++)
          for(k = 0; k < Nmesh; k++)
            {
             coord = (i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k; 

// ****  wrc ****
#ifdef ORTOG_LSS_FNL
              p1p2p3_K12D[coord] = p1p2p3inv[coord]*p1p2p3inv[coord]; // 1/k phi 1/k phi
              p1p2p3_K12E[coord] = p1p2p3nab[coord]*p1p2p3inv[coord];   //k^2 phi 1/k phi
              p1p2p3_K12F[coord] = pot[coord]*p1p2p3inv[coord];   //phi 1/k phi
              p1p2p3_K12G[coord] = partpot[coord]*p1p2p3inv[coord]; // k phi 1/k phk
#endif
// ****  wrc ****
              p1p2p3sym[coord] = partpot[coord]*partpot[coord];
              p1p2p3sca[coord] = pot[coord]*partpot[coord];   
              p1p2p3nab[coord] = pot[coord]*p1p2p3nab[coord];
              p1p2p3tre[coord] = p1p2p3nab[coord]*partpot[coord];
              partpot[coord] = pot[coord]*pot[coord];       /** NOTE: now partpot is potential squared **/

            }

      MPI_Barrier(MPI_COMM_WORLD);
      rfftwnd_mpi(Forward_plan, 1, pot, Workspace, FFTW_NORMAL_ORDER);
      rfftwnd_mpi(Forward_plan, 1, partpot, Workspace, FFTW_NORMAL_ORDER);
      rfftwnd_mpi(Forward_plan, 1, p1p2p3sym, Workspace, FFTW_NORMAL_ORDER);
      rfftwnd_mpi(Forward_plan, 1, p1p2p3sca, Workspace, FFTW_NORMAL_ORDER);
      rfftwnd_mpi(Forward_plan, 1, p1p2p3nab, Workspace, FFTW_NORMAL_ORDER);
      rfftwnd_mpi(Forward_plan, 1, p1p2p3tre, Workspace, FFTW_NORMAL_ORDER);

// ****  wrc ****
#ifdef ORTOG_LSS_FNL
      rfftwnd_mpi(Forward_plan, 1, p1p2p3_K12D, Workspace, FFTW_NORMAL_ORDER);
      rfftwnd_mpi(Forward_plan, 1, p1p2p3_K12E, Workspace, FFTW_NORMAL_ORDER);
      rfftwnd_mpi(Forward_plan, 1, p1p2p3_K12F, Workspace, FFTW_NORMAL_ORDER);
      rfftwnd_mpi(Forward_plan, 1, p1p2p3_K12G, Workspace, FFTW_NORMAL_ORDER);
#endif

// ****  wrc ****
      MPI_Barrier(MPI_COMM_WORLD);

       /* divide by appropiate k's, sum terms according to non-local model */ 
       /* remove the N^3 I got by forwardfurier and put zero to zero mode */

      nmesh3 = ((unsigned int) Nmesh) * ((unsigned int) Nmesh ) * ((unsigned int) Nmesh);    

      for(ii = 0; ii < Local_nx; ii++)
        for(j = 0; j < Nmesh; j++)
          for(k = 0; k <= Nmesh / 2 ; k++)
            {

              coord = (ii * Nmesh + j) * (Nmesh / 2 + 1) + k;
              i = ii + Local_x_start;
    
                   /* if(i == 0 && j == 0 && k == 0); continue; */

                      if(i < Nmesh / 2)
                        kvec[0] = i * 2 * PI / Box;
                      else
                        kvec[0] = -(Nmesh - i) * 2 * PI / Box;

                      if(j < Nmesh / 2)
                        kvec[1] = j * 2 * PI / Box;
                      else
                        kvec[1] = -(Nmesh - j) * 2 * PI / Box;

                      if(k < Nmesh / 2)
                        kvec[2] = k * 2 * PI / Box;
                      else
                        kvec[2] = -(Nmesh - k) * 2 * PI / Box;

                      kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
// ****************************** DSJ *************************
					  kmag_2_over_3 = pow(kmag2, exp_1_over_3);                      
					  kmag_1_over_3 = pow(kmag2, exp_1_over_6);
// ****************************** DSJ *************************
                      
                      if(i == 0 && j == 0 && k == 0)
                        {
                          cpot[0].re=0.;
                          cpot[0].im=0.;
                          continue;
                        }

// ****************************** DSJ *************************
		    if(SphereMode == 1) {
			  kmag = sqrt(kmag2);
			  if(kmag * Box / (2 * PI) > Nsample / 2) {	/* select a sphere in k-space */
				cpot[coord].re = 0.;
				cpot[coord].im = 0.;
				continue;
			  }
			}
// ****************************** DSJ *************************
                     
#ifdef EQUIL_FNL
      /* fnl equilateral */
// ****************************************************************************************** DSJ ********************************************************************************************************
              cpot[coord].re += Fnl * (-3.0 * cpartpot[coord].re - 2.0 * cp1p2p3sym[coord].re / kmag_2_over_3 + 4.0 * cp1p2p3sca[coord].re / kmag_1_over_3 + 2.0 * cp1p2p3nab[coord].re / kmag_2_over_3);
              cpot[coord].im += Fnl * (-3.0 * cpartpot[coord].im - 2.0 * cp1p2p3sym[coord].im / kmag_2_over_3 + 4.0 * cp1p2p3sca[coord].im / kmag_1_over_3 + 2.0 * cp1p2p3nab[coord].im / kmag_2_over_3); 
// ****************************************************************************************** DSJ ********************************************************************************************************
#endif

#ifdef ORTOG_FNL 
// ****************************************************************************************** DSJ ********************************************************************************************************
              cpot[coord].re += Fnl * (-9.0 * cpartpot[coord].re - 8.0 * cp1p2p3sym[coord].re / kmag_2_over_3 + 10.0 * cp1p2p3sca[coord].re / kmag_1_over_3 + 8.0 * cp1p2p3nab[coord].re / kmag_2_over_3);
              cpot[coord].im += Fnl * (-9.0 * cpartpot[coord].im - 8.0 * cp1p2p3sym[coord].im / kmag_2_over_3 + 10.0 * cp1p2p3sca[coord].im / kmag_1_over_3 + 8.0 * cp1p2p3nab[coord].im / kmag_2_over_3);
// ****************************************************************************************** DSJ ********************************************************************************************************
#endif

// ****  wrc ****
#ifdef ORTOG_LSS_FNL
              cpot[coord].re += 3*Fnl * (-(1.0+1.0/3.0*orth_p) * cpartpot[coord].re - 1.0/3.0*(2.0+20./9.*orth_p) * cp1p2p3sym[coord].re / kmag_2_over_3 + 2.0*(1.0-5.0/9.0*orth_p)*(1.0-orth_t) * cp1p2p3sca[coord].re / kmag_1_over_3 + 2.0*(1.0-5.0/9.0*orth_p)*orth_t* cp1p2p3nab[coord].re / kmag_2_over_3);
              cpot[coord].re += 3*Fnl * (orth_p/27.0 * cp1p2p3_K12D[coord].re*kmag_2_over_3 - 20.0*orth_p/27.0 * cp1p2p3_K12E[coord].re / kmag_1_over_3 - 4.0*orth_p/9.0 * cp1p2p3_K12F[coord].re * kmag_1_over_3 + 10.0/9.0*orth_p * cp1p2p3_K12G[coord].re );
              cpot[coord].im += 3*Fnl * (-(1.0+1.0/3.0*orth_p) * cpartpot[coord].im - 1.0/3.0*(2.0+20./9.*orth_p) * cp1p2p3sym[coord].im / kmag_2_over_3 + 2.0*(1.0-5.0/9.0*orth_p)*(1.0-orth_t) * cp1p2p3sca[coord].im / kmag_1_over_3 + 2.0*(1.0-5.0/9.0*orth_p)*orth_t* cp1p2p3nab[coord].im / kmag_2_over_3);
              cpot[coord].im += 3*Fnl * (orth_p/27.0 * cp1p2p3_K12D[coord].im*kmag_2_over_3 - 20.0*orth_p/27.0 * cp1p2p3_K12E[coord].im / kmag_1_over_3 - 4.0*orth_p/9.0 * cp1p2p3_K12F[coord].im * kmag_1_over_3 + 10.0/9.0*orth_p * cp1p2p3_K12G[coord].im );
#endif
// ****  wrc ****
              cpot[coord].re /= (double) nmesh3; 
              cpot[coord].im /= (double) nmesh3; 

            }

      free(cpartpot);
      free(cp1p2p3sym);
      free(cp1p2p3sca);
      free(cp1p2p3nab);
      free(cp1p2p3tre);
// ****  wrc ****
#ifdef ORTOG_LSS_FNL
      free(cp1p2p3inv);
      free(cp1p2p3_K12D);
      free(cp1p2p3_K12E);
      free(cp1p2p3_K12F);
      free(cp1p2p3_K12G);
#endif
// ****  wrc ****


      if (ThisTask == 0 ) print_timed_done(1);
  
#endif

  if (ThisTask == 0) {printf("Computing gradient of non-Gaussian potential..."); fflush(stdout);};

    /****** FINISHED NON LOCAL POTENTIAL OR LOCAL FNL, STILL IN NONGAUSSIAN SECTION ****/   
    /*****  Now 2LPT ****/

  for(axes=0,bytes=0; axes < 3; axes++)
    {
      cdisp[axes] = (fftw_complex *) malloc(bytes += sizeof(fftw_real) * TotalSizePlusAdditional);
      disp[axes] = (fftw_real *) cdisp[axes];
    }

  ASSERT_ALLOC(cdisp[0] && cdisp[1] && cdisp[2]);


#if defined(MULTICOMPONENTGLASSFILE) && defined(DIFFERENT_TRANSFER_FUNC)
  for(Type = MinType; Type <= MaxType; Type++)
#endif
    {

      /* first, clean the array */
      for(i = 0; i < Local_nx; i++)
	for(j = 0; j < Nmesh; j++)
	  for(k = 0; k <= Nmesh / 2; k++)
	    for(axes = 0; axes < 3; axes++)
	      {
		cdisp[axes][(i * Nmesh + j) * (Nmesh / 2 + 1) + k].re = 0;
		cdisp[axes][(i * Nmesh + j) * (Nmesh / 2 + 1) + k].im = 0;
	      }


      for(ii = 0; ii < Local_nx; ii++)
        for(j = 0; j < Nmesh; j++)
          for(k = 0; k <= Nmesh / 2 ; k++)
            {

              coord = (ii * Nmesh + j) * (Nmesh / 2 + 1) + k;
              i = ii + Local_x_start;
    
                   /*   if(i == 0 && j == 0 && k == 0); continue; */

                      if(i < Nmesh / 2)
                        kvec[0] = i * 2 * PI / Box;
                      else
                        kvec[0] = -(Nmesh - i) * 2 * PI / Box;

                      if(j < Nmesh / 2)
                        kvec[1] = j * 2 * PI / Box;
                      else
                        kvec[1] = -(Nmesh - j) * 2 * PI / Box;

                      if(k < Nmesh / 2)
                        kvec[2] = k * 2 * PI / Box;
                      else
                        kvec[2] = -(Nmesh - k) * 2 * PI / Box;

                      kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
                      kmag = sqrt(kmag2);

                     t_of_k = TransferFunc(kmag);

                     twb = t_of_k / Dplus  / Beta;   

                     for(axes = 0; axes < 3; axes++)   
                                        {
                                cdisp[axes][coord].im = kvec[axes] * twb * cpot[coord].re;
                                cdisp[axes][coord].re = - kvec[axes] * twb * cpot[coord].im;
                                        }
            }


      free(cpot);

      if (ThisTask == 0 ) print_timed_done(1);

#endif


       MPI_Barrier(MPI_COMM_WORLD);
 
      /* Compute displacement gradient */

      if (ThisTask == 0) {printf("Computing 2LPT potential..."); fflush(stdout);};

      for(i = 0; i < 6; i++)
	{
	  cdigrad[i] = (fftw_complex *) malloc(bytes = sizeof(fftw_real) * TotalSizePlusAdditional);
	  digrad[i] = (fftw_real *) cdigrad[i];
	  ASSERT_ALLOC(cdigrad[i]);
	}
      
      for(i = 0; i < Local_nx; i++)
	for(j = 0; j < Nmesh; j++)
	  for(k = 0; k <= Nmesh / 2; k++)
	    {
	      coord = (i * Nmesh + j) * (Nmesh / 2 + 1) + k;
	      if((i + Local_x_start) < Nmesh / 2)
		kvec[0] = (i + Local_x_start) * 2 * PI / Box;
	      else
		kvec[0] = -(Nmesh - (i + Local_x_start)) * 2 * PI / Box;
	      
	      if(j < Nmesh / 2)
		kvec[1] = j * 2 * PI / Box;
	      else
		kvec[1] = -(Nmesh - j) * 2 * PI / Box;
	      
	      if(k < Nmesh / 2)
		kvec[2] = k * 2 * PI / Box;
	      else
		kvec[2] = -(Nmesh - k) * 2 * PI / Box;
	      
	      /* Derivatives of ZA displacement  */
	      /* d(dis_i)/d(q_j)  -> sqrt(-1) k_j dis_i */
	      cdigrad[0][coord].re = -cdisp[0][coord].im * kvec[0]; /* disp0,0 */
	      cdigrad[0][coord].im = cdisp[0][coord].re * kvec[0];

	      cdigrad[1][coord].re = -cdisp[0][coord].im * kvec[1]; /* disp0,1 */
	      cdigrad[1][coord].im = cdisp[0][coord].re * kvec[1];

	      cdigrad[2][coord].re = -cdisp[0][coord].im * kvec[2]; /* disp0,2 */
	      cdigrad[2][coord].im = cdisp[0][coord].re * kvec[2];
	      
	      cdigrad[3][coord].re = -cdisp[1][coord].im * kvec[1]; /* disp1,1 */
	      cdigrad[3][coord].im = cdisp[1][coord].re * kvec[1];

	      cdigrad[4][coord].re = -cdisp[1][coord].im * kvec[2]; /* disp1,2 */
	      cdigrad[4][coord].im = cdisp[1][coord].re * kvec[2];

	      cdigrad[5][coord].re = -cdisp[2][coord].im * kvec[2]; /* disp2,2 */
	      cdigrad[5][coord].im = cdisp[2][coord].re * kvec[2];
	    }


      for(i = 0; i < 6; i++) rfftwnd_mpi(Inverse_plan, 1, digrad[i], Workspace, FFTW_NORMAL_ORDER);

      /* Compute second order source and store it in digrad[3]*/

      for(i = 0; i < Local_nx; i++)
	for(j = 0; j < Nmesh; j++)
	  for(k = 0; k < Nmesh; k++)
	    {
	      coord = (i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k;

	      digrad[3][coord] =

		digrad[0][coord]*(digrad[3][coord]+digrad[5][coord])+digrad[3][coord]*digrad[5][coord]
                -digrad[1][coord]*digrad[1][coord]-digrad[2][coord]*digrad[2][coord]-digrad[4][coord]*digrad[4][coord];
	    }

      rfftwnd_mpi(Forward_plan, 1, digrad[3], Workspace, FFTW_NORMAL_ORDER);
      
      /* The memory allocated for cdigrad[0], [1], and [2] will be used for 2nd order displacements */
      /* Freeing the rest. cdigrad[3] still has 2nd order displacement source, free later */

      for(axes = 0; axes < 3; axes++) 
	{
	  cdisp2[axes] = cdigrad[axes]; 
	  disp2[axes] = (fftw_real *) cdisp2[axes];
	}

      free(cdigrad[4]); free(cdigrad[5]); 

      /* Solve Poisson eq. and calculate 2nd order displacements */

      for(i = 0; i < Local_nx; i++)
	for(j = 0; j < Nmesh; j++)
	  for(k = 0; k <= Nmesh / 2; k++)
	    {
	      coord = (i * Nmesh + j) * (Nmesh / 2 + 1) + k;
	      if((i + Local_x_start) < Nmesh / 2)
		kvec[0] = (i + Local_x_start) * 2 * PI / Box;
	      else
		kvec[0] = -(Nmesh - (i + Local_x_start)) * 2 * PI / Box;
	      
	      if(j < Nmesh / 2)
		kvec[1] = j * 2 * PI / Box;
	      else
		kvec[1] = -(Nmesh - j) * 2 * PI / Box;
	      
	      if(k < Nmesh / 2)
		kvec[2] = k * 2 * PI / Box;
	      else
		kvec[2] = -(Nmesh - k) * 2 * PI / Box;

	      kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
#ifdef CORRECT_CIC
	      /* calculate smooth factor for deconvolution of CIC interpolation */
	      fx = fy = fz = 1;
	      if(kvec[0] != 0)
		{
		  fx = (kvec[0] * Box / 2) / Nmesh;
		  fx = sin(fx) / fx;
		}
	      if(kvec[1] != 0)
		{
		  fy = (kvec[1] * Box / 2) / Nmesh;
		  fy = sin(fy) / fy;
		}
	      if(kvec[2] != 0)
		{
		  fz = (kvec[2] * Box / 2) / Nmesh;
		  fz = sin(fz) / fz;
		}
	      ff = 1 / (fx * fy * fz);
	      smth = ff * ff;
	      /*  */
#endif

	      /* cdisp2 = source * k / (sqrt(-1) k^2) */
	      for(axes = 0; axes < 3; axes++)
		{
		  if(kmag2 > 0.0) 
		    {
		      cdisp2[axes][coord].re = cdigrad[3][coord].im * kvec[axes] / kmag2;
		      cdisp2[axes][coord].im = -cdigrad[3][coord].re * kvec[axes] / kmag2;
		    }
		  else cdisp2[axes][coord].re = cdisp2[axes][coord].im = 0.0;
#ifdef CORRECT_CIC
		  cdisp[axes][coord].re *= smth;   cdisp[axes][coord].im *= smth;
		  cdisp2[axes][coord].re *= smth;  cdisp2[axes][coord].im *= smth;
#endif
		}
	    }
      
      /* Free cdigrad[3] */
      free(cdigrad[3]);

      MPI_Barrier(MPI_COMM_WORLD);

      /* Now, both cdisp, and cdisp2 have the ZA and 2nd order displacements */

      for(axes = 0; axes < 3; axes++)
	{
	  rfftwnd_mpi(Inverse_plan, 1, disp[axes], Workspace, FFTW_NORMAL_ORDER);
	  rfftwnd_mpi(Inverse_plan, 1, disp2[axes], Workspace, FFTW_NORMAL_ORDER);

	  /* now get the plane on the right side from neighbour on the right, 
	     and send the left plane */
      
	  recvTask = ThisTask;
	  do
	    {
	      recvTask--;
	      if(recvTask < 0)
		recvTask = NTask - 1;
	    }
	  while(Local_nx_table[recvTask] == 0);
      
	  sendTask = ThisTask;
	  do
	    {
	      sendTask++;
	      if(sendTask >= NTask)
		sendTask = 0;
	    }
	  while(Local_nx_table[sendTask] == 0);
      
	  /* use non-blocking send */
      
	  if(Local_nx > 0)
	    {
	      /* send ZA disp */
	      MPI_Isend(&(disp[axes][0]),
			sizeof(fftw_real) * Nmesh * (2 * (Nmesh / 2 + 1)),
			MPI_BYTE, recvTask, 10, MPI_COMM_WORLD, &request);
	      
	      MPI_Recv(&(disp[axes][(Local_nx * Nmesh) * (2 * (Nmesh / 2 + 1))]),
		       sizeof(fftw_real) * Nmesh * (2 * (Nmesh / 2 + 1)),
		       MPI_BYTE, sendTask, 10, MPI_COMM_WORLD, &status);
	      
	      MPI_Wait(&request, &status);

	      
	      /* send 2nd order disp */
	      MPI_Isend(&(disp2[axes][0]),
			sizeof(fftw_real) * Nmesh * (2 * (Nmesh / 2 + 1)),
			MPI_BYTE, recvTask, 10, MPI_COMM_WORLD, &request);
	      
	      MPI_Recv(&(disp2[axes][(Local_nx * Nmesh) * (2 * (Nmesh / 2 + 1))]),
		       sizeof(fftw_real) * Nmesh * (2 * (Nmesh / 2 + 1)),
		       MPI_BYTE, sendTask, 10, MPI_COMM_WORLD, &status);
	      
	      MPI_Wait(&request, &status);
	    }
	}

      if (ThisTask == 0 ) print_timed_done(21);
      if (ThisTask == 0) {printf("Computing displacements and velocitites..."); fflush(stdout);};

      /* read-out displacements */

      nmesh3 = Nmesh * Nmesh * Nmesh;
      
      for(n = 0; n < NumPart; n++)
	{
#if defined(MULTICOMPONENTGLASSFILE) && defined(DIFFERENT_TRANSFER_FUNC)
	  if(P[n].Type == Type)
#endif
	    {
	      u = P[n].Pos[0] / Box * Nmesh;
	      v = P[n].Pos[1] / Box * Nmesh;
	      w = P[n].Pos[2] / Box * Nmesh;
	      
	      i = (int) u;
	      j = (int) v;
	      k = (int) w;
	      
	      if(i == (Local_x_start + Local_nx))
		i = (Local_x_start + Local_nx) - 1;
	      if(i < Local_x_start)
		i = Local_x_start;
	      if(j == Nmesh)
		j = Nmesh - 1;
	      if(k == Nmesh)
		k = Nmesh - 1;
	      
	      u -= i;
	      v -= j;
	      w -= k;
	      
	      i -= Local_x_start;
	      ii = i + 1;
	      jj = j + 1;
	      kk = k + 1;
	      
	      if(jj >= Nmesh)
		jj -= Nmesh;
	      if(kk >= Nmesh)
		kk -= Nmesh;
	      
	      f1 = (1 - u) * (1 - v) * (1 - w);
	      f2 = (1 - u) * (1 - v) * (w);
	      f3 = (1 - u) * (v) * (1 - w);
	      f4 = (1 - u) * (v) * (w);
	      f5 = (u) * (1 - v) * (1 - w);
	      f6 = (u) * (1 - v) * (w); 
	      f7 = (u) * (v) * (1 - w);
	      f8 = (u) * (v) * (w);
	     
	      for(axes = 0; axes < 3; axes++)
		{
		  dis = disp[axes][(i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k] * f1 +
		    disp[axes][(i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + kk] * f2 +
		    disp[axes][(i * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + k] * f3 +
		    disp[axes][(i * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + kk] * f4 +
		    disp[axes][(ii * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k] * f5 +
		    disp[axes][(ii * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + kk] * f6 +
		    disp[axes][(ii * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + k] * f7 +
		    disp[axes][(ii * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + kk] * f8;

		  dis2 = disp2[axes][(i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k] * f1 +
		    disp2[axes][(i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + kk] * f2 +
		    disp2[axes][(i * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + k] * f3 +
		    disp2[axes][(i * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + kk] * f4 +
		    disp2[axes][(ii * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k] * f5 +
		    disp2[axes][(ii * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + kk] * f6 +
		    disp2[axes][(ii * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + k] * f7 +
		    disp2[axes][(ii * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + kk] * f8;
		  dis2 /= (float) nmesh3;
	      
		  
#ifdef ONLY_ZA
		  P[n].Pos[axes] += dis;
		  P[n].Vel[axes] = dis * vel_prefac;
#else
		  P[n].Pos[axes] += dis - 3./7. * dis2;
		  P[n].Vel[axes] = dis * vel_prefac - 3./7. * dis2 * vel_prefac2;
#endif

		  P[n].Pos[axes] = periodic_wrap(P[n].Pos[axes]);

		  if(fabs(dis - 3./7. * dis2 > maxdisp))
		    maxdisp = fabs(dis - 3./7. * dis2);
		}
	    }
	}
    }
 
      if (ThisTask == 0 ) print_timed_done(6);

  for(axes = 0; axes < 3; axes++) free(cdisp[axes]);
  for(axes = 0; axes < 3; axes++) free(cdisp2[axes]);

  gsl_rng_free(random_generator);

  MPI_Reduce(&maxdisp, &max_disp_glob, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

/*  if(ThisTask == 0)
    {
      printf("\nMaximum displacement (1D): %g kpc/h, in units of the part-spacing= %g\n",
	     max_disp_glob, max_disp_glob / (Box / Nmesh));
    }*/
}

double periodic_wrap(double x)
{
  while(x >= Box)
    x -= Box;

  while(x < 0)
    x += Box;

  return x;
}


void set_units(void)		/* ... set some units */
{
  UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;
  G = GRAVITY / pow(UnitLength_in_cm, 3) * UnitMass_in_g * pow(UnitTime_in_s, 2);
  Hubble = HUBBLE * UnitTime_in_s;
}



void initialize_ffts(void)
{
  int total_size, i, additional;
  int local_ny_after_transpose, local_y_start_after_transpose;
  int *slab_to_task_local;
  size_t bytes;


  Inverse_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD,
					 Nmesh, Nmesh, Nmesh, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);

  Forward_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD,
					 Nmesh, Nmesh, Nmesh, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);

  rfftwnd_mpi_local_sizes(Forward_plan, &Local_nx, &Local_x_start,
			  &local_ny_after_transpose, &local_y_start_after_transpose, &total_size);

  Local_nx_table = malloc(sizeof(int) * NTask);
  MPI_Allgather(&Local_nx, 1, MPI_INT, Local_nx_table, 1, MPI_INT, MPI_COMM_WORLD);

  Slab_to_task = malloc(sizeof(int) * Nmesh);
  slab_to_task_local = malloc(sizeof(int) * Nmesh);

  for(i = 0; i < Nmesh; i++)
    slab_to_task_local[i] = 0;

  for(i = 0; i < Local_nx; i++)
    slab_to_task_local[Local_x_start + i] = ThisTask;

  MPI_Allreduce(slab_to_task_local, Slab_to_task, Nmesh, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  free(slab_to_task_local);

  additional = (Nmesh) * (2 * (Nmesh / 2 + 1));	/* additional plane on the right side */

  TotalSizePlusAdditional = total_size + additional;

  Workspace = (fftw_real *) malloc(bytes = sizeof(fftw_real) * total_size);

  ASSERT_ALLOC(Workspace)

}



void free_ffts(void)
{
  free(Workspace);
  free(Slab_to_task);
  rfftwnd_mpi_destroy_plan(Inverse_plan);
  rfftwnd_mpi_destroy_plan(Forward_plan);
}


int FatalError(int errnum)
{
  printf("FatalError called with number=%d\n", errnum);
  fflush(stdout);
  MPI_Abort(MPI_COMM_WORLD, errnum);
  exit(0);
}




static double A, B, alpha, beta, V, gf;

double fnl(double x)		/* Peacock & Dodds formula */
{
  return x * pow((1 + B * beta * x + pow(A * x, alpha * beta)) /
		 (1 + pow(pow(A * x, alpha) * gf * gf * gf / (V * sqrt(x)), beta)), 1 / beta);
}

void print_spec(void)
{
  double k, knl, po, dl, dnl, neff, kf, kstart, kend, po2, po1, DDD;
  char buf[1000];
  FILE *fd;

  if(ThisTask == 0)
    {
      sprintf(buf, "%s/inputspec_%s.txt", OutputDir, FileBase);

      fd = fopen(buf, "w");

      gf = GrowthFactor(0.001, 1.0) / (1.0 / 0.001);

      DDD = GrowthFactor(1.0 / (Redshift + 1), 1.0);

      fprintf(fd, "%12g %12g\n", Redshift, DDD);	/* print actual starting redshift and 
							   linear growth factor for this cosmology */

      kstart = 2 * PI / (1000.0 * (3.085678e24 / UnitLength_in_cm));	/* 1000 Mpc/h */
      kend = 2 * PI / (0.001 * (3.085678e24 / UnitLength_in_cm));	/* 0.001 Mpc/h */

      for(k = kstart; k < kend; k *= 1.025)
	{
	  po = PowerSpec(k);
	  dl = 4.0 * PI * k * k * k * po;

	  kf = 0.5;

	  po2 = PowerSpec(1.001 * k * kf);
	  po1 = PowerSpec(k * kf);

	  if(po != 0 && po1 != 0 && po2 != 0)
	    {
	      neff = (log(po2) - log(po1)) / (log(1.001 * k * kf) - log(k * kf));

	      if(1 + neff / 3 > 0)
		{
		  A = 0.482 * pow(1 + neff / 3, -0.947);
		  B = 0.226 * pow(1 + neff / 3, -1.778);
		  alpha = 3.310 * pow(1 + neff / 3, -0.244);
		  beta = 0.862 * pow(1 + neff / 3, -0.287);
		  V = 11.55 * pow(1 + neff / 3, -0.423) * 1.2;

		  dnl = fnl(dl);
		  knl = k * pow(1 + dnl, 1.0 / 3);
		}
	      else
		{
		  dnl = 0;
		  knl = 0;
		}
	    }
	  else
	    {
	      dnl = 0;
	      knl = 0;
	    }

	  fprintf(fd, "%12g %12g    %12g %12g\n", k, dl, knl, dnl);
	}
      fclose(fd);
    }
}
