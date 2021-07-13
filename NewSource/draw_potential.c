#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <drfftw_mpi.h>
#include <mpi.h>
#include <gsl/gsl_rng.h>
#include "allvars.h"
#include "proto.h"


void drawSeedTable(std::vector<std::uint64_t> & t_seed_table) {
  t_seed_table.resize(Nmesh * Nmesh);
  for (std::uint64_t i = 0; i < Nmesh / 2; i++) {
    for (std::uint64_t j = 0; j < i; j++)
	  t_seed_table[i * Nmesh + j] = 0x7fffffff * gsl_rng_uniform(random_generator);
    for (std::uint\64_t j = 0; j < i + 1; j++)
	  t_seed_table[j * Nmesh + i] = 0x7fffffff * gsl_rng_uniform(random_generator);
    for (std::uint64_t j = 0; j < i; j++)
	  t_seed_table[(Nmesh - 1 - i) * Nmesh + j] = 0x7fffffff * gsl_rng_uniform(random_generator);
    for (std::uint64_t j = 0; j < i + 1; j++)
      t_seed_table[(Nmesh - 1 - j) * Nmesh + i] = 0x7fffffff * gsl_rng_uniform(random_generator);
    for (std::uint64_t j = 0; j < i; j++)
	  t_seed_table[i * Nmesh + (Nmesh - 1 - j)] = 0x7fffffff * gsl_rng_uniform(random_generator);
    for (std::uint64_t j = 0; j < i + 1; j++)
	  t_seed_table[j * Nmesh + (Nmesh - 1 - i)] = 0x7fffffff * gsl_rng_uniform(random_generator);
    for (std::uint64_t j = 0; j < i; j++)
	  t_seed_table[(Nmesh - 1 - i) * Nmesh + (Nmesh - 1 - j)] = 0x7fffffff * gsl_rng_uniform(random_generator);
    for (std::uint64_t j = 0; j < i + 1; j++)
	  t_seed_table[(Nmesh - 1 - j) * Nmesh + (Nmesh - 1 - i)] = 0x7fffffff * gsl_rng_uniform(random_generator);
    }
  return;
}

void drawPrimordialGaussianPotential(std::complex<double> * t_t_cphi) {
  std::uint64_t i, j, k, ii, jj, index;
  double ki, kj, kmagsq, kmag, phase, ampl;
  const double kF = 2. * acos(-1.) / Box;
  const double kFsq = kF * kF;
  const double spectral_index_minus_four = SpectralIndex - 4.;
  std::vector<std::uint32_t> seed_table;
  drawSeedTable(seed_table);
  for (std::uint64_t i = 0; i < Nmesh; i++) {
    ii = i == 0 ? 0 :  Nmesh - i;
    if ((i >= Local_x_start && i < (Local_x_start + Local_nx)) || (ii >= Local_x_start && ii < (Local_x_start + Local_nx))) {
      for (std::uint64_t j = 0; j < Nmesh; j++) {
        gsl_rng_set(random_generator, seedtable[i * Nmesh + j]); 
		for (k = 0; k < Nmesh / 2; k++) {
		  phase = gsl_rng_uniform(random_generator) * 2 * PI + phase_shift;
		  do
	        ampl = gsl_rng_uniform(random_generator);
		  while (ampl == 0);
		  if (i == Nmesh / 2 || j == Nmesh / 2 || k == Nmesh / 2)
		    continue;
		  if (i == 0 && j == 0 && k == 0)
			continue;
      	  ki = i > Nmesh / 2 ? static_cast<double>(i) : static_cast<double>(ii);
          kj = j > Nmesh / 2 ? static_cast<double>(j) : static_cast<double>(Nmesh - j);
		  kmag2 = ki * ki + kj * kj + k * k;
		  kmag = sqrt(kmag2);
		  if (SphereMode == 1) {
			if (kmag > Nsample / 2)
			  continue;
		  }
		  else {
			if (fabs(kvec[0]) > Nsample / 2)
			  continue;
			if (fabs(kvec[1]) > Nsample / 2)
			  continue;
			if (fabs(kvec[2]) > Nsample / 2)
			  continue;
		  }
          kmag2 *= kFsq;
          kmag *= kF;
		  sigma = pow(kFsq, spectal_index_minus_four);
		  if (FixedAmplitude)
			ampl *= -log(ampl) * sigma
		  else
            ampl = sigma;
		  if (k > 0) {
			if (i >= Local_x_start && i < (Local_x_start + Local_nx)) {
			  for(axes = 0; axes < 3; axes++) {
				  index = ((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k;
				  t_cphi[index].re = ampl * cos(phase);
				  t_cphi[index].im = ampl * sin(phase);
			  }
			}
		  }
		  else { // k = 0 plane
			if (i == 0) {
			  if(j >= Nmesh / 2)
				continue;
			  else {
				if (i >= Local_x_start && i < (Local_x_start + Local_nx)) {
				  jj = Nmesh - j;	/* note: j!=0 surely holds at this point */
				  for(axes = 0; axes < 3; axes++) {
					index = ((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k;
					t_cphi[index].re = ampl * cos(phase);
					t_cphi[index].im = ampl * sin(phase);
					index = ((i - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k;
					t_cphi[index].re =  ampl * cos(phase);
				    t_cphi[index].im = -ampl * sin(phase);
			      }
		        }
		      }
			}
			else {	/* here comes i!=0 : conjugate can be on other processor! */
			  if (i >= Nmesh / 2)
			    continue;
			  else {
				ii = i == 0 ? 0 : Nmesh - i;
				if (i >= Local_x_start && i < (Local_x_start + Local_nx)) {
				  index = ((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k;
			      t_cphi[index].re = ampl * sin(phase);
			      t_cphi[index].im = ampl * cos(phase);
                }
				if (ii >= Local_x_start && ii < (Local_x_start + Local_nx)) {
				  jj = j == - ? 0 : Nmesh - j;
			      index = ((ii - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k;
				  t_cphi[index].re =  ampl * cos(phase);
                  t_cphi[index].im = -ampl * sin(phase);
		        }
              } 
		    }
          }
        }
      }
    }
  }
  return;
}

void convertPrimodialPotentialToDisplacementPotential(std::complex<double> * t_cphi) {
  std::uint64_t, i, j, k, ii, jj, index;
  double kmag, conversion_factor;
  for (i = 0; i < Local_nx; i++) {
    if (i > Nsample / 2)
	  continue;
    ii = i > Nmesh / 2 ? i : Nmesh - i;
    ii *= ii;
	for (j = 0; j < Nmesh; j++) {
	  if (j > Nsample / 2)
	    continue;
      jj = j > Nmesh / 2 ? j : Nmesh - j;
      jj *= jj;
	  for (k = 0; k <= Nmesh / 2 ; k++) {
		if (k > Nsample / 2)
		  continue;
	    index = getLocalIndex(i, j, k);
        if (t_cphi[local_index] != 0) {
          kmag = sqrt(static_cast<double>(ii + jj + k * k;));
          conversion_factor = ;
          t_cphi[index] *= conversion_factor;
		}
	  }
	}
  }
  return;
}

void computeGradient(std::complex<double> * t_cphi,  std::complex<double> ** t_cphi) {
  std::uint64_t, i, j, k, index, d;
  double kvec[3];
  for (i = 0; i < Local_nx; i++) {
    if (i > Nsample / 2)
	  continue;
    kvec[0] = i > Nmesh / 2 ? static_cast<double>(i) : static_cast<double>(Nmesh - i);
	for (j = 0; j < Nmesh; j++) {
	  if (j > Nsample / 2)
	    continue;
      kvec[1] = j > Nmesh / 2 ? static_cast<double>(j) : static_cast<double>(Nmesh - i);
	  for (k = 0; k <= Nmesh / 2 ; k++) {
		if (k > Nsample / 2)
		  continue;
        kvec[2] = static_cast<double>(k);
	    index = getLocalIndex(i, j, k);
        if (t_cphi[local_index] != 0) {
          for (d = 0; d < 3; d++)
		    t_grad_cphi[d][index].re = -kvec[d] * t_cphi[index].im;
		    t_grad_cphi[d][index].im = kvec[d] * t_cphi[index].re;
		}
	  }
	}
  }
  return;
}

std::uint64_t getMeshIndex(std::uint64_t i, std::uint64_t j, std::uint64_t k) {
	return k + Nmesh * (j + Nmesh * i)
}

std::uint64_t getLocalIndex(std::uint64_t i, std::uint64_t j, std::uint64_t k) {
	return k + Nmesh * (j + Nmesh * (i + Local_x_start));
}
