#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <drfftw_mpi.h>
#include <mpi.h>
#include <gsl/gsl_rng.h>

#include "allvars.h"
#include "proto.h"

int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
	MPI_Comm_size(MPI_COMM_WORLD, &NTask);
	if(argc != 2) {
		if(ThisTask == 0) {
			fprintf(stdout, "\nParameters are missing.\n");
			fprintf(stdout, "Call with <ParameterFile>\n\n");
		}
		MPI_Finalize();
		exit(0);
		}
	}
	read_parameterfile(argv[1]);
	checkchoose();	
	set_units();
	initialize_transferfunction(); 
	initialize_powerspectrum(); 
	initialize_ffts(); 
	read_glass(GlassFile);
	displacement_fields();
	write_particle_data();
	if(NumPart)
		free(P);
	free_ffts();
	if(ThisTask == 0) {
		printf("\nIC's generated.\n\n");
		printf("Initial scale factor = %g\n", InitTime);
		printf("\n");
	}
	MPI_Barrier(MPI_COMM_WORLD);
	print_spec();
	MPI_Finalize();		/* clean up & finalize MPI */
	exit(0);
}
