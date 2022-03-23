EXEC   = 2LPTnonlocal

OBJS   = main.o power.o checkchoose.o allvars.o save.o read_param.o  read_glass.o  \
         nrsrc/nrutil.o nrsrc/qromb.o nrsrc/polint.o nrsrc/trapzd.o

INCL   = allvars.h proto.h  nrsrc/nrutil.h  Makefile



#OPT   +=  -DPRODUCEGAS   # Set this to automatically produce gas particles 
                          # for a single DM species in the input file by interleaved by a half a grid spacing
                          # only for Gaussian and ZA

OPT   +=  -DMULTICOMPONENTGLASSFILE  # set this if the initial glass file contains multiple components
                                      # only for Gaussian and ZA

#OPT   +=  -DDIFFERENT_TRANSFER_FUNC  # set this if you want to implement a transfer function that depends on
                                      # particle type
                                      # only for Gaussian and ZA

OPT   +=  -DNO64BITID     # switch this on if you want normal 32-bit IDs

#OPT  +=  -DCORRECT_CIC  # only switch this on if particles start from a glass (as opposed to grid)
                         # only for Gaussian and ZA

#OPT += -DONLY_ZA    # swith this on if you want ZA initial conditions (2LPT otherwise)

#MODE = -DONLY_GAUSSIAN
#MODE = -DLOCAL_FNL
MODE = -DEQUIL_FNL

MODE = -DORTOG_FNL
MODE = -DORTOG_LSS_FNL


ifeq ($(MODE),-DONLY_GAUSSIAN)
	EXEC:=2LPT
else ifeq ($(MODE),-DLOCAL_FNL)
	EXEC:=2LPTNGLC
else ifeq ($(MODE),-DEQUIL_FNL)
	EXEC:=2LPTNGEQ
else ifeq ($(MODE),-DORTOG_FNL)
	EXEC:=2LPTNGOR
else ifeq ($(MODE),-DORTOG_LSS_FNL)
	EXEC:=2LPTNGORLSS
endif
OPT += $(MODE)
OPTIONS =  $(OPT)

# ----- Chichipio/astro values ---- 
#
GSL_INCL= -I/usr/local/gsl_gcc_mpi/include
GSL_LIBS= -L/usr/local/gsl_gcc_mpi/lib
FFTW_INCL= -I/mnt/home/ajamieson/ceph/Software/Source/Test_Roman_2LPT_PNG/fftw2/include
FFTW_LIBS= -L/mnt/home/ajamieson/ceph/Software/Source/Test_Roman_2LPT_PNG/fftw2/lib

CC       = mpicc #-g -Wall -fbounds-check    # sets the C-compiler (default)
MPICHLIB = -L/usr/local/mpich_gcc/lib

#CC       =  /usr/local/mpich_intel/bin/mpicc -g -Wall -fbounds-check    # sets the C-compiler (default)

# -----------------------------------

OPTIMIZE =   -O3 -Wall    # optimization and warning flags (default)

SYSTYPE="flatiron"
#SYSTYPE="GordonS"

ifeq ($(SYSTYPE),"flatiron")
CC       =   mpicc     # sets the C-compiler   
OPT      +=  -DMPICH_IGNORE_CXX_SEEK
OPTIMIZE =   -std=gnu99 -O3 -g -Wall -Wno-unused-but-set-variable -Wno-uninitialized -Wno-unknown-pragmas -Wno-unused-function -march=native
GSL_INCL =  -I$GSL_ROOT/include
GSL_LIBS =  
FFTW_INCL=  -I$FFTW_ROOT/include 
FFTW_LIBS=  -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
MPICHLIB =  -lmpi
HDF5INCL =  -DH5_USE_16_API
HDF5LIB  =  -lhdf5 -lz
ifeq (NUM_THREADS,$(findstring NUM_THREADS,$(CONFIGVARS)))
OPTIMIZE +=  -fopenmp
OPT      += -DIMPOSE_PINNING -DSOCKETS=4 -DMAX_CORES=16
endif
endif

ifeq ($(SYSTYPE),"GordonS")
CC = mpicc -g -O2 #-xW -ipo -Wall
CXX = mpiCC -g -O2 -xW -ipo -Wall
OPTIMIZE =
GMP_INCL = -I/opt/gnu/gmp/include
GMP_LIBS = -L/opt/gnu/gmp/lib
GSL_INCL = -I/opt/gsl/2.1/intel/include
GSL_LIBS = -L/opt/gsl/2.1/intel/lib
FFTW_INCL= -I/opt/fftw/2.1.5/intel/mvapich2_ib/include
FFTW_LIBS= -L/opt/fftw/2.1.5/intel/mvapich2_ib/lib
HDF5INCL = -I/opt/hdf5/intel/mvapich2_ib/include
HDF5LIB = -L/opt/hdf5/intel/mvapich2_ib/lib -lhdf5 -lz
endif


#FFTW_LIB =  $(FFTW_LIBS) -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
FFTW_LIB =  $(FFTW_LIBS)

LIBS   =   -lm  $(MPICHLIB)  $(FFTW_LIB)  $(GSL_LIBS)  -lgsl -lgslcblas

CFLAGS =   $(OPTIONS)  $(OPTIMIZE)  $(FFTW_INCL) $(GSL_INCL)

$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 


.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)



