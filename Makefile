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

#MODE = -DONLY_GAUSSIAN # switch this if you want gaussian initial conditions (fnl otherwise) 
#MODE = -DLOCAL_FNL #switch this if you want only local non-gaussianities
#MODE = -DEQUIL_FNL #switch this if you want equilateral Fnl
MODE = -DORTOG_FNL #switch this if you want ortogonal Fnl
ifeq  ($(MODE), -DONLY_GAUSSIAN))
	EXEC = 2LPT
else ifeq ($(MODE), -DLOCAL_FNL)
	EXEC = 2LPTNGLC
else ifeq ($(MODE), -DEQUIL_FNL)
	EXEC = 2LPTNGEQ
else ifeq ($(MODE), -DORTOG_FNL)
	EXEC = 2LPTNGOR
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


#FFTW_LIB =  $(FFTW_LIBS) -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
FFTW_LIB =  $(FFTW_LIBS) -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw

LIBS   =   -lm  $(MPICHLIB)  $(FFTW_LIB)  $(GSL_LIBS)  -lgsl -lgslcblas

CFLAGS =   $(OPTIONS)  $(OPTIMIZE)  $(FFTW_INCL) $(GSL_INCL)

$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 


.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)



