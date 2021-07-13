#include "allvars.h"

clock_t start_time;
clock_t previous_time;

struct io_header_1 header1, header;

int WhichSpectrum;
int WhichTransfer;  

int SphereMode;
int *Local_nx_table;

FILE *FdTmp, *FdTmpInput;

int Nmesh, Nsample;

long long IDStart;

char GlassFile[500];
char FileWithInputSpectrum[500];
char FileWithInputTransfer[500]; 

int GlassTileFac;

double Box;
int Seed; 

long long TotNumPart;

int NumPart;

int *Slab_to_task;

int NTaskWithN;

struct part_data *P;

int Nglass;

double InitTime;
double Redshift;
double MassTable[6];
double Fnl; 
// *** FAVN/DSJ ***
int FixedAmplitude;
int PhaseFlip;
// *** FAVN/DSJ ***

char OutputDir[100], FileBase[100];
int NumFilesWrittenInParallel;


int ThisTask, NTask;

int Local_nx, Local_x_start;

int IdStart;

rfftwnd_mpi_plan Inverse_plan;
rfftwnd_mpi_plan Forward_plan;
unsigned int TotalSizePlusAdditional;
fftw_real *Workspace;


double UnitTime_in_s, UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
double InputSpectrum_UnitLength_in_cm;
double G, Hubble;
double RhoCrit;

double Omega, OmegaLambda, OmegaDM_2ndSpecies, Sigma8;
double OmegaBaryon, HubbleParam;
double ShapeGamma;
double PrimordialIndex;
double Anorm;
double Dplus;			/* growth factor at output=InitTime and growth factor for initial potential */

#ifdef DIFFERENT_TRANSFER_FUNC
int Type, MinType, MaxType;
#endif

int WDM_On;
int WDM_Vtherm_On;
double WDM_PartMass_in_kev;
