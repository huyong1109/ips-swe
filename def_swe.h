#ifndef SWE_headfile

#define SWE_headfile

/* fixed parameters */
#define DOF        3

#define pi (PETSC_PI)
#define GRAVITY  (9.80616)
#define CORIOLIS (7.292E-5)


#define MESHSIZE  12
#define WIDTH     2
#define XLOW   (-250.)
#define XHI    (250.)
#define YLOW   (-250.)
#define YHI    (250.)
#define EPS    (1.0E-13)

#define max(a,b) ((a)>(b) ? (a) : (b))
#define min(a,b) ((a)>(b) ? (b) : (a))

typedef struct
{
	PetscInt       tscurr;            // the count number of the current time step
	PetscInt       tsmax;             // the total number of time steps we wish to take
	PassiveScalar  tsize;             // the size of (\Delta t)
	PassiveScalar  tstart;            // start time
	PassiveScalar  tfinal;            // final time
	PassiveScalar  tcurr;             // the current time accumulation (in Alfven units)
	PetscInt       torder;            // order of temporal discretization
} TstepCtx; 

typedef struct { PetscScalar eta, U,  V; }   ActiveField;

typedef struct { PetscScalar eta, U,  V; }   BdyField;

typedef struct
{
	BdyField       *L_o, *R_o, *B_o, *T_o;
	BdyField       *L_i, *R_i, *B_i, *T_i;
	PetscBool      onbdy, on_L, on_R, on_B, on_T;
	PetscInt       start, end, length;
} BdyCtx;

typedef struct {
	PassiveScalar  H; 
} MetricField;


typedef struct {
	DM            da;
	Vec           tensor;
} MetricCtx;


typedef struct
{
	MPI_Comm      comm, mycomm;
	PetscMPIInt   rank, size, myrank, mysize, myid;
	DM            da;
} SWE;

typedef struct
{
	PetscLogEvent func;
	PetscLogEvent func_apply,  mpi;
} EventCtx;

typedef struct
{
	PassiveScalar gravity, coriolis, PreLoading;
} ParamCtx; // useful parameters, including physical parameters


typedef struct {
	SWE	       *swe; 
	TstepCtx       *tsctx;
	ParamCtx       *param;
	MetricCtx      *metric;
	BdyCtx         *bdy;
	EventCtx       *event;
	MPI_Comm       comm;
	PetscMPIInt    rank, size;
	PetscInt       meshsize;
	PetscInt       px,py,*lx,*ly;
	Vec	       sol;
	PassiveScalar  dx,dy;
	PetscViewer    viewer;
} UserCtx;

EXTERN_C_BEGIN

extern PetscErrorCode SWEInitialize( UserCtx* user);
extern PetscErrorCode SWEFinalize( UserCtx* user );
extern PetscErrorCode Update(UserCtx *user);
extern PetscErrorCode FormInitialValue(UserCtx *user);

EXTERN_C_END

#endif

