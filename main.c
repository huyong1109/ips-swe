static char help[] = "Solving shallow water Equations using Parallel Splitting Method, ny Yong Hu.\n\n";
/* Part of the code copied from Chao Yang. */ 

#include <petscsnes.h>
#include <petscdmda.h>
#include "def_swe.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
	PetscErrorCode ierr;
	UserCtx        user[1];
	MPI_Comm       comm;
	PetscMPIInt    rank, size;
	SWE            swe[1];
	PetscInt       meshsize;
	ParamCtx       param;
	TstepCtx       tsctx;
	MetricCtx      metric[1];
	BdyCtx         bdy[1];
	EventCtx       event[1];
	PetscViewer    viewer;

	PetscInitialize(&argc,&argv,(char *)0,help);

	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n SET UP \n"); CHKERRQ(ierr);
	PetscPreLoadBegin(PETSC_TRUE,"SetUp");

	param.PreLoading = PetscPreLoading;

	comm = PETSC_COMM_WORLD;
	ierr = MPI_Comm_rank( comm, &rank ); CHKERRQ(ierr);
	ierr = MPI_Comm_size( comm, &size ); CHKERRQ(ierr);

	ierr = PetscViewerCreate(comm, &viewer); CHKERRQ(ierr);
	ierr = PetscViewerSetType(viewer, PETSCVIEWERASCII); CHKERRQ(ierr);
	ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB); CHKERRQ(ierr);

	/* Register a user-defined event for profiling (error checking). */
	PetscPreLoadStage("Initiate");

	/* Problem parameters */
	param.gravity  = GRAVITY;
	param.coriolis = CORIOLIS;

	/* Time-stepping context */
	tsctx.tstart    = 0.0;
	tsctx.tfinal    = 1.0;
	tsctx.tsize     = 0.05;
	tsctx.tsmax     = 5;

	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n Get Options \n"); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(PETSC_NULL,"-px",&(user->px),PETSC_NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(PETSC_NULL,"-py",&(user->py),PETSC_NULL); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n Create DMDA 2D [x,y] = [%d, %d]\n", user->px, user->py); CHKERRQ(ierr);
	ierr = PetscOptionsGetScalar(PETSC_NULL,"-tfinal",&tsctx.tfinal,PETSC_NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetScalar(PETSC_NULL,"-tsize",&tsctx.tsize,PETSC_NULL); CHKERRQ(ierr);
	tsctx.tcurr = tsctx.tstart;
	ierr = PetscOptionsGetInt(PETSC_NULL,"-tsmax",&tsctx.tsmax,PETSC_NULL); CHKERRQ(ierr);

	meshsize = MESHSIZE;
	ierr = PetscOptionsGetInt(PETSC_NULL,"-meshsize",&meshsize,PETSC_NULL); CHKERRQ(ierr);

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   Create user context, set problem data, create vector data structures.
	   Also, compute the initial guess.
	   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	user->viewer   = viewer;
	user->sol      = PETSC_NULL;
	user->lx       = PETSC_NULL;
	user->ly       = PETSC_NULL;
	user->comm     = comm;
	user->rank     = rank;
	user->size     = size;
	user->tsctx    = &tsctx;
	user->event    = event;
	user->param    = &param;
	user->meshsize = meshsize;
	user->dx       = (XHI - XLOW) / (PassiveScalar)(user->meshsize);
	user->dy       = user->dx;
	metric->tensor = PETSC_NULL;
	user->metric   = metric;
	user->bdy      = bdy;
	user->swe      = swe;
	
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n SWE Initialize \n"); CHKERRQ(ierr);
	ierr = SWEInitialize(user); CHKERRQ(ierr);
	
	CHKMEMQ;
	
	ierr = DMGetGlobalVector(user->swe->da, &(user->sol)); CHKERRQ(ierr);
	ierr = DMGetGlobalVector(user->swe->da, &(user->metric->tensor)); CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n PreLoading %G \n", param.PreLoading); CHKERRQ(ierr);
	if (param.PreLoading) {
		ierr = PetscPrintf(PETSC_COMM_WORLD, "\n+++++++++++++++++++++++ Problem parameters +++++++++++++++++++++\n"); CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD," Example: Isolated Mountain.\n"); CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,
		                   " Geometric parameters: Gravity = %G, Coriolis = %G.\n",
		                   param.gravity, param.coriolis); CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD," Problem size %D X %D (six patches), Ncpu = %D, %s %D-point \n",
		                   meshsize, meshsize, size, "Star", 5); CHKERRQ(ierr);

	        ierr = PetscPrintf(PETSC_COMM_WORLD,
			                   " Explicit :  dT = %G (fixed), final time = %G, MaxSteps = %D\n",
			                   tsctx.tsize, tsctx.tfinal, tsctx.tsmax ); CHKERRQ(ierr);
//		}
		ierr = PetscPrintf(PETSC_COMM_WORLD, "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"); CHKERRQ(ierr);
	}

	CHKMEMQ;

	ierr = FormInitialValue(user); CHKERRQ(ierr);
	
	CHKMEMQ;

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   Solve the nonlinear system
	   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	PetscPreLoadStage("Solve");

	//ierr = Update(user); CHKERRQ(ierr);
	//CHKMEMQ;

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   Free work space.  All PETSc objects should be destroyed when they
	   are no longer needed.
	   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	
	ierr = DMRestoreGlobalVector(user->swe->da, &user->sol); CHKERRQ(ierr);

	ierr = SWEFinalize(user); CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
	if (PetscPreLoading) {
		ierr = PetscPrintf(PETSC_COMM_WORLD," PreLoading over!\n"); CHKERRQ(ierr);
	}
	PetscPreLoadEnd();

	ierr = PetscFinalize(); CHKERRQ(ierr);
	return 0;
}

/*******************************************************************************************/

#undef __FUNCT__
#define __FUNCT__ "Update"
PetscErrorCode Update(UserCtx *user)
{
	PetscErrorCode   ierr;
	ParamCtx         *param = user->param;
	TstepCtx        *tsctx = user->tsctx;
	PetscInt         max_steps, size;
	PassiveScalar    error[12];
	PetscLogDouble   time0,time1,totaltime=0.0;
	FILE            *fp;
	char             history[36], filename[256];
	MPI_Comm         comm=PETSC_COMM_WORLD;

	PetscFunctionBegin;

	if (param->PreLoading) {
		max_steps = 1;
	} else {
		max_steps = tsctx->tsmax;
	}

	PetscFunctionReturn(0);
}

/*******************************************************************************************/

/*******************************************************************************************/

#undef __FUNCT__
#define __FUNCT__ "FormInitialValue"
PetscErrorCode FormInitialValue(UserCtx *user)
{
	PetscErrorCode  ierr;
	ParamCtx        *param = user->param;
	TstepCtx        *tsctx = user->tsctx;
	MetricCtx       *metric = user->metric;
	DM              da = user->swe->da;
	PassiveScalar   dx = user->dx, dy = user->dy;
	ActiveField     **x=0;
	MetricField     **tensor;
	PetscInt        i, j, k, xl, yl, nxl, nyl, xg, yg, nxg, nyg;
	char            filename[PETSC_MAX_PATH_LEN-1];
	PetscBool       flg;

	PetscFunctionBegin;

	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n Form Initiate  \n"); CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners( da, &xg, &yg, 0, &nxg, &nyg, 0 ); CHKERRQ(ierr);
	ierr = DMDAGetCorners(      da, &xl, &yl, 0, &nxl, &nyl, 0 ); CHKERRQ(ierr);
	k = user->swe->myid;
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n Get Corners  \n"); CHKERRQ(ierr);

	ierr = DMDAVecGetArray(da,  user->sol, &x); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n Get Corners  \n"); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, metric->tensor, &tensor); CHKERRQ(ierr);

	PassiveScalar grav = param->gravity;
	PassiveScalar xi, yj, h0;
	h0 = 10.;
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n Initail values  \n"); CHKERRQ(ierr);
	for ( j=yl; j < yl+nyl; j++) {
		for ( i=xl; i < xl+nxl; i++) {
			xi = XLOW + ((PassiveScalar)i + 0.5)*dx;
			yj = YLOW + ((PassiveScalar)j + 0.5)*dy;
					
			x[j][i].eta = exp(-abs(xi) -abs(yj));
			x[j][i].U = 0.;
			x[j][i].V = 0.;
			tensor[j][i].H = h0;
		}
	}
	
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n End Initiate  \n"); CHKERRQ(ierr);
	

	ierr = DMDAVecRestoreArray(da, user->sol, &x); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(da, metric->tensor, &tensor); CHKERRQ(ierr);
	PetscFunctionReturn(0);
}


