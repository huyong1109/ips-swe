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
	PetscInt       meshsize;
	ParamCtx       param;
	TstepCtx       tsctx;
	MetricCtx      metric[1];
	BdyCtx         bdy[1];
	EventCtx       event[1];
	PetscLogStage  stage;
	PetscViewer    viewer;

	PetscInitialize(&argc,&argv,(char *)0,help);

	PetscPreLoadBegin(PETSC_TRUE,"SetUp");

	param.PreLoading = PetscPreLoading;

	comm = PETSC_COMM_WORLD;
	ierr = MPI_Comm_rank( comm, &rank ); CHKERRQ(ierr);
	ierr = MPI_Comm_size( comm, &size ); CHKERRQ(ierr);

	ierr = PetscViewerCreate(comm, &viewer); CHKERRQ(ierr);
	ierr = PetscViewerSetType(viewer, PETSCVIEWERASCII); CHKERRQ(ierr);
	ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB); CHKERRQ(ierr);

	/* Register a user-defined event for profiling (error checking). */
	ierr = PetscLogStageRegister("Initiating", &stage);
    	ierr = PetscLogStagePush(stage);

	/* Problem parameters */
	param.gravity  = GRAVITY;
	param.coriolis = CORIOLIS;

	/* Time-stepping context */
	tsctx.torder    = TORDER;
	tsctx.tstart    = 0.0;
	tsctx.tfinal    = 1.0;
	tsctx.tsize     = 0.05;
	tsctx.tsmax     = 5;

	ierr = PetscOptionsGetScalar(PETSC_NULL,"-tstart",&tsctx.tstart,PETSC_NULL); CHKERRQ(ierr);
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
	user->comm     = comm;
	user->rank     = rank;
	user->size     = size;
	user->tsctx    = &tsctx;
	user->event    = event;
	user->param    = &param;
	user->meshsize = meshsize;
	user->dx       = (XHI - XLOW) / (PassiveScalar)(user->meshsize);
	user->dy       = user->dx;
	user->metric   = metric;
	user->bdy      = bdy;
	
	ierr = SWEInitialize(user,0); CHKERRQ(ierr);
	
	CHKMEMQ;
	
	ierr = DMGetGlobalVector(user->swe->da, &(user->sol)); CHKERRQ(ierr);
	ierr = DMGetGlobalVector(user->swe->da, &(user->Q0)); CHKERRQ(ierr);
	ierr = DMGetGlobalVector(user->swe->da, &(user->Q1)); CHKERRQ(ierr);
	ierr = DMGetGlobalVector(user->swe->da, &(user->tmp)); CHKERRQ(ierr);
	ierr = DMGetLocalVector(user->swe->da, &(user->loc)); CHKERRQ(ierr);

	if (!param.PreLoading) {
		ierr = PetscPrintf(PETSC_COMM_WORLD, "\n+++++++++++++++++++++++ Problem parameters +++++++++++++++++++++\n"); CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD," Example: Isolated Mountain.\n"); CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,
		                   " Geometric parameters: Radius = %G, Gravity = %G, Coriolis = %G.\n",
		                   param.radius, param.gravity, param.coriolis); CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD," Problem size %D X %D (six patches), Ncpu = %D, %s %D-point \n",
		                   meshsize, meshsize, size, "Star", 5); CHKERRQ(ierr);

//		if (tsctx.cfl>0.0) {
//			ierr = PetscPrintf(PETSC_COMM_WORLD,
//			                   " Explicit : Torder = %D, CFL = %G, final time = %G, MaxSteps = %D\n",
//			                   tsctx.torder, tsctx.cfl, tsctx.tfinal, tsctx.tsmax ); CHKERRQ(ierr);
//		} else {
			ierr = PetscPrintf(PETSC_COMM_WORLD,
			                   " Explicit : Torder = %D, dT = %G (fixed), final time = %G, MaxSteps = %D\n",
			                   tsctx.torder, tsctx.tsize, tsctx.tfinal, tsctx.tsmax ); CHKERRQ(ierr);
//		}
		ierr = PetscPrintf(PETSC_COMM_WORLD, "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"); CHKERRQ(ierr);
	}

	CHKMEMQ;

	ierr = FormInitialValue(user); CHKERRQ(ierr);
	ierr = VecCopy(user->sol,user->Q1); CHKERRQ(ierr);
	
	CHKMEMQ;
	ierr = PetscLogStagePop();

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   Solve the nonlinear system
	   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	PetscPreLoadStage("Solve");

	ierr = Update(user); CHKERRQ(ierr);
	//CHKMEMQ;

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   Free work space.  All PETSc objects should be destroyed when they
	   are no longer needed.
	   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	
	ierr = DMRestoreGlobalVector(user->cubed->da, &user->sol); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(user->cubed->da, &user->Q0); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(user->cubed->da, &user->Q1); CHKERRQ(ierr);
	ierr = DMRestoreGlobalVector(user->cubed->da, &user->tmp); CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(user->cubed->da, &(user->loc)); CHKERRQ(ierr);

	ierr = SWEdFinalize(user); CHKERRQ(ierr);
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
	ParaCtx         *param = user->param;
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

	if (!param->PreLoading ) {
		ierr = MPI_Comm_size( user->comm, &size ); CHKERRQ(ierr);
		sprintf( history, "history_c%d_n%d.data", user->size, user->meshsize );

		ierr = PetscFOpen(user->comm, history, "a", &fp); CHKERRQ(ierr);
		ierr = PetscFPrintf(comm,fp, "%% step, tstep, tcurr, time, error_1, error_2, error_inf, mass, energy, enstrophy\n" ); CHKERRQ(ierr);
	}
#if use_misc
	if (!param->PreLoading && tsctx->tsback > 0) {
		sprintf( filename, "output.%d", user->tsctx->tscurr );
		ierr = DataSave(user->sol, filename); CHKERRQ(ierr);
//		ierr = Save_Mountain(user); CHKERRQ(ierr);
	}
#endif

	for (tsctx->tscurr = tsctx->tsstart + 1; (tsctx->tscurr <= tsctx->tsstart + max_steps); tsctx->tscurr++) {

		ierr = PetscGetTime(&time0); CHKERRQ(ierr);
		ierr = ExplicitSolve(user); CHKERRQ(ierr);
		ierr = PetscGetTime(&time1); CHKERRQ(ierr);
		tsctx->tcurr += tsctx->tsize;

		if (tsctx->tscurr < tsctx->tsstart + max_steps) {
			ierr = VecCopy(user->sol,user->Q1); CHKERRQ(ierr);
		}
		if (!param->PreLoading ) {
			ierr = PetscPrintf(comm, "\n====================== Step: %D, time: %G ====================\n",
			                   tsctx->tscurr, tsctx->tcurr ); CHKERRQ(ierr);

			ierr = PetscFPrintf(comm,fp, "%D, %G, %G, %G, ",
			                    tsctx->tscurr, tsctx->tsize, tsctx->tcurr, time1-time0 ); CHKERRQ(ierr);
#if use_misc
			if ( tsctx->tscomp > 0 ) {
				if ( (tsctx->tscurr)%(tsctx->tscomp)==0 ) {
					ierr = EstimateConserve(user,error); CHKERRQ(ierr);
					ierr = PetscPrintf(comm,
					                   " Conservation: mass = %G, energy = %G, enstrophy = %G (%20.15e/%20.15e), vorticity = %G/%G, divergence = %G/%G\n",
					                   (error[0]-error[1])/error[1], (error[2]-error[3])/error[3],
					                   (error[4]-error[5])/error[5], error[4], error[5], error[6], error[7], error[8], error[9]); CHKERRQ(ierr);
					ierr = PetscFPrintf(comm,fp, "%G, %G, %G",
					                    (error[0]-error[1])/error[1], (error[2]-error[3])/error[3], (error[4]-error[5])/error[5] ); CHKERRQ(ierr);
				}
			}
#endif
			ierr = PetscFPrintf(comm,fp, "\n");
			totaltime += time1-time0;
			ierr = PetscPrintf(comm, " Time cost: %G, %G\n", time1-time0, totaltime ); CHKERRQ(ierr);
#if use_misc
			if ( tsctx->tsback > 0 ) {
				if ( (tsctx->tscurr-tsctx->tsstart)%(tsctx->tsback)==0 ) {
					sprintf( filename, "output.%d", user->tsctx->tscurr );
					ierr = DataSave(user->sol,filename); CHKERRQ(ierr);
				}
			}
#endif
		}
		ierr = MPI_Barrier(comm); CHKERRQ(ierr);
		if ( tsctx->tcurr > tsctx->tfinal - EPS ) { tsctx->tscurr++; break; }
	}
	tsctx->tscurr--;
#if use_misc
	if ( tsctx->tscomp > 0 ) {
		if ( (tsctx->tscurr)%(tsctx->tscomp)!=0 ) {
			ierr = PetscPrintf(comm, "\n+++++++++++++++++++++++ Summary +++++++++++++++++++++++++\n" ); CHKERRQ(ierr);
			ierr = PetscPrintf(comm, " Final time = %G, Cost time = %G\n", tsctx->tcurr, totaltime ); CHKERRQ(ierr);
			ierr = EstimateConserve(user,error); CHKERRQ(ierr);
			ierr = PetscPrintf(comm,
			                   " Final conservation: mass = %G, energy = %G, enstrophy = %G, vorticity = %G/%G, divergence = %G/%G\n",
			                   (error[0]-error[1])/error[1], (error[2]-error[3])/error[3],
			                   (error[4]-error[5])/error[5], error[6], error[7], error[8], error[9]); CHKERRQ(ierr);
			ierr = PetscPrintf(comm,
			                   " Conservation: mass = %G/%G, energy = %G/%G, enstrophy = %G/%G\n",
			                   error[0],error[1], error[2],error[3], error[4],error[5]); CHKERRQ(ierr);
		}
	}
	if (!param->PreLoading && tsctx->tsback < 0 ) {
		if (tsctx->tscurr > tsctx->tsstart + max_steps)
			tsctx->tscurr--;
		sprintf( filename, "output.%d", user->tsctx->tscurr );
		ierr = DataSave(user->sol,filename); CHKERRQ(ierr);
	}
#endif
	if (!param->PreLoading) {
		ierr = PetscFClose(user->comm, fp); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}

/*******************************************************************************************/

#undef __FUNCT__
#define __FUNCT__ "FormInitialValue"
PetscErrorCode FormInitialValue(UserCtx *user)
{
	PetscErrorCode  ierr;
	ParaCtx         *param = user->param;
	TstepCtx        *tsctx = user->tsctx;
	MetricCtx       *metric = user->metric;
	DM              da = user->cubed->da, da0 = metric->da;
	PassiveScalar   R = param->radius;
	PassiveScalar   dx = user->dx, dy = user->dy;
	ActiveField     **x=0;
	MetricField     **tensor;
	PetscInt        i, j, k, xl, yl, nxl, nyl, xg, yg, nxg, nyg;
	char            filename[PETSC_MAX_PATH_LEN-1];
	PetscBool       flg;

	PetscFunctionBegin;

	ierr = DMDAGetGhostCorners( da, &xg, &yg, 0, &nxg, &nyg, 0 ); CHKERRQ(ierr);
	ierr = DMDAGetCorners(      da, &xl, &yl, 0, &nxl, &nyl, 0 ); CHKERRQ(ierr);
	k = user->cubed->myid;

	ierr = DMDAVecGetArray(da,  user->Q0, &x); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da0, metric->tensor, &tensor); CHKERRQ(ierr);

	PassiveScalar grav = param->gravity, omega = param->coriolis;
	PassiveScalar xi, yj, thetai, lambdai, h0, u0, alpha0=0.0;
	PassiveScalar u_lambda, u_theta, slambda, clambda, stheta, ctheta, uu, vv, tmp, p2c[3][3];
	h0 = 5960.0/Hunit;
	u0 = 20.0/Uunit;
	for ( j=yl; j < yl+nyl; j++) {
		for ( i=xl; i < xl+nxl; i++) {
			xi = XLOW + ((PassiveScalar)i + 0.5)*dx;
			yj = YLOW + ((PassiveScalar)j + 0.5)*dy;
			ierr = Coor_Cube2Polar( k, xi, yj, &thetai, &lambdai ); CHKERRQ(ierr);
			slambda = sin(lambdai);
			clambda = cos(lambdai);
			stheta  = sin(thetai);
			ctheta  = cos(thetai);
			tmp = -sin(alpha0)*clambda*ctheta + cos(alpha0)*stheta;
			x[j][i].H = h0 - ( ( R*omega*u0 + 0.5*u0*u0 )*tmp*tmp )/grav;
			u_lambda  = u0*( cos(alpha0)*ctheta + sin(alpha0)*clambda*stheta )/(R*ctheta); //TODO: theta=+-pi/2?
			u_theta   = -u0*sin(alpha0)*slambda/R;
			ierr = Vel_Polar2Cube_Coe( k, p2c, thetai, lambdai ); CHKERRQ(ierr);
			uu = u_lambda*p2c[1][1] + u_theta*p2c[1][2];
			vv = u_lambda*p2c[2][1] + u_theta*p2c[2][2];
			x[j][i].H -= tensor[j][i].Hs; 
			x[j][i].HU = x[j][i].H * uu;
			x[j][i].HV = x[j][i].H * vv;
		}
	}
	PetscLogFlops((40+5*FLOPS_SIN+4*FLOPS_COS)*nyl*nxl);

	ierr = DMDAVecRestoreArray(da, user->Q0, &x); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(da0, metric->tensor, &tensor); CHKERRQ(ierr);
#if use_misc
	ierr = PetscOptionsGetString(PETSC_NULL,"-f",filename,PETSC_MAX_PATH_LEN-1,&flg); CHKERRQ(ierr);
	if ( flg || tsctx->tsstart>0 ) {
		ierr = PetscPrintf(user->comm, " going to restore file %s\n", filename ); CHKERRQ(ierr);
		ierr = DataLoad(user->Q0,filename); CHKERRQ(ierr);	
	}
#endif
	ierr = VecCopy(user->Q0,user->sol); CHKERRQ(ierr);
	PetscFunctionReturn(0);
}


