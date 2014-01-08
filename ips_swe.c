#include <petscdmda.h>
#include <private/daimpl.h>
#include "def_swe.h"
/*******************************************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SWEInitialize"
PetscErrorCode SWEInitialize( UserCtx* user)
{
	SWE             *swe = user->swe;
	PetscErrorCode  ierr;
	DM		da;

	PetscFunctionBegin;

	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n Create DMDA 2D [x,y] = [%d, %d]\n", user->px, user->py); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n Create DMDA 2D [lx,ly] = [%d, %d]\n", user->lx, user->ly); CHKERRQ(ierr);
	ierr = DMDACreate2d(user->comm,DMDA_BOUNDARY_GHOSTED,DMDA_BOUNDARY_GHOSTED,DMDA_STENCIL_STAR,
		    user->meshsize,user->meshsize,user->px,user->py,
		    DOF,WIDTH,0,0,&da); CHKERRQ(ierr);
	swe->da = da;

	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n END Create DMDA 2D\n"); CHKERRQ(ierr);
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "ExplicitStep"
PetscErrorCode ExplicitStep(UserCtx *user)
{
	PetscErrorCode  ierr;
	ParamCtx        *param = user->param;
	TstepCtx        *tsctx = user->tsctx;
	MetricCtx       *metric = user->metric;
	DM              da = user->swe->da;

	PassiveScalar   dx = user->dx, dy = user->dy;
	PassiveScalar   dt = tsctx->tsize;
	ActiveField     **x=0;
	MetricField     **tensor;
	Vec		loc_X = PETSC_NULL;
	PetscInt        i, j, k, xl, yl, nxl, nyl, xg, yg, nxg, nyg;
	FILE		*fp;
	char            filename[PETSC_MAX_PATH_LEN-1];
	PetscBool       flg;

	PetscFunctionBegin;
    
    
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n Start step %d, times = %G\n", tsctx->tscurr, tsctx->tcurr); CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners( da, &xg, &yg, 0, &nxg, &nyg, 0 ); CHKERRQ(ierr);
	ierr = DMDAGetCorners(      da, &xl, &yl, 0, &nxl, &nyl, 0 ); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n [x, y] = ghost[%d, %d,%d, %d], coner[%d, %d,%d, %d] \n", xg, yg,nxg,nyg,xl, yl,nxl,nyl  ); CHKERRQ(ierr);
	k = user->swe->myid;

	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n create local \n"); CHKERRQ(ierr);
	ierr = DMCreateLocalVector( da, &loc_X  ); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n created locx"); CHKERRQ(ierr);
	ierr = DMGetLocalVector( da, &loc_X  ); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n Get locx"); CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(da, user->sol, INSERT_VALUES, loc_X); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n Get locx"); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(  da, user->sol, INSERT_VALUES, loc_X); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n Get x"); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, loc_X, &x); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n Get tensor"); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, metric->tensor, &tensor); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n Geta  tensor"); CHKERRQ(ierr);

	PassiveScalar grav = param->gravity;
	PetscScalar du, dv;
	
	sprintf( filename, "initinstep%d.dat",user->tsctx->tscurr );
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n Save loc"); CHKERRQ(ierr);
	ierr = DataSave(user->sol, filename); CHKERRQ(ierr);
	
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n Update U, V  \n"); CHKERRQ(ierr);
	/* update u,v  regardless of global boundary */
	for ( j=yl+1; j < yl+nyl-1; j++) {
		for ( i=xl+1; i < xl+nxl-1; i++) {
			//ierr = PetscPrintf(PETSC_COMM_WORLD, "\n id = %d, [x, y] = [%d, %d] \n", k ,i , j ); CHKERRQ(ierr);
			du = -dt*grav*(x[j][i+1].eta-x[j][i].eta)/dx;
			dv = -dt*grav*(x[j+1][i].eta-x[j][i].eta)/dy;
			x[j][i].U += du;
			x[j][i].V += dv;
		}
	}
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n Update Boundary  \n"); CHKERRQ(ierr);
	/* boundary */
	if (xl == user->meshsize - 1 ){
		for ( j=yl; j < yl+nyl; j++) {
			x[j][i].U = 0.;
		}
	}
	if (yl == user->meshsize - 1){
		for ( i=xl; i < xl+nxl; i++) {
			x[j][i].V = 0.;
		}
	}
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n Update eta \n"); CHKERRQ(ierr);
	/* update eta   regardless of global boundary */
	for ( j=yl; j < yl+nyl; j++) {
		for ( i=xl; i < xl+nxl; i++) {
			du  = ((x[j][i+1].H+x[j][i].H)*x[j][i].U-(x[j][i].H+x[j][i-1].H)*x[j][i-1].U)/dx;
			dv  = ((x[j+1][i].H+x[j][i].H)*x[j][i].V-(x[j][i].H+x[j-1][i].H)*x[j-1][i].V)/dy;
			x[j][i].eta -= dt*(du+dv);
		}
	}

	
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n End Step  \n"); CHKERRQ(ierr);
	

	ierr = DMLocalToGlobalBegin(da, loc_X, INSERT_VALUES, user->sol); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n loc to global"); CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(  da,loc_X , INSERT_VALUES, user->sol); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(da, loc_X, &x); CHKERRQ(ierr);
	sprintf( filename, "step%d.dat", user->tsctx->tscurr );
	ierr = DataSave(user->sol, filename); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(da, metric->tensor, &tensor); CHKERRQ(ierr);
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "SWEFinalize"
PetscErrorCode SWEFinalize( UserCtx* user )
{
	PetscErrorCode  ierr;

	PetscFunctionBegin;

	if (user->lx) {ierr = PetscFree(user->lx); CHKERRQ(ierr);}
	if (user->ly) {ierr = PetscFree(user->ly); CHKERRQ(ierr);}

	ierr = DMDestroy(&user->swe->da); CHKERRQ(ierr);
	PetscFunctionReturn(0);
}
