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
	ierr = DMDACreate2d(user->comm,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_STENCIL_STAR,
		    user->meshsize,user->meshsize,user->px,user->py,
		    DOF,WIDTH,0,0,&da); CHKERRQ(ierr);
	swe->da = da;

	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n END Create DMDA 2D\n"); CHKERRQ(ierr);
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
