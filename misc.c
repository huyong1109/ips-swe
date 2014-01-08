#include <petscsnes.h>
#include <petscdmda.h>
#include "def_swe.h"


#undef __FUNCT__
#define __FUNCT__ "DataSave"
PetscErrorCode DataSave(Vec x, char *filename)
{
	PetscErrorCode ierr;
	PetscViewer    dataviewer;

	PetscFunctionBegin;

	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n Save %s", filename); CHKERRQ(ierr);
	ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&dataviewer); CHKERRQ(ierr);
	ierr = PetscViewerSetFormat(dataviewer,PETSC_VIEWER_ASCII_SYMMODU); CHKERRQ(ierr);
	ierr = VecView(x,dataviewer); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n Done with %s", filename); CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&dataviewer); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DataShow"
PetscErrorCode DataShow(Vec x, char *vecname)
{
	PetscErrorCode ierr;
	PetscViewer    showviewer;

	PetscFunctionBegin;

	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n Show %s", vecname); CHKERRQ(ierr);
	ierr = PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL,PETSC_NULL,0, 0 , 1024, 980, &showviewer); CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject)showviewer, vecname); CHKERRQ(ierr);
	ierr = PetscViewerPushFormat(showviewer,PETSC_VIEWER_DRAW_CONTOUR); CHKERRQ(ierr);
	ierr = VecView(x,showviewer); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n Done with %s", vecname); CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&showviewer); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

/*******************************************************************************************/

#undef __FUNCT__
#define __FUNCT__ "DataLoad"
PetscErrorCode DataLoad(Vec x, char *filename)
{
	PetscErrorCode ierr;
	PetscViewer    dataviewer;

	PetscFunctionBegin;

	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&dataviewer); CHKERRQ(ierr);
	ierr = VecLoad(x,dataviewer); CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&dataviewer); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}



