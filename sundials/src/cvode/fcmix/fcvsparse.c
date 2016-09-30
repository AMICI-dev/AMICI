/*
 * -----------------------------------------------------------------
 * $Revision: 4815 $
 * $Date: 2016-07-20 16:51:55 -0700 (Wed, 20 Jul 2016) $
 * -----------------------------------------------------------------
 * Programmer(s): Carol Woodward @ LLNL
 *                Ting Yan and Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2015, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include "fcvode.h"
#include "cvode_impl.h"
#include <cvode/cvode_sparse.h>

/* Prototype of the Fortran routine */
 
#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
 
extern void FCV_SPJAC(realtype *T, realtype *Y, 
		       realtype *FY, int *N, int *NNZ, 
		       realtype *JDATA, int *JRVALS, 
		       int *JCPTRS, realtype *H, 
		       long int *IPAR, realtype *RPAR, 
		       realtype *V1, realtype *V2, 
		       realtype *V3, int *ier);
 
#ifdef __cplusplus
}
#endif
 
/*=============================================================*/

/* Fortran interface to C routine CVSlsSetSparseJacFn; see
   fcvode.h for further information */
void FCV_SPARSESETJAC(int *ier)
{
  *ier = CVSlsSetSparseJacFn(CV_cvodemem, FCVSparseJac);
}

/*=============================================================*/
 
/* C interface to user-supplied Fortran routine FCVSPJAC; see 
   fcvode.h for additional information  */
int FCVSparseJac(realtype t, N_Vector y, N_Vector fy, 
		 SlsMat J, void *user_data, N_Vector vtemp1, 
		 N_Vector vtemp2, N_Vector vtemp3)
{
  int ier;
  realtype *ydata, *fydata, *v1data, *v2data, *v3data;
  realtype h;
  FCVUserData CV_userdata;

  CVodeGetLastStep(CV_cvodemem, &h);
  ydata   = N_VGetArrayPointer(y);
  fydata  = N_VGetArrayPointer(fy);
  v1data  = N_VGetArrayPointer(vtemp1);
  v2data  = N_VGetArrayPointer(vtemp2);
  v3data  = N_VGetArrayPointer(vtemp3);
  CV_userdata = (FCVUserData) user_data;

  FCV_SPJAC(&t, ydata, fydata, &(J->NP), &(J->NNZ),
	    J->data, J->indexvals, J->indexptrs, &h, 
	    CV_userdata->ipar, CV_userdata->rpar, v1data, 
	    v2data, v3data, &ier); 
  return(ier);
}

