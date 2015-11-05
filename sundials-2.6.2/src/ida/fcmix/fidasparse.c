/*
 * -----------------------------------------------------------------
 * $Revision: 4402 $
 * $Date: 2015-02-28 19:35:39 -0800 (Sat, 28 Feb 2015) $
 * -----------------------------------------------------------------
 * Programmer(s): Carol Woodward @ LLNL
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
#include "fida.h"
#include "ida_impl.h"
#include <sundials/sundials_sparse.h>

/* Prototype of the Fortran routine */
 
#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
 
extern void FIDA_SPJAC(realtype *T, realtype *CJ, realtype *Y, 
		       realtype *YP, realtype *R, int *N, int *NNZ, 
		       realtype *JDATA, int *JRVALS, 
		       int *JCPTRS, realtype *H, 
		       long int *IPAR, realtype *RPAR, 
		       realtype *V1, realtype *V2, 
		       realtype *V3, int *ier);
 
#ifdef __cplusplus
}
#endif
 
/*=============================================================*/
 
/* C interface to user-supplied Fortran routine FIDASPJAC; see 
   fida.h for additional information  */
int FIDASparseJac(realtype t, realtype cj, N_Vector y, N_Vector yp,
		  N_Vector fval, SlsMat J, void *user_data, 
		  N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  int ier;
  realtype *ydata, *ypdata, *rdata, *v1data, *v2data, *v3data;
  realtype h;
  FIDAUserData IDA_userdata;

  IDAGetLastStep(IDA_idamem, &h);
  ydata   = N_VGetArrayPointer(y);
  ypdata  = N_VGetArrayPointer(yp);
  rdata  = N_VGetArrayPointer(fval);
  v1data  = N_VGetArrayPointer(vtemp1);
  v2data  = N_VGetArrayPointer(vtemp2);
  v3data  = N_VGetArrayPointer(vtemp3);
  IDA_userdata = (FIDAUserData) user_data;

  FIDA_SPJAC(&t, &cj, ydata, ypdata, rdata, &(J->N), &(J->NNZ),
	    J->data, J->rowvals, J->colptrs, &h, 
	    IDA_userdata->ipar, IDA_userdata->rpar, v1data, 
	    v2data, v3data, &ier); 
  return(ier);
}

