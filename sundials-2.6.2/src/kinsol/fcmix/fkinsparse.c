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
 * This is the implementation file for the Fortran interface to
 * the KINSuperLUMT solver. See fkinsol.h for usage.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include "fkinsol.h"
#include "kinsol_impl.h"
#include <sundials/sundials_sparse.h>

/* Prototype of the Fortran routine */
 
#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
 
extern void FKIN_SPJAC(realtype *Y, realtype *FY, int *N, int *NNZ, 
                       realtype *JDATA, int *JRVALS, int *JCPTRS, 
                       realtype *V1, realtype *V2, int *ier);
 
#ifdef __cplusplus
}
#endif
 
/*=============================================================*/
 
/* C interface to user-supplied Fortran routine FKINSPJAC; see 
   fkinsol.h for additional information  */
int FKINSparseJac(N_Vector y, N_Vector fy, 
                  SlsMat J, void *user_data, N_Vector vtemp1, 
                  N_Vector vtemp2)
{
  int ier;
  realtype *ydata, *fydata, *v1data, *v2data;
 
  ydata   = N_VGetArrayPointer(y);
  fydata  = N_VGetArrayPointer(fy);
  v1data  = N_VGetArrayPointer(vtemp1);
  v2data  = N_VGetArrayPointer(vtemp2);
 
  FKIN_SPJAC(ydata, fydata, &(J->N), &(J->NNZ),
             J->data, J->rowvals, J->colptrs,  
             v1data, v2data, &ier); 
  return(ier);
}

