/*
 * -----------------------------------------------------------------
 * $Revision: 4525 $
 * $Date: 2015-08-12 16:28:00 -0700 (Wed, 12 Aug 2015) $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * Common header file for the direct linear solvers in IDA.
 * -----------------------------------------------------------------
 */

#ifndef _IDADLS_H
#define _IDADLS_H

#include <sundials/sundials_direct.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * =================================================================
 *              I D A D I R E C T     C O N S T A N T S
 * =================================================================
 */

/* 
 * -----------------------------------------------------------------
 * IDADLS return values 
 * -----------------------------------------------------------------
 */

#define IDADLS_SUCCESS           0
#define IDADLS_MEM_NULL         -1
#define IDADLS_LMEM_NULL        -2
#define IDADLS_ILL_INPUT        -3
#define IDADLS_MEM_FAIL         -4

/* Additional last_flag values */

#define IDADLS_JACFUNC_UNRECVR  -5
#define IDADLS_JACFUNC_RECVR    -6

/*
 * =================================================================
 *              F U N C T I O N   T Y P E S
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Types : IDADlsDenseJacFn
 * -----------------------------------------------------------------
 *
 * A dense Jacobian approximation function djac must be of type 
 * IDADlsDenseJacFn.
 * Its parameters are:                     
 *                                                                
 * N   is the problem size, and length of all vector arguments.   
 *                                                                
 * t   is the current value of the independent variable t.        
 *                                                                
 * y   is the current value of the dependent variable vector,     
 *     namely the predicted value of y(t).                     
 *                                                                
 * yp  is the current value of the derivative vector y',          
 *     namely the predicted value of y'(t).                    
 *                                                                
 * r   is the residual vector F(tt,yy,yp).                     
 *                                                                
 * c_j is the scalar in the system Jacobian, proportional to 
 *     the inverse of the step size h.
 *                                                                
 * user_data is a pointer to user Jacobian data - the same as the    
 *     user_data parameter passed to IDASetUserData.                     
 *                                                                
 * Jac is the dense matrix (of type DlsMat) to be loaded by  
 *     an IDADlsDenseJacFn routine with an approximation to the   
 *     system Jacobian matrix                                  
 *            J = dF/dy + c_j *dF/dy'                            
 *     at the given point (t,y,y'), where the ODE system is    
 *     given by F(t,y,y') = 0.
 *     Note that Jac is NOT preset to zero!
 *                                                                
 * tmp1, tmp2, tmp3 are pointers to memory allocated for          
 *     N_Vectors which can be used by an IDADlsDenseJacFn routine 
 *     as temporary storage or work space.                     
 *                                                                
 * A IDADlsDenseJacFn should return                                
 *     0 if successful,                                           
 *     a positive int if a recoverable error occurred, or         
 *     a negative int if a nonrecoverable error occurred.         
 * In the case of a recoverable error return, the integrator will 
 * attempt to recover by reducing the stepsize (which changes c_j).
 *
 * -----------------------------------------------------------------
 *
 * NOTE: The following are two efficient ways to load a dense Jac:         
 * (1) (with macros - no explicit data structure references)      
 *     for (j=0; j < Neq; j++) {                                  
 *       col_j = LAPACK_DENSE_COL(Jac,j);                                 
 *       for (i=0; i < Neq; i++) {                                
 *         generate J_ij = the (i,j)th Jacobian element           
 *         col_j[i] = J_ij;                                       
 *       }                                                        
 *     }                                                          
 * (2) (without macros - explicit data structure references)      
 *     for (j=0; j < Neq; j++) {                                  
 *       col_j = (Jac->data)[j];                                   
 *       for (i=0; i < Neq; i++) {                                
 *         generate J_ij = the (i,j)th Jacobian element           
 *         col_j[i] = J_ij;                                       
 *       }                                                        
 *     }                                                          
 * A third way, using the LAPACK_DENSE_ELEM(A,i,j) macro, is much less   
 * efficient in general.  It is only appropriate for use in small 
 * problems in which efficiency of access is NOT a major concern. 
 *                                                                
 * NOTE: If the user's Jacobian routine needs other quantities,   
 *     they are accessible as follows: hcur (the current stepsize)
 *     and ewt (the error weight vector) are accessible through   
 *     IDAGetCurrentStep and IDAGetErrWeights, respectively, but this
 *     requires including in user_data a pointer to the solver memory.
 *     The unit roundoff is available as UNIT_ROUNDOFF defined in
 *     sundials_types.h.
 *
 * -----------------------------------------------------------------
 */
  
  
typedef int (*IDADlsDenseJacFn)(long int N, realtype t, realtype c_j,
				N_Vector y, N_Vector yp, N_Vector r, 
				DlsMat Jac, void *user_data,
				N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/*
 * -----------------------------------------------------------------
 * Types : IDADlsBandJacFn
 * -----------------------------------------------------------------
 * A banded Jacobian approximation function bjac must have the    
 * prototype given below. Its parameters are:                     
 *                                                                
 * Neq is the problem size, and length of all vector arguments.   
 *                                                                
 * mupper is the upper bandwidth of the banded Jacobian matrix.   
 *                                                                
 * mlower is the lower bandwidth of the banded Jacobian matrix.   
 *                                                                
 * tt is the current value of the independent variable t.        
 *                                                                
 * yy is the current value of the dependent variable vector,     
 *    namely the predicted value of y(t).                     
 *                                                                
 * yp is the current value of the derivative vector y',          
 *    namely the predicted value of y'(t).                    
 *                                                                
 * rr is the residual vector F(tt,yy,yp).                     
 *                                                                
 * c_j is the scalar in the system Jacobian, proportional to 1/hh.
 *                                                                
 * user_data  is a pointer to user Jacobian data - the same as the    
 *    user_data parameter passed to IDASetUserData.                      
 *                                                                
 * Jac is the band matrix (of type BandMat) to be loaded by    
 *     an IDADlsBandJacFn routine with an approximation to the    
 *     system Jacobian matrix                                  
 *            J = dF/dy + c_j *dF/dy'                             
 *     at the given point (t,y,y'), where the DAE system is    
 *     given by F(t,y,y') = 0.  Jac is preset to zero, so only 
 *     the nonzero elements need to be loaded.  See note below.
 *                                                                
 * tmp1, tmp2, tmp3 are pointers to memory allocated for          
 *     N_Vectors which can be used by an IDADlsBandJacFn routine  
 *     as temporary storage or work space.                     
 *                                                                
 * An IDADlsBandJacFn function should return                                 
 *     0 if successful,                                           
 *     a positive int if a recoverable error occurred, or         
 *     a negative int if a nonrecoverable error occurred.         
 * In the case of a recoverable error return, the integrator will 
 * attempt to recover by reducing the stepsize (which changes c_j).
 *
 * -----------------------------------------------------------------
 *
 * NOTE: The following are two efficient ways to load Jac:
 *                                                                
 * (1) (with macros - no explicit data structure references)      
 *    for (j=0; j < Neq; j++) {                                   
 *       col_j = BAND_COL(Jac,j);                                  
 *       for (i=j-mupper; i <= j+mlower; i++) {                   
 *         generate J_ij = the (i,j)th Jacobian element           
 *         BAND_COL_ELEM(col_j,i,j) = J_ij;                       
 *       }                                                        
 *     }                                                          
 *                                                                
 * (2) (with BAND_COL macro, but without BAND_COL_ELEM macro)     
 *    for (j=0; j < Neq; j++) {                                   
 *       col_j = BAND_COL(Jac,j);                                  
 *       for (k=-mupper; k <= mlower; k++) {                      
 *         generate J_ij = the (i,j)th Jacobian element, i=j+k    
 *         col_j[k] = J_ij;                                       
 *       }                                                        
 *     }                                                          
 *                                                                
 * A third way, using the BAND_ELEM(A,i,j) macro, is much less    
 * efficient in general.  It is only appropriate for use in small 
 * problems in which efficiency of access is NOT a major concern. 
 *                                                                
 * NOTE: If the user's Jacobian routine needs other quantities,   
 *       they are accessible as follows: hcur (the current stepsize)
 *       and ewt (the error weight vector) are accessible through   
 *       IDAGetCurrentStep and IDAGetErrWeights, respectively, but this
 *       requires including in user_data a pointer to the solver memory.
 *       The unit roundoff is available as UNIT_ROUNDOFF defined in
 *       sundials_types.h                   
 *                                                                
 * -----------------------------------------------------------------
 */

typedef int (*IDADlsBandJacFn)(long int N, long int mupper, long int mlower,
			       realtype t, realtype c_j, 
			       N_Vector y, N_Vector yp, N_Vector r,
			       DlsMat Jac, void *user_data,
			       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  
/*
 * =================================================================
 *            E X P O R T E D    F U N C T I O N S 
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Optional inputs to the IDADLS linear solver
 * -----------------------------------------------------------------
 * IDADlsSetDenseJacFn specifies the dense Jacobian approximation
 * routine to be used for a direct dense linear solver.
 *
 * IDADlsSetBandJacFn specifies the band Jacobian approximation
 * routine to be used for a direct band linear solver.
 *
 * By default, a difference quotient approximation, supplied with
 * the solver is used.
 *
 * The return value is one of:
 *    IDADLS_SUCCESS   if successful
 *    IDADLS_MEM_NULL  if the IDA memory was NULL
 *    IDADLS_LMEM_NULL if the linear solver memory was NULL
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDADlsSetDenseJacFn(void *ida_mem, IDADlsDenseJacFn jac);
SUNDIALS_EXPORT int IDADlsSetBandJacFn(void *ida_mem, IDADlsBandJacFn jac);

/*
 * -----------------------------------------------------------------
 * Optional outputs from the IDADLS linear solver
 * -----------------------------------------------------------------
 *
 * IDADlsGetWorkSpace   returns the real and integer workspace used
 *                      by the direct linear solver.
 * IDADlsGetNumJacEvals returns the number of calls made to the
 *                      Jacobian evaluation routine jac.
 * IDADlsGetNumResEvals returns the number of calls to the user
 *                      F routine due to finite difference Jacobian
 *                      evaluation.
 * IDADlsGetLastFlag    returns the last error flag set by any of
 *                      the IDADLS interface functions.
 *
 * The return value of IDADlsGet* is one of:
 *    IDADLS_SUCCESS   if successful
 *    IDADLS_MEM_NULL  if the IDA memory was NULL
 *    IDADLS_LMEM_NULL if the linear solver memory was NULL
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDADlsGetWorkSpace(void *ida_mem, long int *lenrwLS, long int *leniwLS);
SUNDIALS_EXPORT int IDADlsGetNumJacEvals(void *ida_mem, long int *njevals);
SUNDIALS_EXPORT int IDADlsGetNumResEvals(void *ida_mem, long int *nfevalsLS);
SUNDIALS_EXPORT int IDADlsGetLastFlag(void *ida_mem, long int *flag);

/*
 * -----------------------------------------------------------------
 * The following function returns the name of the constant 
 * associated with a IDADLS return flag
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT char *IDADlsGetReturnFlagName(long int flag);


#ifdef __cplusplus
}
#endif

#endif
