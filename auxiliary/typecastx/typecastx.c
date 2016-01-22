/*************************************************************************************
 *
 * MATLAB (R) is a trademark of The Mathworks (R) Corporation
 *
 * Function:    typecast
 * Filename:    typecast.c
 * Programmer:  James Tursa
 * Version:     3.00
 * Date:        March 17, 2011
 * Copyright:   (c) 2009, 2011 by James Tursa, All Rights Reserved
 *
 %  This code uses the BSD License:
 %
 %  Redistribution and use in source and binary forms, with or without 
 %  modification, are permitted provided that the following conditions are 
 %  met:
 %
 %     * Redistributions of source code must retain the above copyright 
 %       notice, this list of conditions and the following disclaimer.
 %     * Redistributions in binary form must reproduce the above copyright 
 %       notice, this list of conditions and the following disclaimer in 
 %       the documentation and/or other materials provided with the distribution
 %      
 %  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 %  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 %  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
 %  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
 %  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 %  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
 %  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
 %  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
 %  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 %  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 %  POSSIBILITY OF SUCH DAMAGE.
 *
 * typecast is a mex function intended to mimic the MATLAB intrinsic typecast function
 * for those users with older versions of MATLAB that do not have this intrinsic. Users
 * of newer versions of MATLAB may be interested in this C-mex version to take advantage
 * of the several extensions offered. This C-mex version of typecast differs from the
 * intrinsic typecast in the following important aspects: 
 *
 *                      Intrinsic typecast    C-mex typecast
 *                      ------------------    --------------
 * Type of copy:        Deep Data Copy        Shared Data Copy
 * Allows complex:      No                    Yes
 * Allows logical:      No                    Yes (cannot convert from complex)
 * Allows char:         No                    Yes (cannot convert from complex)
 * Allows non-vector:   No                    Yes
 *
 * Since this C-mex typecast produces a shared data copy of the original, it
 * is more efficient than the MATLAB intrinsic typecast, which may be important
 * if you are working with large variables. For non-vector inputs, the first
 * non-singleton dimension must be compatible for the conversion.
 *
 * Building:
 *
 * typecast requires that a mex routine be built (one time only). This
 * process is typically self-building the first time you call the function
 * as long as you have the files typecast.m and typecast.c in the same
 * directory somewhere on the MATLAB path. If you need to manually build
 * the mex function, here are the commands:
 *
 * >> mex -setup
 *   (then follow instructions to select a C or C++ compiler of your choice)
 * >> mex typecast.c
 *
 * The usage is as follows (from The Mathworks website documentation):
 *
 * Syntax
 *
 * Y = typecast(X, type)
 *
 * Description
 *
 * Y = typecast(X, type) converts a value in X to the data type specified by type.
 * Input X must be a full numeric or char or logical variable. The type input is a string
 * set to one of the following: 'uint8', 'int8', 'uint16', 'int16', 'uint32', 'int32',
 * 'uint64', 'int64', 'single', 'double', 'char' or 'logical'. typecast is different
 * from the MATLAB cast function in that it does not alter the input data. typecast
 * always returns the same number of bytes in the output Y as were in the input X.
 * For example, casting the 16-bit integer 1000 to uint8 with typecast returns the full
 * 16 bits in two 8-bit segments (3 and 232) thus keeping its original value
 * (3*256 + 232 = 1000). The cast function, on the other hand, truncates the input value
 * to 255.
 *
 * The output of typecast can be formatted differently depending on what system you use it on.
 * Some computer systems store data starting with its most significant byte (an ordering
 * called big-endian), while others start with the least significant byte (called little-endian). 
 *
 * typecast issues an error if X contains fewer values than are needed to make an output value. 
 *
 */

#include "mex.h"
#include <string.h>

/* Needed for older versions of MATLAB that do not have the mwSize typedef */

#ifndef  MWSIZE_MAX
#define  mwIndex        int
#define  mwSignedIndex  int
#define  mwSize         int
#endif

/* Used for setting the varible type */

#define VariableType_Temporary  4

/* mxArray struct for hacking into the fields */

struct mxArray_Tag {
	char *name; /* Name of variable in workspace, NULL for R2009a and later */
	mxClassID ClassID; /*  0 = unknown
					       1 = cell
						   2 = struct
						   3 = logical
						   4 = char
						   5 = void
						   6 = double
						   7 = single
						   8 = int8
						   9 = uint8
					      10 = int16
						  11 = uint16
						  12 = int32
						  13 = uint32
						  14 = int64
						  15 = uint64
						  16 = function_handle
						  17 = opaque (User Defined Class indicator new classdef style)
						  18 = object (User Defined Class indicator old @directory style)
						  19 = index (deprecated)
						  20 = sparse (deprecated)
					   */
	int VariableType;  /*  0 = normal
					       1 = persistent
						   2 = global
						   3 = sub-element (field or cell)
						   4 = temporary
						   5 = (unknown)
						   6 = property of opaque class object
					   */
	mxArray *CrossLink;  /* Address of next shared-data variable in linked list */
	size_t ndim;
	unsigned int RefCount; /* Number of sub-elements identical to this one */
	unsigned int flags;  /*  bit  0 = is scalar double full
						     bit  2 = is empty double full
							 bit  4 = is temporary
							 bit  5 = is sparse
							 bit  8 = is constant from parsed m-file
							 bit  9 = is numeric
							 bits 24 - 31 = User Bits
						 */
	union {
		size_t M;        /* Row size for 2D matrices, or            */
		size_t *dims;    /* Pointer to dims array for nD > 2 arrays */
	} Mdims;
	size_t N;            /* Column size for 2D matrices */
	void *pr;            /* Pointer to real data (or pointer to cell or field elements */
	void *pi;            /* Pointer to imag data (or pointer to field information */
	union {
		mwIndex *ir;        /* Pointer to row values for sparse arrays, or           */
		mxClassID ClassID;  /* User Defined Class ID (opaque new style), or          */
		char *ClassName;    /* Pointer to User Defined Class Name (object old style) */
	} irClassNameID;
	union {
		mwIndex *jc;        /* Pointer to column values for sparse arrays, or */
		mxClassID ClassID;  /* User Defined Class ID (object old style)       */
	} jcClassID;
	size_t nzmax;           /* Number of elements allocated for sparse matrix */
	unsigned int reserved;  /* (unknown) */
};

/* Prototype */

mxArray *mxCreateSharedDataCopy(const mxArray *mx);  /* Undocumented function */

/* Gateway function */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	struct mxArray_Tag *mx;
    mwSize inbytes, outbytes, i, k, idim, ndim;
	mwSize *dims_old, *dims_new;
    char *outstring;
    mxClassID outclass;
	int out_numeric;

/* Check input arguments */
    
	if( nrhs > 2 ) {
        mexErrMsgTxt("Too many input arguments.");
	}
	if( nrhs < 2 ) {
        mexErrMsgTxt("Not enough input arguments.");
	}
	if( nlhs > 1 ) {
        mexErrMsgTxt("Too many output arguments.");
	}
	if( mxIsSparse(prhs[0]) || (!mxIsNumeric(prhs[0]) && !mxIsChar(prhs[0]) && !mxIsLogical(prhs[0])) ) {
        mexErrMsgTxt("The first input argument must be a full numeric value, or char, or logical.");
	}
	if( !mxIsChar(prhs[1]) ) {
        mexErrMsgTxt("The second input argument must be a character array.");
	}

/* Get input argument byte length */

   inbytes = mxGetElementSize(prhs[0]);
   
/* Check second input argument for desired output type */
   
   outstring = mxArrayToString(prhs[1]);

   out_numeric = 1;
   if(        strcmp(outstring,"int8") == 0 ) {
       outclass = mxINT8_CLASS;
       outbytes = 1;
   } else if( strcmp(outstring,"uint8") == 0 ) {
       outclass = mxUINT8_CLASS;
       outbytes = 1;
   } else if( strcmp(outstring,"int16") == 0 ) {
       outclass = mxINT16_CLASS;
       outbytes = 2;
   } else if( strcmp(outstring,"uint16") == 0 ) {
       outclass = mxUINT16_CLASS;
       outbytes = 2;
   } else if( strcmp(outstring,"int32") == 0 ) {
       outclass = mxINT32_CLASS;
       outbytes = 4;
   } else if( strcmp(outstring,"uint32") == 0 ) {
       outclass = mxUINT32_CLASS;
       outbytes = 4;
   } else if( strcmp(outstring,"int64") == 0 ) {
       outclass = mxINT64_CLASS;
       outbytes = 8;
   } else if( strcmp(outstring,"uint64") == 0 ) {
       outclass = mxUINT64_CLASS;
       outbytes = 8;
   } else if( strcmp(outstring,"double") == 0 ) {
       outclass = mxDOUBLE_CLASS;
       outbytes = 8;
   } else if( strcmp(outstring,"single") == 0 ) {
       outclass = mxSINGLE_CLASS;
       outbytes = 4;
   } else if( strcmp(outstring,"char") == 0 ) {
       outclass = mxCHAR_CLASS;
       outbytes = 2;
	   out_numeric = 0;
   } else if( strcmp(outstring,"logical") == 0 ) {
       outclass = mxLOGICAL_CLASS;
       outbytes = 1;
	   out_numeric = 0;
   } else {
       mxFree(outstring);
       mexErrMsgTxt("Unsupported class.\n");
   }
   mxFree(outstring);

/* Check for complex coversion to non-numeric */

   if( mxIsComplex(prhs[0]) && !out_numeric ) {
       mexErrMsgTxt("Cannot typecast a complex input to a non-numeric class.\n");
   }

/* Check for empty input. No data to share, so simply create a new empty variable */

    if( mxIsEmpty(prhs[0]) ) {
	    if( out_numeric ) {
		    plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]),
			                               mxGetDimensions(prhs[0]),
							               outclass, mxREAL);
	    } else if( outclass == mxCHAR_CLASS ) {
		    plhs[0] = mxCreateCharArray(mxGetNumberOfDimensions(prhs[0]),
		 	                            mxGetDimensions(prhs[0]));
	    } else {
		    plhs[0] = mxCreateLogicalArray(mxGetNumberOfDimensions(prhs[0]),
			                              mxGetDimensions(prhs[0]));
	    }
	    return;
    }

/* Check the old & new sizes for compatibility */
   
    ndim = mxGetNumberOfDimensions(prhs[0]);
    dims_old = mxGetDimensions(prhs[0]);
    for( i=0; i<ndim; i++ ) {
	    if( dims_old[i] != 1 || i == ndim-1 ) {
		    k = (dims_old[i] * inbytes) / outbytes;
			if( k * outbytes != dims_old[i] * inbytes ) {
                mexErrMsgTxt("Too few input values to make output type.\n");
			}
			idim = i;
		    break;
	    }
    }
	dims_new = mxMalloc(ndim * sizeof(*dims_new));
    for( i=0; i<ndim; i++ ) {
		dims_new[i] = dims_old[i];
    }
	dims_new[idim] = k;
   
/* Create the output array as a shared data copy, then manually set the class
 * and size parameters by accessing the structure fields directly. Note that
 * this is using undocumented MATLAB API functions and hacking of the
 * mxArray_tag structure, so this may not work in future versions of MATLAB.
 */
   
   plhs[0] = mxCreateSharedDataCopy(prhs[0]);
   mx = (struct mxArray_Tag *) plhs[0];
   mx->ClassID = outclass;
   mx->VariableType = VariableType_Temporary;
   mxSetDimensions(plhs[0],dims_new,ndim);
   mxFree(dims_new);

/* Also need to fix up the flags */

   if( outclass == mxDOUBLE_CLASS ) {
	   if( mxGetNumberOfElements(plhs[0]) == 1 ) {
		   mx->flags = 0x0211;  /* Set numeric, temporary, and full double scalar bits */
	   } else {
		   mx->flags = 0x0210;  /* Set numeric and temporary bits */
	   }
   } else if( out_numeric ) {
	   mx->flags = 0x0210;  /* Set numeric and temporary bits */
   } else {
	   mx->flags = 0x0010;  /* Set the temporary bit */
   }
}
