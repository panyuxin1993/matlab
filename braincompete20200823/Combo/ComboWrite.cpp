/*=================================================================
 * ComboWrite.cpp
 * ��ָ���ĵ�����һ����豸д�����ݡ�
 * �﷨��ʽ
 *		count = ComboWrite(hDevice, sizeA, precision);
 * ���������	hDevice	- �豸���
 *				sizeA	- д������ݸ���
 *				precision - д�����ݵ����ݾ���
 * ���������	count	- ʵ��д������ݸ���
 
 * This is a MEX-file for MATLAB.
 * Copyright 1949-2015 by Bluesky.
 * All rights reserved.
 *=================================================================*/
#define _AFXDLL
#include "mex.h"
#include <string.h>
#include "stdarg.h"
#include "tchar.h"
#include "afx.h"
#include "afxstr.h"
#include "mmsystem.h"
#include "common.h"

void mexFunction(
	int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
    int        i;
       
    /* Examine input (right-hand-side) arguments. */
    mexPrintf("\nThere are %d right-hand-side argument(s).", nrhs);
    for (i=0; i<nrhs; i++)  {
        mexPrintf("\n\tInput Arg %i is of type:\t%s ",i,mxGetClassName(prhs[i]));
    }
    
    /* Examine output (left-hand-side) arguments. */
    mexPrintf("\n\nThere are %d left-hand-side argument(s).\n", nlhs);
    if (nlhs > nrhs)
      mexErrMsgIdAndTxt( "MATLAB:mexfunction:inputOutputMismatch",
              "Cannot specify more outputs than inputs.\n");
    
    for (i=0; i<nlhs; i++)  {
        plhs[i]=mxCreateDoubleMatrix(1,1,mxREAL);
        *mxGetPr(plhs[i])=(double)mxGetNumberOfElements(prhs[i]);
    }
}

