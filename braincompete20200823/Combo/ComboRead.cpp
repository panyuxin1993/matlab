/*=================================================================
 * ComboRead.cpp 
 * ��ָ���ĵ�����һ����豸��ȡ���ݡ�
 * �﷨��ʽ
 *		[A, count] = ComboRead(hDevice, sizeA, precision);
 * ���������	hDevice	- �豸���
 *				sizeA	- ��ȡ�����ݸ���
 *				precision - ÿ�����ݵ�Ԫռ�ݵ��ֽڿռ䣬��INT32ռ��4�ֽڣ��ò�������4��
 * ���������	A		- ��ȡ�����ݾ��󣬰��������ŷš�
 *				count	- ʵ�ʶ�ȡ�����ݸ���
 
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

#define StrLen 20

void mexFunction(
	int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	UINT8 precision = 1;	// һ������ͨѶ��Ԫռ�ݵ��ֽڿռ䡣Ĭ��Ϊuint8�ͣ�ÿ�����ݵ�ռ��1���ֽڡ�
	UINT32 count = 0;		
	char* pPrecisionStr = (char*) mxMalloc(StrLen);
	//
	// chech the input argument 
	//
	if (nlhs < 2)
		mexErrMsgTxt("�����������\n���÷�����\t[A, count] = ComboRead(hDevice, sizeA, precision);\n\tprecisionĬ��Ϊuint8����");
	if (nrhs == 3)
	{
		mxGetString(prhs[2], pPrecisionStr, StrLen);
		// int nLen = WideCharToMultiByte(CP_ACP, 0, (LPCWSTR)pPrecisionWStr, -1, NULL, 0, NULL, NULL);	// ��һ�ε�����ȡ����Ҫת���ռ�Ĵ�С
		// WideCharToMultiByte(CP_ACP, 0, (LPCWSTR)pPrecisionWStr, -1, pPrecisionStr, nLen, NULL, NULL);	// �ڶ��ε������ת������
		CString strPrecision(pPrecisionStr);
		if (strPrecision.Compare("uchar") == 0)			precision = sizeof(UCHAR);
		else if (strPrecision.Compare("schar") == 0)	precision = sizeof(CHAR);
		else if (strPrecision.Compare("char") == 0)		precision = sizeof(CHAR);
		else if (strPrecision.Compare("int8") == 0)		precision = sizeof(INT8);
		else if (strPrecision.Compare("int16") == 0)	precision = sizeof(INT16);
		else if (strPrecision.Compare("int32") == 0)	precision = sizeof(INT32);
		else if (strPrecision.Compare("uint8") == 0)	precision = sizeof(UINT8);
		else if (strPrecision.Compare("uint16") == 0)	precision = sizeof(UINT16);
		else if (strPrecision.Compare("uint32") == 0)	precision = sizeof(UINT32);
		else if (strPrecision.Compare("short") == 0)	precision = sizeof(SHORT);
		else if (strPrecision.Compare("int") == 0)		precision = sizeof(INT);
		else if (strPrecision.Compare("long") == 0)		precision = sizeof(LONG);
		else if (strPrecision.Compare("ushort") == 0)	precision = sizeof(USHORT);
		else if (strPrecision.Compare("uint") == 0)		precision = sizeof(UINT);
		else if (strPrecision.Compare("ulong") == 0)	precision = sizeof(ULONG);
		else if (strPrecision.Compare("single") == 0)	precision = sizeof(FLOAT);
		else if (strPrecision.Compare("float32") == 0)	precision = sizeof(FLOAT);
		else if (strPrecision.Compare("float") == 0)	precision = sizeof(FLOAT);
		else if (strPrecision.Compare("double") == 0)	precision = sizeof(DOUBLE);
		else if (strPrecision.Compare("float64") == 0)	precision = sizeof(DOUBLE);
		else	mexErrMsgTxt("�����������precision�������󣬲μ�fread�е�precision��Ч������\n");
	}


	//
	// ��ȡ�������
	//
	PHANDLE pDevice = (PHANDLE) mxGetData(prhs[0]);
	HANDLE hDevice = *pDevice;
	PINT32 pSizeA = (PINT32) mxGetData(prhs[1]);
	INT32 size = *pSizeA;


/*
	//
	// ��ȡ����ӿڣ����������������
	//
	plhs[0] = mxCreateNumericMatrix(1, nDataCount, mxINT32_CLASS, mxREAL);	// ���ɻ���������ָ�ɸ���һ���������
	INT32 *pData = (INT32 *)mxGetData(plhs[0]);								// ȡ�û�����ʵ�ʴ�����ݵĵ�ַָ��

	
	int        i;
       
    / * Examine input (right-hand-side) arguments. * /
    mexPrintf("\nThere are %d right-hand-side argument(s).", nrhs);
    for (i=0; i<nrhs; i++)  {
        mexPrintf("\n\tInput Arg %i is of type:\t%s ",i,mxGetClassName(prhs[i]));
    }
    
    / * Examine output (left-hand-side) arguments. * /
    mexPrintf("\n\nThere are %d left-hand-side argument(s).\n", nlhs);
    if (nlhs > nrhs)
      mexErrMsgIdAndTxt( "MATLAB:mexfunction:inputOutputMismatch",
              "Cannot specify more outputs than inputs.\n");
    
    for (i=0; i<nlhs; i++)  {
        plhs[i]=mxCreateDoubleMatrix(1,1,mxREAL);
        *mxGetPr(plhs[i])=(double)mxGetNumberOfElements(prhs[i]);
    }
*/




	if (nlhs == 2)
	{
		plhs[1] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);	// ���ɻ���������ָ�ɸ���һ���������
		PINT32 pDataCount = (PINT32) mxGetData(plhs[1]);				// ȡ�û�����ʵ�ʴ�����ݵĵ�ַָ�룬�ں�ʵ�ʶ�ȡ��������
		*pDataCount = count;
	}

	mxFree(pPrecisionStr);
}

