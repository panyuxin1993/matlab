/*=================================================================
 * ComboClose.cpp
 * �ر�ָ���ĵ�����һ����豸��
 * �������(Matlab��)��ComboOpen���ص�һ���豸���
 * �������(Matlab��)����
 * 
 * ���뷽����
 *		��Matlab�����д�����ִ�� mex -g ComboClose.cpp

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
	int nlhs, mxArray *plhs[],			// ��ֵ�����������
	int nrhs, const mxArray *prhs[])	// ��ֵ�����������
{
	//
	// chech the input argument 
	//
	if (nrhs < 1)
		mexErrMsgTxt("Has not input. Input argument is the handle of ComboAmplifier device.");
	if (!mxIsInt32(prhs[0]))
		mexErrMsgTxt("Input argument is not type of HANDLE.");

	//
	// �رմ���
	//
	PHANDLE pDevice = (PHANDLE) mxGetData(prhs[0]);
	CloseHandle(*pDevice);
}
