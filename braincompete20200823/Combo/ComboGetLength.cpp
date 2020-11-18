/*=================================================================
 * ComboGetLength.cpp 
 * ��ѯ�������Ѳɼ���������
 * �������(Matlab��)��ComboOpen���ص�һ���豸���
 * �������(Matlab��)��һ������Ѿ��ɼ��������ݳ���
 * 
 * ���뷽����
 *		��Matlab�����д�����ִ�� mex -g ComboGetLength.cpp

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
	// ��ȡ���ھ��
	//
	PHANDLE pDevice = (PHANDLE) mxGetData(prhs[0]);
	HANDLE hDevice = *pDevice;
	DWORD dwRet;
	COMMAND_MESSAGE Cmd = {0};
	CString str;

	//
	// ��ȡ�豸�������ڵ�������
	// 
	Cmd.SyncMark = 0xAA5500FF;
	Cmd.nCommand = ADC24_LOCK_BUFFER;
	if (!WriteFile(hDevice, &Cmd, sizeof(COMMAND_MESSAGE), &dwRet, NULL))
	{
		str.Format("ComboGetLength���������豸д������ʱ��������Handle = %d\n", hDevice);
		mexWarnMsgTxt(str);
	}
	Sleep(1);
	if (!ReadFile(hDevice, &Cmd, sizeof(COMMAND_MESSAGE), &dwRet, NULL))
	{
		str.Format("ComboGetLength��������ȡ�豸����ʱ��������Handle = %d\n", hDevice);
		mexWarnMsgTxt(str);
	}
	if (Cmd.SyncMark != 0xAA5500FF)
	{
		str.Format("ComboGetLength�������豸����ͨѶʱ�������䣬��������豸���������ݡ�Handle = %d\n", hDevice);
		mexWarnMsgTxt(str);
		INT32 Err[1000];
		for (int m=0; m<10; m++)
			ReadFile(hDevice, &Err, sizeof(INT32)*1000, &dwRet, NULL);
	}

	plhs[0] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
	INT32 *pHandle = (INT32 *)mxGetData(plhs[0]);
	*pHandle = (INT32)Cmd.Param0;
}
