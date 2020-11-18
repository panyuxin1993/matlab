/*=================================================================
 * ComboOpen.cpp 
 *	Version 1.10
 * ��ָ���ĵ�����һ����豸
 * �������(Matlab��)��ComboQuery���ص�һ���豸���ں�
 * �������(Matlab��)���豸���
 * 
 * ���뷽����
 *		��Matlab�����д�����ִ�� mex -g ComboOpen.cpp

 * This is a MEX-file for MATLAB.
 * Copyright 1949-2017 by Bluesky.
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

HANDLE OpenComm(CString strCom);

void mexFunction(
	int nlhs, mxArray *plhs[],			// ��ֵ�����������
	int nrhs, const mxArray *prhs[])	// ��ֵ�����������
{
	//
	// chech the input argument
	//
	if (nrhs < 1)
		mexErrMsgTxt("Has not input. Input argument is the name of ComboAmplifier device.");
	if (!mxIsChar(prhs[0]))
		mexErrMsgTxt("Input argument is not type of CHAR.");

	//
	// open the device
	//
	size_t nBufLen = mxGetM(prhs[0]) * mxGetN(prhs[0]) * sizeof(mxChar) + 1;
	char *strCom = (char *)calloc(nBufLen, sizeof(char));
	if (strCom==NULL)	mexErrMsgTxt("Not enough memory in heap.");
	mxGetString(prhs[0], strCom, (int)nBufLen);
	HANDLE hDevice = OpenComm(strCom);
	free(strCom);

	plhs[0] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
	INT32 *pHandle = (INT32 *)mxGetData(plhs[0]);
	*pHandle = (INT32)hDevice;
}

// ��ָ���Ĵ���
HANDLE OpenComm(CString strCom)
{
	HANDLE hCom = NULL;
	strCom = "\\\\.\\" + strCom;		//// ΢��Ԥ����ı�׼�豸��ֻ���С�COM1��-��COM9����COM10�����ϴ��ж˿�Ӧ��ʹ��"\\.\COM10"�����Ĵ���ר������
	//// 1. �򿪴���
	hCom = CreateFile(strCom,		// COM������
		GENERIC_READ|GENERIC_WRITE, // �������д
		0,							// ��ռ��ʽ
		NULL,						// ���ð�ȫ�����Խṹ��ȱʡֵΪNULL��
		OPEN_EXISTING,				// �򿪶����Ǵ���
		0,							// ͬ����ʽ
		NULL);
	if (hCom == (HANDLE)-1)
	{
		// ErrorMessage(strCom+_T(": �򿪴���ʧ��!"));
		return hCom = NULL;
	}
	//// 2. ���ô���
	//	�������ʱ���ã����ٶ�ȡ���⴮������ʱ�����
	//	���������ʱ�����迹�����豸����Ĵ��ڶ�ȡ�����ַ�ʱ������ͬ���ȴ����޷����ء�
	//	����ڲ����迹�����豸ʱ��������ʱ�ķ�ʽ�򿪴��ڣ�
	//	������ʹ���迹�����豸ʱ���Խ�ֹ��ʱ�ķ�ʽ�򿪴��ڡ�
	//	�˴�Ϊ����ʹ���豸����˽�ֹ��ʱ���á�
	SetupComm(hCom, 102400, 524288);		//���(PC-->USB)��������С100KB������(USB-->PC)�������Ĵ�С��512KB
	// COMMTIMEOUTS TimeOuts;					//�趨����ʱ����λ��Ϊms��
	// TimeOuts.ReadIntervalTimeout = 0;			// ���������ֽ�֮��ĳ�ʱʱ�䡣=0����ʾ��ʹ�øó�ʱ����
	// TimeOuts.ReadTotalTimeoutMultiplier = 0;	// ���������ܳ�ʱʱ�䡣ʵ��ִ�еĳ�ʱʱ��Ҫ���Զ�ȡ���ֽ�����
	// TimeOuts.ReadTotalTimeoutConstant = 50;		// �������ĳ�ʱ�������ܳ�ʱ = ReadTotalTimeoutConstant + ReadTotalTimeoutMultiplier * �ֽ��� = 50 ms��
	// TimeOuts.WriteTotalTimeoutMultiplier = 10;
	// TimeOuts.WriteTotalTimeoutConstant = 10;
	// SetCommTimeouts(hCom, &TimeOuts);			// ���ó�ʱ
	DCB dcb;
	GetCommState(hCom, &dcb);
	dcb.BaudRate = 115200;						//������Ϊ115200
	dcb.ByteSize = 8;							//ÿ���ֽ���8λ
	dcb.Parity = NOPARITY;						//����żУ��λ
	dcb.StopBits = ONESTOPBIT;					//1��ֹͣλ
	SetCommState(hCom, &dcb);
	PurgeComm( hCom, PURGE_TXABORT | PURGE_RXABORT | PURGE_TXCLEAR | PURGE_RXCLEAR ); //��ɾ����롢���������
	return hCom;
}
