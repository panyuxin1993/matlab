/*=================================================================
 * ComboGetData.cpp
 * ��ȡ�����������ݡ����ñ�����֮ǰ�������ȵ���һ��ComboGetLength��������ȡ�豸�����������ݳ��ȡ�
 * �������(Matlab��)��	prhs[0]��ComboOpen���ص�һ���豸�����
 *						prhs[1]����Ҫ��ȡ�����ݳ��ȣ������ݵ�Ϊ��λ��
 * �������(Matlab��)��	plhs[0]��һ����豸���Ѿ�������ADC24���ݣ�
 *						plhs[1]��һ����豸���Ѿ�������Port���ݡ�
 * 
 * ���뷽����
 *		��Matlab�����д�����ִ�� mex -g ComboGetData.cpp -lWinmm.lib
 
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
	if (nrhs < 2)
		mexErrMsgTxt("�����������\n���÷�����\t[Adc24Data PortData] = ComboGetData(hDeviceHandle, nDataCount);");
	if (!mxIsInt32(prhs[0]))
		mexErrMsgTxt("Input argument is not type of HANDLE.");
	//
	// ��ȡ�������
	//
	PHANDLE pDevice = (PHANDLE) mxGetData(prhs[0]);
	HANDLE hDevice = *pDevice;
	PINT32 pDataCount = (PINT32) mxGetData(prhs[1]);
	INT32 nDataCount = *pDataCount;

	//
	// ��ȡ����ӿڣ����������������
	//
	plhs[0] = mxCreateNumericMatrix(1, nDataCount, mxINT32_CLASS, mxREAL);	// ����Adc24���ݻ���������ָ�ɸ���һ���������
	INT32 *pData = (INT32 *)mxGetData(plhs[0]);								// ȡ�û�����ʵ�ʴ�����ݵĵ�ַָ��
	plhs[1] = mxCreateNumericMatrix(2, nDataCount, mxINT16_CLASS, mxREAL);	// ����Adc1���ݻ���������ָ�ɸ��ڶ����������
	UINT8 *pAdc1Data = (UINT8 *)mxGetData(plhs[1]);							// ȡ�û�����ʵ�ʴ�����ݵĵ�ַָ��
	plhs[2] = mxCreateNumericMatrix(1, nDataCount, mxUINT8_CLASS, mxREAL);	// ����Port���ݻ���������ָ�ɸ��������������
	UINT8 *pPortData = (UINT8 *)mxGetData(plhs[2]);							// ȡ�û�����ʵ�ʴ�����ݵĵ�ַָ��


	DWORD dwRet;
	COMMAND_MESSAGE Cmd = {0};
	CString str;
	DWORD dwTransNum = 0;

	//
	// Ҫ���ȡ�豸�������ڵ�ADC��Port����
	// 
	Cmd.SyncMark = 0xAA5500FF;
	Cmd.nCommand = ADC24_READ_BUFFER;
	Cmd.Param0 = nDataCount;
	if (!WriteFile(hDevice, &Cmd, sizeof(COMMAND_MESSAGE), &dwRet, NULL))
	{
		str.Format("ComboGetData���������豸д������ʱ��������Handle = %d\n", hDevice);
		mexWarnMsgTxt(str);
	}
	Sleep(0);
	//
	// ��ȡPort�˿�����
	// 
	COMMAND_MESSAGE Vendor;
	Vendor.SyncMark = 0xAA5500FF;
	Vendor.nCommand = CMD_GET_VENDOR;		// ʶ���豸�����
	INT32 nNeedCount = nDataCount;
	PUINT8 pReadBuf = pPortData;
	INT32 nExtDataNum = 0;					// �������������
	do 
	{
		DWORD dwTime = timeGetTime();
		ReadFile(hDevice, pPortData, nNeedCount, &dwRet, NULL);
		if (nNeedCount != dwRet)
		{
			DWORD dwDelayTime = timeGetTime();
			pReadBuf += dwRet;
			nNeedCount -= dwRet;
			//// ����USBͨ������ <-- д��һ������󣬶����ö���USBlyzer�������ӵ������ݴ�MCU��������ڻ�������
			//// �ɴ�֤ʵ��IN���ݿ����������ǳ�����STM32�ϣ�STM32û�������ݳ�����
			//// ���������д��ʲô������USB���еģ�����Ҫ��STM32�����µ����ݣ�ԭ����USB�����������ݲ���һ�������
			WriteFile(hDevice, &Vendor, sizeof(COMMAND_MESSAGE), &dwRet, NULL);
			nExtDataNum++;
//			str.Format("STM32��USB IN���������һ�Σ���ʱ%d ms���Ѿ��ָ��ֳ�\n", dwDelayTime - dwTime);
			str.Format("\n", dwDelayTime - dwTime);
			mexWarnMsgTxt(str);
		}
		else
			break;
	} while (nNeedCount);
	while (nExtDataNum--)
		ReadFile(hDevice, &Vendor, sizeof(COMMAND_MESSAGE), &dwRet, NULL);	// ����ȫ���Ŀ�������Ϣ
	//
	// ��ȡADC1����
	// 
	nNeedCount = nDataCount * 2 * sizeof(UINT16);
	pReadBuf = (PUINT8) pAdc1Data;
	nExtDataNum = 0;					// �������������
	do 
	{
		DWORD dwTime = timeGetTime();
		ReadFile(hDevice, pAdc1Data, nNeedCount, &dwRet, NULL);
		if (nNeedCount != dwRet)
		{
			DWORD dwDelayTime = timeGetTime();
			pReadBuf += dwRet;
			nNeedCount -= dwRet;
			//// ����USBͨ������ <-- д��һ������󣬶����ö���USBlyzer�������ӵ������ݴ�MCU��������ڻ�������
			//// �ɴ�֤ʵ��IN���ݿ����������ǳ�����STM32�ϣ�STM32û�������ݳ�����
			//// ���������д��ʲô������USB���еģ�����Ҫ��STM32�����µ����ݣ�ԭ����USB�����������ݲ���һ�������
			WriteFile(hDevice, &Vendor, sizeof(COMMAND_MESSAGE), &dwRet, NULL);
			nExtDataNum++;
//			str.Format("STM32��USB IN���������һ�Σ���ʱ%d ms���Ѿ��ָ��ֳ�\n", dwDelayTime - dwTime);
			str.Format("\n", dwDelayTime - dwTime);
			mexWarnMsgTxt(str);
		}
		else
			break;
	} while (nNeedCount);
	while (nExtDataNum--)
		ReadFile(hDevice, &Vendor, sizeof(COMMAND_MESSAGE), &dwRet, NULL);	// ����ȫ���Ŀ�������Ϣ
	//
	// ��ȡAdc24����
	// 
	nNeedCount = nDataCount * sizeof(INT32);
	pReadBuf = (PUINT8) pData;
	nExtDataNum = 0;
	do 
	{
		DWORD dwTime = timeGetTime();
		ReadFile(hDevice, pReadBuf, nNeedCount, &dwRet, NULL);
		if (nNeedCount != dwRet)
		{
			DWORD dwDelayTime = timeGetTime();
			pReadBuf += dwRet;
			nNeedCount -= dwRet;
			//// ����USBͨ������ <-- д��һ������󣬶����ö���USBlyzer�������ӵ������ݴ�MCU��������ڻ�������
			//// �ɴ�֤ʵ��IN���ݿ����������ǳ�����STM32�ϣ�STM32û�������ݳ�����
			//// ���������д��ʲô������USB���еģ�����Ҫ��STM32�����µ����ݣ�ԭ����USB�����������ݲ���һ�������
			Vendor.SyncMark = 0xAA5500FF;
			Vendor.nCommand = CMD_GET_VENDOR;		// ʶ���豸�����
			WriteFile(hDevice, &Vendor, sizeof(COMMAND_MESSAGE), &dwRet, NULL);
			nExtDataNum++;
//			str.Format("STM32��USB IN���������һ�Σ���ʱ%d ms���Ѿ��ָ��ֳ�\n", dwDelayTime - dwTime);
			str.Format("\n", dwDelayTime - dwTime);
			mexWarnMsgTxt(str);
		}
		else
			break;
	} while (nNeedCount);
	while (nExtDataNum--)
		ReadFile(hDevice, &Vendor, sizeof(COMMAND_MESSAGE), &dwRet, NULL);	// ����ȫ���Ŀ�������Ϣ
}
