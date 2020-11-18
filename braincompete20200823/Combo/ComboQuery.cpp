/*=================================================================
* ComboQuery.cpp 
*	Version 1.10
* ��ѯ������һ����豸
* �������(Matlab��)����
* �������(Matlab��)��Device���豸���ں��б���ά�ַ������飬ÿһ����һ�������豸���ơ�
* 
* ���뷽����
*		��Matlab�����д�����ִ�� mex -g ComboQuery.cpp

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

#define StrLen 20					// �豸���Ƶĳ��ȣ����20���ַ�

HANDLE OpenComm(CString strCom);
int SearchDevice(void);
BOOL IdentifyComboAmp(CString tchValue);

TCHAR (*ptcDevice)[StrLen];			// ��ά����ָ�룬�������в��ҵ����豸���ơ�ÿһ�д��һ���豸��������C���԰ڷŷ�ʽ��


void mexFunction(
	int nlhs, mxArray *plhs[],			// ��ֵ�����������
	int nrhs, const mxArray *prhs[])	// ��ֵ�����������
{
	int nDeviceNum = 0;
	ptcDevice = NULL;
	mxArray *pDeviceName;

	nDeviceNum = SearchDevice();
	// mexPrintf("\nFound %d device(s).\n", nDeviceNum);

	if (nDeviceNum > 0)
	{
		// ���ɷ���ֵmxArray����
		mwSize ndim = 2;	// ָ�����������ά������ά
		mwSize dims[2];	// ���� dims[0] * dims[1] ����
		dims[0] = nDeviceNum; dims[1] = StrLen;
		pDeviceName = mxCreateCharArray (ndim, dims);
		mxChar *charData = (mxChar *)mxGetData(pDeviceName);

		// Matlab�ڲ����ַ��������ǰ����ŷŵģ�����Ҫת��ptcDevice�����ݵ�����˳��
		TCHAR* pMatlabDeviceName = new TCHAR[StrLen*nDeviceNum];
		TCHAR* pTemp = pMatlabDeviceName;
		for (int i=0; i<StrLen; i++)
			for (int k=0; k<nDeviceNum; k++)
				*pTemp++ = ptcDevice[k][i];
		// Matlab���ַ���Ϊwchar_t���ַ��ͣ���Ҫ��char_t��ת��Ϊwchar_t��
		int wLen = MultiByteToWideChar(CP_ACP, 0, (LPCSTR)pMatlabDeviceName, nDeviceNum*StrLen, NULL, 0);	// ��һ�ε�����ȡ����Ҫת���ռ�Ĵ�С
		MultiByteToWideChar(CP_ACP, 0, (LPCSTR)pMatlabDeviceName, nDeviceNum*StrLen, charData, wLen);		// �ڶ��ε������ת������


		delete [] ptcDevice;
		delete [] pMatlabDeviceName;
	}
	else
		pDeviceName = mxCreateCharArray (0, 0);

	// ��Ĭ�Ϸ���ֵans�����豸����
	plhs[0] = pDeviceName;
}


// ��������Ӳ���豸������ֵ���豸����
int SearchDevice(void)
{
	delete [] ptcDevice;
	ptcDevice = NULL;
	int nDeviceNum = 0;
	HKEY hkey;
	LONG lRes = RegOpenKeyEx(HKEY_LOCAL_MACHINE, _T( "HARDWARE\\DEVICEMAP\\SERIALCOMM"), NULL, KEY_QUERY_VALUE | KEY_ENUMERATE_SUB_KEYS | KEY_READ, &hkey);
	if (lRes == ERROR_SUCCESS)
	{ 
		TCHAR tchKey[MAX_PATH];
		TCHAR tchValue[StrLen];
		DWORD dwIndex = 0;
		DWORD dwType = REG_SZ;
		while(lRes == ERROR_SUCCESS)
		{
			DWORD dwCount = MAX_PATH;
			DWORD dwVCount = StrLen;
			for(int k=0; k<StrLen; k++)	tchValue[k] = 0;				// ���ַ�������0
			lRes = RegEnumValue(hkey, dwIndex++, tchKey, &dwCount, NULL, &dwType, (LPBYTE)tchValue, &dwVCount);
			if(lRes == ERROR_SUCCESS)
			{
				if((dwVCount > 0) && (dwCount > 0))
				{
					if (IdentifyComboAmp(tchValue))
					{
						/*	// ��������ҵ����豸
						mexPrintf("Device %d: ", nDeviceNum);
						for (int k=0; k<StrLen; k++)
							mexPrintf("%c(%d) ", tchValue[k], tchValue[k]);
						mexPrintf("\n");
						*/

						nDeviceNum++;
						// m_DeviceComBox.AddString(tchValue);
						TCHAR (*DeviceName)[StrLen] = new TCHAR[nDeviceNum][StrLen];
						if (DeviceName == NULL)	mexErrMsgTxt("Not enough memory in heap.\n");
						else
						{
							for (int i=0; i<nDeviceNum-1; i++)
							{
								for(int k=0; k<StrLen; k++)
									DeviceName[i][k] = ptcDevice[i][k];	// ��ԭ�е��豸������������
							}
							for(int k=0; k<StrLen; k++)
								DeviceName[nDeviceNum-1][k] = tchValue[k];			// �����»�ȡ���豸������������
							delete [] ptcDevice;
							ptcDevice = DeviceName;
						}
					}
				}
			}
		}
	}
	RegCloseKey(hkey);

	return nDeviceNum;
}


// �����Ƿ���һ�廯�Ŵ����豸
BOOL IdentifyComboAmp(CString tchValue)
{
	HANDLE hCom = OpenComm(tchValue);	// ���ھ��
	DWORD dwRet;
	COMMAND_MESSAGE Cmd = {0};

	Cmd.SyncMark = 0xAA5500FF;
	Cmd.nCommand = CMD_GET_DEVICE;		// ʶ���豸
	BOOL bRet = WriteFile(hCom, &Cmd, sizeof(COMMAND_MESSAGE), &dwRet, NULL);		// if (!bRet)	ErrorMessage("WriteFile Error");
	Sleep(20);
	bRet = ReadFile(hCom, &Cmd, sizeof(COMMAND_MESSAGE), &dwRet, NULL);				// if (!bRet)	ErrorMessage("ReadFile Error");
	if (strncmp((char*)&Cmd.Param0, "ComoboAmp", 9) == 0)							//// һ�廯�Ŵ�����ش�ComoboAmp��
		bRet = TRUE;
	else
	{
		//// �رմ���
		CloseHandle(hCom);
		return FALSE;
	}

	Cmd.SyncMark = 0xAA5500FF;
	Cmd.nCommand = CMD_GET_VENDOR;		// ʶ���豸�����
	bRet = WriteFile(hCom, &Cmd, sizeof(COMMAND_MESSAGE), &dwRet, NULL);		// if (!bRet)	ErrorMessage("WriteFile Error");
	Sleep(2);
	bRet = ReadFile(hCom, &Cmd, sizeof(COMMAND_MESSAGE), &dwRet, NULL);				// if (!bRet)	ErrorMessage("ReadFile Error");
	if (strncmp((char*)&Cmd.Param0, "Bluesky", 9) == 0)							//// һ�廯�Ŵ�����ش�Bluesky��
		bRet = TRUE;
	else
		bRet = FALSE;

	//// �رմ���
	CloseHandle(hCom);
	return bRet;
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
	//	�ڴ˴�Ϊ�����豸���������ʱ���á�
	SetupComm(hCom, 102400, 524288);		//���(PC-->USB)��������С100KB������(USB-->PC)�������Ĵ�С��512KB
	COMMTIMEOUTS TimeOuts;					//�趨����ʱ����λ��Ϊms��
	TimeOuts.ReadIntervalTimeout = 0;			// ���������ֽ�֮��ĳ�ʱʱ�䡣=0����ʾ��ʹ�øó�ʱ����
	TimeOuts.ReadTotalTimeoutMultiplier = 0;	// ���������ܳ�ʱʱ�䡣ʵ��ִ�еĳ�ʱʱ��Ҫ���Զ�ȡ���ֽ�����
	TimeOuts.ReadTotalTimeoutConstant = 50;		// �������ĳ�ʱ�������ܳ�ʱ = ReadTotalTimeoutConstant + ReadTotalTimeoutMultiplier * �ֽ��� = 50 ms��
	TimeOuts.WriteTotalTimeoutMultiplier = 10;
	TimeOuts.WriteTotalTimeoutConstant = 10;
	SetCommTimeouts(hCom, &TimeOuts);			// ���ó�ʱ
	DCB dcb; GetCommState(hCom, &dcb);
	dcb.BaudRate = 115200;						//������Ϊ115200
	dcb.ByteSize = 8;							//ÿ���ֽ���8λ
	dcb.Parity = NOPARITY;						//����żУ��λ
	dcb.StopBits = ONESTOPBIT;					//1��ֹͣλ
	SetCommState(hCom, &dcb);
	PurgeComm( hCom, PURGE_TXABORT | PURGE_RXABORT | PURGE_TXCLEAR | PURGE_RXCLEAR ); //��ɾ����롢���������
	return hCom;
}
