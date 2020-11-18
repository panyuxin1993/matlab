#pragma once

typedef enum __MESSAGE_TYPE__		// USB������Ϣ����ö��
{
	CMD_GET_DEVICE = 0,				// PC-->MCU��ѯ���豸���ơ�����Ӧ��ComoboAmp��
	CMD_GET_VENDOR,					// PC-->MCU��ѯ���豸���ʦ������Ӧ��Bluesky��
	ADC24_LOCK_BUFFER,				// PC-->MCU��
	ADC24_READ_BUFFER,				// PC-->MCU��
	ADC16_CONFIG,					// PC-->MCU��
	ADC16_LOCK_BUFFER,				// PC-->MCU��
	ADC16_READ_BUFFER,				// PC-->MCU��
	LED_CMD,						// PC-->MCU������LEDָ��
	TYPE_STATE,						// MCU-->PC��״̬��Ϣ�������ִ�����Ϣ
} USB_Message_Data_Type;

typedef struct						// ָ����Ϣ��
{
	UINT32 SyncMark;				// ��ϢУ���������Ϊ0xAA5500FF
	UINT32 nCommand;				// ������
	UINT32 Param0;					// ��Ϣ����0
	UINT32 Param1;					// ��Ϣ����1
	UINT32 Param2;					// ��Ϣ����2
} COMMAND_MESSAGE, *PCOMMAND_MESSAGE;
