% ������һ����豸��ComboXxxͨѶ�������������Matlab���ݲɼ����̡�
% ComboXxxͨѶ������������5���������ɣ�
%	DeviceNames = ComboQuery;			��ѯ������PC���Ѿ����ӵ����е�����һ����豸��
%	handle = ComboOpen(DeviceName);		��ָ�����豸�������ظ��豸�Ĳ��������
%	len = ComboGetLength(handle);		��ѯָ���豸��ǰ�Ѿ��ɼ�������������
%	data = ComboGetData(handle, len);	�����豸�Ѿ��ɼ��������ݣ�len��������Ҫ��ComboGetLength()�������صĳ���һ�£�
%	ComboClose(handle);					�ر��Ѿ��򿪵��豸���ͷ��豸��
%
%
% ʾ������
%
% % 1. ���豸
% com = ComboQuery;						% comΪ�ַ������飬�������еĵ�����һ����豸����
% h = ComboOpen(com(1,:));				% �򿪵�һ���豸
%										% com(1,:)�ں���һ���豸��com(2,:)�ں��ڶ����豸�������豸������������
% % 2. �ɼ����ݡ���ο����ظ���ν��С�
% figure(9)
% data = [];
% for i=1:20
%	len = ComboGetLength(h);				% ȡ��һ������Ѿ��ɼ������ݳ���
%	[A E P] = ComboGetData(h, len);			% �����������ݳ��ȣ���������
%	data = [data, double(A).*9.9341e-09];	% ADC24���ݡ�1 LSB = ��2.5V/30/2^24
%	ExtAdc = [ExtAdc; E'];					% ��չADC����
%	port = [port, P];						% ���ֶ˿�����
%	data = [data, double(A).*9.9341e-09];	% 1 LSB = ��2.5V/30/2^24
%	figure(9), plot(data)					% ��ʾ��ЩADC����
%	figure(10), plot(bitget(port, 1));		% ��ʾPB9����
% end
% % 3. �ر��豸���˳�
% ComboClose(h)								% �ر��豸

help ComboHelp

% 1. ���豸
delete(instrfindall);						% �ȹر����д���
com = ComboQuery;							% comΪ�ַ������飬�������еĵ�����һ����豸����
h = ComboOpen(com(1,:));					% �򿪵�һ���豸

% 2. �ɼ����ݡ���ο����ظ���ν��С�
figure(9), clf
data = [];
port = [];
ExtAdc = [];
t1 = tic;
while 1
	len = ComboGetLength(h);				% ȡ��һ������Ѿ��ɼ������ݳ���
	[A, E, P] = ComboGetData(h, len);	% �����������ݳ��ȣ���������
	data = [data, double(A).*9.9341e-09];	% 1 LSB = ��2.5V/30/2^24
	port = [port, P];
	ExtAdc = [ExtAdc; double((2048 - E') .* 10) / 4095];	% ���ȡ�10V������Ŵ�����
	t2 = toc(t1);
	figure(9)
	ax(1) = subplot(3,1,1); plot(data)		% ��ʾ��ЩADC24����
	ax(2) = subplot(3,1,2); plot(ExtAdc)	% ��ʾ��ЩExtADC����
	PB9 = bitget(port, 1) + 0;				% ȡ��PB9(port�˿ڵ�0ͨ��)
	PC6 = bitget(port, 2) + 1;				% ȡ��PC6(port�˿ڵ�1ͨ��)
	PC7 = bitget(port, 3) + 2;				% ȡ��PC7(port�˿ڵ�2ͨ��)
	PC8 = bitget(port, 4) + 3;				% ȡ��PC8(port�˿ڵ�3ͨ��)
	ax(3) = subplot(3,1,3); hold on
	plot(PB9, 'b')							% ���ƶ˿�����
	plot(PC6, 'r')
	plot(PC7, 'k')
	plot(PC8, 'g')
	hold off
	if (t2 >= 15)							% �ɼ�5��
		break;
	end
end
linkaxes(ax,'x')							% ��������

% 3. �ر��豸���˳�
ComboClose(h)								% �ر��豸
