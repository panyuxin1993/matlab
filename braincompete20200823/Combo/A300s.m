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

% �� deploytool --> Application Compiler ���ɶ�����ִ�г���
% A300s.prj�����úõ�deploytool������Ŀ������ֱ�Ӵ򿪱��롣

help ComboHelp
clear all
fileName = sprintf('%4d%02d%02d_%02d%02d%02d.mat', year(now),month(now),day(now),hour(now),minute(now),int32(second(now)));

% 1. ���豸
delete(instrfindall);						% �ȹر����д���
com = ComboQuery;							% comΪ�ַ������飬�������еĵ�����һ����豸����
h = ComboOpen(com(1,:));					% �򿪵�һ���豸

% 2. �ɼ����ݡ���ο����ظ���ν��С�
figure(9), clf
nDuration = 300;							% �趨�ɼ�30��ʱ��
nLength = (nDuration+1) * 30000;
data = zeros(nLength, 1);
port = zeros(nLength, 1, 'uint8');
ExtAdc = zeros(nLength, 2);
t1 = tic;
startPoint = 1;
endPoint = 1;
while 1
	len = ComboGetLength(h);				% ȡ��һ������Ѿ��ɼ������ݳ���
	[A, E, P] = ComboGetData(h, len);	% �����������ݳ��ȣ���������
	endPoint = startPoint + len - 1;
	% data = [data, double(A).*9.9341e-09];	% 1 LSB = ��2.5V/30/2^24
	%ExtAdc = [ExtAdc; double((2048 - E') .* 10) / 4095];	% ���ȡ�10V������Ŵ�����
	data(startPoint:endPoint) = double(A).*9.9341e-09;	% 1 LSB = ��2.5V/30/2^24
	port(startPoint:endPoint) = P;
	ExtAdc(startPoint:endPoint,:) = double((2048 - E') .* 10) / 4095;	% ���ȡ�10V������Ŵ�����
	t2 = toc(t1);
	% if endPoint < 30000*3
	figure(9)
	subplot(3,1,1), clf											% �������figure�е����ݣ����򵽺�����Խ��Խ��
	subplot(3,1,2), clf
	subplot(3,1,3), clf
	leadPoint = max(1, endPoint-2*30000);
	Ncomp = int32(floor(single(endPoint - leadPoint) / 30));	% 30������ѹ�����������
	dataRange = endPoint-Ncomp*30+1:endPoint;					% ��ʾ�����ݷ�Χ
	viewEEG = mean(reshape(data(dataRange), 30, Ncomp), 1);		% EEG����ѹ��30��
	viewExt = squeeze(mean(reshape(ExtAdc(dataRange,:), 30, Ncomp, 2), 1));	% ͬ���ź�ѹ��30��
	ax(1) = subplot(3,1,1); plot(viewEEG)						% ��ʾADC24����
	ax(2) = subplot(3,1,2); plot(viewExt)						% ��ʾͬ���ź�����
	PB9 = bitget(port(dataRange), 1) + 0;						% ȡ��PB9(port�˿ڵ�0ͨ��)
	PC6 = bitget(port(dataRange), 2) + 1;						% ȡ��PC6(port�˿ڵ�1ͨ��)
	PC7 = bitget(port(dataRange), 3) + 2;						% ȡ��PC7(port�˿ڵ�2ͨ��)
	PC8 = bitget(port(dataRange), 4) + 3;						% ȡ��PC8(port�˿ڵ�3ͨ��)
	ax(3) = subplot(3,1,3); hold on
	plot(mean(reshape(PB9, 30, Ncomp), 1), 'b')		% ���ƶ˿�����
	plot(mean(reshape(PC6, 30, Ncomp), 1), 'r')
	plot(mean(reshape(PC7, 30, Ncomp), 1), 'k')
	plot(mean(reshape(PC8, 30, Ncomp), 1), 'g')
	hold off
	% end
	if (t2 >= nDuration)					% �ɼ�nDuration��
		break;
	end
	startPoint = endPoint + 1;
	fprintf('.');
	pause(0.001)							% ��������
end
linkaxes(ax,'x')							% ��������

% 3. �ر��豸���˳�
ComboClose(h)								% �ر��豸
save(fileName)
