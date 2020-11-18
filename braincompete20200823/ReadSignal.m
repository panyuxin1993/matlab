function [ data, Npoint] = ReadSignal( h )
%READSIGNAL  Read signal from Combo Amplifier
%   Input, obj, the COM port connected to Combo Amplifier
%   Output, data��  signal
%           Npoint, number of signal point
% Citing from YQ Wen, Modified by Zhiwei Wang
% 2015/06/27

    len = ComboGetLength(h);				% ȡ��һ������Ѿ��ɼ������ݳ���
	[data Adc1 DIO] = ComboGetData(h, len);	% �����������ݳ��ȣ���������
    data = double(data).*9.9341e-09;
    Npoint = length(data);
end

