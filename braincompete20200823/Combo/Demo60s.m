% �ɼ�60���ʾ������

% 1. ���豸
com = ComboQuery;						% comΪ�ַ������飬�������еĵ�����һ����豸����
h = ComboOpen(com(1,:));				% �򿪵�һ���豸
% com(1,:)�ں���һ���豸��com(2,:)�ں��ڶ����豸�������豸������������
% 2. �ɼ����ݡ���ο����ظ���ν��С�
data = [];
tic
while(length(data)<60*30000)			% �ɼ�60��
len = ComboGetLength(h);				% ȡ��һ������Ѿ��ɼ������ݳ���
A = ComboGetData(h, len);				% �����������ݳ��ȣ���������
data = [data, double(A).*9.9341e-09];	% 1 LSB = ��2.5V/30/2^24
end
toc
figure(9), plot(data)					% ��ʾ��Щ����
% 3. �ر��豸���˳�
ComboClose(h)								% �ر��豸
clear h