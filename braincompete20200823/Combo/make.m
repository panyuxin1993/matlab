% ����һ����豸��ComboXxxͨѶ������
% ����ǵ�һ��ʹ��mex���룬��Ҫ������ָ�����ñ��뻷��
%	mex -setup
% 

% for Debug version
fprintf('���� ComboQuery.cpp ...\n');
mex -g ComboQuery.cpp
fprintf('���� ComboOpen.cpp ...\n');
mex -g ComboOpen.cpp
fprintf('���� ComboGetLength.cpp ...\n');
mex -g ComboGetLength.cpp
fprintf('���� ComboGetData.cpp ...\n');
mex -g ComboGetData.cpp -lWinmm.lib
fprintf('���� ComboClose.cpp ...\n');
mex -g ComboClose.cpp

% % for Release version
% fprintf('���� ComboQuery.cpp ...\n');
% mex ComboQuery.cpp
% fprintf('���� ComboOpen.cpp ...\n');
% mex ComboOpen.cpp
% fprintf('���� ComboGetLength.cpp ...\n');
% mex ComboGetLength.cpp
% fprintf('���� ComboGetData.cpp ...\n');
% mex ComboGetData.cpp -lWinmm.lib
% fprintf('���� ComboClose.cpp ...\n');
% mex ComboClose.cpp
% delete *.pdb
