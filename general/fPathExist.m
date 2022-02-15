function bool = fPathExist(pathstr)
% �ж�·���ַ���pathstr�Ƿ��Ѿ�������MATLAB Search Path ��
% bool=1�����Ѵ��ڣ�0 ��������

ps=';'; % path seperate
MPC = regexp([matlabpath ps],['.[^' ps ']*' ps],'match')'; %MATLAB Path Cell,����Ǵ�MATLAB�Դ��ĺ��� rmpath ��ѧ���ġ�
bool = any(strcmp(pathstr,MPC));
end

