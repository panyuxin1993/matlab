  function BrainCompeteMain
global  obj1 obj2 stop com Closetag isstop;
 Closetag = 0; isstop = 1;
% 1. �򿪴����豸
    if ~isempty(instrfind)
        delete(instrfindall);
    end
    com = ComboQuery							% comΪ�ַ������飬�������еĵ�����һ����豸����
    if isempty(com)
        msgbox('û���豸���ӣ���ô�氡��');
        return;
    elseif length(com(:,1))<2
        msgbox('ֻ��һ���豸���ҶԿ�ɶ����');
        return;
    end
    obj1 = 0;
    obj2 = 0;
    BrainCompete;
    
end
