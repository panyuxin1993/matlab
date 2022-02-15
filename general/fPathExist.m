function bool = fPathExist(pathstr)
% 判断路径字符串pathstr是否已经存在于MATLAB Search Path 中
% bool=1代表已存在，0 代表不存在

ps=';'; % path seperate
MPC = regexp([matlabpath ps],['.[^' ps ']*' ps],'match')'; %MATLAB Path Cell,这句是从MATLAB自带的函数 rmpath 中学到的。
bool = any(strcmp(pathstr,MPC));
end

