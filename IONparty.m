path='D:\党员\脑智中心第一届两委选举';
file=[path,filesep,'酝酿提名统计.xlsx'];
T=readtable(file,'ReadVariableNames',true,'ReadRowNames',true);
strvar={'党委','纪委'};
titlestr={'姓名','票数'};
for ncol=1:2
    varname='';
    for nrow=1:size(T,1)
        namescell=table2cell(T(nrow,ncol));
        names=strsplit(namescell{1},{' ','，','、',' ','  '},'CollapseDelimiters',1);
        if length(names)==1
            names=strsplit(names{1});
        end
        names=strtrim(names);
        varname=union(varname,names);
    end
%     xlsrange=['A1:A',num2str(length(varname))];
%     xlswrite(file,varname',ncol+1,xlsrange);
end
%%
%手动将人名调整好再程序计数
path='D:\党员\脑智中心第一届两委选举';
file=[path,filesep,'酝酿提名统计.xlsx'];
T=readtable(file,'ReadVariableNames',true,'ReadRowNames',true);
for ncol=1:2
    [~,~,varname]=xlsread(file,ncol+1);
    varnum=cell(size(varname));
    for i=1:length(varname)
        if isnan(varname{i})
            varname{i}='';
        end
    end
    data=table2cell(T(:,ncol));
    for iname=1:length(varname)
        temp=cellfun(@(x) contains(x,varname{iname}),data);
        varnum{iname}=sum(temp);
    end
    cellout=cell(length(varnum),2);
    [cellout(:,1)]=deal(varname);
    [cellout(:,2)]=deal(varnum);
    xlsrange=['A2:B',num2str(length(varname)+1)];
    xlswrite(file,cellout,ncol+1,xlsrange);
end
    
        
        