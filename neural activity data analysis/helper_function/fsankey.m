function sankeyHdl = fsankey(varargin)
%FSANKEY sankey plot ref sankey2 from https://blog.csdn.net/slandarer/article/details/118084972
%   input data in table form, with variable 'left','n','right'
%   Detailed explanation goes here
if strcmp(get(varargin{1},'type'),'axes' )
    ax=varargin{1};
else
    ax=gca;
end
hold(ax,'on')

%��δ���ã���ͼ��ĳ�ʼֵ==================================================
prop.Color=[0,0,0];
prop.FontSize=10;
prop.FontColor=[0,0,0];
prop.Xlim=[0,1];
prop.YLim=[0,1];
prop.PieceWidth=0.15;
prop.List=[];
prop.Margin=0.05;
prop.Sep=1/8;
prop.EdgeColor=[0 0 0];
prop.Table=[];
prop.layer_str=[];

%�ӿɱ䳤�ȱ�������ȡ������Ϣ==============================================
for i=1:length(varargin)
    tempVar=varargin{i};
    if ischar(tempVar)&&length(tempVar)>1
        prop.(tempVar)=varargin{i+1};
    end
end

%�������󹹽�==============================================================
nameList=unique([prop.Table.left;prop.Table.right],'sorted');
blockMat=zeros(length(nameList));
for i=1:size(prop.Table,1)
    s=(nameList==prop.Table{i,1});
    e=(nameList==prop.Table{i,3});
    blockMat(s,e)=prop.Table{i,2};
end
totalFlow=max([sum(blockMat,1);sum(blockMat,2)'],[],1);


%����ɣ��ͼ���============================================================
List_L=prop.Table.left;
List_R=prop.Table.right;
prop.layer=[];layerRoot=[];n=1;
for i=length(List_R):-1:1%�����һ��
    if ~any(List_L==List_R(i))
        layerRoot=[layerRoot;find(nameList==List_R(i))];
    end
end
layerRoot=unique(layerRoot,'stable');
while ~isempty(List_L)
    layer_n=[];
    for i=length(List_L):-1:1%�ӵ�һ�п�ʼ�����
        if ~any(List_R ==List_L(i))
            layer_n=[layer_n;find(nameList==List_L(i))];
            List_L(i)=[];
            List_R(i)=[];
        end
    end
    layer_n=unique(layer_n,'stable');
    prop.layer(length(layer_n),n)=0;
    prop.layer(1:length(layer_n),n)=layer_n;
    n=n+1;
end
prop.layer(length(layerRoot),n)=0;
prop.layer(1:length(layerRoot),n)=layerRoot;
prop.layerNum=size(prop.layer,2);




%���Ʒ���==================================================================
baseBlockX=[0,1,1,0];
baseBlockY=[0,0,1,1];
bnul=max(sum(prop.layer~=0,1));   %block number upper limit
baseLenY=(diff(prop.YLim)-2*prop.Margin)/(bnul+(bnul-1)*prop.Sep)*bnul;
baseLenX=1;%(diff(prop.XLim)-2*prop.Margin)/(prop.layerNum-0.5);

for i=1:prop.layerNum
    tempY=prop.Margin;
    elemSet=prop.layer(prop.layer(:,i)~=0,i);
    flowSet=totalFlow(elemSet);
    offSet=(diff(prop.YLim)-2*prop.Margin-baseLenY/length(elemSet)*((length(elemSet)+(length(elemSet)-1)*prop.Sep)))/2;
    for j=1:length(elemSet)
        tempLenY=baseLenY./sum(flowSet).*flowSet(j);
        colorIndex=elemSet(j);
        sankeyHdl.block(prop.layer(j,i))=...
        fill(baseBlockX.*prop.PieceWidth+i*baseLenX-mean(baseBlockX).*prop.PieceWidth,...
            baseBlockY.*tempLenY+tempY+offSet,...
            prop.Color(colorIndex,:),'EdgeColor',prop.EdgeColor);
        tempY=tempY+tempLenY+baseLenY/length(elemSet)*prop.Sep;
    end
end

%ÿ��layer��ע��һ����epoch��
set(gca,'XTick',1:prop.layerNum,'XTickLabel',prop.layer_str);

%��������
layerList=prop.layer(:);
for i=1:length(nameList)
    for j=i:length(nameList)
        if blockMat(i,j)~=0
            Hdl_L=sankeyHdl.block(i);
            Hdl_R=sankeyHdl.block(j);
            list_L=find(blockMat(i,:)~=0);
            list_R=find(blockMat(:,j)~=0);
            [~,pl,~]=intersect(layerList,list_L(:));
            [~,pr,~]=intersect(layerList,list_R(:));
            list_L=layerList(sort(pl));
            list_R=layerList(sort(pr));
            flow_L=blockMat(i,list_L);
            flow_R=blockMat(list_R,j);
            XData_L=Hdl_L.XData;YData_L=Hdl_L.YData;
            XData_R=Hdl_R.XData;YData_R=Hdl_R.YData;
            xx=[XData_L(1:2);XData_R(1:2)]';
            k_L=find(list_L==j);
            k_R=find(list_R==i);
            yy=[YData_L(1:2)+(YData_L(3:4)-YData_L(1:2))./sum(flow_L).*sum(flow_L(1:k_L-1));
                YData_R(1:2)+(YData_R(3:4)-YData_R(1:2))./sum(flow_R).*sum(flow_R(1:k_R-1))]';
            xxq=XData_L(2):0.01:XData_R(1);
            yyq=interp1(xx,yy,xxq,'pchip');%������ֵ�����������˵���״��ע����Ҫ4���㣬����aabb������������ֵ�����γ�aa-cdef...-bb��������������״�Ĳ�ֵ
            tempColor=Hdl_L.FaceColor;
            width=(YData_R(3)-YData_R(1))./sum(flow_R).*flow_R(k_R);
             sankeyHdl.connect(i,k_L)=...
            fill([xxq,xxq(end:-1:1)],[yyq,yyq(end:-1:1)+width],tempColor,'EdgeColor','none','FaceAlpha',0.3);
        end    
    end
end

%�����ı�
for i=1:prop.layerNum
    tempY=prop.Margin;
    elemSet=prop.layer(prop.layer(:,i)~=0,i);
    flowSet=totalFlow(elemSet);
    offSet=(diff(prop.YLim)-2*prop.Margin-baseLenY/length(elemSet)*((length(elemSet)+(length(elemSet)-1)*prop.Sep)))/2;
    for j=1:length(elemSet)
        tempLenY=baseLenY./sum(flowSet).*flowSet(j);
        sankeyHdl.txt(prop.layer(j,i))=...
        text(prop.PieceWidth+prop.Margin+i*baseLenX,tempLenY/2+tempY+offSet,[string(nameList(elemSet(j)))],...
            'FontSize',prop.FontSize,'Color',prop.FontColor);
        
        tempY=tempY+tempLenY+baseLenY/length(elemSet)*prop.Sep;
    end
end
sankeyHdl.nameList=nameList';
end

