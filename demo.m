a=1:10;
b=[2,4,6];
for i=1:length(a)
    if sum(a(i)==b)>0
        a(i)=0;
    end
end
a(a==0)=[];

for j=1:length(b)
    a(a==b(i))=[];
end