lim=0:100;
[x,y]=meshgrid(lim,lim);
z1=zeros(size(x));
z1(mean(lim),mean(lim))=1;
z2=(50-(x-mean(lim)).^2-(y-mean(lim)).^2);
z2(z2<0)=0;
z2=z2.^0.5;
z2=z2/(max(max(z2)));
figure;
subplot(1,2,1);
surf(x,y,z1);
subplot(1,2,2);
surf(x,y,z2);

