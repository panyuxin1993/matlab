[y_before,fs]=audioread('D:\CloudMusic\Matteo - Panama.mp3');
a=mean(y_before,2);
freq=fft(a);
m2=abs(freq);
n=length(freq);
m1=m2(n/2+1:end);%每个频率成分的大小
f = (0:n/2-1)*fs/length(freq)/1000;%频率值, 单位kHz
plot(f,m1);
xlabel('kHz');