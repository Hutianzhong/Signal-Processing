clear all;
close all;
clc;

str = '/Users/hutianzhong/Desktop/��������Ӧ��Ȩ��߶������̬�˲�����й���������ȡ�о�/��������/213.mat';
c = load(str);
x = c.X213_FE_time(1:48000)';
N = 48000; %��������
l = 0:N-1; 
fs = 48000; %����Ƶ��
dt = 1/fs; %ʱ����С�������ʱ��ֱ���
t = l*dt; %����ʱ�䳤��
high = 8;%���߶�
fxy	 = 800; %Ƶ��ͼ��ʾ��Χ
Es = 800; %������Ƶ��ȡֵ��Χ

y1 = x;
fprintf('Rate(y1)��%s\n',RATE(y1,N,Es));

figure(1); 
plot(t,y1);
title('ԭʼ�ź� '); 
xlabel('ʱ��/s');
ylabel('��ֵA/(m/s^2)');

figure(2);
NFFT = 2^nextpow2(N);% ����N�������2���������ݵĵ���
Y = fft(y1,NFFT)/N;
f = fs/2*linspace(0,1,NFFT/2+1);
plot(f(4:fxy),2*abs(Y(4:fxy)));
title('Ƶ����� '); 
xlabel('Ƶ��f/Hz'); 
ylabel('��ֵA/(m/s^2)');

% WMAVG
for length = 2 : 2
    %C(+-)
    Fc = y1;
    for i = 2 : high
        Fc =[Fc; y1];
    end
    Fc = operator(Fc,'+',length,high);
    Fc = operator(Fc,'-',length,high);
    
    %O(-+)
    Fo = y1;
    for i = 2 : high
        Fo =[Fo; y1];
    end
    Fo = operator(Fo,'-',length,high);
    Fo = operator(Fo,'+',length,high);
end
w = weight(y1,(Fc + Fo) / 2,high);
WMAVG = w * (Fc + Fo) / 2;
figure(8); plot(t,WMAVG); title('WMAVG �ź� '); xlabel('ʱ��/s'); ylabel('��ֵA/(m/s^2)');
WMAVG = hilbert(WMAVG);%ϣ�����ر任 
WMAVG = abs(WMAVG);%�����ź� 
figure(9); NFFT = 2^nextpow2(N); Y = fft(WMAVG,NFFT)/N; f = fs/2*linspace(0,1,NFFT/2+1); plot(f(4:fxy),2*abs(Y(4:fxy))); xlabel('Ƶ��f/Hz'); ylabel('��ֵA/(m/s^2)');
fprintf('Rate(WMAVG)��%s\n',RATE(WMAVG,N,Es));

% WMDIF
for length = 2 : 2
    %C(+-)
    Fc = y1;
    for i = 2 : high
        Fc =[Fc; y1];
    end
    Fc = operator(Fc,'+',length,high);
    Fc = operator(Fc,'-',length,high);
    
    %O(-+)
    Fo = y1;
    for i = 2 : high
        Fo =[Fo; y1];
    end
    Fo = operator(Fo,'-',length,high);
    Fo = operator(Fo,'+',length,high);
end
w = weight(y1,(Fc - Fo) / 2,high);
WMDIF = w * (Fc - Fo) / 2;
figure(10); plot(t,WMDIF); title('WMDIF �ź� '); xlabel('ʱ��/s'); ylabel('��ֵA/(m/s^2)');
WMDIF = hilbert(WMDIF);%ϣ�����ر任 
WMDIF = abs(WMDIF);%�����ź� 
figure(10); plot(t,WMDIF); title('WMDIF �ź� '); xlabel('ʱ��/s'); ylabel('��ֵA/(m/s^2)');
figure(11); NFFT = 2^nextpow2(N); Y = fft(WMDIF,NFFT)/N; f = fs/2*linspace(0,1,NFFT/2+1); plot(f(4:fxy),2*abs(Y(4:fxy))); xlabel('Ƶ��f/Hz'); ylabel('��ֵA/(m/s^2)');
fprintf('Rate(WMDIF)��%s\n',RATE(WMDIF,N,Es));

%WMOCCO
for length = 2 : 2
    %OC(-++-)
    Foc = y1;
    for i = 2 : high
        Foc =[Foc; y1];
    end
    Foc = operator(Foc,'-',length,high);
    Foc = operator(Foc,'+',length,high);
    Foc = operator(Foc,'+',length,high);
    Foc = operator(Foc,'-',length,high);
    
    %CO(+--+)
    Fco = y1;
    for i = 2 : high
        Fco =[Fco; y1];
    end
    Fco = operator(Fco,'+',length,high);
    Fco = operator(Fco,'-',length,high);
    Fco = operator(Fco,'-',length,high);
    Fco = operator(Fco,'+',length,high);
end
w = weight(y1,(Foc + Fco) / 2,high);
WMOCCO = w * (Foc + Fco) / 2;
figure(12); plot(t,WMOCCO); title('WMOCCO �ź� '); xlabel('ʱ��/s'); ylabel('��ֵA/(m/s^2)');
WMOCCO = hilbert(WMOCCO);%ϣ�����ر任 
WMOCCO = abs(WMOCCO);%�����ź� 
figure(13); NFFT = 2^nextpow2(N); Y = fft(WMOCCO,NFFT)/N; f = fs/2*linspace(0,1,NFFT/2+1); plot(f(4:fxy),2*abs(Y(4:fxy))); xlabel('Ƶ��f/Hz'); ylabel('��ֵA/(m/s^2)');
fprintf('Rate(WMOCCO)��%s\n',RATE(WMOCCO,N,Es));

%WMMG
for length = 2 : 2
    %D(+)
    Fd = y1;
    for i = 2 : high
        Fd =[Fd; y1];
    end
    Fd = operator(Fd,'+',length,high);
    
    %E(-+)
    Fe = y1;
    for i = 2 : high
        Fe =[Fe; y1];
    end
    Fe = operator(Fe,'-',length,high);
end
w = weight(y1,Fd - Fe,high);
WMMG = w * (Fd - Fe);
figure(14); plot(t,WMMG); title('WMMG �ź� '); xlabel('ʱ��/s'); ylabel('��ֵA/(m/s^2)');
WMMG = hilbert(WMMG);%ϣ�����ر任 
WMMG = abs(WMMG);%�����ź� 
figure(15); NFFT = 2^nextpow2(N); Y = fft(WMMG,NFFT)/N; f = fs/2*linspace(0,1,NFFT/2+1); plot(f(4:fxy),2*abs(Y(4:fxy))); xlabel('Ƶ��f/Hz'); ylabel('��ֵA/(m/s^2)');
fprintf('Rate(WMMG)��%s\n',RATE(WMMG,N,Es));

%WMEDMF
tek_WMEDMF = 0;
for length = 2 : 2
    
    %CDE(+-+-)
    Fcde = y1;
    for i = 2 : high
        Fcde =[Fcde; y1];
    end
    
    Fcde = operator(Fcde,'+',length,high);
    Fcde = operator(Fcde,'-',length,high);
    Fcde = operator(Fcde,'+',length,high);
    Fcde = operator(Fcde,'-',length,high);
    
    %DCE(++--)
    Fdce = y1;
    for i = 2 : high
        Fdce =[Fdce; y1];
    end
    
    Fdce = operator(Fdce,'+',length,high);
    Fdce = operator(Fdce,'+',length,high);
    Fdce = operator(Fdce,'-',length,high);
    Fdce = operator(Fdce,'-',length,high);
    
    %EOD(--++)
    Feod = y1;
    for i = 2 : high
        Feod =[Feod; y1];
    end
    
    Feod = operator(Feod,'-',length,high);
    Feod = operator(Feod,'-',length,high);
    Feod = operator(Feod,'+',length,high);
    Feod = operator(Feod,'+',length,high);
    
    %OED(-+-+)
    Foed = y1;
    for i = 2 : high
        Foed =[Foed; y1];
    end
    
    Foed = operator(Foed,'-',length,high);
    Foed = operator(Foed,'+',length,high);
    Foed = operator(Foed,'-',length,high);
    Foed = operator(Foed,'+',length,high);
    
    if (tek_WMEDMF < kurtosis(weight(y1,(Fcde - Feod) + (Fdce - Foed) / 2,high) * ((Fcde - Feod) + (Fdce - Foed)) / 2))
    WMEDMF = ((Fcde - Feod) + (Fdce - Foed)) / 2;
    tek_WMEDMF = kurtosis(weight(y1,(Fcde - Feod) + (Fdce - Foed) / 2,high) * ((Fcde - Feod) + (Fdce - Foed)) / 2);
    end
end
wmedmf = pso(high,WMEDMF) * WMEDMF;
figure(20); plot(t,wmedmf); title('WMEDMF �ź�'); xlabel('ʱ��/s'); ylabel('��ֵA/(m/s^2)');
wmedmf = hilbert(wmedmf);%ϣ�����ر任 
wmedmf = abs(wmedmf);%�����ź� 
figure(21); NFFT = 2^nextpow2(N); Y = fft(wmedmf,NFFT)/N; f = fs/2*linspace(0,1,NFFT/2+1); plot(f(4:fxy),2*abs(Y(4:fxy))); title('Ƶ����� ');  xlabel('Ƶ��f/Hz');  ylabel('��ֵA/(m/s^2)');
fprintf('Rate(wmedmf)��%s\n',RATE(wmedmf,N,Es));

hil = hilbert(y1); %ϣ�����ر任 
hil = abs(hil); %�����ź� 
figure(22); plot(t,hil); title('hilbert �����ź�'); xlabel('ʱ��/s'); ylabel('��ֵA/(m/s^2)');
figure(23); NFFT = 2^nextpow2(N); Y = fft(hil,NFFT)/N; f = fs/2*linspace(0,1,NFFT/2+1); plot(f(4:fxy),2*abs(Y(4:fxy))); xlabel('Ƶ��f/Hz'); ylabel('��ֵA/(m/s^2)');
fprintf('Rate(hilbert)��%s\n',RATE(hil,N,Es));
