clear all;
close all;
clc;

str = '/Users/hutianzhong/Desktop/基于自适应加权多尺度组合形态滤波的轴承故障特征提取研究/西楚数据/213.mat';
c = load(str);
x = c.X213_FE_time(1:48000)';
N = 48000; %采样点数
l = 0:N-1; 
fs = 48000; %采样频率
dt = 1/fs; %时域最小间隔，即时域分辨率
t = l*dt; %采样时间长度
high = 8;%最大尺度
fxy	 = 800; %频谱图显示范围
Es = 800; %能量比频率取值范围

y1 = x;
fprintf('Rate(y1)：%s\n',RATE(y1,N,Es));

figure(1); 
plot(t,y1);
title('原始信号 '); 
xlabel('时间/s');
ylabel('幅值A/(m/s^2)');

figure(2);
NFFT = 2^nextpow2(N);% 计算N点最近的2的整数次幂的点数
Y = fft(y1,NFFT)/N;
f = fs/2*linspace(0,1,NFFT/2+1);
plot(f(4:fxy),2*abs(Y(4:fxy)));
title('频域分析 '); 
xlabel('频率f/Hz'); 
ylabel('幅值A/(m/s^2)');

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
figure(8); plot(t,WMAVG); title('WMAVG 信号 '); xlabel('时间/s'); ylabel('幅值A/(m/s^2)');
WMAVG = hilbert(WMAVG);%希尔伯特变换 
WMAVG = abs(WMAVG);%包络信号 
figure(9); NFFT = 2^nextpow2(N); Y = fft(WMAVG,NFFT)/N; f = fs/2*linspace(0,1,NFFT/2+1); plot(f(4:fxy),2*abs(Y(4:fxy))); xlabel('频率f/Hz'); ylabel('幅值A/(m/s^2)');
fprintf('Rate(WMAVG)：%s\n',RATE(WMAVG,N,Es));

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
figure(10); plot(t,WMDIF); title('WMDIF 信号 '); xlabel('时间/s'); ylabel('幅值A/(m/s^2)');
WMDIF = hilbert(WMDIF);%希尔伯特变换 
WMDIF = abs(WMDIF);%包络信号 
figure(10); plot(t,WMDIF); title('WMDIF 信号 '); xlabel('时间/s'); ylabel('幅值A/(m/s^2)');
figure(11); NFFT = 2^nextpow2(N); Y = fft(WMDIF,NFFT)/N; f = fs/2*linspace(0,1,NFFT/2+1); plot(f(4:fxy),2*abs(Y(4:fxy))); xlabel('频率f/Hz'); ylabel('幅值A/(m/s^2)');
fprintf('Rate(WMDIF)：%s\n',RATE(WMDIF,N,Es));

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
figure(12); plot(t,WMOCCO); title('WMOCCO 信号 '); xlabel('时间/s'); ylabel('幅值A/(m/s^2)');
WMOCCO = hilbert(WMOCCO);%希尔伯特变换 
WMOCCO = abs(WMOCCO);%包络信号 
figure(13); NFFT = 2^nextpow2(N); Y = fft(WMOCCO,NFFT)/N; f = fs/2*linspace(0,1,NFFT/2+1); plot(f(4:fxy),2*abs(Y(4:fxy))); xlabel('频率f/Hz'); ylabel('幅值A/(m/s^2)');
fprintf('Rate(WMOCCO)：%s\n',RATE(WMOCCO,N,Es));

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
figure(14); plot(t,WMMG); title('WMMG 信号 '); xlabel('时间/s'); ylabel('幅值A/(m/s^2)');
WMMG = hilbert(WMMG);%希尔伯特变换 
WMMG = abs(WMMG);%包络信号 
figure(15); NFFT = 2^nextpow2(N); Y = fft(WMMG,NFFT)/N; f = fs/2*linspace(0,1,NFFT/2+1); plot(f(4:fxy),2*abs(Y(4:fxy))); xlabel('频率f/Hz'); ylabel('幅值A/(m/s^2)');
fprintf('Rate(WMMG)：%s\n',RATE(WMMG,N,Es));

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
figure(20); plot(t,wmedmf); title('WMEDMF 信号'); xlabel('时间/s'); ylabel('幅值A/(m/s^2)');
wmedmf = hilbert(wmedmf);%希尔伯特变换 
wmedmf = abs(wmedmf);%包络信号 
figure(21); NFFT = 2^nextpow2(N); Y = fft(wmedmf,NFFT)/N; f = fs/2*linspace(0,1,NFFT/2+1); plot(f(4:fxy),2*abs(Y(4:fxy))); title('频域分析 ');  xlabel('频率f/Hz');  ylabel('幅值A/(m/s^2)');
fprintf('Rate(wmedmf)：%s\n',RATE(wmedmf,N,Es));

hil = hilbert(y1); %希尔伯特变换 
hil = abs(hil); %包络信号 
figure(22); plot(t,hil); title('hilbert 包络信号'); xlabel('时间/s'); ylabel('幅值A/(m/s^2)');
figure(23); NFFT = 2^nextpow2(N); Y = fft(hil,NFFT)/N; f = fs/2*linspace(0,1,NFFT/2+1); plot(f(4:fxy),2*abs(Y(4:fxy))); xlabel('频率f/Hz'); ylabel('幅值A/(m/s^2)');
fprintf('Rate(hilbert)：%s\n',RATE(hil,N,Es));
