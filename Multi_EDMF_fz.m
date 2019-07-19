clear all;
close all;
clc;

str = 'y202.txt';
c = load(str);
x = c.';
N = 1024; %采样点数
l = 1:N; 
fs = 1024; %采样频率
dt = 1/fs; %时域最小间隔，即时域分辨率
t = l*dt;  %采样时间长度
high = 8; %最大尺度
fxy	 = 200; %频谱图显示范围
Es = 200; %能量比频率取值范围

y1 = x;

% fprintf('Rate(y1)：%s\n',RATE(y1,N,Es));

figure(1); 
plot(y1);
title('原始信号'); 
xlabel('Sample Number');
ylabel('Amplitude(m/s^2)');

figure(2);
NFFT = 2^nextpow2(N);% 计算N点最近的2的整数次幂的点数
Y = fft(y1,NFFT)/N;
f = fs/2*linspace(0,1,NFFT/2+1);
plot(f(2:fxy),2*abs(Y(2:fxy)));
title('频域分析'); 
xlabel('Frequency(Hz)'); 
ylabel('Amplitude(m/s^2)');

refSE = strel('line',8,0);

figure(111);
plot(imdilate(y1,refSE));

figure(222);
plot(imerode(y1,refSE));

figure(333);
plot(imclose(y1,refSE));

figure(444);
plot(imopen(y1,refSE));

kur_edmf = 0;
for i = 2 : 16
refSE = strel('line',i,0);
Fcde = imerode(imdilate(imerode(imdilate(y1,refSE),refSE),refSE),refSE);
Fdce = imerode(imerode(imdilate(imdilate(y1,refSE),refSE),refSE),refSE);
Feod = imdilate(imdilate(imerode(imerode(y1,refSE),refSE),refSE),refSE);
Foed = imdilate(imerode(imdilate(imerode(y1,refSE),refSE),refSE),refSE);

    if (kur_edmf < kurtosis((Fdce - Foed) + (Fcde - Feod) / 2))
    edmf = ((Fdce - Foed) + (Fcde - Feod)) / 2;
    kur_edmf = kurtosis((Fdce - Foed) + (Fcde - Feod) / 2);
    end
end
%fprintf('Rate(edmf)：%s\n',RATE(edmf,N,Es)); 
figure(6); plot(edmf); title('edmf 信号 '); xlabel('Frequency(Hz)'); ylabel('Amplitude(m/s^2)');
figure(7); NFFT = 2^nextpow2(N); Y = fft(edmf,NFFT)/N; f = fs/2*linspace(0,1,NFFT/2+1); plot(f(2:fxy),2*abs(Y(2:fxy))); xlabel('频率f/Hz'); ylabel('幅值A/(m/s^2)');
qu1=2*abs(Y(1:NFFT/2+1));
single = [qu1(13) qu1(25) qu1(37) qu1(49) qu1(61) qu1(73)];

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
    
    if (tek_WMEDMF < kurtosis(weight(y1,(Fcde - Feod) + (Fcde - Foed) / 2,high) * ((Fcde - Feod) + (Fcde - Foed)) / 2))
    WMEDMF = ((Fcde - Feod) + (Fcde - Foed)) / 2;
    tek_WMEDMF = kurtosis(weight(y1,(Fcde - Feod) + (Fcde - Foed) / 2,high) * ((Fcde - Feod) + (Fcde - Foed)) / 2);
    end
    
    %fprintf('kur(wmedmf)：%s\n',kurtosis(weight(y1,(Fcde - Feod) + (Fcde - Foed) / 2,high) * ((Fcde - Feod) + (Fcde - Foed)) / 2));
end

figure(101);
plot(Fcde(6,:));

figure(102);
plot(Feod(6,:));

figure(103);
plot(Foed(6,:));

figure(104);
plot(WMEDMF(6,:));

% figure(101);
% NFFT = 2^nextpow2(N);% 计算N点最近的2的整数次幂的点数
% [Y,X] = meshgrid(1:1024,1:12);
% Y = fs/2*linspace(0,1,NFFT/2+1);
% for i = 1 : high
% Z(i,:) = fft(WMEDMF(i,:),NFFT)/N;
% end
% plot3(X(1,2:fxy),Y(2:fxy),2*abs(Z(1,2:fxy)),'b-');
% hold on;
% plot3(X(2,2:fxy),Y(2:fxy),2*abs(Z(2,2:fxy)),'c-');
% hold on;
% plot3(X(3,2:fxy),Y(2:fxy),2*abs(Z(3,2:fxy)),'k-');
% hold on;
% plot3(X(4,2:fxy),Y(2:fxy),2*abs(Z(4,2:fxy)),'m-');
% hold on;
% plot3(X(5,2:fxy),Y(2:fxy),2*abs(Z(5,2:fxy)),'r-');
% hold on;
% plot3(X(6,2:fxy),Y(2:fxy),2*abs(Z(6,2:fxy)),'y-');
% hold on;
% plot3(X(7,2:fxy),Y(2:fxy),2*abs(Z(7,2:fxy)),'b-');
% hold on;
% plot3(X(8,2:fxy),Y(2:fxy),2*abs(Z(8,2:fxy)),'c-');
% hold on;
% plot3(X(9,2:fxy),Y(2:fxy),2*abs(Z(9,2:fxy)),'k-');
% hold on;
% plot3(X(10,2:fxy),Y(2:fxy),2*abs(Z(10,2:fxy)),'m-');
% hold on;
% plot3(X(11,2:fxy),Y(2:fxy),2*abs(Z(11,2:fxy)),'r-');
% hold on;
% plot3(X(12,2:fxy),Y(2:fxy),2*abs(Z(12,2:fxy)),'y-');
% hold on;
% title('频域分析 '); 
% xlabel('SE尺度');
% ylabel('频率f/Hz'); 
% zlabel('幅值A/(m/s^2)');

wmedmf = pso(high,WMEDMF) * WMEDMF;
figure(20); plot(wmedmf); title('WMEDMF 信号'); xlabel('Frequency(Hz)'); ylabel('Amplitude(m/s^2)');
figure(21); NFFT = 2^nextpow2(N); Y = fft(wmedmf,NFFT)/N; f = fs/2*linspace(0,1,NFFT/2+1); plot(f(2:fxy),2*abs(Y(2:fxy))); title('频域分析 ');  xlabel('频率f/Hz');  ylabel('幅值A/(m/s^2)');
% fprintf('Rate(wmedmf)：%s\n',RATE(wmedmf,N,Es));
% FFT 幅值
qu2=2*abs(Y(1:NFFT/2+1));
multi = [qu2(13) qu2(25) qu2(37) qu2(49) qu2(61) qu2(73)];


% insnr=[13 25 37 49 61 73];
% figure(25); 
%     h1=plot(insnr,single,'m-*','LineWidth',1.5);
%     hold on;
%     h2=plot(insnr,multi,'b-+','LineWidth',1.5);
%     hold on;
%     xlabel('频率f/Hz');ylabel('幅值A/(m/s^2)');
%     legend([h2,h1],'WMEDMF','EDMF');
