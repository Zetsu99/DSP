% Figure 11.17: Analog Lowpass Elliptic Filter Design
%               using Matlab
% Fp = 40 Hz; Fs = 50 Hz; Ap = 0.1; As = 50;

clc; close all; %echo on;

% Given Design Parameters
Fp = 40; Fs = 50; Ap = 0.1; As = 50;
Omegap = 2*pi*Fp; Omegas = 2*pi*Fs;
% Analog Design Parameters (Eq. 10.9)
epsilon = sqrt(10^(0.1*Ap)-1); A = 10^(0.05*As);
Rp = 1/sqrt(1+epsilon^2);

% Design using SP Toolbox functions
[N, Omegac] = ellipord(Omegap, Omegas, Ap, As, 's');
N
Omegac/(2*pi)
[C,D] = ellip(N,Ap,As,Omegac,'s');
zk = roots(C); pk = roots(D);
% D1 = real(poly(pk(1:2)));
% D2 = real(poly(pk(3)));

%% Design Plots
Hf11_05 = figure('units','inches','position',[1,1,5.8,3.6],...
    'paperunits','inches','paperposition',[0,0,5.8,3.6]);
set(Hf11_05,'number','off','name','Ex11.5');

Fmax = 100; Ommax = 2*pi*Fmax; F = linspace(0,Fmax,101);
Om = 2*pi*F; H = freqs(C,D,Om);
Hmag = abs(H); Hpha = angle(H); Hdb = 20*log10(Hmag);
Hgdl = -diff(unwrap(Hpha))./diff(Om); Hgdl = [Hgdl,Hgdl(end)];
Hgdl = medfilt1(Hgdl,3);

subplot(2,2,1); % Magnitude Response
plot(F,Hmag,'b','linewidth',1); axis([0,Fmax,0,1.1]);
xlabel('Frequency in Hz'); ylabel('Magnitude');
title('Magnitude Response');
set(gca,'xtick',[0,Fp,Fs,Fmax]); set(gca,'ytick',[0,1]); grid; 

subplot(2,2,2); % Log-Magnitude Response in dB
plot(F,Hdb,'b','linewidth',1); axis([0,Fmax,-80,10]);
xlabel('Frequency in Hz'); ylabel('Decibels');
title('Log-Magnitude Response');
set(gca,'xtick',[0,Fp,Fs,Fmax]);
set(gca,'ytick',[-80,-As,0]); grid; 

subplot(2,2,3); % Group-Delay Response in Samples
plot(F,Hgdl,'b','linewidth',1); axis([0,Fmax,0,0.15]); 
xlabel('Frequency in Hz'); ylabel('Samples');
title('Group Delay Response');
set(gca,'xtick',[0,Fp,Fs,Fmax]);
set(gca,'ytick',(0:0.05:0.15)); grid;

subplot(2,2,4); % Zoomed Passband Magnitude Response
Fmax = Fp; F = linspace(0,Fmax,201);
H = freqs(C,D,2*pi*F); Hmag = abs(H);
plot(F,Hmag,'b','linewidth',1); 
axis([0,Fp,0.988,1.001]);
xlabel('Frequency in Hz'); ylabel('Magnitude');
title('Zoomed Magnitude Response');
set(gca,'xtick',(0:20:Fp)); 
set(gca,'ytick',[Rp,1]);

% Print Plot
print -depsc2 ../artfiles/1117_elliptic3.eps;