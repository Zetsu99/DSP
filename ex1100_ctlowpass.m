% Script for Analog Lowpass Butterworth Filter Design
% Fp = 40 Hz; Fs = 50 Hz; Ap = 1; As = 30;

clc; close all; %echo on;

% Given Design Parameters
Fp = 40; Fs = 50; Ap = 1; As = 30;
Omegap = 2*pi*Fp; Omegas = 2*pi*Fs;
% Analog Design Parameters (Eq. 10.9)
epsilon = sqrt(10^(0.1*Ap)-1); A = 10^(0.05*As);
Rp = 1/sqrt(1+epsilon^2);

% Plotting Parameters
Fmax = 100; Ommax = 2*pi*Fmax; 
F = linspace(0,Fmax,501); Om = 2*pi*F; 

%% Designs using SP Toolbox functions

% Butterworth Approximation
disp('** Butterworth **');
[N, Omegac] = buttord(Omegap, Omegas, Ap, As, 's');
N
[C,D] = butter(N,Omegac,'s'); H = freqs(C,D,Om);
Hmag = abs(H); Hpha = angle(H); Hdb = 20*log10(Hmag);
Hgdl = -diff(unwrap(Hpha))./diff(Om); Hgdl = [Hgdl,Hgdl(end)];
Hgdl = medfilt1(Hgdl,3);
pk = roots(D); sigmak = real(pk); Omegak = imag(pk);
%% Design Plots
Hf11_05a = figure('units','inches','position',[1,1,5.8,3.6],...
    'paperunits','inches','paperposition',[0,0,5.8,3.6]);
set(Hf11_05a,'number','off','name','Ex11.5: Butterworth');

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

subplot(2,2,4); % Pole-Zero Plot
plot(sigmak,Omegak,'bx','linewidth',1.5); hold on;
plot([-500,500],[0,0],'k','linewidth',0.75);
plot([0,0],[-500,500],'k','linewidth',0.75);
%plot(Omegac*cos(0.5*pi*(1:0.01:3)),Omegac*sin(0.5*pi*(1:0.01:3)),'k:');
axis([-500,500,-500,500]); axis square;
xlabel('Real-axis'); ylabel('Imaginary-axis');
title('Pole Locations');
set(gca,'xtick',[-500,0,500],'ytick',[-500,0,500]);
% Print Plot
%print -depsc2 ../artfiles/1116_ex1105a.eps;

% Chebyshev-I Approximation
disp('** Chebyshev-I **');
[N, Omegac] = cheb1ord(Omegap, Omegas, Ap, As, 's');
N
[C,D] = cheby1(N,Ap,Omegac,'s'); H = freqs(C,D,Om);
Hmag = abs(H); Hpha = angle(H); Hdb = 20*log10(Hmag);
Hgdl = -diff(unwrap(Hpha))./diff(Om); Hgdl = [Hgdl,Hgdl(end)];
Hgdl = medfilt1(Hgdl,3);
pk = roots(D); sigmak = real(pk); Omegak = imag(pk);
%% Design Plots
Hf11_05b = figure('units','inches','position',[7,1,5.8,3.6],...
    'paperunits','inches','paperposition',[0,0,5.8,3.6]);
set(Hf11_05b,'number','off','name','Ex11.5: Chebyshev-I');

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

subplot(2,2,4); % Pole-Zero Plot
plot(sigmak,Omegak,'bx','linewidth',1.5); hold on;
plot([-500,500],[0,0],'k','linewidth',0.75);
plot([0,0],[-500,500],'k','linewidth',0.75);
%plot(Omegac*cos(0.5*pi*(1:0.01:3)),Omegac*sin(0.5*pi*(1:0.01:3)),'k:');
axis([-500,500,-500,500]); axis square;
xlabel('Real-axis'); ylabel('Imaginary-axis');
title('Pole Locations');
set(gca,'xtick',[-500,0,500],'ytick',[-500,0,500]);
% Print Plot
%print -depsc2 ../artfiles/1116_ex1105b.eps;

% Chebyshev-II Approximation
disp('** Chebyshev-II **');
[N, Omegac] = cheb2ord(Omegap, Omegas, Ap, As, 's');
N
[C,D] = cheby2(N,As,Omegac,'s'); H = freqs(C,D,Om);
Hmag = abs(H); Hpha = angle(H); Hdb = 20*log10(Hmag);
Hgdl = -diff(unwrap(Hpha))./diff(Om); Hgdl = [Hgdl,Hgdl(end)];
Hgdl = medfilt1(Hgdl,3);
pk = roots(D); sigmapk = real(pk); Omegapk = imag(pk);
zk = roots(C); sigmazk = real(zk); Omegazk = imag(zk);
%% Design Plots
Hf11_05c = figure('units','inches','position',[1,5.5,5.8,3.6],...
    'paperunits','inches','paperposition',[0,0,5.8,3.6]);
set(Hf11_05c,'number','off','name','Ex11.5: Chebyshev-II');

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

subplot(2,2,4); % Pole-Zero Plot
plot(sigmapk,Omegapk,'bx','linewidth',1.5); hold on;
plot(sigmazk,Omegazk,'ro','linewidth',1.5); 
plot([-500,500],[0,0],'k','linewidth',0.75);
plot([0,0],[-500,500],'k','linewidth',0.75);
%plot(Omegac*cos(0.5*pi*(1:0.01:3)),Omegac*sin(0.5*pi*(1:0.01:3)),'k:');
axis([-500,500,-500,500]); axis square;
xlabel('Real-axis'); ylabel('Imaginary-axis');
title('Pole Locations');
set(gca,'xtick',[-500,0,500],'ytick',[-500,0,500]);
% Print Plot
%print -depsc2 ../artfiles/1116_ex1105c.eps;

% Elliptic Approximation
disp('** Elliptic **');
[N, Omegac] = ellipord(Omegap, Omegas, Ap, As, 's');
N
[C,D] = ellip(N,Ap,As,Omegac,'s'); H = freqs(C,D,Om);
Hmag = abs(H); Hpha = angle(H); Hdb = 20*log10(Hmag);
Hgdl = -diff(unwrap(Hpha))./diff(Om); Hgdl = [Hgdl,Hgdl(end)];
Hgdl = medfilt1(Hgdl,3);
pk = roots(D); sigmapk = real(pk); Omegapk = imag(pk);
zk = roots(C); sigmazk = real(zk); Omegazk = imag(zk);
%% Design Plots
Hf11_05d = figure('units','inches','position',[7,5.5,5.8,3.6],...
    'paperunits','inches','paperposition',[0,0,5.8,3.6]);
set(Hf11_05d,'number','off','name','Ex11.5: Elliptic');

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

subplot(2,2,4); % Pole-Zero Plot
plot(sigmapk,Omegapk,'bx','linewidth',1.5); hold on;
plot(sigmazk,Omegazk,'ro','linewidth',1.5); 
plot([-500,500],[0,0],'k','linewidth',0.75);
plot([0,0],[-500,500],'k','linewidth',0.75);
%plot(Omegac*cos(0.5*pi*(1:0.01:3)),Omegac*sin(0.5*pi*(1:0.01:3)),'k:');
axis([-500,500,-500,500]); axis square;
xlabel('Real-axis'); ylabel('Imaginary-axis');
title('Pole Locations');
set(gca,'xtick',[-500,0,500],'ytick',[-500,0,500]);
% Print Plot
%print -depsc2 ../artfiles/1116_ex1105d.eps;