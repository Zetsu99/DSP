% Example 11.20: Analog Bandpass Elliptic Filter Design
%               using Matlab function lp2bp
% Fp1 = 20 Hz, Fp2 = 60 Hz; Ap = 1; As = 40;

clc; close all; %echo on;

% Given Design Parameters
Fp1 = 20; Fp2 = 60; Ap = 1; As = 40;
Omegap1 = 2*pi*Fp1; Omegap2 = 2*pi*Fp2; %Omegas = 2*pi*Fs;

% Design a Unity Bandwidth Lowpass Elliptic Filter
wp = 1; ws = 1.5; [N,wc] = ellipord(wp,ws,Ap,As,'s'); N
[C,D] = ellip(N,Ap,As,wc,'s');
wmax = 5; w = linspace(0,wmax,501); 
H = freqs(C,D,w); Hmag = abs(H); Hdb = 20*log10(Hmag);

Hf11_20 = figure('units','inches','position',[1,1,5.8,2],...
    'paperunits','inches','paperposition',[0,0,5.8,2]);
set(Hf11_20,'number','off','name','Ex11.20: LP2BP_Elliptic');

% Log-Magnitude Response in dB (Lowpass)
subplot('position',[.08,.2,.4,.7]); 
plot(w,Hdb,'b','linewidth',1); axis([0,wmax,-60,5]);
xlabel('Frequency in rad/s','verticalalignment','middle'); 
ylabel('Decibels');
title('Analog Lowpass Elliptic Prototype');
set(gca,'xtick',[0,wp,ws,wmax]); 
set(gca,'ytick',[-60,-As,0],'yticklabel','60|40|0'); 
grid; box off;

% Lowpass to Bandpass Transformation using lp2bp
Omega0 = sqrt(Omegap1*Omegap2); BW = Omegap2-Omegap1;
[B,A] = lp2bp(C,D,Omega0,BW); Nbp = length(A)-1
Fmax = 100; F = linspace(0,Fmax,501); Om = 2*pi*F;
H = freqs(B,A,Om); Hmag = abs(H); Hdb = 20*log10(Hmag);

% Log-Magnitude Response in dB (Bandpass)
subplot('position',[.58,.2,.4,.7]); 
plot(F,Hdb,'b','linewidth',1); axis([0,Fmax,-60,5]);
xlabel('Frequency in Hz','verticalalignment','middle'); 
ylabel('Decibels');
title('Analog Bandpass Elliptic Design');
set(gca,'xtick',[0,Fp1,Fp2,Fmax]);
set(gca,'ytick',[-60,-As,0],'yticklabel','60|40|0'); 
grid; box off;

% Print Plot
%print -depsc2 ../artfiles/1132_ex1120.eps;