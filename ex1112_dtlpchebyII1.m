% Example 11.12: Digital Lowpass Chebyshev II Filter Design
%                using Impulse-Invariance (Matlab)
%                wp = 0.25*pi; ws = 0.35*pi; Ap = 1; As = 30;

clc; close all; %echo on;

% Given Design Parameters
omegap = 0.25*pi; omegas = 0.4*pi; Ap = 1; As = 30;
% Analog Design Parameters (Eq. 10.9)
epsilon = sqrt(10^(0.1*Ap)-1); A = 10^(0.05*As);
Rp = 1/sqrt(1+epsilon^2);

%% Step by Step Impulse Invariance Design
% Step-1: Choose Td
Td = 0.1;
% Step-2: Compute Analog Edge Frequencies
Omegap = omegap/Td; Omegas = omegas/Td;
% Step-3: Design Analog Chebyshev II Approximation
[N,Omegac] = cheb2ord(Omegap,Omegas,Ap,As,'s'); N
[C,D] = cheby2(N,As,Omegac,'s');
% Step-4: Obtain Digital Chebyshev II Filter
[B,A] = impinvar(C,D,1/Td);

% Plotting Parameters and Filter Responses
om = linspace(0,1,501)*pi; H = freqz(B,A,om);
Hmag = abs(H); Hpha = angle(H); Hdb = 20*log10(Hmag);
Hgdl = -diff(unwrap(Hpha))./diff(om); Hgdl = [Hgdl,Hgdl(end)];
Hgdl = medfilt1(Hgdl,3);
N = 50; n = 0:N; x = (n==0); h = filter(B,A,x);
t = linspace(0,N*Td,501); hc = impulse(C,D,t);

%% Design Plots
Hf11_12 = figure('units','inches','position',[1,1,5.8,3.6],...
    'paperunits','inches','paperposition',[0,0,5.8,3.6]);
set(Hf11_12,'number','off','name','Ex11.12: Digital Chebyshev II');

subplot(2,2,1); % Magnitude Response
plot(om/pi,Hmag,'b','linewidth',1); axis([0,1,0,1.2]);
xlabel('\omega/\pi'); ylabel('Magnitude');
title('Magnitude Response');
set(gca,'xtick',[0,omegap,omegas,pi]/pi); 
set(gca,'ytick',[0,Rp,1],'yticklabel','0|0.9|1'); 
grid; box off; 

subplot(2,2,2); % Log-Magnitude Response in dB
plot(om/pi,Hdb,'b','linewidth',1); axis([0,1,-80,10]);
xlabel('\omega/\pi'); ylabel('Decibels');
title('Log-Magnitude Response');
set(gca,'xtick',[0,omegap,omegas,pi]/pi);
set(gca,'ytick',[-80,-As,0]); grid; box off; 

subplot(2,2,3); % Group-Delay Response in Samples
plot(om/pi,Hgdl,'b','linewidth',1); axis([0,1,-2,15]); 
xlabel('\omega/\pi'); ylabel('Samples');
title('Group Delay Response');
set(gca,'xtick',[0,omegap,omegas,pi]/pi);
set(gca,'ytick',(0:5:15)); grid; box off;

subplot(2,2,4); % Impulse Response Plots
stem(n*Td,h/Td,'filled','markersize',3); hold on;
plot(t,hc,'b','linewidth',1); axis([0,N*Td,-1.5,3.5])
xlabel('Time in seconds'); ylabel('Amplitude');
title('Impulse Responses');
set(gca,'xtick',[0:1:5],'ytick',[-1:1:3]); box off;

% Print Plot
print -depsc2 ../artfiles/1122_ex1112.eps;