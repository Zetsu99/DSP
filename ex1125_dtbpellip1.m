% Example 11.25: Digital Bandpass Elliptic Filter Design
%                using Matlab's digital design functions
%                wp = [0.3,0.6]*pi; ws = [0.2,0.7]*pi; 
%                Ap = 0.1; As = 60;

clc; close all;

% Given Design Parameters
omegap = [0.3,0.6]; omegas = [0.2,0.7]; Ap = 0.1; As = 60;

% Design using Matlab Functions
[N,omegac] = ellipord(omegap,omegas,Ap,As)
[B,A] = ellip(N,Ap,As,omegac);

% Plotting Parameters and Filter Responses
om = linspace(0,1,501)*pi; H = freqz(B,A,om);
Hmag = abs(H); Hpha = angle(H); Hdb = 20*log10(Hmag);
Hgdl = -diff(unwrap(Hpha))./diff(om); Hgdl = [Hgdl,Hgdl(end)];
Hgdl = medfilt1(Hgdl,3); Hgrd = grpdelay(B,A,om);
N = 100; n = 0:N; x = (n==0); h = filter(B,A,x);

%% Design Plots
Hf11_25 = figure('units','inches','position',[1,1,5.8,3.6],...
    'paperunits','inches','paperposition',[0,0,5.8,3.6]);
set(Hf11_25,'number','off','name','Ex11.25: Digital Elliptic BPF');

subplot(2,2,1); % Magnitude Response
plot(om/pi,Hmag,'b','linewidth',1); axis([0,1,0,1.1]);
xlabel('\omega/\pi'); ylabel('Magnitude');
title('Magnitude Response');
set(gca,'xtick',[0,omegas(1),omegap(1),omegap(2),omegas(2),1],...
    'xticklabel','0|.2|.3|.6|.7|1'); 
set(gca,'ytick',[0,1],'yticklabel','0|1'); 
grid; box off; 

subplot(2,2,2); % Log-Magnitude Response in dB
plot(om/pi,Hdb,'b','linewidth',1); axis([0,1,-80,10]);
xlabel('\omega/\pi'); ylabel('Decibels');
title('Log-Magnitude Response');
set(gca,'xtick',[0,omegas(1),omegap(1),omegap(2),omegas(2),1],...
    'xticklabel','0|.2|.3|.6|.7|1'); 
set(gca,'ytick',[-80,-As,0]); grid; box off; 

subplot(2,2,3); % Group-Delay Response in Samples
plot(om/pi,Hgrd,'b','linewidth',1); axis([0,1,0,50]); 
xlabel('\omega/\pi'); ylabel('Samples');
title('Group Delay Response');
set(gca,'xtick',[0,omegas(1),omegap(1),omegap(2),omegas(2),1],...
    'xticklabel','0|.2|.3|.6|.7|1'); 
set(gca,'ytick',(0:10:50)); grid; box off;

subplot(2,2,4); % Impulse Response Plots
stem(n,h,'filled','markersize',3); axis([0,N,-.4,0.4])
xlabel('Time in seconds'); ylabel('Amplitude');
title('Impulse Response');
set(gca,'xtick',[0:20:100],'ytick',[-.4:.2:.4]); box off;

% Print Plot
%print -depsc2 ../artfiles/1138_ex1125.eps;