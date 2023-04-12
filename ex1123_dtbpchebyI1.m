% Example 11.23: Digital Bandpass Chebyshev I Filter Design
%                using Matlab's digital design functions
%                wp = [0.3,0.5]*pi; ws = [0.25,0.6]*pi; 
%                Ap = 0.5; As = 50;

clc; close all;

% Given Design Parameters
omegap = [0.3,0.5]; omegas = [0.25,0.6]; Ap = 0.5; As = 50;

% Design using Matlab Functions
[N,omegac] = cheb1ord(omegap,omegas,Ap,As)
[B,A] = cheby1(N,Ap,omegac);

% Plotting Parameters and Filter Responses
om = linspace(0,1,501)*pi; H = freqz(B,A,om);
Hmag = abs(H); Hpha = angle(H); Hdb = 20*log10(Hmag);
Hgdl = -diff(unwrap(Hpha))./diff(om); Hgdl = [Hgdl,Hgdl(end)];
Hgdl = medfilt1(Hgdl,3); Hgrd = grpdelay(B,A,om);
N = 100; n = 0:N; x = (n==0); h = filter(B,A,x);

%% Design Plots
Hf11_23 = figure('units','inches','position',[1,1,5.8,3.6],...
    'paperunits','inches','paperposition',[0,0,5.8,3.6]);
set(Hf11_23,'number','off','name','Ex11.23: Digital Chenyshev I BPF');

subplot(2,2,1); % Magnitude Response
plot(om/pi,Hmag,'b','linewidth',1); axis([0,1,0,1.1]);
xlabel('\omega/\pi'); ylabel('Magnitude');
title('Magnitude Response');
set(gca,'xtick',[0,omegas(1),omegap(1),omegap(2),omegas(2),1],...
    'xticklabel','0|.25  | .3|.5|.6|1'); 
set(gca,'ytick',[0,1],'yticklabel','0|1'); 
grid; box off; 

subplot(2,2,2); % Log-Magnitude Response in dB
plot(om/pi,Hdb,'b','linewidth',1); axis([0,1,-80,10]);
xlabel('\omega/\pi'); ylabel('Decibels');
title('Log-Magnitude Response');
set(gca,'xtick',[0,omegas(1),omegap(1),omegap(2),omegas(2),1],...
    'xticklabel','0|.25  | .3|.5|.6|1'); 
set(gca,'ytick',[-80,-As,0]); grid; box off; 

subplot(2,2,3); % Group-Delay Response in Samples
plot(om/pi,Hgrd,'b','linewidth',1); axis([0,1,0,80]); 
xlabel('\omega/\pi'); ylabel('Samples');
title('Group Delay Response');
set(gca,'xtick',[0,omegas(1),omegap(1),omegap(2),omegas(2),1],...
    'xticklabel','0|.25  | .3|.5|.6|1'); 
set(gca,'ytick',(0:20:80)); grid; box off;

subplot(2,2,4); % Impulse Response Plots
stem(n,h,'filled','markersize',3); axis([0,N,-.2,0.2])
xlabel('Time in seconds'); ylabel('Amplitude');
title('Impulse Response');
set(gca,'xtick',[0:20:100],'ytick',[-.2:.1:.2]); box off;

% Print Plot
%print -depsc2 ../artfiles/1136_ex1123.eps;