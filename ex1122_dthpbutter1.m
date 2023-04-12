% Example 11.22: Digital Highpass Butterworth Filter Design
%                using Matlab's digital design functions
%                wp = 0.7*pi; ws = 0.5*pi; Ap = 1; As = 40;

clc; close all;

% Given Design Parameters
omegap = 0.7; Ap = 1; omegas = 0.5; As = 40;

% Design using Matlab Functions
[N,omegac] = buttord(omegap,omegas,Ap,As)
[B,A] = butter(N,omegac,'high');

% Plotting Parameters and Filter Responses
om = linspace(0,1,501)*pi; H = freqz(B,A,om);
Hmag = abs(H); Hpha = angle(H); Hdb = 20*log10(Hmag);
Hgdl = -diff(unwrap(Hpha))./diff(om); Hgdl = [Hgdl,Hgdl(end)];
Hgdl = medfilt1(Hgdl,3); Hgrd = grpdelay(B,A,om);
N = 50; n = 0:N; x = (n==0); h = filter(B,A,x);

%% Design Plots
Hf11_22 = figure('units','inches','position',[1,1,5.8,3.6],...
    'paperunits','inches','paperposition',[0,0,5.8,3.6]);
set(Hf11_22,'number','off','name','Ex11.2: Digital Butterworth HPF');

subplot(2,2,1); % Magnitude Response
plot(om/pi,Hmag,'b','linewidth',1); axis([0,1,0,1.1]);
xlabel('\omega/\pi'); ylabel('Magnitude');
title('Magnitude Response');
set(gca,'xtick',[0,omegas,omegac,1],...
    'xticklabel','0|.5|.67|1'); 
set(gca,'ytick',[0,1/sqrt(2),1],'yticklabel','0|0.7|1'); 
grid; box off; 

subplot(2,2,2); % Log-Magnitude Response in dB
plot(om/pi,Hdb,'b','linewidth',1); axis([0,1,-80,10]);
xlabel('\omega/\pi'); ylabel('Decibels');
title('Log-Magnitude Response');
set(gca,'xtick',[0,omegas,omegap,1],...
    'xticklabel','0|.5|.7|1');
set(gca,'ytick',[-80,-As,0]); grid; box off; 

subplot(2,2,3); % Group-Delay Response in Samples
plot(om/pi,Hgrd,'b','linewidth',1); axis([0,1,0,12]); 
xlabel('\omega/\pi'); ylabel('Samples');
title('Group Delay Response');
set(gca,'xtick',[0,omegas,omegap,1],...
    'xticklabel','0|.5|.7|1');
set(gca,'ytick',(0:2:12)); grid; box off;

subplot(2,2,4); % Impulse Response Plots
stem(n,h,'filled','markersize',3); axis([0,N,-.4,0.4])
xlabel('Time in seconds'); ylabel('Amplitude');
title('Impulse Response');
set(gca,'xtick',[0:10:50],'ytick',[-.4:.2:.4]); box off;

% Print Plot
%print -depsc2 ../artfiles/1135_ex1122.eps;