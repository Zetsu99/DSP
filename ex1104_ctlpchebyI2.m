% Example 11.4: Analog Lowpass Chebyshev I Filter Design
%               using Matlab
%               Fp = 40 Hz; Fs = 50 Hz; Ap = 1; As = 30;

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

% Chebyshev-I Approximation
disp('** Chebyshev-I **');
[N, Omegac] = cheb1ord(Omegap, Omegas, Ap, As, 's'); N
Fc = Omegac/(2*pi)
[C,D] = cheby1(N,Ap,Omegac,'s'); H = freqs(C,D,Om);
Hmag = abs(H); Hpha = angle(H); Hdb = 20*log10(Hmag);
Hgdl = -diff(unwrap(Hpha))./diff(Om); Hgdl = [Hgdl,Hgdl(end)];
Hgdl = medfilt1(Hgdl,3);
pk = roots(D); sigmak = real(pk); Omegak = imag(pk);

% Butterworth Approximation for group-delay Comparison
[NB, OmegacB] = buttord(Omegap, Omegas, Ap, As, 's');
[CB,DB] = butter(NB,OmegacB,'s'); HB = freqs(CB,DB,Om);
HphaB = angle(HB); HgdlB = -diff(unwrap(HphaB))./diff(Om); 
HgdlB = [HgdlB,HgdlB(end)]; HgdlB = medfilt1(HgdlB,3);

%% Design Plots
Hf11_04 = figure('units','inches','position',[7,1,5.8,3.6],...
    'paperunits','inches','paperposition',[0,0,5.8,3.6]);
set(Hf11_04,'number','off','name','Ex11.4: Chebyshev-I');

subplot(2,2,1); % Magnitude Response
plot(F,Hmag,'b','linewidth',1); axis([0,Fmax,0,1.1]);
xlabel('Frequency in Hz'); ylabel('Magnitude');
title('Magnitude Response');
set(gca,'xtick',[0,Fp,Fs,Fmax]); 
set(gca,'ytick',[0,Rp,1],'yticklabel','0|0.9|1'); 
grid; box off;

subplot(2,2,2); % Log-Magnitude Response in dB
plot(F,Hdb,'b','linewidth',1); axis([0,Fmax,-80,10]);
xlabel('Frequency in Hz'); ylabel('Decibels');
title('Log-Magnitude Response');
set(gca,'xtick',[0,Fp,Fs,Fmax]);
set(gca,'ytick',[-80,-As,0]); grid; box off;

subplot(2,2,3); % Group-Delay Response in Samples
plot(F,Hgdl,'b','linewidth',1); axis([0,Fmax,0,0.15]); hold on;
plot(F,HgdlB,'b--','linewidth',1);
xlabel('Frequency in Hz'); ylabel('Samples');
title('Group Delay Response');
set(gca,'xtick',[0,Fp,Fs,Fmax]);
set(gca,'ytick',(0:0.05:0.15)); grid; box off;
text(50,0.05,'Butterworth','HorizontalAlignment','left',...
    'VerticalAlignment','bottom');
text(44,0.015,'Chebyshev I','HorizontalAlignment','right'); hold off;

subplot(2,2,4); % Pole-Zero Plot
plot(sigmak,Omegak,'bx','linewidth',1.5); hold on;
plot([-300,300],[0,0],'k','linewidth',0.75);
plot([0,0],[-300,300],'k','linewidth',0.75);
%plot(Omegac*cos(0.5*pi*(1:0.01:3)),Omegac*sin(0.5*pi*(1:0.01:3)),'k:');
axis([-300,300,-300,300]); axis square;
xlabel('Real-axis'); ylabel('Imaginary-axis');
title('Pole Locations');
set(gca,'xtick',[-300,0,300],'ytick',[-300,0,300]);

% Print Plot
%print -depsc2 ../artfiles/1111_ex1104.eps;