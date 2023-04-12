% Example 11.4: Analog Lowpass Chebyshev II Filter Design
%               using Matlab
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

% Chebyshev-II Approximation
disp('** Chebyshev-II **');
[N, Omegac] = cheb2ord(Omegap, Omegas, Ap, As, 's'); N
Fc = Omegac/(2*pi)
[C,D] = cheby2(N,As,Omegac,'s'); H = freqs(C,D,Om);
Hmag = abs(H); Hpha = angle(H); Hdb = 20*log10(Hmag);
Hgdl = -diff(unwrap(Hpha))./diff(Om); Hgdl = [Hgdl,Hgdl(end)];
Hgdl = medfilt1(Hgdl,3);
pk = roots(D); sigmapk = real(pk); Omegapk = imag(pk);
zk = roots(C); sigmazk = real(zk); Omegazk = imag(zk);

% Cheby-I Approximation for group-delay Comparison
[NI, OmegacI] = cheb1ord(Omegap, Omegas, Ap, As, 's');
[CI,DI] = cheby1(NI,Ap,OmegacI,'s'); HI = freqs(CI,DI,Om);
HphaI = angle(HI); HgdlI = -diff(unwrap(HphaI))./diff(Om); 
HgdlI = [HgdlI,HgdlI(end)]; HgdlI = medfilt1(HgdlI,3);

%% Design Plots
Hf11_06 = figure('units','inches','position',[7,1,5.8,3.6],...
    'paperunits','inches','paperposition',[0,0,5.8,3.6]);
set(Hf11_06,'number','off','name','Ex11.6: Chebyshev-I');

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
plot(F,HgdlI,'b--','linewidth',1);
xlabel('Frequency in Hz'); ylabel('Samples');
title('Group Delay Response');
set(gca,'xtick',[0,Fp,Fs,Fmax]);
set(gca,'ytick',(0:0.05:0.15)); grid; box off;
text(40,0.1,'Chebyshev I','HorizontalAlignment','left',...
    'VerticalAlignment','bottom');
text(50,0.025,'Chebyshev II','HorizontalAlignment','left',...
    'VerticalAlignment','bottom'); hold off;

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
%print -depsc2 ../artfiles/1113_ex1106.eps;