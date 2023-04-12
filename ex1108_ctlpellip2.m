% Example 11.8: Analog Lowpass Elliptic Filter Design
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


% % Cheby-I Approximation for group-delay Comparison
% [NII, OmegacII] = cheb1ord(Omegap, Omegas, Ap, As, 's');
% [CII,DII] = cheby1(NII,Ap,OmegacII,'s'); HII = freqs(CII,DII,Om);
% HphaII = angle(HII); HgdlII = -diff(unwrap(HphaII))./diff(Om); 
% HgdlII = [HgdlII,HgdlII(end)]; HgdlII = medfilt1(HgdlII,3);

% Cheby-II Approximation for group-delay Comparison
[NII, OmegacII] = cheb2ord(Omegap, Omegas, Ap, As, 's');
[CII,DII] = cheby2(NII,As,OmegacII,'s'); HII = freqs(CII,DII,Om);
HphaII = angle(HII); HgdlII = -diff(unwrap(HphaII))./diff(Om); 
HgdlII = [HgdlII,HgdlII(end)]; HgdlII = medfilt1(HgdlII,3);

%% Design Plots
Hf11_08 = figure('units','inches','position',[1,1,5.8,3.6],...
    'paperunits','inches','paperposition',[0,0,5.8,3.6]);
set(Hf11_08,'number','off','name','Ex11.8: Elliptic');

subplot(2,2,1); % Magnitude Response
plot(F,Hmag,'b','linewidth',1); axis([0,Fmax,0,1.1]);
xlabel('Frequency in Hz'); ylabel('Magnitude');
title('Magnitude Response');
set(gca,'xtick',[0,Fp,Fs,Fmax]); 
set(gca,'ytick',[0,Rp,1],'yticklabel','0|0.9|1'); grid; box off;

subplot(2,2,2); % Log-Magnitude Response in dB
plot(F,Hdb,'b','linewidth',1); axis([0,Fmax,-80,10]);
xlabel('Frequency in Hz'); ylabel('Decibels');
title('Log-Magnitude Response');
set(gca,'xtick',[0,Fp,Fs,Fmax]);
set(gca,'ytick',[-80,-As,0]); grid; box off;

subplot(2,2,3); % Group-Delay Response in Samples
plot(F,Hgdl,'b','linewidth',1); axis([0,Fmax,0,0.15]); hold on;
plot(F,HgdlII,'b--','linewidth',1);
xlabel('Frequency in Hz'); ylabel('Samples');
title('Group Delay Response');
set(gca,'xtick',[0,Fp,Fs,Fmax]);
set(gca,'ytick',(0:0.05:0.15)); grid; box off; 
text(42,0.12,'Elliptic','HorizontalAlignment','left',...
    'VerticalAlignment','bottom');
text(50,0.025,'Chebyshev II','HorizontalAlignment','left',...
    'VerticalAlignment','bottom'); hold off;

subplot(2,2,4); % Pole-Zero Plot
plot(sigmapk,Omegapk,'bx','linewidth',1.5); hold on;
plot(sigmazk,Omegazk,'ro','linewidth',1.5); 
plot([-400,400],[0,0],'k','linewidth',0.75);
plot([0,0],[-400,400],'k','linewidth',0.75);
%plot(Omegac*cos(0.5*pi*(1:0.01:3)),Omegac*sin(0.5*pi*(1:0.01:3)),'k:');
axis([-400,400,-400,400]); axis square;
xlabel('Real-axis'); ylabel('Imaginary-axis');
title('Pole Locations');
set(gca,'xtick',[-400,0,400],'ytick',[-400,0,400]);

% Print Plot
print -depsc2 ../artfiles/1116_ex1108.eps;