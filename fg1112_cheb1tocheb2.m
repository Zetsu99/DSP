% Script for creating Chebyshev-I to Chebyshev-II transition

clc; close all;
Hf = figure('units','inches','position',[1,1,5.8,1.5],...
    'paperunits','inches','paperposition',[0,0,5.8,1.5]);
set(Hf,'number','off','name','Cheb1toCheb2');

% Chebyshev-I Lowpass
Omegap = 2; Omegas = 3; Ap = 1; As = 10;
epsilon = sqrt(10^(0.1*Ap)-1); A = 10^(0.05*As);
Rp = 1/sqrt(1+epsilon^2);
[N, Wp] = cheb1ord(Omegap, Omegas, Ap, As, 's');
[C1,D1] = cheby1(N,Ap, Wp,'s');
Ommax = 5; Om = linspace(0,Ommax,101); H1 = freqs(C1,D1,Om);
Hmag1 = abs(H1); Hmag1sq = Hmag1.*Hmag1;

subplot(1,3,1);
plot(Om,Hmag1sq,'b','linewidth',1); axis([0,Ommax,0,1.1]);
%xlabel('Frequency in rad/sec.'); ylabel('Magnitude');
%title('Magnitude Response');
set(gca,'xtick',[0,Omegap,Ommax]);
set(gca,'ytick',[0,1/A^2,Rp^2,1],'yticklabel',' | | | '); grid;box off

% Chebyshev-I Highpass
Omegap = 2; Omegas = 1.5; Ap = 1; As = 10;
epsilon = sqrt(10^(0.1*Ap)-1); A = 10^(0.05*As);
Rp = 1/sqrt(1+epsilon^2);
[N, Wp] = cheb1ord(Omegap, Omegas, Ap, As, 's');
[C2,D2] = cheby1(N,Ap, Wp,'high','s');
Ommax = 5; Om = linspace(0,Ommax,101); H2 = freqs(C2,D2,Om);
Hmag2 = abs(H2); Hmag2sq = Hmag2.*Hmag2;

subplot(1,3,2);
plot(Om,Hmag2sq,'b','linewidth',1); axis([0,Ommax,0,1.1]);
%xlabel('Frequency in rad/sec.'); ylabel('Magnitude');
%title('Magnitude Response');
set(gca,'xtick',[0,Omegap,Ommax]);
set(gca,'ytick',[0,1/A^2,Rp^2,1],'yticklabel',' | | | '); grid;box off;

% Chebyshev-II Lowpass
Omegap = 2; Omegas = 1.5; Ap = 1; As = 10;
epsilon = sqrt(10^(0.1*Ap)-1); A = 10^(0.05*As);
Rp = 1/sqrt(1+epsilon^2);
[N, Wp] = cheb1ord(Omegap, Omegas, Ap, As, 's');
[C2,D2] = cheby1(N,Ap, Wp,'high','s');
Ommax = 5; Om = linspace(0,Ommax,101); H2 = freqs(C2,D2,Om);
Hmag2 = abs(H2); Hmag2sq = Hmag2.*Hmag2;
Hmag3sq = 1-Hmag2sq;

subplot(1,3,3);
plot(Om,Hmag3sq,'b','linewidth',1); axis([0,Ommax,0,1.1]);
%xlabel('Frequency in rad/sec.'); ylabel('Magnitude');
%title('Magnitude Response');
set(gca,'xtick',[0,Omegap,Ommax]);
set(gca,'ytick',[0,1-Rp^2,1-1/A^2,1],'yticklabel',' | | | '); grid;
box off;

% Print Plot
print -depsc2 ../artfiles/1112_cheb1tocheb2.eps;