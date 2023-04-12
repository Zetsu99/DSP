% Example 11.7: Analog Lowpass Elliptic Filter Design
%               using hand calculations
%               Omegap = 20; Omegas = 30; Ap = 6; As = 20;

clc; close all; %echo on;

% Given Design Parameters
Omegap = 2; Omegas = 3; Ap = 6; As = 20;
% Analog Design Parameters (Eq. 10.9)
epsilon = sqrt(10^(0.1*Ap)-1); A = 10^(0.05*As);
Rp = 1/sqrt(1+epsilon^2);

% Design using SP Toolbox functions
[N, Omegap] = ellipord(Omegap, Omegas, Ap, As, 's')
[C,D] = ellip(N,Ap,As,Omegap,'s');
zk = roots(C); pk = roots(D);
% D1 = real(poly(pk(1:2)));
% D2 = real(poly(pk(3)));

%% Design of Chebyshev-II lowpass filter for Group-Delay Plot
[NI, WsI] = cheb2ord(Omegap, Omegas, Ap, As, 's');
[cI,dI] = cheby2(NI,As, WsI,'s');
Ommax = 5; Om = linspace(0,Ommax,101); HI = freqs(cI,dI,Om);
HIpha = angle(HI); HIgdl = -diff(unwrap(HIpha))./diff(Om); 
HIgdl = [HIgdl,HIgdl(end)]; HIgdl = medfilt1(HIgdl,3);

%% Design Plots
Hf11_07 = figure('units','inches','position',[1,1,5.8,3.6],...
    'paperunits','inches','paperposition',[0,0,5.8,3.6]);
set(Hf11_07,'number','off','name','Ex11.7');

Ommax = 5; Om = linspace(0,Ommax,101); H = freqs(C,D,Om);
Hmag = abs(H); Hpha = angle(H); Hdb = 20*log10(Hmag);
Hgdl = -diff(unwrap(Hpha))./diff(Om); Hgdl = [Hgdl,Hgdl(end)];
Hgdl = medfilt1(Hgdl,3);

subplot(2,2,1); % Magnitude Plot
plot(Om,Hmag,'b','linewidth',1); axis([0,5,0,1.1]);
xlabel('Frequency in rad/sec.'); ylabel('Magnitude');
title('Magnitude Response');
set(gca,'xtick',[0,Omegap,Omegas,Ommax]);
set(gca,'ytick',[0,1/A,Rp,1]); grid;

subplot(2,2,2); % Log-Magnitude Plot in dB
plot(Om,Hdb,'b','linewidth',1); axis([0,5,-40,1]);
xlabel('Frequency in rad/sec.'); ylabel('Decibels');
title('Log-Magnitude Response');
set(gca,'xtick',[0,Omegap,Omegas,Ommax]);
set(gca,'ytick',[-20,-6,0]); grid;

subplot(2,2,3); % Group-Delay Plot in samples
plot(Om,Hgdl,'b','linewidth',1); axis([0,5,0,4]); hold on;
plot(Om,HIgdl,'b--','linewidth',1);
%plot(Om,Hpha/pi,'b','linewidth',1); axis([0,5,-1.1,1.1]);
xlabel('Frequency in rad/sec.'); ylabel('Samples');
title('Group Delay Response');
set(gca,'xtick',[0,Omegap,Omegas,Ommax]);
set(gca,'ytick',(0:2:6)); grid;
legend('Ellip','ChbII','location','best');
%set(gca,'yticklabel','-p|0|p','fontname','symbol');

subplot(2,2,4); % Ploe-Zero Plot
plot(real(zk),imag(zk),'bo','linewidth',1.5); hold on;
plot([-4,4],[0,0],'k','linewidth',0.75);
plot([0,0],[-4,4],'k','linewidth',0.75);
plot(real(pk),imag(pk),'bx','linewidth',1.5);
axis([-4,4,-4,4]); axis square;
xlabel('Real-axis'); ylabel('Imaginary-axis');
title('Pole-Zero Plot');
set(gca,'xtick',(-4:2:4)); set(gca,'ytick',(-4:2:4)); hold off;