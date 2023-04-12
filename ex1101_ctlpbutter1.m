% Example 11.1: Analog Lowpass Butterworth Filter Design
%               using hand calculations
%               Omegap = 2; Omegas = 3; Ap = 6; As = 20;

clc; close all;

% Given Design Parameters
Omegap = 2; Omegas = 3; Ap = 6; As = 20;
%Omegap = 0.2*pi; Omegas = 0.3*pi; Ap = 7; As = 16;
% Analog Design Parameters (Eq. 10.9)
epsilon = sqrt(10^(0.1*Ap)-1), A = 10^(0.05*As)

%% Design Steps
% Step-1: Calulation of N
alpha = Omegas/Omegap, beta = (1/epsilon)*sqrt(A^2-1)
N = log(beta)/log(alpha)
N = ceil(N)
% Step-2: Calculation of Omegac
OmegacL = Omegap/(10^(0.1*Ap)-1)^(1/(2*N))
OmegacH = Omegas/(10^(0.1*As)-1)^(1/(2*N))
Omegac = OmegacH;
% Step-3: Calculations of Poles
k = 1:N; thetak = pi/2+(2*k-1)*pi/(2*N);
sigmak = Omegac*cos(thetak); Omegak = Omegac*sin(thetak);
sk = cplxpair(sigmak + 1j*Omegak); 
% Step-4: Calculation of the system function
C = Omegac^N, D = real(poly(sk)), % Direct Form
D1 = real(poly(sk(1:2)))
D2 = real(poly(sk(3:4)))
D3 = real(poly(sk(5)))

% Design using SP Toolbox functions
[N, Wn] = buttord(Omegap, Omegas, Ap, As, 's');
[C,D] = butter(N,Wn,'s')

%% Design Plots
Hf11_01 = figure('units','inches','position',[1,1,5.8,3.6],...
    'paperunits','inches','paperposition',[0,0,5.8,3.6]);
set(Hf11_01,'number','off','name','Ex11.1');

Ommax = 5; Om = linspace(0,Ommax,101); H = freqs(c,d,Om);
Hmag = abs(H); Hpha = angle(H); Hdb = 20*log10(Hmag);
Hgdl = -diff(unwrap(Hpha))./diff(Om); Hgdl = [Hgdl,Hgdl(end)];

subplot(2,2,1); % Magnitude Plot
plot(Om,Hmag,'b','linewidth',1); axis([0,5,0,1.1]);
xlabel('Frequency in rad/sec.'); ylabel('Magnitude');
title('Magnitude Response');
set(gca,'xtick',[0,Omegac,Omegas,Ommax]);
set(gca,'ytick',[0,0.707,1]); grid;

subplot(2,2,2); % Log-Magnitude Plot in dB
plot(Om,Hdb,'b','linewidth',1); axis([0,5,-30,1]);
xlabel('Frequency in rad/sec.'); ylabel('Decibels');
title('Log-Magnitude Response');
set(gca,'xtick',[0,Omegap,Omegas,Ommax]);
set(gca,'ytick',[-20,-6,0]); grid;

subplot(2,2,3); % Group-Delay Plot in Samples
plot(Om,Hgdl,'b','linewidth',1); axis([0,5,0,3]);
%plot(Om,Hpha/pi,'b','linewidth',1); axis([0,5,-1.1,1.1]);
xlabel('Frequency in rad/sec.'); ylabel('Samples');
title('Group Delay Response');
set(gca,'xtick',[0,Omegap,Omegas,Ommax]);
set(gca,'ytick',[0:3]); grid;
%set(gca,'yticklabel','-p|0|p','fontname','symbol');

subplot(2,2,4); % Pole-zero Plot
plot(sigmak,Omegak,'bx','linewidth',1.5); hold on;
plot([-2.5,2.5],[0,0],'k','linewidth',0.75);
plot([0,0],[-2.5,2.5],'k','linewidth',0.75);
plot(Omegac*cos(0.5*pi*[1:0.01:3]),Omegac*sin(0.5*pi*[1:0.01:3]),'k:');
axis([-2.5,2.5,-2.5,2.5]); axis square;
xlabel('Real-axis'); ylabel('Imaginary-axis');
title('Pole Locations');
set(gca,'xtick',[-2:2:2]); set(gca,'ytick',[-2:2:2]);