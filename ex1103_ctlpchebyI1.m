% Example 11.3: Analog Lowpass Chebyshev-1 Filter Design
%               using hand calculations
%               Omegap = 20; Omegas = 30; Ap = 6; As = 20;

clc; close all;

% Given Design Parameters
Omegap = 2; Omegas = 3; Ap = 6; As = 20;
%Omegap = 0.2*pi; Omegas = 0.3*pi; Ap = 7; As = 16;
% Analog Design Parameters (Eq. 10.9)
epsilon = sqrt(10^(0.1*Ap)-1); A = 10^(0.05*As);
Rp = 1/sqrt(1+epsilon^2);

%% Design Steps
% Step-1: Compute alpha and beta
alpha = Omegas/Omegap; beta = (1/epsilon)*sqrt(A^2-1);
% Step-2: Calulation of N
N = log(beta+sqrt(beta^2-1))/log(alpha+sqrt(alpha^2-1));
N = ceil(N);
% Step-3: Calculation of a and b
Omegac = Omegap;
gamma = (1/epsilon+sqrt(1+1/epsilon^2))^(1/N);
a = 0.5*(gamma-1/gamma); b = 0.5*(gamma+1/gamma);
% Step-4: Calculations of Poles
k = 1:N; thetak = pi/2+(2*k-1)*pi/(2*N);
sigmak = (a*Omegac)*cos(thetak); Omegak = (b*Omegac)*sin(thetak);
sk = cplxpair(sigmak + 1j*Omegak);
% Step-5: Calculation of the system function
D = real(poly(sk)); % Direct Form
if iseven(N)
    G = D(end)*Rp; %1/sqrt(1+epsilon^2);
else
    G = D(end);
end
C = G; 

% Design using SP Toolbox functions
[N, Wp] = cheb1ord(Omegap, Omegas, Ap, As, 's');
[c,d] = cheby1(N,Ap, Wp,'s');
D1 = real(poly(sk(1:2)));
D2 = real(poly(sk(3)));

%% Design of Butterworth lowpass filter for Group-Delay Plot
[NB, WnB] = buttord(Omegap, Omegas, Ap, As, 's');
[CB,DB] = butter(NB,WnB,'s')
Ommax = 5; Om = linspace(0,Ommax,101); HB = freqs(CB,DB,Om);
HBpha = angle(HB); HBgdl = -diff(unwrap(HBpha))./diff(Om); 
HBgdl = [HBgdl,HBgdl(end)];

%% Design Plots
Hf11_03 = figure('units','inches','position',[1,1,5.8,3.6],...
    'paperunits','inches','paperposition',[0,0,5.8,3.6]);
set(Hf11_03,'number','off','name','Ex11.3: Chebyshev-I');

Ommax = 5; Om = linspace(0,Ommax,101); H = freqs(c,d,Om);
Hmag = abs(H); Hpha = angle(H); Hdb = 20*log10(Hmag);
Hgdl = -diff(unwrap(Hpha))./diff(Om); Hgdl = [Hgdl,Hgdl(end)];

subplot(2,2,1);
plot(Om,Hmag,'b','linewidth',1); axis([0,5,0,1.1]);
xlabel('Frequency in rad/sec.'); ylabel('Magnitude');
title('Magnitude Response');
set(gca,'xtick',[0,Omegac,Omegas,Ommax]);
set(gca,'ytick',[0,Rp,1]); grid;


subplot(2,2,2);
plot(Om,Hdb,'b','linewidth',1); axis([0,5,-30,1]);
xlabel('Frequency in rad/sec.'); ylabel('Decibels');
title('Log-Magnitude Response');
set(gca,'xtick',[0,Omegap,Omegas,Ommax]);
set(gca,'ytick',[-20,-6,0]); grid;


subplot(2,2,3);
plot(Om,Hgdl,'b','linewidth',1); axis([0,5,0,6]);hold on;
plot(Om,HBgdl,'b--','linewidth',1);
%plot(Om,Hpha/pi,'b','linewidth',1); axis([0,5,-1.1,1.1]);
xlabel('Frequency in rad/sec.'); ylabel('Samples');
title('Group Delay Response');
set(gca,'xtick',[0,Omegap,Omegas,Ommax]);
set(gca,'ytick',[0:2:6]); grid;
legend('Cheby','Butter','location','best');
%set(gca,'yticklabel','-p|0|p','fontname','symbol');


subplot(2,2,4);
plot(sigmak,Omegak,'bx','linewidth',1.5); hold on;
plot([-2.5,2.5],[0,0],'k','linewidth',0.75);
plot([0,0],[-2.5,2.5],'k','linewidth',0.75);
plot(a*Omegac*cos(0.5*pi*[1:0.01:3]),b*Omegac*sin(0.5*pi*[1:0.01:3]),'k:');
axis([-1.5,1.5,-2.5,2.5]); %axis square;
xlabel('Real-axis'); ylabel('Imaginary-axis');
title('Pole Locations');
set(gca,'xtick',[-1:1:1]); set(gca,'ytick',[-2:2:2]);