% Example 11.5: Analog Lowpass Chebyshev-2 Filter Design
%               using hand calculations
% Omegap = 20; Omegas = 30; Ap = 6; As = 20;

clc; close all; echo on;

% Given Design Parameters
Omegap = 2; Omegas = 3; Ap = 6; As = 20;
%Omegap = 0.2*pi; Omegas = 0.3*pi; Ap = 7; As = 16;
% Analog Design Parameters (Eq. 10.9)
epsilon = sqrt(10^(0.1*Ap)-1); A = 10^(0.05*As);
Rp = 1/sqrt(1+epsilon^2);

%% Design Steps
% Step-1: Compute alpha and beta
alpha = Omegas/Omegap;
beta = (1/epsilon)*sqrt(A^2-1);
% Step-2: Calulation of N
N = log(beta+sqrt(beta^2-1))/log(alpha+sqrt(alpha^2-1));
N = ceil(N);
% Step-3; Calculation of exact Omegas
Ws = Omegap*cosh(acosh(beta)/N); Omegac = Ws; 
% Step-3: Calculation of a and b
%gamma = (1/epsilon+sqrt(1+1/epsilon^2))^(1/N);
gamma = (A+sqrt(A^2-1))^(1/N);
a = 0.5*(gamma-1/gamma);
b = 0.5*(gamma+1/gamma);
%ep = 1/sqrt(A^2-1); v0 = asinh(1/ep)/N;
%a = sinh(v0), b = cosh(v0),
% Step-4: Calculations of Poles
k = 1:N; thetak = pi/2+(2*k-1)*pi/(2*N);
sigmak = (a*Omegac)*cos(thetak); Omegak = (b*Omegac)*sin(thetak);
pk = (sigmak + 1j*Omegak);
pk = (Omegac^2)./pk;
% Step-5: Calculation of Zeros
thetak = (2*k-1)*pi/(2*N);
if mod(N,2) == 1
    thetak = [thetak(1:(N-1)/2),thetak((N+3)/2:N)];
end
zk = 1j*Omegac./cos(thetak);

% Step-5: Calculation of the system function
c = real(poly(zk)); % Direct Form
d = real(poly(pk));
c = [0,c]; c = d(end)*c/c(end);

%echo off;

% Design using SP Toolbox functions
[N, Ws] = cheb2ord(Omegap, Omegas, Ap, As, 's');
[C,D] = cheby2(N,As, Ws,'s');
% zk = roots(c); pk = roots(d);
% D1 = real(poly(pk(1:2)));
% D2 = real(poly(pk(3)));

%% Design of Chebyshev-I lowpass filter for Group-Delay Plot
[NI, WpI] = cheb1ord(Omegap, Omegas, Ap, As, 's');
[cI,dI] = cheby1(NI,Ap, WpI,'s');
Ommax = 5; Om = linspace(0,Ommax,101); HI = freqs(cI,dI,Om);
HIpha = angle(HI); HIgdl = -diff(unwrap(HIpha))./diff(Om); 
HIgdl = [HIgdl,HIgdl(end)];

%% Design Plots
Hf11_05 = figure('units','inches','position',[1,1,5.8,3.6],...
    'paperunits','inches','paperposition',[0,0,5.8,3.6]);
set(Hf11_05,'number','off','name','Ex11.5');

Ommax = 5; Om = linspace(0,Ommax,101); H = freqs(c,d,Om);
Hmag = abs(H); Hpha = angle(H); Hdb = 20*log10(Hmag);
Hgdl = -diff(unwrap(Hpha))./diff(Om); Hgdl = [Hgdl,Hgdl(end)];
Hgdl = medfilt1(Hgdl,3);

subplot(2,2,1);
plot(Om,Hmag,'b','linewidth',1); axis([0,5,0,1.1]);
xlabel('Frequency in rad/sec.'); ylabel('Magnitude');
title('Magnitude Response');
set(gca,'xtick',[0,Omegap,Omegas,Ommax]);
set(gca,'ytick',[0,1/A,Rp,1]); grid;


subplot(2,2,2);
plot(Om,Hdb,'b','linewidth',1); axis([0,5,-40,1]);
xlabel('Frequency in rad/sec.'); ylabel('Decibels');
title('Log-Magnitude Response');
set(gca,'xtick',[0,Omegap,Omegas,Ommax]);
set(gca,'ytick',[-20,-6,0]); grid;


subplot(2,2,3);
plot(Om,Hgdl,'b','linewidth',1); axis([0,5,0,6]); hold on;
plot(Om,HIgdl,'b--','linewidth',1);
%plot(Om,Hpha/pi,'b','linewidth',1); axis([0,5,-1.1,1.1]);
xlabel('Frequency in rad/sec.'); ylabel('Samples');
title('Group Delay Response');
set(gca,'xtick',[0,Omegap,Omegas,Ommax]);
set(gca,'ytick',(0:2:6)); grid;
legend('Cheb-II','Cheb-I','location','best');
%set(gca,'yticklabel','-p|0|p','fontname','symbol');


subplot(2,2,4);
plot(real(zk),imag(zk),'bo','linewidth',1.5); hold on;
plot([-4,4],[0,0],'k','linewidth',0.75);
plot([0,0],[-4,4],'k','linewidth',0.75);
plot(real(pk),imag(pk),'bx','linewidth',1.5);
axis([-4,4,-4,4]); axis square;
xlabel('Real-axis'); ylabel('Imaginary-axis');
title('Pole-Zero Plot');
set(gca,'xtick',(-4:2:4)); set(gca,'ytick',(-4:2:4)); hold off;