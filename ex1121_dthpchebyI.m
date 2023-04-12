% Example 11.21: Digital High Chebyshev I Filter Design
%                using Transformation Function
%                ws = 0.45*pi; wp = 0.6*pi; Ap = 1; As = 15;

clc; close all;

% Given Design Parameters
omegas = 0.45*pi; omegap = 0.6*pi; Ap = 1; As = 15;
epsilon = sqrt(10^(0.1*Ap)-1); A = 10^(0.05*As);
Rp = 1/sqrt(1+epsilon^2);

thetap = 1; % Lowpass prototype passband cutoff
al = -(cos((thetap+omegap)/2))/(cos((thetap-omegap)/2)) % alpha
thetas = atan(-(1-al^2)*sin(omegas)/(2*al+(1+al^2)*cos(omegas)))
% Lowpass Prototype Filter Design
[N,omegac] = cheb1ord(thetap/pi,thetas/pi,Ap,As)
[Bp,Ap] = cheby1(N,Ap,omegac)

% Lowpass to Highpass transformation
[Bh,Ah] = zmapping(Bp,Ap,-[al,1],[1,al]);
Ah0 = Ah(1); Bh = Bh/Ah0, Ah = Ah/Ah0

% Plotting Parameters and Filter Responses
om = linspace(0,1,501)*pi; Hp = freqz(Bp,Ap,om);
Hmagp = abs(Hp); Hphap = angle(Hp); Hdbp = 20*log10(Hmagp);
Hh = freqz(Bh,Ah,om);
Hmagh = abs(Hh); Hphah = angle(Hh); Hdbh = 20*log10(Hmagh);

N = 100; n = 0:N; x = (n==0); h = filter(Bh,Ah,x);

%% Design Plots
Hf11_21 = figure('units','inches','position',[1,1,5.8,3.6],...
    'paperunits','inches','paperposition',[0,0,5.8,3.6]);
set(Hf11_21,'number','off','name','Ex11.21: Digital Chebyshev I');

subplot(2,2,1); % Magnitude Response (Lowpass)
plot(om,Hmagp,'b','linewidth',1); axis([0,pi,0,1.1]);
xlabel('\omega'); ylabel('Magnitude');
title('Lowpass Prototype Filter'); 
set(gca,'xtick',[0,thetap,thetas,pi],...
    'xticklabel','0|1|1.44|p','fontname','symbol'); 
set(gca,'ytick',[0,1/A,Rp,1],'yticklabel','0|0.18|0.9|1'); 
grid; box off; 

% subplot(2,2,3); % Log-Magnitude Response in dB (Lowpass)
% plot(om,Hdbp,'b','linewidth',1); axis([0,pi,-50,5]);
% xlabel('\omega'); ylabel('Decibels');
% title('Lowpass Prototype Filter');
% set(gca,'xtick',[0,thetap,thetas,pi],...
%     'xticklabel','0|1|1.44|p','fontname','symbol');
% set(gca,'ytick',[-50,-30,-As,0]); grid; box off; 

subplot(2,2,2); % Magnitude Response (Highpass)
plot(om/pi,Hmagh,'b','linewidth',1); axis([0,1,0,1.1]);
xlabel('\omega/\pi'); ylabel('Magnitude');
title('Digital Highpass Filter'); 
set(gca,'xtick',[0,omegas/pi,omegap/pi,1],...
    'xticklabel','0|0.45|0.6|1'); 
set(gca,'ytick',[0,1/A,Rp,1],'yticklabel','0|0.18|0.9|1'); 
grid; box off; 

% subplot(2,2,4); % Log-Magnitude Response in dB (Highpass)
% plot(om/pi,Hdbh,'b','linewidth',1); axis([0,1,-50,5]);
% xlabel('\omega/\pi'); ylabel('Decibels');
% title('Digital Highpass Filter');
% set(gca,'xtick',[0,omegas/pi,omegap/pi,1],...
%     'xticklabel','0|0.45|0.6|1'); 
% set(gca,'ytick',[-50,-30,-As,0]); grid; box off; return


% Print Plot
%print -depsc2 ../artfiles/1134_ex1121.eps;