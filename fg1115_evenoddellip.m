% Script for Figure 11.15 
%  Odd and Even Magnitude Responses of Elliptic Filters

clc; close all;
epsilon = 0.75; Ap = 10*log10(1+epsilon^2); Rp = 1/sqrt(1+epsilon^2);
A = 10; As = 20*log10(A);
Omegap = 3; Omegas = 4;

% Plotting Parameters
Ommax = 10; Om = linspace(0,Ommax,501);

% Odd N
[N, Omegac] = ellipord(Omegap, Omegas, Ap, As, 's');
[C,D] = ellip(N,Ap,As,Omegac,'s'); H = freqs(C,D,Om); Hmag = abs(H);
Hf11_15 = figure('units','inches','position',[1,1,4,2],...
    'paperunits','inches','paperposition',[0,0,4,2]);
set(Hf11_15,'number','off','name','Figure 11.15');
subplot('position',[0.15,0.15,0.8,0.7]); % Magnitude Response
plot(Om,Hmag,'b','linewidth',1); axis([0,Ommax,0,1.1]);
%xlabel('Frequency in rad/sec'); 
ylabel('Magnitude'); %title('{\itN} Odd');
set(gca,'xtick',[0,Omegac,Ommax],'ytick',[0,1/A,Rp,1]); 
set(gca,'yticklabel','0|0.1|0.8|1');

% Even N
Omegas = 3.5;
[N, Omegac] = ellipord(Omegap, Omegas, Ap, As, 's');
[C,D] = ellip(N,Ap,As,Omegac,'s'); H = freqs(C,D,Om); Hmag = abs(H);
%subplot('position',[0.62,0.15,0.35,0.7]); % Magnitude Response
hold on;
plot(Om,Hmag,'b--','linewidth',1); 
grid; box off;

print -depsc2 ../artfiles/1115_evenoddellip.eps;