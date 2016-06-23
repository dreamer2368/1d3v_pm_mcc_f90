close all
clear all
clc

N = 100000; E = 100.0;
fileID = fopen('output.bin');
vp = fread(fileID,N*3,'double');
vp = reshape(vp, [N,3]);

% scatter3(vp(:,1),vp(:,2),vp(:,3));
chi = acos(vp(:,2));

dchi = linspace(0,pi,N);
f = E*sin(dchi)/2/log(1+E)./( 1+E*sin(dchi/2).^2 );
figure(1)
nbins = 60;
histogram(chi,nbins);
hold on
plot(dchi, f*pi/nbins*N, 'linewidth', 4);
xlabel('$\chi(0\sim\pi)$','interpreter','latex');
ylabel('$\sigma(E,\chi)\sin(\chi)$','interpreter','latex');
title('$10^5$ sample of $\sigma(E,\chi)$ at $E=100.0eV$','interpreter','latex');
legend('sample','expected distribution');
set(gca,'fontsize',25);