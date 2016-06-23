close all
clear all
clc

N = 10000; E = 100.0;
fileID = fopen('output_Ar.bin');
vp = fread(fileID,N*3,'double');
vp = reshape(vp, [N,3]);

vel = sqrt( sum( vp.^2, 2 ) );
chi = acos(vp(:,2)./vel);

figure(1)
scatter3(vp(:,1),vp(:,2),vp(:,3),'.');

figure(2)
nbins = 60;
histogram(chi, nbins);
hold on
dchi = linspace(0,pi/2,N);
plot(dchi,sin(2*dchi)*pi/2/nbins*N,'linewidth',2);

figure(3)
plot(chi,vel,'.',dchi,cos(dchi),'-r');