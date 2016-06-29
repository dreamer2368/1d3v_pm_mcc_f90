close all
clear all
clc

q_e = 1.602E-19;
m_e = 9.10938356e-31;
m_Ar = 6.6335209e-26;

Nsample = 10000;
fileID = fopen('coll_sample.bin');
ncoll = fread(fileID,Nsample*4,'int32');
ncoll = reshape(ncoll,[Nsample,4]);
fileID = fopen('prob.bin');
prob = fread(fileID,4,'double');

fileID = fopen('before/np.bin');
np = fread(fileID,1,'int32');

fileID = fopen('before/xp_e.bin');
xp_e1 = fread(fileID, np, 'double');
fileID = fopen('before/vp_e.bin');
vp_e1 = fread(fileID, np*3, 'double');
vp_e1 = reshape(vp_e1,[np,3]);
fileID = fopen('before/xp_Ar.bin');
xp_Ar1 = fread(fileID, np, 'double');
fileID = fopen('before/vp_Ar.bin');
vp_Ar1 = fread(fileID, np*3, 'double');
vp_Ar1 = reshape(vp_Ar1,[np,3]);

E_e1 = 0.5*m_e*sum( vp_e1.^2, 2 )/q_e;
E_Ar1 = 0.5*m_Ar*sum( vp_Ar1.^2, 2 )/q_e;

fileID = fopen('after/np_e.bin');
np_e = fread(fileID,1,'int32');
fileID = fopen('after/np_Ar.bin');
np_Ar = fread(fileID,1,'int32');

fileID = fopen('after/xp_e.bin');
xp_e2 = fread(fileID, np_e, 'double');
fileID = fopen('after/vp_e.bin');
vp_e2 = fread(fileID, np_e*3, 'double');
vp_e2 = reshape(vp_e2,[np_e,3]);
fileID = fopen('after/xp_Ar.bin');
xp_Ar2 = fread(fileID, np_Ar, 'double');
fileID = fopen('after/vp_Ar.bin');
vp_Ar2 = fread(fileID, np_Ar*3, 'double');
vp_Ar2 = reshape(vp_Ar2,[np_Ar,3]);

E_e2 = 0.5*m_e*sum( vp_e2.^2, 2 )/q_e;
E_Ar2 = 0.5*m_Ar*sum( vp_Ar2.^2, 2 )/q_e;

exp = np*prob;
std = sqrt(np*prob.*(1-prob/prob(1)));

%%
close all
figure(1)
semilogy(xp_e1,E_e1,'.k',xp_Ar1,E_Ar1,'.r');

figure(2)
semilogy(xp_e2,E_e2,'.k',xp_Ar2,E_Ar2,'.r');

figure(3)
plot(ncoll(218,2:4));
hold on
errorbar(exp(2:4),std(2:4));
title('collision distribution for $N=10^5$, $E=20eV$','interpreter','latex');
xlabel('collision type');
ylabel('number of collisions');
set(gca,'fontsize',25);