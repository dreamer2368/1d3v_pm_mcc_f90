close all
clear all
clc

fileID = fopen('npk.bin');
npk = fread(fileID,4,'int32');
fileID = fopen('dtk.bin');
dtk = fread(fileID,4,'double');

Nncoll = zeros(4,4);
Nprob = zeros(4,4);
for i=1:4
    fileID = fopen(strcat('ncoll',num2str(i),'.bin'));
    Nncoll(i,:) = fread(fileID,4,'int32');
    fileID = fopen(strcat('Nprob_',num2str(i),'.bin'));
    Nprob(i,:) = fread(fileID,4,'double');
end

DTncoll = zeros(4,4);
DTprob = zeros(4,4);
for i=1:4
    fileID = fopen(strcat('ncoll',num2str(i+4),'.bin'));
    DTncoll(i,:) = fread(fileID,4,'int32');
    fileID = fopen(strcat('DTprob_',num2str(i),'.bin'));
    DTprob(i,:) = fread(fileID,4,'double');
end

%%
close all

Nexp = diag(npk)*Nprob;
Nstd = sqrt( diag(npk)*(Nprob.*(1-Nprob)) );
Nerr = abs(Nncoll - Nexp);

figure(1)
plot(Nncoll(1,:));
hold on
errorbar(Nexp(1,:),Nstd(1,:));

figure(2)
% loglog(npk,Nerr(:,1),'.-');
% hold on
loglog(npk,Nerr(:,2),'.-');
hold on
loglog(npk,Nerr(:,3),'.-');
loglog(npk,Nerr(:,4),'.-');
loglog(npk,Nstd(:,1),'--');
loglog(npk,Nstd(:,2),'--');
loglog(npk,Nstd(:,3),'--');
loglog(npk,Nstd(:,4),'--');

%%
close all

DTexp = npk(4)*DTprob;
DTstd = sqrt( npk(4)*(DTprob.*(1-DTprob)) );
DTerr = abs(DTncoll - DTexp);

figure(1)
plot(DTncoll(1,:));
hold on
errorbar(DTexp(1,:),DTstd(1,:));

figure(2)
loglog(dtk,DTerr(:,2),'.-');
hold on
loglog(dtk,DTerr(:,3),'.-');
loglog(dtk,DTerr(:,4),'.-');
loglog(dtk,DTstd(:,1),'--');
loglog(dtk,DTstd(:,2),'--');
loglog(dtk,DTstd(:,3),'--');
loglog(dtk,DTstd(:,4),'--');
loglog(dtk,5e6*dtk.^2,'-r');