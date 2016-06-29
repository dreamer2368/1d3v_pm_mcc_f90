clear all
close all
clc

spec = importdata('record');
N = spec(1); Ng = spec(2); Nt = spec(3); L = spec(4); mod = spec(5);
Nt = floor(Nt/mod);

fileID = fopen('Np.bin');
Np = fread(fileID,N*Nt,'int32');
Np = reshape(Np, [N,Nt]);

% fileID = fopen('xp_1.bin');
% xp1 = fread(fileID,sum(Np(1,:)),'double');
% 
% fileID = fopen('vp_1.bin');
% vp1 = fread(fileID,sum(Np(1,:)),'double');
% 
% fileID = fopen('xp_2.bin');
% xp2 = fread(fileID,sum(Np(2,:)),'double');
% 
% fileID = fopen('vp_2.bin');
% vp2 = fread(fileID,sum(Np(2,:)),'double');

fileID = fopen('E.bin');
E = fread(fileID,Ng*Nt,'double');
E = reshape(E,[Ng,Nt]);

fileID = fopen('PE.bin');
PE = fread(fileID,Nt,'double');

fileID = fopen('phi.bin');
phi = fread(fileID,Ng*Nt,'double');
phi = reshape(phi,[Ng,Nt]);

fileID = fopen('rho.bin');
rho = fread(fileID,Ng*Nt,'double');
rho = reshape(rho,[Ng,Nt]);


EV_TO_K = 11604.52;
mTorr_to_Pa = 0.13332237;
Kb = 1.38065E-23;
qe = 1.602e-19; me = 9.10938215E-31;
PN = 50.0; TN = 0.026;
T0 = 1.0; Np0 = 1e6; L = 0.02;
ionengy = 15.76;

gden = PN*mTorr_to_Pa/(qe*TN);
f0 = 2*exp(-ionengy/T0);
n0 = gden*f0;
spwt = n0*L/Np0;

%%
close all
clc
npsum1 = 1; npsum2 = 1;
dx = L/(Ng-1);
xg = dx*(0:Ng-1);

% %video clip
% writerObj = VideoWriter('phi.avi');
% writerObj.FrameRate = 20;
% open(writerObj);

for i=1:Nt
%     figure(1)
%     plot(xp1(npsum1:npsum1-1+Np(1,i)),vp1(npsum1:npsum1-1+Np(1,i)),'.k',xp2(npsum2:npsum2-1+Np(2,i)),vp2(npsum2:npsum2-1+Np(2,i)),'.r');
%     npsum1 = npsum1 + Np(1,i);
%     npsum2 = npsum2 + Np(2,i);
    
    figure(2)
    plot(xg,phi(:,i) - phi(floor(Ng/2),i),'-k');
%     axis([0 L -1e0 1]);
%     axis([0 L 0 4]);
    title('potential');
    xlabel('$x$(m)','interpreter','latex');
    ylabel('$\phi$(V)','interpreter','latex');
    set(gca,'fontsize',25);
    
    
%     %videoclip
%     frame = getframe(gcf);
%     writeVideo(writerObj,frame);
    pause(.0001);
end

% % videoclip close
% close(writerObj);

%%
close all

figure(1)
plot(phi(200,:));

%%
close all
clc
dx = L/(Ng-1);
xg = dx*(0:Ng-1);

EV_TO_K = 11604.52;
Te = 50.0*EV_TO_K;
tau = 100.0;
Ti = Te/tau;
me = 9.10938215E-31;
mu = 1836;
mi = mu*me;
K = 1.38065E-23;
vB = sqrt(K*(Te+3*Ti)/mi);

ve = sqrt(K*Te/me); vi = sqrt(K*Ti/mi);

% %video clip
% writerObj1 = VideoWriter('electron.avi');
% writerObj1.FrameRate = 60;
% open(writerObj1);
% 
% %video clip
% writerObj2 = VideoWriter('ion.avi');
% writerObj2.FrameRate = 60;
% open(writerObj2);
% 
% %video clip
% writerObj3 = VideoWriter('phi.avi');
% writerObj3.FrameRate = 60;
% open(writerObj3);

for i=1:Nt
    fileID = fopen(strcat('xp/',num2str(i),'_1.bin'));
    xp_e = fread(fileID,Np(1,i),'double');
    fileID = fopen(strcat('vp/',num2str(i),'_1.bin'));
    vp_e = fread(fileID,Np(1,i)*3,'double');
    vp_e = reshape(vp_e,[Np(1,i), 3]);
    
    fileID = fopen(strcat('xp/',num2str(i),'_2.bin'));
    xp_i = fread(fileID,Np(2,i),'double');
    fileID = fopen(strcat('vp/',num2str(i),'_2.bin'));
    vp_i = fread(fileID,Np(2,i)*3,'double');
    vp_i = reshape(vp_i,[Np(2,i), 3]);
    
    f1=figure(1);
    plot(xp_e,vp_e(:,1),'.k');
    axis([0 L -3*ve 3*ve]);
    title('Electron distribution');
    xlabel('$x$(m)','interpreter','latex');
    ylabel('$v$(m/s)','interpreter','latex');
    set(gca,'fontsize',25);
    
    f2=figure(2);
    plot(xp_i,vp_i(:,1),'.r',[0 L], [vB vB], '-b', [0 L], [-vB -vB],'-b');
    axis([0 L -5*vi 5*vi]);
    title('Ion distribution');
    xlabel('$x$(m)','interpreter','latex');
    ylabel('$v$(m/s)','interpreter','latex');
    set(gca,'fontsize',25);
    
    f3=figure(3);
    plot(xg,phi(:,i),'-k');
    axis([0 L -400 400]);
    title('potential');
    xlabel('$x$(m)','interpreter','latex');
    ylabel('$\phi$(V)','interpreter','latex');
    set(gca,'fontsize',25);
    
    EEDF = 0.5*me*sum( vp_e.^2, 2 )/qe;
    [n,x] = hist(EEDF,50);
    figure(4)
    semilogy(x,n/Np(1,Nt)/diff(x(1:2)));
    
%     %videoclip
%     frame = getframe(f1);
%     writeVideo(writerObj1,frame);
%     
%     %videoclip
%     frame = getframe(f2);
%     writeVideo(writerObj2,frame);
% 
%     %videoclip
%     frame = getframe(f3);
%     writeVideo(writerObj3,frame);

    fclose('all');
    pause(.000001);
end

% % videoclip close
% close(writerObj1);
% close(writerObj2);
% close(writerObj3);

%%
close all

figure(1)
plot(1:Nt,Np(1,:),'-k',1:Nt,Np(2,:),'-r','linewidth',2);
xlabel('Time step');
ylabel('Particles');
title('Number of ion/electron');
legend('Electron','Ion');
set(gca,'fontsize',25);
figure(2)
plot(1:Nt,abs(Np(1,:)-Np(2,:)),'-k');
figure(3)
plot(1:Nt,phi(Ng,:),'-k');
% axis([0 Nt 0 25]);

%%
close all

figure(1)
sample = ( abs(xp_i-0.5*L)<0.001 );
histogram( xp_i );

%%
close all

t = 8e-11*20*(1:Nt);
figure(1)
plot(t,phi(300,:),t,min(phi(300,:))*sin(2*pi*13.56e6*t));

%%
close all

i=312;
fileID = fopen(strcat('xp/',num2str(i),'_1.bin'));
xp_e = fread(fileID,Np(1,i),'double');
fileID = fopen(strcat('vp/',num2str(i),'_1.bin'));
vp_e = fread(fileID,Np(1,i)*3,'double');
vp_e = reshape(vp_e,[Np(1,i), 3]);

figure(1);
plot(xp_e,vp_e(:,1),'.k');
axis([0 L -3*ve 3*ve]);
title('Electron distribution');
xlabel('$x$(m)','interpreter','latex');
ylabel('$v$(m/s)','interpreter','latex');
set(gca,'fontsize',25);

vt = mean(vp_e,1); vt_e = zeros(Np(1,i),3);
vt_e(:,1) = vp_e(:,1) - vt(1);
vt_e(:,2) = vp_e(:,2) - vt(2);
vt_e(:,3) = vp_e(:,3) - vt(3);

figure(2)
histogram(vt_e(:,1));

EEDF = 0.5*me*sum( vt_e.^2, 2 )/qe;
[n,x] = hist(EEDF,500);

b1 = (gamma(2.5))^1.5*(gamma(1.5))^(-2.5);
b2 = gamma(2.5)*(gamma(1.5))^(-1);
T1 = 4.0; T2 = 8.5;
EEDF1 = T1^(-1.5)*b1*exp( -x*b2/T1 );
EEDF2 = T2^(-1.5)*b1*exp( -x*b2/T2 );

figure(3)
semilogy(x,n/Np(1,i)/diff(x(1:2)),x,EEDF1,x,2*EEDF2);
axis([0 15 1e-2 3e-1]);