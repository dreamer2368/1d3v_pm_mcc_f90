clear all
close all
clc

spec = importdata('record');
N = spec(1); Ng = spec(2); Nt = spec(3); L = spec(4); mod = spec(5);
Nt = floor(Nt/mod);

fileID = fopen('Np.bin');
Np = fread(fileID,N*(Nt+1),'int32');
Np = reshape(Np, [N,Nt+1]);

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
    plot(xg,phi(:,i),'-k');
    axis([0 L -140 30]);
%     axis([0 L 0 4]);
    title('potential');
    xlabel('$x$(m)','interpreter','latex');
    ylabel('$\phi$(V)','interpreter','latex');
    set(gca,'fontsize',25,'ticklabelinterpreter','latex');
    
    
%     %videoclip
%     frame = getframe(gcf);
%     writeVideo(writerObj,frame);
%     pause();
    drawnow;
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

n0 = 2e14; wt = n0*L/5E4;

% %video clip
% writerObj1 = VideoWriter('electron.avi');
% writerObj1.FrameRate = 70;
% open(writerObj1);
% 
% %video clip
% writerObj2 = VideoWriter('ion.avi');
% writerObj2.FrameRate = 70;
% open(writerObj2);
% 
% %video clip
% writerObj3 = VideoWriter('phi.avi');
% writerObj3.FrameRate = 70;
% open(writerObj3);
% 
% %video clip
% writerObj4 = VideoWriter('number.avi');
% writerObj4.FrameRate = 70;
% open(writerObj4);

% for i=1:Nt
    i=Nt
    fileID = fopen(strcat('xp/',num2str(i),'_1.bin'));
    xp_e = fread(fileID,Np(1,i+1),'double');
    fileID = fopen(strcat('vp/',num2str(i),'_1.bin'));
    vp_e = fread(fileID,Np(1,i+1)*3,'double');
    vp_e = reshape(vp_e,[Np(1,i+1), 3]);
    
    fileID = fopen(strcat('xp/',num2str(i),'_2.bin'));
    xp_i = fread(fileID,Np(2,i+1),'double');
    fileID = fopen(strcat('vp/',num2str(i),'_2.bin'));
    vp_i = fread(fileID,Np(2,i+1)*3,'double');
    vp_i = reshape(vp_i,[Np(2,i+1), 3]);
    
    f1=figure(1);
    plot(xp_e,vp_e(:,1),'.k');
    axis([0 L -5*ve 5*ve]);
    title('Electron distribution','interpreter','latex');
    xlabel('$x$(m)','interpreter','latex');
    ylabel('$v$(m/s)','interpreter','latex');
    set(gca,'fontsize',25,'ticklabelinterpreter','latex');
    
    f2=figure(2);
    plot([0 L], [-vB -vB],'-b', xp_i,vp_i(:,1),'.r',[0 L], [vB vB], '-b');
    axis([0 L -35*vi 35*vi]);
    title('Ion distribution','interpreter','latex');
    xlabel('$x$(m)','interpreter','latex');
    ylabel('$v$(m/s)','interpreter','latex');
    h=legend('$v_B$');
    set(h,'interpreter','latex','location','southwest');
    set(gca,'fontsize',25,'ticklabelinterpreter','latex');
    
    f3=figure(3);
    plot(xg,phi(:,i),'-k');
    axis([0 L -140 30]);
    title('Potential $\phi(x)$','interpreter','latex');
    xlabel('$x$(m)','interpreter','latex');
    ylabel('$\phi$(V)','interpreter','latex');
    set(gca,'fontsize',25,'ticklabelinterpreter','latex');
    
    f4=figure(4);
    n_e = histcounts(xp_e,Ng); n_e = n_e*wt;
    n_i = histcounts(xp_i,Ng); n_i = n_i*wt;
    semilogy(xg,n_e,'-k',xg,n_i,'-r');
    axis([0 L 1e9 1e12]);
    title('Number density','interpreter','latex');
    xlabel('$x$(m)','interpreter','latex');
    ylabel('$n$(m$^{-1}$)','interpreter','latex');
    h=legend('Electron','Ion');
    set(h,'interpreter','latex','location','southwest');
    set(gca,'fontsize',25,'ticklabelinterpreter','latex');
    
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
%     
%     %videoclip
%     frame = getframe(f4);
%     writeVideo(writerObj4,frame);

    fclose('all');
%     pause();
    drawnow;
% end

% % videoclip close
% close(writerObj1);
% close(writerObj2);
% close(writerObj3);
% close(writerObj4);

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
plot(1:Nt,phi(20,:),'-k');
axis([0 Nt 0 25]);

%%
close all

figure(1)
sample = ( abs(xp_i-0.5*L)<0.001 );
histogram( xp_i );