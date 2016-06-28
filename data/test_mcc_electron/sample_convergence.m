close all
clear all
clc

Nsample = 10000; Np = 100000;

fileID = fopen('prob.bin');
prob = fread(fileID,4,'double');
fileID = fopen('coll_sample.bin');
Ncoll = fread(fileID,Nsample*4,'int32');
Ncoll = reshape(Ncoll,[Nsample,4]);

%%
close all

Nexp = Np*prob;
Nstd = sqrt( Np*prob.*(1-prob) );
r = zeros(Nsample,4);
r(:,1) = Ncoll(:,1)-Nexp(1);
r(:,2) = Ncoll(:,2)-Nexp(2);
r(:,3) = Ncoll(:,3)-Nexp(3);
r(:,4) = Ncoll(:,4)-Nexp(4);
m = zeros(4,1); o=m;
for i=1:4
    m(i) = mean(r(:,i));
    o(i) = std(r(:,i));
end

Ng = 10000;
for i=1:4
    figure(i)
    h=histogram(r(:,i));
    hold on
    L = max(abs(r(:,i)));
    x = linspace(-L,L,Ng);
    f = Nsample*2*L/h.NumBins/sqrt(2*pi)/Nstd(i)*exp( -x.^2/2/Nstd(i)/Nstd(i) );
    xm = linspace(-L,L,Ng)+m(i);
    fm = Nsample*2*L/h.NumBins/sqrt(2*pi)/o(i)*exp( -(xm-m(i)).^2/2/o(i)/o(i) );
    plot(xm,fm,'-k');
    plot(x,f,'-r');
    xlabel(strcat('$r_',num2str(i-1),'$'),'interpreter','latex');
    ylabel('population');
    title(strcat(num2str(Nsample),' sample for $N=',num2str(Np),'$, $P_',num2str(i-1),'=',num2str(prob(i)),'$'),'interpreter','latex');
    set(gca,'fontsize',25);
end