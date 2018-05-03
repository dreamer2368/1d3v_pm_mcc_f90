
%% PIC - base parameters at N=1e1

spec = importdata('../data/N1E2/Debye/record');
N = spec(1); Ng = spec(2); Nt = spec(3); L = spec(4); mod = spec(5); Nt = floor(Nt/mod);

fileID = fopen('../data/N1E2/Debye/Np.bin');
Np = fread(fileID,N*(Nt+1),'int32');
Np = reshape(Np, [N,Nt+1]);

fileID = fopen('../data/N1E2/Debye/vp/0_1.bin');
vp0 = fread(fileID,Np(1)*3,'double');
vp0 = reshape(vp0,[Np(1),3]);

fileID = fopen('../data/N1E2/Debye/xp/0_1.bin');
xp0 = fread(fileID,Np(1),'double');

fileID = fopen('../data/N1E2/Debye/PE.bin');
PE = fread(fileID,Nt,'double');
fileID = fopen('../data/N1E2/Debye/KE_1.bin');
KE = fread(fileID,Nt,'double');

dx = L/Ng; P = dx*ones(Ng,1);
xg = dx*(1:Ng) - 0.5*dx;
t = .5*(1:Nt);

w = 0.1*L; Q = 2.0;
rho_back = -(-1) + Q*( -1/L + 1/sqrt(2*pi)/w*exp( -(xg-0.5*L).^2/2/w/w ) );

%% PIC - field variables

dir0 = {'../data/N1E2/Debye/','../data/N1E3/Debye/',	...
		'../data/N1E4/Debye/','../data/N1E5/Debye/','../data/N1E6/Debye/'};
dir1 = {'../data/N1E2/Debye_perturbed/','../data/N1E3/Debye_perturbed/',	...
		'../data/N1E4/Debye_perturbed/','../data/N1E5/Debye_perturbed/','../data/N1E6/Debye_perturbed/'};
%dir0 = {'../data/Ng64ppn160/Debye/','../data/Ng256ppn160/Debye/',	...
%		'../data/Ng1024ppn160/Debye/','../data/Ng4096ppn160/Debye/'};
%dir1 = {'../data/Ng64ppn160/Debye_perturbed/','../data/Ng256ppn160/Debye_perturbed/',	...
%		'../data/Ng1024ppn160/Debye_perturbed/','../data/Ng4096ppn160/Debye_perturbed/'};
Ng = 64;
Np = 10.^(2:6);

%% wasserstein_rho

dW_rho1 = zeros(Nt,4); dW_pp = zeros(Nt,1);
j=1;
%Nw = Np(j);
%Ridx = zeros(Nw,4);
%for j=1:4
%    Ridx(:,j) = randi(Np(j),[Nw,1]);
%end
parfor (i=1:100,20)
%     for j=1:4
        fileID = fopen(strcat(dir0{j},'xp/',num2str(i),'_1.bin'));
        xp0 = fread(fileID,Np(j),'double');
        fileID = fopen(strcat(dir1{j},'xp/',num2str(i),'_1.bin'));
        xp1 = fread(fileID,Np(j),'double');

%     rho_t = assign(xp0,1e5);
    
%         dW_rho1(i,j) = dBL(-rho0{j}(:,i)/L,-rho1{j}(:,i)/L,Ng(j),dx(j));
		temp = dW_pp(i);
        temp = dBL_ptc(xp0,xp1,Np(j),L);
        dW_pp(i) = temp;
%     dW_rho2(i) = dBL(-rho0(:,i)/L,-rho2(:,i)/L,Ng,dx);
%     dW_pm(i) = dBL_pm(xp0(Ridx),-rho0(:,i+1)/L,Nw,Ng,L);

        fclose('all');
%     end
    i
end
dW_pp(1:100)

save N1E6/1.dat dW_pp -ASCII;

exit;
