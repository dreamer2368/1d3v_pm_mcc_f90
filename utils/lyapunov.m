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

dir0 = {'../data/N1E2/Debye/','../data/N1E3/Debye/',  ...
        '../data/N1E4/Debye/','../data/N1E5/Debye/','../data/N1E6/Debye/'};
dir1 = {'../data/N1E2/Debye_perturbed/','../data/N1E3/Debye_perturbed/',    ...
        '../data/N1E4/Debye_perturbed/','../data/N1E5/Debye_perturbed/','../data/N1E6/Debye_perturbed/'};

for i=1:5
    fileID = fopen(strcat(dir0{i},'E.bin'));
    E0{i} = fread(fileID,Ng*(Nt+1),'double');
    E0{i} = reshape(E0{i},[Ng,Nt+1]);

    fileID = fopen(strcat(dir0{i},'rho.bin'));
    rho0{i} = fread(fileID,Ng*(Nt+1),'double');
    rho0{i} = reshape(rho0{i},[Ng,Nt+1]);

    fileID = fopen(strcat(dir0{i},'phi.bin'));
    phi0{i} = fread(fileID,Ng*(Nt+1),'double');
    phi0{i} = reshape(phi0{i},[Ng,Nt+1]);

    fileID = fopen(strcat(dir1{i},'E.bin'));
    E1{i} = fread(fileID,Ng*(Nt+1),'double');
    E1{i} = reshape(E1{i},[Ng,Nt+1]);

    fileID = fopen(strcat(dir1{i},'rho.bin'));
    rho1{i} = fread(fileID,Ng*(Nt+1),'double');
    rho1{i} = reshape(rho1{i},[Ng,Nt+1]);

    fileID = fopen(strcat(dir1{i},'phi.bin'));
    phi1{i} = fread(fileID,Ng*(Nt+1),'double');
    phi1{i} = reshape(phi1{i},[Ng,Nt+1]);
end

%% PIC - particles

dL1 = zeros(Nt,5); dLx = dL1; dW1 = dL1;
Ridx = 1:1e2:1e5; Nw = length(Ridx);
Np = 10.^(2:6);

for i=1:Nt
% i = Nt-19;
    for j=1:5
        fileID = fopen(strcat(dir0{j},'xp/',num2str(i),'_1.bin'));
        xp0 = fread(fileID,Np(j),'double');
        fileID = fopen(strcat(dir0{j},'spwt/',num2str(i),'_1.bin'));
        spwt0 = fread(fileID,Np(j),'double');
        fileID = fopen(strcat(dir0{j},'vp/',num2str(i),'_1.bin'));
        vp0 = fread(fileID,Np(j)*3,'double');
        vp0 = reshape(vp0,[Np(j), 3]);

        fileID = fopen(strcat(dir1{j},'xp/',num2str(i),'_1.bin'));
        xp1 = fread(fileID,Np(j),'double');
        fileID = fopen(strcat(dir1{j},'spwt/',num2str(i),'_1.bin'));
        spwt1 = fread(fileID,Np(j),'double');
        fileID = fopen(strcat(dir1{j},'vp/',num2str(i),'_1.bin'));
        vp1 = fread(fileID,Np(j)*3,'double');
        vp1 = reshape(vp1,[Np(j), 3]);

    %     figure(2)
    %     plot(xp0,vp0(:,1),'.k');
    %     drawnow();

        dxp1 = min( [ abs(xp1-xp0), L-abs(xp1-xp0) ], [], 2);
        dL1(i,j) = spwt0'*sqrt( dxp1.^2 + (vp1(:,1) - vp0(:,1)).^2 );
        dLx1(i,j) = spwt0'*dxp1;
    %     dW(i) = dBL_ptc(xp0(Ridx),xp1(Ridx),Nw,L);

        fclose('all');
    end
	i
end

save PIC_dL.dat dL1 -ASCII;
save PIC_dLx.dat dLx1 -ASCII;

