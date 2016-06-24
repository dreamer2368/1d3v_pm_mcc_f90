clear all
close all
clc

N=64;
fileID = fopen('xg.bin');
xg = fread(fileID,N,'double');
fileID = fopen('phi.bin');
phi = fread(fileID,N,'double');
fileID = fopen('phi_sol.bin');
sol = fread(fileID,N,'double');

%%

plot(xg,phi,'.k',xg,sol,'-r');