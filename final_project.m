clc
clear all
close all


% initialize variables
% refer to the google doc for a description of each parameter & their sources:
% https://docs.google.com/spreadsheets/d/13LxiwIeWfVTxi7gsmRAGJGKGRv3-lwYtxaVlOxKLEW8/edit#gid=0
k1 = 2.31E-11;
k2 = 3.70E-11;
th1 = 5.79E+03;
th2 = 5.79E+03;
a1 = 5.79E+03;
a2 = 5.79E+03;
Dn1 = 2.31E-15;
Dn2 = 2.31E-15;
d1 = 4.63E-11;
d2 = 7.41E-11;
H1opt = 3.98E-08;
H1width = 3.98E-12;
H2opt = 1.58E-07;
H2width = 1.60E-11;
Dh = 5.79E-11;
kacid = 2.55E-22;
dh = 3.47E-07;
Ho = 3.98E-08;
kneut = 1; %estimate
Db = 1.27E-10;
Dp = 1.10E-13;
kp = 0.01; %estimate
dp = 0.001; %estimate

Bpulse = 1; %estimate
Mo = 1; %estimate
N1o = th1*0.1; %estimated as 1/10th of carrying capacity 

dx = 0.1; %cm
dy = 0.1; %cm
dt = 0.02; %day
sx = 10; %cm
sy = 10; %cm
tfinal = 20; %days

N1 = zeros(sx/dx, sy/dx, tfinal/dt);
N2 = zeros(sx/dx, sy/dx, tfinal/dt);
H = zeros(sx/dx, sy/dx, tfinal/dt);
B = zeros(sx/dx, sy/dx, tfinal/dt);
M = zeros(sx/dx, sy/dx, tfinal/dt);
P = zeros(sx/dx, sy/dx, tfinal/dt);

load('cellmaps.mat');

N1(:, :, 1) = n1init*N1o; 
N2(:, :, 1) = n2init*th2; %assume tumor cells are at carrying capacity
H(:, :, 1) = Ho;
B(:, :, 1) = Bpulse; %If treatment is immediatelly administered at t=1
M(:, :, 1) = Mo; 
P(:, :, 1) = 0; %probably dont need P initial condition, but maybe?


% run simulation
for t = 1:tfinal/dt
    for x = 1:sx/dx
        for y = 1:sy/dy
        end
    end
end



% print plots

subplot(2, 3, 1)
imagesc(N1(:, :, 1))
colorbar
subplot(2, 3, 4)
imagesc(N2(:, :, 1))
colorbar
subplot(2, 3, 2)
imagesc(H(:, :, 1))
colorbar
subplot(2, 3, 5)
imagesc(B(:, :, 1))
colorbar
subplot(2, 3, 3)
imagesc(M(:, :, 1))
colorbar
subplot(2, 3, 6)
imagesc(P(:, :, 1))
colorbar