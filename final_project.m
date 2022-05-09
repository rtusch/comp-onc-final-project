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

fign = 1;

% run simulation
for t = 2:tfinal/dt
    for x = 2:sx/dx-1     %im skipping the boundaries for now because i really
        for y = 2:sy/dy-1 %dont want to deal with the goddamn boudary conditions rn
            %z = 1-(N1(x,y,t)/(th1-a1*M(x,y,t)))-(N2(x,y,t)/(th2-a2*M(x,y,t)))
            N1(x,y,t) = ...
                (k1*N1(x,y,t-1)*(1-(N1(x,y,t-1)/(th1-(a1*M(x,y,t-1))))-(N2(x,y,t-1)/(th2-(a2*M(x,y,t-1)))))... %proliferation
                +Dn1*((1-(N1(x,y,t-1)/(th1-a1*M(x,y,t-1)))-(N2(x,y,t-1)/(th2-a2*M(x,y,t-1))))*(((N1(x-1,y,t-1)-2*N1(x,y,t-1)+N1(x+1,y,t-1))/dx^2) ...
                +(1/dx)*((N1(x+1,y+1,t-1)-N1(x+1,y-1,t-1))/(2*dy)-(N1(x-1,y+1,t-1)-N1(x-1,y-1,t-1))/(2*dy)) ...
                +((N1(x,y-1,t-1)-2*N1(x,y,t-1)+N1(x,y+1,t-1))/dy^2)) ...
                +((1-(N1(x+1,y,t)/(th1-a1*M(x+1,y,t)))-(N2(x+1,y,t)/(th2-a2*M(x+1,y,t)))-1+(N1(x-1,y,t)/(th1-a1*M(x-1,y,t)))+(N2(x-1,y,t)/(th2-a2*M(x-1,y,t))))/dx)...
                *((N1(x+1,y,t-1)-N1(x-1,y,t-1))/dx+(N1(x,y+1,t-1)-N1(x,y-1,t-1))/dy) ...
                +((1-(N1(x,y+1,t)/(th1-a1*M(x,y+1,t)))-(N2(x,y+1,t)/(th2-a2*M(x,y+1,t)))-1+(N1(x,y-1,t)/(th1-a1*M(x,y-1,t)))+(N2(x,y-1,t)/(th2-a2*M(x,y-1,t))))/dx)...
                *((N1(x+1,y,t-1)-N1(x-1,y,t-1))/dx+(N1(x,y+1,t-1)-N1(x,y-1,t-1))/dy))... %diffusion
                -d1*(1-exp(((H(x,y,t-1)-H1opt)/H1width)^2))... %pH-dependence
                )*dt+N1(x,y,t-1); %finite difference stuff

            N2(x,y,t) = ...
                (k2*N2(x,y,t-1)*(1-(N1(x,y,t-1)/(th1-(a1*M(x,y,t-1))))-(N2(x,y,t-1)/(th2-(a2*M(x,y,t-1)))))... %proliferation
                +Dn2*((1-(N1(x,y,t-1)/(th1-a1*M(x,y,t-1)))-(N2(x,y,t-1)/(th2-a2*M(x,y,t-1))))*(((N2(x-1,y,t-1)-2*N2(x,y,t-1)+N2(x+1,y,t-1))/dx^2) ...
                +(1/dx)*((N2(x+1,y+1,t-1)-N2(x+1,y-1,t-1))/(2*dy)-(N2(x-1,y+1,t-1)-N2(x-1,y-1,t-1))/(2*dy)) ...
                +((N2(x,y-1,t-1)-2*N2(x,y,t-1)+N2(x,y+1,t-1))/dy^2)) ...
                +((1-(N1(x+1,y,t)/(th1-a1*M(x+1,y,t)))-(N2(x+1,y,t)/(th2-a2*M(x+1,y,t)))-1+(N1(x-1,y,t)/(th1-a1*M(x-1,y,t)))+(N2(x-1,y,t)/(th2-a2*M(x-1,y,t))))/dx)...
                *((N2(x+1,y,t-1)-N2(x-1,y,t-1))/dx+(N2(x,y+1,t-1)-N2(x,y-1,t-1))/dy) ...
                +((1-(N1(x,y+1,t)/(th1-a1*M(x,y+1,t)))-(N2(x,y+1,t)/(th2-a2*M(x,y+1,t)))-1+(N1(x,y-1,t)/(th1-a1*M(x,y-1,t)))+(N2(x,y-1,t)/(th2-a2*M(x,y-1,t))))/dx)...
                *((N2(x+1,y,t-1)-N2(x-1,y,t-1))/dx+(N2(x,y+1,t-1)-N2(x,y-1,t-1))/dy))... %diffusion
                -d2*(1-exp(((H(x,y,t-1)-H2opt)/H2width)^2))... %pH-dependence
                )*dt+N2(x,y,t-1); %finite difference stuff

            H(x,y,t) = (Dh*((H(x-1,y,t-1)-2*H(x,y,t-1)+H(x+1,y,t-1))/dx^2 ...
                +(H(x,y-1,t-1)-2*H(x,y,t-1)+H(x,y+1,t-1))/dy^2)... %diffusion
                +kacid*N2(x,y,t-1)...%acid production by tumor cells
                -dh*(H(x,y,t-1)-Ho)...  %hydrogen uptake
                -kneut*H(x,y,t-1)*B(x,y,t-1)...   %acid-base neutralization
                )*dt+H(x,y,t-1); %finite difference

            B(x,y,t) = (Db*((B(x-1,y,t-1)-2*B(x,y,t-1)+B(x+1,y,t-1))/dx^2 ...
                +(B(x,y-1,t-1)-2*B(x,y,t-1)+B(x,y+1,t-1))/dy^2)... %diffusion
                -kneut*H(x,y,t-1)*B(x,y,t-1)...   %acid-base neutralization
                )*dt+B(x,y,t-1); %finite difference

            M(x,y,t) = (-(3.4*exp(-1*(-1*log10(H(x,y,t-1))-5.9)^2)+0.5)*P(x,y,t-1)*M(x,y,t-1)...
                )*dt+M(x,y,t-1);

            P(x,y,t) = (Dp*((P(x-1,y,t-1)-2*P(x,y,t-1)+P(x+1,y,t-1))/dx^2 ...
                +(P(x,y-1,t-1)-2*P(x,y,t-1)+P(x,y+1,t-1))/dy^2)... %diffusion
                +kp*N2(x,y,t-1)... %MMP production by tumor cells
                -dp*P(x,y,t-1)...   %MMP degredation
                )*dt+P(x,y,t-1);
        end
    end
    if mod(t, 100) == 0
        % print plots
        figure(fign)
        disp(t);
        subplot(2, 3, 1)
        imagesc(N1(:, :, t))
        colorbar
        subplot(2, 3, 4)
        imagesc(N2(:, :, t))
        colorbar
        subplot(2, 3, 2)
        imagesc(H(:, :, t))
        colorbar
        subplot(2, 3, 5)
        imagesc(B(:, :, t))
        colorbar
        subplot(2, 3, 3)
        imagesc(M(:, :, t))
        colorbar
        subplot(2, 3, 6)
        imagesc(P(:, :, t))
        colorbar

        fign = fign+1;
    end
end
