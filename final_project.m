clc
clear all
close all


%% initialize variables
% refer to the google doc for a description of each parameter & their sources:
% https://docs.google.com/spreadsheets/d/13LxiwIeWfVTxi7gsmRAGJGKGRv3-lwYtxaVlOxKLEW8/edit#gid=0
k1 = 0.1728;
k2 = 2.76E-01;
th1 = 5.00E+08;
th2 = 5.00E+08;
a1 = 5.00E+08;
a2 = 5.00E+08;
Dn1 = 1.73E-05;
Dn2 = 1.73E-05;
d1 = 3.46E-01;
d2 = 5.53E-01;
H1opt = 3.98E-08;
H1width = 3.98E-12;
H2opt = 1.58E-07;
H2width = 1.60E-11;
Dh = 4.32E-01;
kacid = 2.55E-22;
dh = 2.59E+03;
Ho = 3.98E-08;
kneut = 1; %estimate
Db = 9.50E-01;
Dp = 8.21E-04;
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

%% Initial Conditions

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

%% run simulation
for t = 2:tfinal/dt
    for x = 1:sx/dx
        for y = 1:sy/dy
            %z = 1-(N1(x,y,t)/(th1-a1*M(x,y,t)))-(N2(x,y,t)/(th2-a2*M(x,y,t)))

            if x == 1
                N1_xx = (1/dx^2)*(2*N1(x+1,y,t-1)-2*N1(x,y,t-1)); %second derivative of N1 wrt x
                N2_xx = (1/dx^2)*(2*N2(x+1,y,t-1)-2*N2(x,y,t-1)); %second derivative of N2 wrt x
                H_xx = (1/dx^2)*(2*H(x+1,y,t-1)-2*H(x,y,t-1)); %second derivative of H wrt x
                B_xx = (1/dx^2)*(2*B(x+1,y,t-1)-2*B(x,y,t-1)); %second derivative of B wrt x
                P_xx = (1/dx^2)*(2*P(x+1,y,t-1)-2*P(x,y,t-1)); %second derivative of P wrt x
            elseif x == sx/dx
                N1_xx = (1/dx^2)*(2*N1(x-1,y,t-1)-2*N1(x,y,t-1));
                N2_xx = (1/dx^2)*(2*N2(x-1,y,t-1)-2*N2(x,y,t-1));
                H_xx = (1/dx^2)*(2*H(x-1,y,t-1)-2*H(x,y,t-1));
                B_xx = (1/dx^2)*(2*B(x-1,y,t-1)-2*B(x,y,t-1));
                P_xx = (1/dx^2)*(2*P(x-1,y,t-1)-2*P(x,y,t-1));
            else
                N1_xx = (1/dx^2)*(N1(x+1,y,t-1)+N1(x-1,y,t-1)-2*N1(x,y,t-1));
                N2_xx = (1/dx^2)*(N2(x+1,y,t-1)+N2(x-1,y,t-1)-2*N2(x,y,t-1));
                H_xx = (1/dx^2)*(H(x+1,y,t-1)+H(x-1,y,t-1)-2*H(x,y,t-1));
                B_xx = (1/dx^2)*(B(x+1,y,t-1)+B(x-1,y,t-1)-2*B(x,y,t-1));
                P_xx = (1/dx^2)*(P(x+1,y,t-1)+P(x-1,y,t-1)-2*P(x,y,t-1));
            end
            
            if y == 1
                N1_yy = (1/dy^2)*(2*N1(x,y+1,t-1)-2*N1(x,y,t-1));
                N2_yy = (1/dy^2)*(2*N2(x,y+1,t-1)-2*N2(x,y,t-1));
                H_yy = (1/dy^2)*(2*H(x,y+1,t-1)-2*H(x,y,t-1));
                B_yy = (1/dy^2)*(2*B(x,y+1,t-1)-2*B(x,y,t-1));
                P_yy = (1/dy^2)*(2*P(x,y+1,t-1)-2*P(x,y,t-1));
            elseif y == sy/dy
                N1_yy = (1/dy^2)*(2*N1(x,y-1,t-1)-2*N1(x,y,t-1));
                N2_yy = (1/dy^2)*(2*N2(x,y-1,t-1)-2*N2(x,y,t-1));
                H_yy = (1/dy^2)*(2*H(x,y-1,t-1)-2*H(x,y,t-1));
                B_yy = (1/dy^2)*(2*B(x,y-1,t-1)-2*B(x,y,t-1));
                P_yy = (1/dy^2)*(2*P(x,y-1,t-1)-2*P(x,y,t-1));
            else
                N1_yy = (1/dy^2)*(N1(x,y+1,t-1)+N1(x,y-1,t-1)-2*N1(x,y,t-1));
                N2_yy = (1/dy^2)*(N2(x,y+1,t-1)+N2(x,y-1,t-1)-2*N2(x,y,t-1));
                H_yy = (1/dy^2)*(H(x,y+1,t-1)+H(x,y-1,t-1)-2*H(x,y,t-1));
                B_yy = (1/dy^2)*(B(x,y+1,t-1)+B(x,y-1,t-1)-2*B(x,y,t-1));
                P_yy = (1/dy^2)*(P(x,y+1,t-1)+P(x,y-1,t-1)-2*P(x,y,t-1));
            end

            CARCAP_LIM = 1-(N1(x,y,t-1)/th1)-(N2(x,y,t-1)/th2); %limiting term for cell # due to carcap
            MAT_LIM = 1-(M(x,y,t-1)/Mo); %limiting term due to ECM
            LIM_TOT = CARCAP_LIM*MAT_LIM;

            N1_PLF = k1*N1(x,y,t-1)*LIM_TOT; %proliferative term
            N1_DIF = Dn1*LIM_TOT*(N1_xx + N1_yy); %diffusion term (this isn't right because I didn't do del • lim term)
%             N1_PH = -d1*(1-exp(((H(x,y,t-1)-H1opt)/H1width)^2))*N1(x,y,t-1); %pH-dependence
            N1(x,y,t) = N1(x,y,t-1) + dt*(N1_PLF + N1_DIF);

            N2_PLF = k2*N2(x,y,t-1)*LIM_TOT; %proliferative term
            N2_DIF = Dn2*LIM_TOT*(N2_xx + N2_yy); %diffusion term (this isn't right because I didn't do del • lim term)
%             N2_PH = -d2*(1-exp(((H(x,y,t-1)-H2opt)/H2width)^2))*N2(x,y,t-1); %pH-dependence
            N2(x,y,t) = N2(x,y,t-1) + dt*(N2_PLF + N2_DIF);

            H_DIF = Dh*(H_xx + H_yy); %diffusion
            H_PROD = kacid*N2(x,y,t-1); %acid production by tumor cells
            H_UPT = -dh*(H(x,y,t-1)-Ho);  %hydrogen uptake
            H_NEU = -kneut*H(x,y,t-1)*B(x,y,t-1);   %acid-base neutralization
            H(x,y,t) = H(x,y,t-1) + dt*(H_DIF + H_PROD + H_UPT + H_NEU);

            B_DIF = Db*(B_xx + B_yy); %diffusion
            B_NEU = -kneut*H(x,y,t-1)*B(x,y,t-1);   %acid-base neutralization
            B(x,y,t) = B(x,y,t-1) + dt*(B_DIF + B_NEU);

            M(x,y,t) = (-(3.4*exp(-1*(-1*log10(H(x,y,t-1))-5.9)^2)+0.5)*P(x,y,t-1)*M(x,y,t-1)...
                )*dt+M(x,y,t-1);

            P_DIF = Dp*(P_xx + P_yy);
            P_PROD = kp*N2(x,y,t-1); %MMP production by tumor cells
            P_DEG = -dp*P(x,y,t-1);   %MMP degredation
            P(x,y,t) = P(x,y,t-1) + dt*(P_DIF + P_PROD + P_DEG);
        end
    end
    if mod(t, 100) == 0
        % print plots
%         disp(t);
%         figure(fign)
%         subplot(2, 3, 1)
%         imagesc(N1(:, :, t))
%         colorbar
%         subplot(2, 3, 4)
%         imagesc(N2(:, :, t))
%         colorbar
%         subplot(2, 3, 2)
%         imagesc(H(:, :, t))
%         colorbar
%         subplot(2, 3, 5)
%         imagesc(B(:, :, t))
%         colorbar
%         subplot(2, 3, 3)
%         imagesc(M(:, :, t))
%         colorbar
%         subplot(2, 3, 6)
%         imagesc(P(:, :, t))
%         colorbar
% 
%         fign = fign+1;
%         drawnow;

        if t == 100
            subplot(2, 1, 1)
            imagesc(N1(:, :, 1))
            colorbar
            subplot(2, 1, 2)
            imagesc(N2(:, :, 1))
            colorbar
            drawnow
        end
        subplot(2, 1, 1)
        imagesc(N1(:, :, t))
        colorbar
        subplot(2, 1, 2)
        imagesc(N2(:, :, t))
        colorbar
        drawnow
        

    end
end
