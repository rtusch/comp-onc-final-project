clc
clear all
close all
%% Plotting Preamble
set(gca,'fontname','Times')
set(groot,'defaultAxesTitleFontSizeMultiplier',1)
set(groot,'defaultErrorbarLineWidth', 0.5)
set(groot,'defaultAxesFontName','Times')
set(groot,'defaultAxesPlotBoxAspectRatioMode','manual')
set(groot,'defaultAxesPlotBoxAspectRatio',[1 1 1])
set(groot,'defaultAxesDataAspectRatioMode','auto')
set(groot,'defaultAxesDataAspectRatio',[1 1 1])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% set(groot,'defaultAxesDataAspectRatio','default')

set(groot,'defaultAxesFontWeight','Normal')
set(0,'DefaultAxesTitleFontWeight','normal');
set(groot,'defaultAxesFontSizeMode','manual')
set(groot,'defaultAxesFontSize',25)
% set(groot,'defaultAxesFontSize',25)
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesLabelFontSizeMultiplier',1)
set(groot,'defaultAxesLineWidth',2)
set(groot,'defaultScatterLineWidth',2)
% set(groot,'defaultScatterMarkerFaceColor','k')
% set(groot,'defaultScatterMarkerEdgeColor','k')
set(groot,'defaultScatterMarkerFaceColor','default')
set(groot,'defaultScatterMarkerEdgeColor','default')
set(groot,'defaultLineColor','k')
set(groot,'defaultLineLineWidth',2)
set(groot,'defaultLineMarkerSize',2)
co = [0,0,0;0 0 1;0 0.5 0;1 0 0;0 0.75 0.75;0.75 0 0.75;0.75 0.75 0;0.25 0.25 0.25];
set(groot,'defaultAxesColorOrder',co)

set(0, 'DefaultFigureWindowState', 'normal');
set(groot, 'DefaultFigureVisible', 'on')

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
Dn2 = 1.73E-05*1000; %multiplied because cells werent diffusing in time
d1 = 3.46E-01;
d2 = 5.53E-01;
%H1opt = 3.98E-08;
%H1width = 3.98E-12;
%H2opt = 1.58E-07;
%H2width = 1.60E-11;
H1opt = 7.4;
H1width = 11.4;
H2opt = 6.8;
H2width = 2.5; %10.8
Dh = 4.32E-01/10; %had to divide this by 10 because it was causing difficulties with the finite difference, the behavior still appears essentially the same 
kacid = 1.9E-12/1000; %changed this because uptake/production was dominating diffusion
dh = 2.59E+03/500; %changed this because uptake/production was dominating diffusion
Ho = 3.98E-08;  %pH 7.4 (10^-7.4)
Htumor = 1.58E-07;  %pH 6.8 (10^-6.8)
kneut = 1E6; %estimate
Db = 9.50E-01/10; %had to lower this for the same reason as Dh
Dp = 8.21E-04;
kp = 1.3E-12;
dp = 0.013; %a lot slower: P was being degraded too fast and wasn't changing M
km = 8.64E-3*1000; %make this lower if M() keeps going negative

Btot = 5E-5; %mmol/cm^3 (molar, 50E-6=50uM)
% treattimes = []; %change this to assess treatment efficacy
%treattimes = [5000];
treattimes = [2500 5000 7500];
%treattimes = [1500 3000 4500 6000 7500 9000];
Bpulse = Btot/size(treattimes,2);
Mo = 1.33E-2;
N1o = th1*0.1; %estimated as 1/10th of carrying capacity 

dx = 0.1; %cm
dy = 0.1; %cm
dt = 0.02; %day
sx = 10; %cm
sy = 10; %cm
tfinal = 200; %days


%% Fringe matrix
%build a matrix so that we can determine where the edge of the tumor is
isFringe = zeros(100,100);
radius_squared = 80;
for x = 1:100
    for y = 1:50
        if (y <= 50.5 - sqrt(radius_squared-(x-50.5)^2))
            isFringe(x,y) = 1;
        end
    end
    for y = 51:100
        if (y >= 50.5 + sqrt(radius_squared-(x-50.5)^2))
            isFringe(x,y) = 1;
        end
    end
end


%% Initial Conditions

N1 = zeros(sx/dx, sy/dx, tfinal/dt);
N2 = zeros(sx/dx, sy/dx, tfinal/dt);
H = zeros(sx/dx, sy/dx, tfinal/dt);
B = zeros(sx/dx, sy/dx, tfinal/dt);
M = zeros(sx/dx, sy/dx, tfinal/dt);
P = zeros(sx/dx, sy/dx, tfinal/dt);
P_initial = zeros(sx/dx, sy/dx, 1);

load('cellmaps.mat');

%N1(:, :, 1) = n1init*N1o; 
N2(:, :, 1) = n2init*th2; %assume tumor cells are at carrying capacity

H(:, :, 1) = Ho*n1init + Htumor*n2init; %initialize tumor and healthy tissue with respective pH values
%B(:, :, 1) = Bpulse; %If treatment is immediatelly administered at t=1
M(:, :, 1) = Mo-n2init*0.5*Mo; %assume matrix is half degraded where tumor is 

radius_squared = 100;
for x = linspace(1,100,1000)
    for y = 1:50
        if (y == floor(50 - sqrt(radius_squared-(x-50)^2))) || (y == ceil(50 - sqrt(radius_squared-(x-50)^2)))
            P_initial(floor(x),y,1) = 0.01;
            P_initial(ceil(x),y,1) = 0.01;
        end
    end
    for y = 51:100
        if (y == floor(50 + sqrt(radius_squared-(x-50)^2))) || (y == ceil(50 + sqrt(radius_squared-(x-50)^2)))
            P_initial(floor(x),y,1) = 0.01;
            P_initial(ceil(x),y,1) = 0.01;
        end
    end
end
P(:, :, 1) = P_initial; %ring of MMPs initially around tumor

%make video animation
% myVideo = VideoWriter('Treatment'); %open video file
% myVideo.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
% open(myVideo);

fign = 1;
figure(fign)
subplot(2, 3, 1)
imagesc(N2(:, :, 1))
title("N")
colorbar
caxis([0, th2])
subplot(2, 3, 4)
tlabel = text(0,0,"0 days",'FontSize',20);
set(subplot(2, 3, 4),'visible','off')
subplot(2, 3, 2)
imagesc(-log10(H(:, :, 1)))
title("pH")
colorbar
caxis([6.5, 7.5])
subplot(2, 3, 5)
imagesc(B(:, :, 1))
title("B")
colorbar
caxis([0, Bpulse])
subplot(2, 3, 3)
imagesc(M(:, :, 1))
title("M")
colorbar
caxis([0, Mo])
subplot(2, 3, 6)
imagesc(P(:, :, 1))
title("P")
colorbar
caxis([0, 0.015])
fign = fign+1;

%% run simulation
for t = 2:tfinal/dt
    if mod(t,10) == 0
        disp(t)
    end
    for x = 1:sx/dx
        for y = 1:sy/dy

            if x == 1
                N1_xx = (1/dx^2)*(2*N1(x+1,y,t-1)-2*N1(x,y,t-1)); %second derivative of N1 wrt x
                N2_xx = (1/dx^2)*(2*N2(x+1,y,t-1)-2*N2(x,y,t-1)); %second derivative of N2 wrt x
                H_xx = (1/dx^2)*(2*H(x+1,y,t-1)-2*H(x,y,t-1)); %second derivative of H wrt x
                B_xx = (1/dx^2)*(2*B(x+1,y,t-1)-2*B(x,y,t-1)); %second derivative of B wrt x
                P_xx = (1/dx^2)*(2*P(x+1,y,t-1)-2*P(x,y,t-1)); %second derivative of P wrt x
                ADJ_LIM_X1 = (1-(N1(x+1,y,t-1)/th1)-(N2(x+1,y,t-1)/th2))*(1-(M(x+1,y,t-1)/Mo));
                ADJ_LIM_X2 = 0;
            elseif x == sx/dx
                N1_xx = (1/dx^2)*(2*N1(x-1,y,t-1)-2*N1(x,y,t-1));
                N2_xx = (1/dx^2)*(2*N2(x-1,y,t-1)-2*N2(x,y,t-1));
                H_xx = (1/dx^2)*(2*H(x-1,y,t-1)-2*H(x,y,t-1));
                B_xx = (1/dx^2)*(2*B(x-1,y,t-1)-2*B(x,y,t-1));
                P_xx = (1/dx^2)*(2*P(x-1,y,t-1)-2*P(x,y,t-1));
                ADJ_LIM_X1 = 0;
                ADJ_LIM_X2 = (1-(N1(x-1,y,t-1)/th1)-(N2(x-1,y,t-1)/th2))*(1-(M(x-1,y,t-1)/Mo));
            else
                N1_xx = (1/dx^2)*(N1(x+1,y,t-1)+N1(x-1,y,t-1)-2*N1(x,y,t-1));
                N2_xx = (1/dx^2)*(N2(x+1,y,t-1)+N2(x-1,y,t-1)-2*N2(x,y,t-1));
                H_xx = (1/dx^2)*(H(x+1,y,t-1)+H(x-1,y,t-1)-2*H(x,y,t-1));
                B_xx = (1/dx^2)*(B(x+1,y,t-1)+B(x-1,y,t-1)-2*B(x,y,t-1));
                P_xx = (1/dx^2)*(P(x+1,y,t-1)+P(x-1,y,t-1)-2*P(x,y,t-1));
                ADJ_LIM_X1 = (1-(N1(x+1,y,t-1)/th1)-(N2(x+1,y,t-1)/th2))*(1-(M(x+1,y,t-1)/Mo));
                ADJ_LIM_X2 = (1-(N1(x-1,y,t-1)/th1)-(N2(x-1,y,t-1)/th2))*(1-(M(x-1,y,t-1)/Mo));
            end
            
            if y == 1
                N1_yy = (1/dy^2)*(2*N1(x,y+1,t-1)-2*N1(x,y,t-1));
                N2_yy = (1/dy^2)*(2*N2(x,y+1,t-1)-2*N2(x,y,t-1));
                H_yy = (1/dy^2)*(2*H(x,y+1,t-1)-2*H(x,y,t-1));
                B_yy = (1/dy^2)*(2*B(x,y+1,t-1)-2*B(x,y,t-1));
                P_yy = (1/dy^2)*(2*P(x,y+1,t-1)-2*P(x,y,t-1));
                ADJ_LIM_Y1 = (1-(N1(x,y+1,t-1)/th1)-(N2(x,y+1,t-1)/th2))*(1-(M(x,y+1,t-1)/Mo));
                ADJ_LIM_Y2 = 0;
            elseif y == sy/dy
                N1_yy = (1/dy^2)*(2*N1(x,y-1,t-1)-2*N1(x,y,t-1));
                N2_yy = (1/dy^2)*(2*N2(x,y-1,t-1)-2*N2(x,y,t-1));
                H_yy = (1/dy^2)*(2*H(x,y-1,t-1)-2*H(x,y,t-1));
                B_yy = (1/dy^2)*(2*B(x,y-1,t-1)-2*B(x,y,t-1));
                P_yy = (1/dy^2)*(2*P(x,y-1,t-1)-2*P(x,y,t-1));
                ADJ_LIM_Y1 = 0;
                ADJ_LIM_Y2 = (1-(N1(x,y-1,t-1)/th1)-(N2(x,y-1,t-1)/th2))*(1-(M(x,y-1,t-1)/Mo));
            else
                N1_yy = (1/dy^2)*(N1(x,y+1,t-1)+N1(x,y-1,t-1)-2*N1(x,y,t-1));
                N2_yy = (1/dy^2)*(N2(x,y+1,t-1)+N2(x,y-1,t-1)-2*N2(x,y,t-1));
                H_yy = (1/dy^2)*(H(x,y+1,t-1)+H(x,y-1,t-1)-2*H(x,y,t-1));
                B_yy = (1/dy^2)*(B(x,y+1,t-1)+B(x,y-1,t-1)-2*B(x,y,t-1));
                P_yy = (1/dy^2)*(P(x,y+1,t-1)+P(x,y-1,t-1)-2*P(x,y,t-1));
                ADJ_LIM_Y1 = (1-(N1(x,y+1,t-1)/th1)-(N2(x,y+1,t-1)/th2))*(1-(M(x,y+1,t-1)/Mo));
                ADJ_LIM_Y2 = (1-(N1(x,y-1,t-1)/th1)-(N2(x,y-1,t-1)/th2))*(1-(M(x,y-1,t-1)/Mo));
            end

            CARCAP_LIM = 1-(N1(x,y,t-1)/th1)-(N2(x,y,t-1)/th2); %limiting term for cell # due to carcap
            MAT_LIM = 1-(M(x,y,t-1)/Mo); %limiting term due to ECM
            LIM_TOT = CARCAP_LIM*MAT_LIM;

            %diffusion limited by the matrix concentration surrounding the cells
            %the limiting term ADJ_LIM is the average of the surrounding voxels
            ADJ_LIM = 0.2*(LIM_TOT+ADJ_LIM_X1+ADJ_LIM_X2+ADJ_LIM_Y1+ADJ_LIM_Y2);

            N1_PLF = k1*N1(x,y,t-1)*LIM_TOT; %proliferative term
            N1_DIF = Dn1*ADJ_LIM*(N1_xx + N1_yy); %diffusion term
            %N1_PH = -d1*(1-exp(-1*(((H(x,y,t-1)-H1opt)/H1width)^2)))*N1(x,y,t-1); %pH-dependence
            N1_PH = -d1*(1-exp(-1*((((-1*log10(H(x,y,t-1)))-H1opt)/H1width)^2)))*N1(x,y,t-1); %pH-dependence
            N1(x,y,t) = N1(x,y,t-1) + dt*(N1_PLF + N1_DIF + N1_PH);

            N2_PLF = k2*N2(x,y,t-1)*LIM_TOT; %proliferative term
            N2_DIF = Dn2*ADJ_LIM*(N2_xx + N2_yy); %diffusion term
            N2_PH = -d2*(1-exp(-1*((((-1*log10(H(x,y,t-1)))-H2opt)/H2width)^2)))*N2(x,y,t-1); %pH-dependence
            N2(x,y,t) = N2(x,y,t-1) + dt*(N2_PLF + N2_DIF + N2_PH);

            H_DIF = Dh*(H_xx + H_yy); %diffusion
            H_PROD = kacid*N2(x,y,t-1); %acid production by tumor cells
            H_UPT = -dh*(H(x,y,t-1)-Ho);  %hydrogen uptake
            H_NEU = -kneut*H(x,y,t-1)*B(x,y,t-1);   %acid-base neutralization
            H(x,y,t) = H(x,y,t-1) + dt*(H_DIF + H_PROD + H_UPT + H_NEU);

            B_DIF = Db*(B_xx + B_yy); %diffusion
            B_NEU = -kneut*H(x,y,t-1)*B(x,y,t-1);   %acid-base neutralization
            B(x,y,t) = B(x,y,t-1) + dt*(B_DIF + B_NEU);

            M_DEG = -1*km*(0.9*exp(-1*(-1*log10(H(x,y,t-1))-5.9)^2)+0.1)*P(x,y,t-1)*M(x,y,t-1);
            M(x,y,t) = M(x,y,t-1) + dt*M_DEG;

            P_DIF = Dp*(P_xx + P_yy);
            if isFringe(x,y)
                P_PROD = 1*kp*N2(x,y,t-1); %MMP production by tumor cells on outside of tumor
            else
                P_PROD = 0.01*kp*N2(x,y,t-1); %MMP production by tumor cells at core of tumor
            end
            P_DEG = -dp*P(x,y,t-1);   %MMP degredation
            P(x,y,t) = P(x,y,t-1) + dt*(P_DIF + P_PROD + P_DEG);


            %DRUG TREATMENT: change this part to try different treatments
            if ismember(t, treattimes)
                B(:, :, t) = Bpulse;
            end

            %im going to do a check here for any negative values
            %things keep being negative and screwing all the other
            %equations up
            N1(x,y,t) = checkbounds(N1(x,y,t),0,th1,'N1');
            N2(x,y,t) = checkbounds(N2(x,y,t),0,th2,'N2');
            H(x,y,t) = checkbounds(H(x,y,t),0,1,'H'); %1mmol/cm^3 is pH 0
            B(x,y,t) = checkbounds(B(x,y,t),0,1,'B');
            M(x,y,t) = checkbounds(M(x,y,t),0,Mo,'M'); %should not ever increase over starting concentration
            P(x,y,t) = checkbounds(P(x,y,t),0,1,'P');
        end
    end
    if mod(t, 100) == 0
        % print plots
        disp(t);
        figure(1)
        subplot(2, 3, 1)
        imagesc(N2(:, :, t))
        title("N")
        colorbar
        caxis([0, th2])
        subplot(2, 3, 4)
        delete(tlabel)
        tlabel = text(0,0,string(t*dt)+" days",'FontSize',20);
        set(subplot(2, 3, 4),'visible','off')
        subplot(2, 3, 2)
        imagesc(-log10(H(:, :, t)))
        title("pH")
        colorbar
        caxis([6, 8])
        subplot(2, 3, 5)
        imagesc(B(:, :, t))
        title("B")
        colorbar
        caxis([0, Bpulse])
        subplot(2, 3, 3)
        imagesc(M(:, :, t))
        title("M")
        colorbar
        caxis([0, Mo])
        subplot(2, 3, 6)
        imagesc(P(:, :, t))
        title("P")
        colorbar
        caxis([0, 0.015])

        fign = fign+1;
        drawnow;
%         frame = getframe(gcf); %get frame
%         writeVideo(myVideo, frame);
%         if t == 100
%             subplot(2, 1, 1)
%             imagesc(N1(:, :, 1))
%             colorbar
%             subplot(2, 1, 2)
%             imagesc(N2(:, :, 1))
%             colorbar
%             drawnow
%         end
%         subplot(2, 1, 1)
%         imagesc(N1(:, :, t))
%         colorbar
%         subplot(2, 1, 2)
%         imagesc(N2(:, :, t))
%         colorbar
%         drawnow
        

    end
end

% close(myVideo)
N2_fringe = zeros(100,100,tfinal/dt);
N2_core = zeros(100,100,tfinal/dt);
for x = 1:100
    for y = 1:100
        if isFringe(x,y)
            N2_fringe(x,y,:) = N2(x,y,:);
        else
            N2_core(x,y,:) = N2(x,y,:);
        end
    end
end
% 
% figure(2)
% plot((1:200/dt)*dt, squeeze(sum(sum(N2)))/(10^10))
% hold on;
% plot((1:200/dt)*dt, squeeze(sum(sum(N2_core(:,:,:))))/(10^10))
% plot((1:200/dt)*dt, squeeze(sum(sum(N2_fringe(:,:,:))))/(10^10))
% % ylim([0,60])
% xlabel("Time (Days)")
% ylabel("Number of Cells ($10^{10}$ Cells)")
% legend("Total Number of Tumor Cells", "Cells at the Core", "Cells at the Fringe")
% hold off;
% 
% figure(3)
% plot((1:200/dt)*dt, squeeze(-log10(H(50, 50, :))))
% xlabel("Time (Days)")
% ylabel("pH at center of tumor")
% 
% figure(4)
% plot(1:100, 1000*squeeze(P(:,50,1)))
% hold on;
% plot(1:100, 1000*squeeze(P(:,50,1250)))
% plot(1:100, 1000*squeeze(P(:,50,5000)))
% plot(1:100, 1000*squeeze(P(:,50,10000)))
% title("MMP Concentration Through Tumor at y = 50 mm")
% xlabel("Distance (mm)")
% ylabel("MMP Concentration (mM)")
% legend("0 Days","25 Days", "100 Days", "200 Days")
% hold off;
% 
% figure(5)
% plot(1:100, 1000*squeeze(M(:,50,1)))
% hold on;
% plot(1:100, 1000*squeeze(M(:,50,1250)))
% plot(1:100, 1000*squeeze(M(:,50,5000)))
% plot(1:100, 1000*squeeze(M(:,50,10000)))
% title("ECM Concentration Through Tumor at y = 50 mm")
% xlabel("Distance (mm)")
% ylabel("ECM Concentration (mM)")
% legend("0 Days","25 Days", "100 Days", "200 Days")
% hold off;

%Use this code for saving variables for plotting later:
% Nt0 = N2;
% pHt0 = -log10(H);
% Bt0 = B;
% Mt0 = M;
% Pt0 = P;
% save('t0.mat','Nt0','pHt0','Bt0','Mt0','Pt0')

function r = checkbounds(f, lb, ub, var)
    if f < lb
        %disp(["Lower bound reached: ", var, f])
        f = lb;
    elseif f > ub
        %disp(["Upper bound reached: ", var, f])
        f = ub;
    end
    r = f;
end
