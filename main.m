%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paper title: Utility maximization for large-scale cell-free massive MIMO downlink
% Journal: IEEE Transactions on Communications
% Authors: Muhammad Farooq, Hien Quoc Ngo, and Le Nam Tran
% Written by: Muhammad Farooq
% Email: Muhammad.Farooq@ucdconnect.ie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
warning('off','all')
rng('default')

% Consider a rectangular area with DxD m^2 nAPs distributed APs serves nUsers terminals
% they all randomly located in the area
% profile on
tic
nAPs = 20; %number of access points (APs)
nUsers=5; %number of users
nTx=1; %number of antennas/AP
B=20; %bandwidth in Mhz

Tp=20; %uplink training phase in samples
Tc=200; %coherence time in samples
D=1; %area length in kilometer

Hb = 15; %base station height in m
Hm = 1.65; %mobile height in m
f = 1900; %frequency in MHz
aL = (1.1*log10(f)-0.7)*Hm-(1.56*log10(f)-0.8);
L = 46.3+33.9*log10(f)-13.82*log10(Hb)-aL; %path loss constant
d0=0.01;%reference distance in km
d1=0.05;%reference distance in km

power_f=nTx*1; %downlink power: 1W
noise_p = 10^((-203.975+10*log10(B*10^6)+9)/10); %noise power
zeta_d = power_f/noise_p;%nomalized downlink power
zeta_p=0.2/noise_p;%normalized pilot power
sigma_sh=8; %log-normal shadowing in dB

[U,S,V11]=svd(randn(Tp,Tp));%U includes tau orthogonal sequences
%U=ones(tau,tau);

maxIter = 1000; %maximum iterations
alphay = 0.01; %stepsize
alphax = 0.01; %stepsize

% Large-scale fading matrix
myBeta=get_slow_fading(nAPs,nUsers,L,D,d0,d1,sigma_sh);

% Pilot Asignment: (random choice)
myPsi=zeros(Tp,nUsers);
if Tp<nUsers
    myPsi(:,[1:1:Tp])=U;
    for k=(Tp+1):nUsers
        Point=randi([1,Tp]);
        myPsi(:,k)=U(:,Point);
    end
else
    myPsi=U(:,[1:1:nUsers]);
end

% Create Nu matrix
mau=zeros(nAPs,nUsers);
for m=1:nAPs
    for k=1:nUsers
        mau(m,k)=norm( (myBeta(m,:).^(1/2)).*(myPsi(:,k)'*myPsi))^2;
    end
end

myNu=Tp*zeta_p*(myBeta.^2)./(Tp*zeta_p*mau + 1);

%% Sum-Rate Maximization

x_prev = rand(nAPs*nUsers,1); %previous x
t_now=1; %current value of t
t_prev=1; %previous value of t
z=x_prev; % previous projected extrapolated point
x_now=x_prev; %current x
bestobj = zeros(maxIter,1);
for nIter=1:maxIter
    y = x_now+(t_prev/t_now)*(z-x_now)+((t_prev-1)/t_now)*(x_now-x_prev); %extrapolated point
    Dy = gradSumRate_v2(nAPs,nTx,nUsers,myNu,myBeta,myPsi,zeta_d,Tc,Tp,y); %gradient of extrapolated point
    y_next= y+alphay*Dy; %take a step in the direction of gradient to find the new y

    Dx = gradSumRate_v2(nAPs,nTx,nUsers,myNu,myBeta,myPsi,zeta_d,Tc,Tp,x_now); %gradient of current x
    x_next = x_now+alphax*Dx; %take a step in the direction of gradient to find the new x
    
    z = proj_Euclidean_ball_posorthant(y_next,nTx,nUsers,nAPs); %projection of new y
    v = proj_Euclidean_ball_posorthant(x_next,nTx,nUsers,nAPs); %projection of new x
    
    Fz = computeSumRate(nAPs,nTx,nUsers,myNu,myBeta,myPsi,zeta_d,Tc,Tp,z); %caclulate objective at z
    Fv = computeSumRate(nAPs,nTx,nUsers,myNu,myBeta,myPsi,zeta_d,Tc,Tp,v); %calculate objective at v

    x_prev = x_now; %assing current x to previous x
    if(Fz>Fv)
        x_now = z; %if objective at z is higher then z is the next x
    else
        x_now = v; %if objective at v is higher then v is the next x
    end
    
    bestobj(nIter)=max(Fz,Fv); %the objective is set to the maximum objective among the two
    t_prev = t_now; %set current t to previous t
    t_now = (sqrt(4*t_now^2)+1)/2; %find current t for next iteration
end
SE_FOM(1) = bestobj(end)/nUsers; %find final value of average sum rate
figure
plot(1:maxIter,bestobj/nUsers,'k','DisplayName','Average Sum Rate') %plot convergence
hold on

%% Proportional Fariness

x_prev = rand(nAPs*nUsers,1); %previous x
t_now=1; %current value of t
t_prev=1; %previous value of t
z=x_prev; % previous projected extrapolated point
x_now=x_prev; %current x
bestobj = zeros(maxIter,1);
for nIter=1:maxIter
    y = x_now+(t_prev/t_now)*(z-x_now)+((t_prev-1)/t_now)*(x_now-x_prev); %extrapolated point
    Dy = gradFairnessRate(nAPs,nTx,nUsers,myNu,myBeta,myPsi,zeta_d,Tc,Tp,y); %gradient at y
    y_next= y+alphay*Dy; %new y

    Dx = gradFairnessRate(nAPs,nTx,nUsers,myNu,myBeta,myPsi,zeta_d,Tc,Tp,x_now); %gradient at x
    x_next = x_now+alphax*Dx; %new x
    
    z = proj_Euclidean_ball_posorthant(y_next,nTx,nUsers,nAPs); %find z
    v = proj_Euclidean_ball_posorthant(x_next,nTx,nUsers,nAPs); %find v
    
    Fz = propFairnessRate(nAPs,nTx,nUsers,myNu,myBeta,myPsi,zeta_d,Tc,Tp,z); %calculate objective at z
    Fv = propFairnessRate(nAPs,nTx,nUsers,myNu,myBeta,myPsi,zeta_d,Tc,Tp,v); %calcualte objective at v

    x_prev = x_now; %assing current x to previous x
    if(Fz>Fv)
        x_now = z; %if objective at z is higher then z is the next x
    else
        x_now = v; %if objective at v is higher then v is the next x
    end
    
    bestobj(nIter)=max(Fz,Fv); %the objective is set to the maximum objective among the two
    t_prev = t_now; %set current t to previous t
    t_now = (sqrt(4*t_now^2)+1)/2; %find current t for next iteration
end
SE_FOM(2)= bestobj(end); %find value of proportional fariness rate
plot(1:maxIter,bestobj,'r','DisplayName','Fairness Rate') %plot convergence
%% Harmonic Rate Maximization

t_now=1; %current value of t
t_prev=1; %previous value of t
z=x_prev; % previous projected extrapolated point
x_now=x_prev; %current x
bestobj = zeros(maxIter,1);
for nIter=1:maxIter
    y = x_now+(t_prev/t_now)*(z-x_now)+((t_prev-1)/t_now)*(x_now-x_prev); %extrapolated point 
    Dy = gradHarmonicRate_v2(nAPs,nTx,nUsers,myNu,myBeta,myPsi,zeta_d,Tc,Tp,y); %gradient at y
    y_next= y+alphay*Dy; %new y

    Dx = gradHarmonicRate_v2(nAPs,nTx,nUsers,myNu,myBeta,myPsi,zeta_d,Tc,Tp,x_now); %gradient at x
    x_next = x_now+alphax*Dx; %new x
    
    z = proj_Euclidean_ball_posorthant(y_next,nTx,nUsers,nAPs); %find z
    v = proj_Euclidean_ball_posorthant(x_next,nTx,nUsers,nAPs); %find v
    
    Fz = harmonicSumRate_v2(nAPs,nTx,nUsers,myNu,myBeta,myPsi,zeta_d,Tc,Tp,z); %objective at z   
    Fv = harmonicSumRate_v2(nAPs,nTx,nUsers,myNu,myBeta,myPsi,zeta_d,Tc,Tp,v); %objective at v

    x_prev = x_now; %assing current x to previous x
    if(Fz>Fv)
        x_now = z; %if objective at z is higher then z is the next x
    else
        x_now = v; %if objective at v is higher then v is the next x
    end
    
    bestobj(nIter)=max(Fz,Fv); %the objective is set to the maximum objective among the two
    t_prev = t_now; %set current t to previous t
    t_now = (sqrt(4*t_now^2)+1)/2; %find current t for next iteration
end
SE_FOM(3)= bestobj(end); %final value of harmonic rate
plot(1:maxIter,bestobj,'b','DisplayName','Harmonic Rate') %plot convergence

%% Min Rate Maximization

t_now=1; %current value of t
t_prev=1; %previous value of t
z=x_prev; % previous projected extrapolated point
x_now=x_prev; %current x
mygamma = 100; %smoothing parameter
bestobj = zeros(maxIter,1);
for nIter=1:maxIter
    y = x_now+(t_prev/t_now)*(z-x_now)+((t_prev-1)/t_now)*(x_now-x_prev); %extrapolated point
    Dy = gradMinRate_v2(nAPs,nTx,nUsers,myNu,myBeta,myPsi,zeta_d,Tc,Tp,mygamma,y); %gradient at y
    y_next= y+alphay*Dy; %new y

    Dx = gradMinRate_v2(nAPs,nTx,nUsers,myNu,myBeta,myPsi,zeta_d,Tc,Tp,mygamma,x_now); %gradient at x
    x_next = x_now+alphax*Dx; %new x
    
    z = proj_Euclidean_ball_posorthant(y_next,nTx,nUsers,nAPs); %find z
    v = proj_Euclidean_ball_posorthant(x_next,nTx,nUsers,nAPs); %find v
    
    Fz = minRateUtility_v2(nAPs,nTx,nUsers,myNu,myBeta,myPsi,zeta_d,Tc,Tp,mygamma,z); %objective at z
    Fv = minRateUtility_v2(nAPs,nTx,nUsers,myNu,myBeta,myPsi,zeta_d,Tc,Tp,mygamma,v); %objective at v
    x_prev = x_now; %assing current x to previous x
    if(Fz>Fv)
        x_now = z; %if objective at z is higher then z is the next x
    else
        x_now = v; %if objective at v is higher then v is the next x
    end
    
    bestobj(nIter)=max(Fz,Fv); %the objective is set to the maximum objective among the two
    t_prev = t_now; %set current t to previous t
    t_now = (sqrt(4*t_now^2)+1)/2; %find current t for next iteration
end
SE_FOM(4)= bestobj(end); %final value of minimum rate
plot(1:maxIter,bestobj,'g','DisplayName','Minimum Rate') %plot convergence

%%
SE_FOM