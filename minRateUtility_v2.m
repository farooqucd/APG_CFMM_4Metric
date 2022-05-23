%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paper title: Utility maximization for large-scale cell-free massive MIMO downlink
% Journal: IEEE Transactions on Communications
% Authors: Muhammad Farooq, Hien Quoc Ngo, and Le Nam Tran
% Written by: Muhammad Farooq
% Email: Muhammad.Farooq@ucdconnect.ie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [minRate,minRateApprox,sumRate,userRate] = minRateUtility_v2(nAPs,nTx,nUsers,myNu,myBeta,myPsi,zeta_d,Tc,Tp,mygamma,x)
%minRateUtility_v2: To calculate the minimum rate from a power control vector
userRate = zeros(nUsers,1);
x = reshape(x,nUsers,[])';
for iUser =1:nUsers
    %signal
    sig = (zeta_d)*(sqrt(myNu(:,iUser))'*x(:,iUser))^2; 
    interference = computeInterference(nAPs,nTx,nUsers,x,myNu,myBeta,myPsi,iUser,zeta_d);
    %interference
    userRate(iUser) = (1-Tp/Tc)*log2(1+sig/(interference+1/nTx^2));
end
sumRate = sum(userRate); %sum rate
minRate = min(userRate); %minimum rate
minRateApprox = -log(sum(exp(-mygamma*userRate)))/mygamma; %approximate minimum rate
end

