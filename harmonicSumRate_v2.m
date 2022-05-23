%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paper title: Utility maximization for large-scale cell-free massive MIMO downlink
% Journal: IEEE Transactions on Communications
% Authors: Muhammad Farooq, Hien Quoc Ngo, and Le Nam Tran
% Written by: Muhammad Farooq
% Email: Muhammad.Farooq@ucdconnect.ie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [harmonicRate,sumRate,userRate] = harmonicSumRate_v2(nAPs,nTx,nUsers,myNu,myBeta,myPsi,zeta_d,Tc,Tp,x)
%harmonicSumRate_v2: To calculate the harmonic rate from a power control vector
userRate = zeros(nUsers,1);
x = reshape(x,nUsers,[])';
for iUser =1:nUsers
    %signal
    sig = (zeta_d)*(sqrt(myNu(:,iUser))'*x(:,iUser))^2; 
    interference = computeInterference(nAPs,nTx,nUsers,x,myNu,myBeta,myPsi,iUser,zeta_d);
    %interference
    userRate(iUser) = (1-Tp/Tc)*log2(1+sig/(interference+1/nTx^2));
end
sumRate = sum(userRate);
harmonicRate = nUsers/(sum(1./userRate));
end

