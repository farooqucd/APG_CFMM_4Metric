%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paper title: Utility maximization for large-scale cell-free massive MIMO downlink
% Journal: IEEE Transactions on Communications
% Authors: Muhammad Farooq, Hien Quoc Ngo, and Le Nam Tran
% Written by: Muhammad Farooq
% Email: Muhammad.Farooq@ucdconnect.ie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fairnessRate,sumRate,userRate] = propFairnessRate(nAPs,nTx,nUsers,myNu,myBeta,myPsi,zeta_d,tau_c,tau_p,x)
%propFairnessRate: To calculate the proportional fairness rate from a power control vector
userRate = zeros(nUsers,1);
x = reshape(x,nUsers,[])';
for iUser =1:nUsers
    %signal
    sig = (zeta_d)*(sqrt(myNu(:,iUser))'*x(:,iUser))^2; 
    interference = computeInterference(nAPs,nTx,nUsers,x,myNu,myBeta,myPsi,iUser,zeta_d);
    %interference
    userRate(iUser) = (1-tau_p/tau_c)*log2(1+sig/(interference+1/nTx^2));
end
sumRate = sum(userRate);
fairnessRate = exp((1/nUsers)*sum(log(userRate)));
end

