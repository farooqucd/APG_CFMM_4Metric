%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paper title: Utility maximization for large-scale cell-free massive MIMO downlink
% Journal: IEEE Transactions on Communications
% Authors: Muhammad Farooq, Hien Quoc Ngo, and Le Nam Tran
% Written by: Muhammad Farooq
% Email: Muhammad.Farooq@ucdconnect.ie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out,SR] = computeSumRate(nAPs,nTx,nUsers,myNu,myBeta,myPsi,zeta_d,Tc,Tp,x)
%computeSumRate: To calculate the sum rate from a power control vector
SR = zeros(nUsers,1);
x = reshape(x,nUsers,[])';
for iUser =1:nUsers
    sig = (zeta_d)*(sqrt(myNu(:,iUser))'*x(:,iUser))^2; 
    interference = computeInterference(nAPs,nTx,nUsers,x,myNu,myBeta,myPsi,iUser,zeta_d);
    SR(iUser) = (1-Tp/Tc)*log2(1+sig/(interference+1/nTx^2));
end
out = sum(SR);
end

