%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paper title: Utility maximization for large-scale cell-free massive MIMO downlink
% Journal: IEEE Transactions on Communications
% Authors: Muhammad Farooq, Hien Quoc Ngo, and Le Nam Tran
% Written by: Muhammad Farooq
% Email: Muhammad.Farooq@ucdconnect.ie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out] = proj_Euclidean_ball_posorthant(x,nTx,nUsers,nAPs)
%Function to calculate projection onto Euclidean ball
out = zeros(nUsers*nAPs,1);
for iAP=1:nAPs
    out((iAP-1)*nUsers+1:iAP*nUsers) = 1/sqrt(nTx)/max(norm(max(x((iAP-1)*nUsers+1:iAP*nUsers),0)),...
        1/sqrt(nTx))*max(x((iAP-1)*nUsers+1:iAP*nUsers),0);
end
end

