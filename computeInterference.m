%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paper title: Utility maximization for large-scale cell-free massive MIMO downlink
% Journal: IEEE Transactions on Communications
% Authors: Muhammad Farooq, Hien Quoc Ngo, and Le Nam Tran
% Written by: Muhammad Farooq
% Email: Muhammad.Farooq@ucdconnect.ie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = computeInterference(nAPs,nTx,nUsers,x,myNu,myBeta,myPsi,iUser,zeta_d)
out=0;
inteferencevector =zeros(nUsers,1);
for jUser =1:nUsers
    if jUser~=iUser
        gamma_tilde = abs(myPsi(:,iUser)'*myPsi(:,jUser))*...
            sqrt(myNu(:,jUser))./myBeta(:,jUser).*((myBeta(:,iUser)));
        inteferencevector(jUser) = x(:,jUser)'*gamma_tilde;
    end
end
out = out+(zeta_d)*norm(inteferencevector)^2;
term2=x.*(repmat(sqrt(myBeta(:,iUser)/nTx),1,nUsers));

out=out+(zeta_d)*norm(term2,'fro')^2;

end

