%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paper title: Utility maximization for large-scale cell-free massive MIMO downlink
% Journal: IEEE Transactions on Communications
% Authors: Muhammad Farooq, Hien Quoc Ngo, and Le Nam Tran
% Written by: Muhammad Farooq
% Email: Muhammad.Farooq@ucdconnect.ie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out] = gradFairnessRate(nAPs,nTx,nUsers,myNu,myBeta,myPsi,zeta_d,Tc,Tp,x)
%gradFairnessRate: Function to calculate the gradient of the fariness rate
x_new = reshape(x,nUsers,[])';
out = zeros(nAPs*nUsers,1);
[fairnessRate,sumRate,userRate] = propFairnessRate(nAPs,nTx,nUsers,myNu,myBeta,myPsi,zeta_d,Tc,Tp,x);
for iUser =1:nUsers
    sig = (zeta_d)*((sqrt(myNu(:,iUser))'*x(iUser:nUsers:end)))^2; 
    
    % gradient of the signal part
    gradsig=zeros(nAPs*nUsers,1);
    gradsig(iUser:nUsers:end)=sqrt(myNu(:,iUser));
    gradsig=gradsig*(sqrt(myNu(:,iUser))'*x(iUser:nUsers:end));
   
    % compute the gradient of the interference
    gradinterference = zeros(nAPs*nUsers,1);
    for jUser = 1:nUsers
        % the first sum in the intererence term
         
        gradinterference(jUser:nUsers:end) =gradinterference(jUser:nUsers:end)+ (myBeta(:,iUser)).*...
            x(jUser:nUsers:end)/nTx;
       if(jUser~=iUser)
            mygammatilde = abs((myPsi(:,iUser)'*myPsi(:,jUser)))...
                *sqrt(myNu(:,jUser))./myBeta(:,jUser)...
                .*myBeta(:,iUser);            
            
            gradinterference(jUser:nUsers:end)=gradinterference(jUser:nUsers:end)+...
                mygammatilde*(mygammatilde'*x(jUser:nUsers:end));
       end        
    end
    %interference
    interference=computeInterference(nAPs,nTx,nUsers,x_new,myNu,myBeta,myPsi,iUser,zeta_d);
    mu = sig+interference+1/nTx^2;
    mubar = interference+1/nTx^2;
    
    gradSEk =  2*(1-Tp/Tc)*zeta_d*log2(exp(1))*(gradsig+gradinterference)/mu-2*zeta_d*log2(exp(1))*(1-Tp/Tc)*gradinterference/mubar;
    gradFairness = (1/nUsers)*gradSEk*fairnessRate/userRate(iUser);
    out = out + gradFairness;
end
end

