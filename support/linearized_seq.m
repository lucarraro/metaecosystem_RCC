function [sol] = linearized_seq(parameters,OCNparam)
%LINEARIZED_RSGC Summary of this function goes here
%   Detailed explanation goes here

[alphaR,alphaG,alphaS,alphaC,alphaF,alphaP,...
    betaR,betaG,betaS,betaC,betaF,betaP,...
    muR,muG,muS,muC,muF,muP,epsilonR,epsilonG,epsilonS,epsilonC,epsilonF,epsilonP,...
    deltaMC,deltaMF,deltaMD,deltaN,aF,aN,kMD,kN,MC0,MF0,MD0,N0,...
    lambdaMC,lambdaMF,lambdaMD,lambdaN,zeta,Kd,Bbar]=v2struct(parameters);
[AS,B,L,Q,V,W,downNode,A,z]=v2struct(OCNparam);


% recalculate
Beff = min(B,Bbar);
kMC = 1000/sum(aF.*Beff.*L); % m^-2 day^-1ind_S
light = (1 - exp(-Kd*z))./(Kd*z).*(1-exp(-1./aF.*B/Bbar)); 


N_reach=length(L);
N_var=10;
N_maxIter=1e5;

ind_R=1;  % PREDATORS
ind_G=2;  % GRAZERS
ind_S=3;  % SHREDDERS
ind_C=4;  % COLLECTORS
ind_F=5;  % FILTER FEEDERS
ind_P=6;  % PRODUCERS
ind_MC=7;  % CPOM
ind_MF=8;  % FPOM
ind_MD=9;  % DOM
ind_N=10;  % NUTRIENTS

sol=zeros(N_var,N_reach);

%TransportMatrix=86400*W'.*repmat(Q',N_reach,1); % only upstream transport
%TransportMatrix=sparse(TransportMatrix);
[aa,ii]=sort(A); % B increases monotonically as A
for ind=ii' % solve nodes by increasing A (so that a node is solved after its upstream nodes have been already solved)
    
    phiMC=kMC*aF(ind)*Beff(ind)*L(ind)/V(ind);
    phiMD=kMD*AS(ind)/V(ind);
    phiN=kN*aN(ind)*AS(ind)/V(ind);
    
    VEC=sparse(N_var,1);
    MAT=sparse(N_var,N_var);
    
    ups_ind=find(downNode==ind);
    
    if length(ups_ind)>0
        VEC(ind_MC) = -phiMC - sol(ind_MC,ups_ind)*deltaMC*86400*Q(ups_ind)/V(ind);
        VEC(ind_MF) = -sol(ind_MF,ups_ind)*deltaMF*86400*Q(ups_ind)/V(ind);
        VEC(ind_MD) = -phiMD - sol(ind_MD,ups_ind)*deltaMD*86400*Q(ups_ind)/V(ind);
        VEC(ind_N) = -phiN - sol(ind_N,ups_ind)*deltaN*86400*Q(ups_ind)/V(ind);
    else
        VEC(ind_MC) = -phiMC;
        VEC(ind_MF) = 0;
        VEC(ind_MD) = -phiMD;
        VEC(ind_N) = -phiN;
    end
    
    MAT(ind_R,ind_R) = -betaR; % effect of R on R
    MAT(ind_G,ind_G) = -betaG; % effect of G on G
    MAT(ind_S,ind_S) = -betaS; % effect of S on S
    MAT(ind_C,ind_C) = -betaC; % effect of C on C
    MAT(ind_F,ind_F) = -betaF; % effect of F on F
    MAT(ind_P,ind_P) = -betaP; % effect of P on P
    
    
    X=1e-7*ones(N_var,1); Xmat=ones(N_var,1000);
    iter=1; maxabstol=1e-3; maxreltol=1;
    while  iter <= N_maxIter && maxreltol > 1e-6 %abstol > 1e-10
        Xold=X;
        Xmat(:,iter)=Xold;
        flag=ones(N_var,1);
        [MAT,VEC] = eval_MAT_VEC_seq(MAT,VEC,parameters,OCNparam,X,flag,ind);
        
        xnew = full(MAT\VEC);
        flag(xnew<0)=0;
        
        while sum(xnew<0)>0
            [MAT,VEC] = eval_MAT_VEC_seq(MAT,VEC,parameters,OCNparam,X,flag,ind);
            xnew = full(MAT\VEC);
            flag=flag.*(xnew>=0);
        end
        X=xnew;
 
        if iter > 1
            abstol=(abs(X-Xold));
            abstol(abstol==Inf)=0; abstol(abstol==-Inf)=0; abstol(isnan(abstol))=NaN;
            maxabstol=max(abstol);
            reltol=abs(X-Xold)./X;
            reltol(reltol==Inf)=0; reltol(reltol==-Inf)=0; reltol(isnan(reltol))=NaN;
            maxreltol=max(reltol);
        end
        iter=iter+1;
        if mod(iter,100)==0 %&& iter > 2
            pos_abs=find(abstol==max(abstol)); pos_abs=pos_abs(1);
            pos_rel=find(reltol==max(reltol)); pos_rel=pos_rel(1);
            disp(sprintf('Iter = %d  -  MaxAbstol = %.2e  at node %d  -  MaxReltol = %.2e at node %d',...
                iter,maxabstol,pos_abs,maxreltol,pos_rel))
            %disp(sprintf('Iter = %d  -  MaxReltol = %.2e at node %d',iter,maxreltol,pos))
        end
    end
%     if iter <= N_maxIter
%         flag_converge=1;
%     else
%         flag_converge=0;
%     end
    sol(:,ind)=X;
    Niter=iter-1;
    disp(sprintf('node: %d  -  Niter: %d',ind,Niter))
end

