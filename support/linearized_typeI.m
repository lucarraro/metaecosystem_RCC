function [Xfinal,flag_converge,maxreltol,Niter,maxabstol,Xmat] = linearized_typeI(parameters,OCNparam)
%LINEARIZED_RSGC Summary of this function goes here
%   Detailed explanation goes here

[alphaR,alphaG,alphaS,alphaC,alphaF,alphaP,...
    betaR,betaG,betaS,betaC,betaF,betaP,...
    muR,muG,muS,muC,muF,muP,epsilonR,epsilonG,epsilonS,epsilonC,epsilonF,epsilonP,...
    deltaMC,deltaMF,deltaMD,deltaN,aF,aN,Beff,kMC,kMD,kN,light,MC0,MF0,MD0,N0,...
    lambdaMC,lambdaMF,lambdaMD,lambdaN,zeta]=v2struct(parameters);
[AS,B,L,Q,V,W]=v2struct(OCNparam);

N_reach=length(L);
N_var=10;
N_maxIter=1e5;

ind_R=1:N_var:N_reach*N_var;  % PREDATORS
ind_G=2:N_var:N_reach*N_var;  % GRAZERS
ind_S=3:N_var:N_reach*N_var;  % SHREDDERS
ind_C=4:N_var:N_reach*N_var;  % COLLECTORS
ind_F=5:N_var:N_reach*N_var;  % FILTER FEEDERS
ind_P=6:N_var:N_reach*N_var;  % PRODUCERS
ind_MC=7:N_var:N_reach*N_var;  % CPOM
ind_MF=8:N_var:N_reach*N_var;  % FPOM
ind_MD=9:N_var:N_reach*N_var;  % DOM
ind_N=10:N_var:N_reach*N_var;  % NUTRIENTS


X=1e-7*ones(N_reach*N_var,1);
Xmat=ones(N_reach*N_var,1000);

TransportMatrix=86400*W'.*repmat(Q',N_reach,1); % only upstream transport
TransportMatrix=sparse(TransportMatrix);

phiMC=kMC*aF.*Beff.*L./V;
phiMD=kMD*AS./V;
phiN=kN*aN.*AS./V;

VEC=sparse(N_reach*N_var,1);
MAT=sparse(N_reach*N_var,N_reach*N_var);

VEC(ind_MC) = -phiMC;
VEC(ind_MF) = 0;
VEC(ind_MD) = -phiMD;
VEC(ind_N) = -phiN;

MAT(ind_R + (ind_R-1)*N_reach*N_var) = -betaR; % effect of R on R
MAT(ind_G + (ind_G-1)*N_reach*N_var) = -betaG; % effect of G on G
MAT(ind_S + (ind_S-1)*N_reach*N_var) = -betaS; % effect of S on S
MAT(ind_C + (ind_C-1)*N_reach*N_var) = -betaC; % effect of C on C
MAT(ind_F + (ind_F-1)*N_reach*N_var) = -betaF; % effect of F on F
MAT(ind_P + (ind_P-1)*N_reach*N_var) = -betaP; % effect of P on P

MAT(ind_MC,ind_MC) = zeta*deltaMC./V.*TransportMatrix; % effect of upstream MC on MC
MAT(ind_MF,ind_MF) = zeta*deltaMF./V.*TransportMatrix; % effect of upstream MF on MF
MAT(ind_MD,ind_MD) = zeta*deltaMD./V.*TransportMatrix; % effect of upstream MD on MD
MAT(ind_N,ind_N)   = zeta*deltaN./V.*TransportMatrix; % effect of upstream N on N

iter=1; maxabstol=1e-3; maxreltol=1;
while  iter <= N_maxIter && maxreltol > 1e-6 %abstol > 1e-10
    Xold=X;
    Xmat(:,iter)=Xold;
    flag=ones(N_reach*N_var,1);
    [MAT,VEC] = eval_MAT_VEC(MAT,VEC,parameters,OCNparam,X,flag);
    
    xnew = full(MAT\VEC);
    flag(xnew<0)=0;
    
    if sum(flag==0)>0
        [MAT,VEC] = eval_MAT_VEC(MAT,VEC,parameters,OCNparam,X,flag);
        X = full(MAT\VEC);
    else
        X=xnew;
    end
    
    %xnew(xnew<0)=0;
    %xnew(xnew==Inf)=0; xnew(isnan(xnew))=0;
    %X=xnew;
    
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
if iter <= N_maxIter
    flag_converge=1;
else
    flag_converge=0;
end
Xfinal=X;
Niter=iter-1;
end

