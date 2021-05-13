function [MAT,VEC] = eval_MAT_VEC(MAT,VEC,parameters,OCNparam,X,flag)
%EVALMATVEC Summary of this function goes here
%   Detailed explanation goes here

[alphaR,alphaG,alphaS,alphaC,alphaF,alphaP,...
    betaR,betaG,betaS,betaC,betaF,betaP,...
    muR,muG,muS,muC,muF,muP,epsilonR,epsilonG,epsilonS,epsilonC,epsilonF,epsilonP,...
    deltaMC,deltaMF,deltaMD,deltaN,aF,aN,Beff,kMC,kMD,kN,light,MC0,MF0,MD0,N0,...
    lambdaMC,lambdaMF,lambdaMD,lambdaN,zeta]=v2struct(parameters);
[AS,B,L,Q,V,W]=v2struct(OCNparam);

N_reach=length(L);
N_var=10;

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

VEC(ind_R) = flag(ind_R)*muR;
VEC(ind_G) = flag(ind_G)*muG;
VEC(ind_S) = flag(ind_S)*muS;
VEC(ind_C) = flag(ind_C)*muC;
VEC(ind_F) = flag(ind_F)*muF;
VEC(ind_P) = flag(ind_P)*muP;


MAT(ind_R + (ind_G-1)*N_reach*N_var) = flag(ind_R)*epsilonR*alphaR; % effect of G on R
MAT(ind_R + (ind_S-1)*N_reach*N_var) = flag(ind_R)*epsilonR*alphaR; % effect of S on R
MAT(ind_R + (ind_C-1)*N_reach*N_var) = flag(ind_R)*epsilonR*alphaR; % effect of C on R
MAT(ind_R + (ind_F-1)*N_reach*N_var) = flag(ind_R)*epsilonR*alphaR; % effect of F on R

MAT(ind_G + (ind_P-1)*N_reach*N_var) = flag(ind_G)*epsilonG*alphaG; % effect of P on G
MAT(ind_G + (ind_R-1)*N_reach*N_var) = -flag(ind_G)*alphaR; % effect of R on G

MAT(ind_S + (ind_R-1)*N_reach*N_var) = -flag(ind_S)*alphaR; % effect of R on S

MAT(ind_C + (ind_R-1)*N_reach*N_var) = -flag(ind_C)*alphaR; % effect of R on C

MAT(ind_F + (ind_R-1)*N_reach*N_var) = -flag(ind_F)*alphaR; % effect of R on F

MAT(ind_P + (ind_G-1)*N_reach*N_var) = -flag(ind_P)*alphaG; % effect of G on P


MAT(ind_S + (ind_MC-1)*N_reach*N_var) = flag(ind_S)*epsilonS*alphaS; % effect of MC on S
MAT(ind_C + (ind_MF-1)*N_reach*N_var) = flag(ind_C)*epsilonC*alphaC; % effect of MF on C
MAT(ind_F + (ind_MD-1)*N_reach*N_var) = flag(ind_F)*epsilonF*alphaF; % effect of MD on F


% vary with X
MAT(ind_MC + (ind_MC-1)*N_reach*N_var) = -deltaMC*86400*Q./V - (deltaMC~=0)*alphaS*X(ind_S) - lambdaMC; % effect of local MC on MC
MAT(ind_MF + (ind_MF-1)*N_reach*N_var) = -deltaMF*86400*Q./V - (deltaMF~=0)*alphaC*X(ind_C) - lambdaMF; % effect of local MF on MF
MAT(ind_MF + (ind_MC-1)*N_reach*N_var) = (1-epsilonS)*alphaS*X(ind_S); % effect of of MC on MF
MAT(ind_MD + (ind_MD-1)*N_reach*N_var) = -deltaMD*86400*Q./V  - lambdaMD  - (deltaMD~=0)*alphaF*X(ind_F); % effect of local MC on MC
MAT(ind_N + (ind_N-1)*N_reach*N_var) = -deltaN*86400*Q./V - lambdaN   - (deltaN~=0).*alphaP.*light.*X(ind_P); % effect of local N on N

MAT(ind_P + (ind_N-1)*N_reach*N_var) =  flag(ind_P)*epsilonP*alphaP.*light./(X(ind_N) + N0); %effect of N on P

% only active when delta==0 (more stable solution apparently)
MAT(ind_MC + (ind_S-1)*N_reach*N_var) = - (deltaMC==0)*alphaS*X(ind_MC);
MAT(ind_MF + (ind_C-1)*N_reach*N_var) = - (deltaMF==0)*alphaC*X(ind_MF);
MAT(ind_MD + (ind_F-1)*N_reach*N_var) = - (deltaMD==0)*alphaF*X(ind_MD);
MAT(ind_N + (ind_P-1)*N_reach*N_var) = - (deltaN==0).*alphaP.*light.*X(ind_N)./(X(ind_N) + N0);



end

