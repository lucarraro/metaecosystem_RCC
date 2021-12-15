function [MAT,VEC] = eval_MAT_VEC_seq(MAT,VEC,parameters,OCNparam,X,flag,ind)
%EVALMATVEC Summary of this function goes here
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

N_var=10;

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

VEC(ind_R) = flag(ind_R)*muR;
VEC(ind_G) = flag(ind_G)*muG;
VEC(ind_S) = flag(ind_S)*muS;
VEC(ind_C) = flag(ind_C)*muC;
VEC(ind_F) = flag(ind_F)*muF;
VEC(ind_P) = flag(ind_P)*muP;

MAT(ind_R, ind_G) = flag(ind_R)*epsilonR*alphaR; % effect of G on R
MAT(ind_R, ind_S) = flag(ind_R)*epsilonR*alphaR; % effect of S on R
MAT(ind_R, ind_C) = flag(ind_R)*epsilonR*alphaR; % effect of C on R
MAT(ind_R, ind_F) = flag(ind_R)*epsilonR*alphaR; % effect of F on R
MAT(ind_G, ind_P) = flag(ind_G)*epsilonG*alphaG; % effect of P on G
MAT(ind_G, ind_R) = -flag(ind_G)*alphaR; % effect of R on G
MAT(ind_S, ind_R) = -flag(ind_S)*alphaR; % effect of R on S
MAT(ind_C, ind_R) = -flag(ind_C)*alphaR; % effect of R on C
MAT(ind_F, ind_R) = -flag(ind_F)*alphaR; % effect of R on F
MAT(ind_P, ind_G) = -flag(ind_P)*alphaG; % effect of G on P


MAT(ind_S, ind_MC) = flag(ind_S)*epsilonS*alphaS; % effect of MC on S
MAT(ind_C, ind_MF) = flag(ind_C)*epsilonC*alphaC; % effect of MF on C
MAT(ind_F, ind_MD) = flag(ind_F)*epsilonF*alphaF; % effect of MD on F
MAT(ind_P, ind_N) =  flag(ind_P)*epsilonP*alphaP*light(ind); %effect of N on P

% vary with X
MAT(ind_MC, ind_MC) = -deltaMC*86400*Q(ind)/V(ind) - (deltaMC~=0)*alphaS*X(ind_S) - lambdaMC; % effect of local MC on MC
MAT(ind_MF, ind_MF) = -deltaMF*86400*Q(ind)/V(ind) - (deltaMF~=0)*alphaC*X(ind_C) - lambdaMF; % effect of local MF on MF
MAT(ind_MF, ind_MC) = (1-epsilonS)*alphaS*X(ind_S); % effect of of MC on MF
MAT(ind_MD, ind_MD) = -deltaMD*86400*Q(ind)/V(ind)  - lambdaMD  - (deltaMD~=0)*alphaF*X(ind_F); % effect of local MC on MC
MAT(ind_N, ind_N) = -deltaN*86400*Q(ind)/V(ind) - lambdaN   - (deltaN~=0).*alphaP*light(ind)*X(ind_P); % effect of local N on N

% only active when delta==0 (more stable solution apparently)
MAT(ind_MC, ind_S) = - (deltaMC==0)*alphaS*X(ind_MC);
MAT(ind_MF, ind_C) = - (deltaMF==0)*alphaC*X(ind_MF);
MAT(ind_MD, ind_F) = - (deltaMD==0)*alphaF*X(ind_MD);
MAT(ind_N, ind_P) = - (deltaN==0).*alphaP*light(ind)*X(ind_N);



end

