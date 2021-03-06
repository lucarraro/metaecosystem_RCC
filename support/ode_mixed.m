function [t,y] = ode_typeI(parameters,OCNparam,tspan,y0)

[alphaR,alphaG,alphaS,alphaC,alphaF,alphaP,...
    betaR,betaG,betaS,betaC,betaF,betaP,...
    muR,muG,muS,muC,muF,muP,epsilonR,epsilonG,epsilonS,epsilonC,epsilonF,epsilonP,...
    deltaMC,deltaMF,deltaMD,deltaN,aF,aN,Beff,kMC,kMD,kN,light,...
MC0,MF0,MD0,N0,lambdaMC,lambdaMF,lambdaMD,lambdaN,zeta]=v2struct(parameters);
[AS,B,L,Q,V,W]=v2struct(OCNparam);

%L(L<1000)=0;
%AS(AS<500000)=0;
Wt = W';

N_reach=length(L);
N_var=10;

ind_R=1:N_var:N_reach*N_var;  % Predators
ind_G=2:N_var:N_reach*N_var;  % GRAZERS
ind_S=3:N_var:N_reach*N_var;  % SHREDDERS
ind_C=4:N_var:N_reach*N_var;  % CONSUMERS
ind_F=5:N_var:N_reach*N_var;  % FILTER FEEDERS
ind_P=6:N_var:N_reach*N_var;  % PRODUCERS
ind_MC=7:N_var:N_reach*N_var;  % CPOM
ind_MF=8:N_var:N_reach*N_var;  % FPOM
ind_MD=9:N_var:N_reach*N_var;  % DOM
ind_N=10:N_var:N_reach*N_var;  % NUTRIENTS

[t,y]=ode23(@(t,y)eqs(t,y),tspan,y0,odeset('RelTol',1e-6,'AbsTol',1e-12));  %'NonNegative',1:(N_reach*N_var),  ) ,odeset('NonNegative',1:(N_reach*6) 

function dy=eqs(t,y)
    
    dy=zeros(N_reach*N_var,1);
    dy(ind_R) = epsilonR*alphaR*y(ind_R).*(y(ind_S) + y(ind_C) + y(ind_G) + y(ind_F)) - muR*y(ind_R) - betaR*y(ind_R).*y(ind_R);
    dy(ind_G) = epsilonG*alphaG*y(ind_G).*y(ind_P) - alphaR*y(ind_R).*y(ind_G) - muG*y(ind_G) - betaG*y(ind_G).*y(ind_G);
    dy(ind_S) = epsilonS*alphaS*y(ind_MC).*y(ind_S) - alphaR*y(ind_R).*y(ind_S) - muS*y(ind_S) - betaS*y(ind_S).*y(ind_S);
    dy(ind_C) = epsilonC*alphaC*y(ind_MF).*y(ind_C) - alphaR*y(ind_R).*y(ind_C) - muC*y(ind_C) - betaC*y(ind_C).*y(ind_C);
    dy(ind_F) = epsilonF*alphaF*y(ind_MD).*y(ind_F) - alphaR*y(ind_R).*y(ind_F) - muF*y(ind_F) - betaF*y(ind_F).*y(ind_F);
    dy(ind_P) = epsilonP*alphaP*light.*y(ind_P).*y(ind_N)./(y(ind_N)+N0) - alphaG*y(ind_P) .*y(ind_G) - muP*y(ind_P) - betaP*y(ind_P).*y(ind_P); 
    dy(ind_MC) = 86400*deltaMC*(Wt*(y(ind_MC).*Q) - y(ind_MC).*Q)./V +...
        kMC*aF.*Beff.*L./V - alphaS*y(ind_MC).*y(ind_S) - lambdaMC*y(ind_MC); 
    dy(ind_MF) = 86400*deltaMF*(Wt*(y(ind_MF).*Q) - y(ind_MF).*Q)./V +...
        (1-epsilonS)*alphaS*y(ind_MC).*y(ind_S) - alphaC*y(ind_MF).*y(ind_C) - lambdaMF*y(ind_MF);
    dy(ind_MD) = 86400*deltaMD*(Wt*(y(ind_MD).*Q) - y(ind_MD).*Q)./V +...
        kMD*AS./V  - alphaF*y(ind_MD).*y(ind_F) - lambdaMD*y(ind_MD);
    dy(ind_N) = 86400*deltaN*(Wt*(y(ind_N).*Q) - y(ind_N).*Q)./V +...
        kN*aN.*AS./V - alphaP*light.*y(ind_N).*y(ind_P)./(y(ind_N)+N0) - lambdaN*y(ind_N);
       
end



end


