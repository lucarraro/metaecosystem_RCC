clear all; close all; clc

addpath('support')
load('support/OCNparam.mat')

scenario='flat';
run_ODE=0; % if =1, run also ODE system; if =0, run only linearized system

z=depth;

N_var=10;
N_reach=length(downNode);
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

W=zeros(N_reach);
Ldown=zeros(N_reach,1);
for i=1:N_reach
    if downNode(i)~=0
        W(i,downNode(i))=1;
    end
    j=i;
    while j~=0
        Ldown(i)=Ldown(i)+L(j);
        j=downNode(j);
    end
end

muR=0.001; muG=0.01; muS=0.01; muC=0.01; muF=0.01; muP=0.01; % net growth rates
epsilonR=0.01; epsilonG=0.1; epsilonS=0.01; epsilonC=0.1; epsilonF=0.1; epsilonP=1; % nondimensional
MC0=1e-2; MF0=1e-2; MD0=1e-3; N0=5e-4; % m^-3
deltaMC=0.01; deltaMF=0.5; deltaMD=1; deltaN=1; % relative velocity of CPOM/FPOM/DOM/N as compared to water
alphaR=1e6; alphaG=2.5e6; alphaS=1e5; alphaC=1e5; alphaF=1e5; alphaP=2.5e5; %
betaR=1e6; betaG=1e6; betaS=1e6; betaC=1e6; betaF=1e6; betaP=1e6; % 
lambdaMC=0.01; lambdaMF=0.01; lambdaMD=0.01; lambdaN=0.01;
zeta=1; % dummy variable controlling input from upstream (if zeta=0, only downstream input but no input from upstream)

Bbar=5;
aF=0.5*ones(N_reach,1); aN=1-aF;

Beff = min(B,Bbar);

kMC = 1000/sum(aF.*Beff.*L); % m^-2 day^-1ind_S
kMD = 1000/sum(AS);         % m^-2 day^-1
kN = 1000/sum(aN.*AS); % m^-2 day^-1

Kd = 1; % m^-1
light = (1 - exp(-Kd*z))./(Kd*z).*(1-exp(-1./aF.*B/Bbar)); 

parameters=v2struct(alphaR,alphaG,alphaS,alphaC,alphaF,alphaP,...
    betaR,betaG,betaS,betaC,betaF,betaP,...
    muR,muG,muS,muC,muF,muP,epsilonR,epsilonG,epsilonS,epsilonC,epsilonF,epsilonP,...
    deltaMC,deltaMF,deltaMD,deltaN,aF,aN,Beff,kMC,kMD,kN,light,...
    MC0,MF0,MD0,N0,lambdaMC,lambdaMF,lambdaMD,lambdaN,zeta);
OCNparam=v2struct(AS,B,L,Q,V,W);

%% Default scenarios

w=warning('off','MATLAB:singularMatrix');
disp(sprintf('Default simulation'))
try load('results/y.mat')
catch
    [X,flag_converge,maxreltol,Niter] = linearized_typeI(parameters,OCNparam);
    if run_ODE
        rng(1)
        y0= X + rand(length(X),1)*1e-14 + 2e-3*(rand(length(X),1)-0.5).*X; y=y0';
        for i=1:100
            tic;
            [t,y_tmp] = ode_typeI(parameters,OCNparam,[0:25:50],y0);
            y=[y; y_tmp(2:end,:)];
            y0=y_tmp(end,:);
            flagsave=1;
            try save('y.mat','y','X')
            catch
                flagsave=0;
            end
            disp(sprintf('%d - Elapsed time: %.2f s  -  Saved: %d',i,toc,flagsave));
        end
    else save('y.mat','X')
    end
end

disp(sprintf('Default (mixed) simulation'))
parameters_mixed=parameters;
parameters_mixed.alphaP=alphaP*(N0 + mean(X(ind_N)));
try load('results/y_mixed.mat')
catch
    [X_mixed,flag_converge,maxreltol,Niter] = linearized_mixed(parameters_mixed,OCNparam);
    if run_ODE
        rng(2)
        y0= X_mixed + rand(length(X_mixed),1)*1e-14 + 2e-3*(rand(length(X_mixed),1)-0.5).*X_mixed; y_mixed=y0';
        for i=1:100
            tic;
            [t,y_tmp] = ode_mixed(parameters_mixed,OCNparam,[0:25:50],y0);
            y_mixed=[y_mixed; y_tmp(2:end,:)];
            y0=y_tmp(end,:);
            flagsave=1;
            try save('results/y_mixed.mat','y_mixed','X_mixed')
            catch
                flagsave=0;
            end
            disp(sprintf('%d - Elapsed time: %.2f s  -  Saved: %d',i,toc,flagsave));
        end
    else save('results/y_mixed.mat','X_mixed')
    end
end


%% Run scenarios

disp(sprintf('noLife simulation'))
try load('results/y_noLife.mat')
catch
    parameters_noLife=parameters;
    parameters_noLife.alphaR=0; parameters_noLife.alphaS=0; parameters_noLife.alphaC=0;
    parameters_noLife.alphaG=0; parameters_noLife.alphaP=0; parameters_noLife.alphaF=0;
    [XnoLife,flag_converge,maxreltol,Niter_noLife] = linearized_typeI(parameters_noLife,OCNparam);
    if run_ODE
        rng(3)
        y0= XnoLife + rand(length(XnoLife),1)*1e-14 + 2e-3*(rand(length(XnoLife),1)-0.5).*X; y_noLife=y0';
        for i=1:100
            tic;
            [t,y_tmp] = ode_typeI(parameters_noLife,OCNparam,[0:25:50],y0);
            y_noLife=[y_noLife; y_tmp(2:end,:)];
            y0=y_tmp(end,:);
            flagsave=1;
            try save('results/y_noLife.mat','y_noLife','XnoLife')
            catch
                flagsave=0;
            end
            disp(sprintf('%d - Elapsed time: %.2f s  -  Saved: %d',i,toc,flagsave));
        end
    else save('results/y_noLife.mat','XnoLife')
    end
end

disp(sprintf('noFlow simulation'))
try load('results/y_noFlow.mat')
catch
    parameters_noFlow=parameters;
    parameters_noFlow.deltaN=0; parameters_noFlow.deltaMD=0; parameters_noFlow.deltaMC=0; parameters_noFlow.deltaMF=0;
    [X_noFlow,flag_converge,maxreltol,Niter_noFlow] = linearized_typeI(parameters_noFlow,OCNparam);
    if run_ODE
        rng(4)
        y0 = X_noFlow + rand(length(X_noFlow),1)*1e-14 + 2e-3*(rand(length(X_noFlow),1)-0.5).*X_noFlow; y_noFlow=y0';
        for i=1:100
            tic;
            [t,y_tmp] = ode_typeI(parameters_noFlow,OCNparam,[0:25:50],y0);
            y_noFlow=[y_noFlow; y_tmp(2:end,:)];
            y0=y_tmp(end,:);
            flagsave=1;
            try save('results/y_noFlow.mat','y_noFlow','X_noFlow')
            catch
                flagsave=0;
            end
            disp(sprintf('%d - Elapsed time: %.2f s  -  Saved: %d',i,toc,flagsave));
        end
    else save('results/y_noFlow.mat','X_noFlow')
    end
end

disp(sprintf('noFlow (mixed) simulation'))
try load('results/y_noFlow_mixed.mat')
catch
    parameters_noFlow_mixed=parameters;
    parameters_noFlow_mixed.deltaN=0; parameters_noFlow_mixed.deltaMD=0; parameters_noFlow_mixed.deltaMC=0; parameters_noFlow_mixed.deltaMF=0;
    parameters_noFlow_mixed.alphaP=alphaP*(N0 + mean(X(ind_N)));
    [X_noFlow_mixed,flag_converge,maxreltol,Niter_noFlow_mixed] = linearized_mixed(parameters_noFlow_mixed,OCNparam);
    if run_ODE
        rng(5)
        y0 = X_noFlow_mixed + rand(length(X_noFlow_mixed),1)*1e-14 + 2e-3*(rand(length(X_noFlow_mixed),1)-0.5).*X_noFlow_mixed; y_noFlow_mixed=y0';
        for i=1:100
            tic;
            [t,y_tmp] = ode_mixed(parameters_noFlow_mixed,OCNparam,[0:25:50],y0);
            y_noFlow_mixed=[y_noFlow_mixed; y_tmp(2:end,:)];
            y0=y_tmp(end,:);
            flagsave=1;
            try save('results/y_noFlow_mixed.mat','y_noFlow_mixed','X_noFlow_mixed')
            catch
                flagsave=0;
            end
            disp(sprintf('%d - Elapsed time: %.2f s  -  Saved: %d',i,toc,flagsave));
        end
    else save('results/y_noFlow_mixed.mat','X_noFlow_mixed')
    end
end

disp(sprintf('fast simulation'))
try load('results/y_fast.mat')
catch
    parameters_fast=parameters;
    parameters_fast.deltaMC=1; parameters_fast.deltaMF=1;
    [X_fast,flag_converge,maxreltol,Niter_fast] = linearized_typeI(parameters_fast,OCNparam);
    if run_ODE
        rng(6)
        y0 = X_fast + rand(length(X_fast),1)*1e-14 + 2e-3*(rand(length(X_fast),1)-0.5).*X_fast; y_fast=y0';
        for i=1:100
            tic;
            [t,y_tmp] = ode_typeI(parameters_fast,OCNparam,[0:25:50],y0);
            y_fast=[y_fast; y_tmp(2:end,:)];
            y0=y_tmp(end,:);
            flagsave=1;
            try save('results/y_fast.mat','y_fast','X_fast')
            catch
                flagsave=0;
            end
            disp(sprintf('%d - Elapsed time: %.2f s  -  Saved: %d',i,toc,flagsave));
        end
    else save('results/y_fast.mat','X_fast')
    end
end

%% save results for export in R
save('support/resultsR.mat','X','XnoLife','X_noFlow','X_fast','parameters')

%% Fig. 3 + Fig. S1
type='Area';
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,4,1); binplot(A,X(ind_R),'k',type,3); title('Predators'); ylabel('Density [m^{-3}]')
text(2e6,0.92e-7,sprintf('%.3e',sum(X(ind_R).*V)))
set(gca,'xlim',[1e6 1e10],'ylim',[0 1e-7],'ytick',[0:5e-8:1e-7]);
subplot(2,4,2); l1=binplot(A,X(ind_G),[0 0.5 0],type,3); title('Preys');
set(gca,'xlim',[1e6 1e10],'ylim',[0 3e-6],'ytick',[0:1e-6:3e-6])
hold on; l2=binplot(A,X(ind_S),[1 0 0],type,3); l3=binplot(A,X(ind_C),[0.5 0.5 0],type,3); l4=binplot(A,X(ind_F),[0 0 1],type,3);
text(2e6,2.8e-6,sprintf('%.3e',sum(X(ind_G).*V)),'Color',[0 0.5 0]); text(2e6,2.7e-6,sprintf('%.3e',sum(X(ind_S).*V)),'Color',[1 0 0]);
text(2e6,2.6e-6,sprintf('%.3e',sum(X(ind_C).*V)),'Color',[0.5 0.5 0]); text(2e6,2.5e-6,sprintf('%.3e',sum(X(ind_F).*V)),'Color',[0 0 1])
legend([l1 l2 l3 l4],'G','S','C','F','location','northeast')
subplot(2,4,3); binplot(A,X(ind_P),'k',type,3); title('Producers');
text(2e6,0.9e-5,sprintf('%.3e',sum(X(ind_P).*V)))
set(gca,'xlim',[1e6 1e10],'ylim',[0 1e-5],'ytick',[0:0.5e-5:1e-5]);
subplot(2,4,4); l1=binplot(A,kMC*aF.*Beff.*L./V,[0.5 0 0.25],type,3); title('Resource input concentrations');
set(gca,'xlim',[1e6 1e10],'ylim',[1e-6 1e-1],'yscale','log');
hold on; l2=binplot(A,kMD*AS./V,[0 0.5 0.25],type,3); l3=binplot(A,kN*aN.*AS./V,[0.25 0 0.5],type,3);
legend([l1 l2 l3],'CPOM','DOM','N','location','northeast')
text(2e6,7e-2,sprintf('%.3e',sum(kMC*aF.*Beff.*L./V.*V)),'Color',[0.5 0 0.25]);
subplot(2,4,5); binplot(A,X(ind_MC),'k',type,3); title('CPOM (Log scale)'); xlabel('Drainage Area [m^2]'); ylabel('Density [m^{-3}]')
set(gca,'xlim',[1e6 1e10],'yscale','log','ylim',[1e-5 1e-2],'ytick',[1e-6 1e-5 1e-4 1e-3 1e-2])
hold on; binplot(A,XnoLife(ind_MC),[1 0.5 0],type,3);
text(2e6,8e-3,sprintf('%.3e',sum(X(ind_MC).*V)));
text(2e6,6e-3,sprintf('%.3e (%.1f%%)',sum(XnoLife(ind_MC).*V),100*sum(XnoLife(ind_MC).*V)/sum(X(ind_MC).*V)),'Color',[1 0.5 0])
plot([1e6 1e10],[MC0 MC0],'--r')
subplot(2,4,6); binplot(A,X(ind_MF),'k',type,3); title('FPOM');
plot([1e6 1e10],[MF0 MF0],'--r')
text(2e6,0.8e-2,sprintf('%.3e',sum(X(ind_MF).*V)))
text(2e6,6e-3,sprintf('%.3e (%.1f%%)',sum(XnoLife(ind_MF).*V),100*sum(XnoLife(ind_MF).*V)/sum(X(ind_MF).*V)),'Color',[1 0.5 0])
set(gca,'xlim',[1e6 1e10],'ylim',[9.99e-6 1e-2],'yscale','log');
subplot(2,4,7); binplot(A,X(ind_MD),'k',type,3); title('DOM (Log scale)'); xlabel('Drainage Area [m^2]');
set(gca,'xlim',[1e6 1e10],'yscale','log','ylim',[1e-5 1e-2],'ytick',[1e-6 1e-5 1e-4 1e-3 1e-2]);
hold on; binplot(A,XnoLife(ind_MD),[1 0.5 0],type,3);
text(2e6,8e-3,sprintf('%.3e',sum(X(ind_MD).*V)));
text(2e6,6e-3,sprintf('%.3e (%.1f%%)',sum(XnoLife(ind_MD).*V),100*sum(XnoLife(ind_MD).*V)/sum(X(ind_MD).*V)),'Color',[1 0.5 0])
plot([1e6 1e10],[MD0 MD0],'--r')
subplot(2,4,8); l1=binplot(A,X(ind_N),'k',type,3); title('Nutrients (Log scale)');
set(gca,'xlim',[1e6 1e10],'yscale','log','ylim',[1e-5 1e-2],'ytick',[1e-6 1e-5 1e-4 1e-3 1e-2]);
hold on;  xlabel('Drainage Area [m^2]');  l2=binplot(A,XnoLife(ind_N),[1 0.5 0],type,3);
l3=plot([1e6 1e10],[N0 N0],'--r');
text(2e6,8e-3,sprintf('%.3e',sum(X(ind_N).*V)));
text(2e6,6e-3,sprintf('%.3e (%.1f%%)',sum(XnoLife(ind_N).*V),100*sum(XnoLife(ind_N).*V)/sum(X(ind_N).*V)),'Color',[1 0.5 0])
legend([l1 l2 l3],'Conc. in water','Conc. w/o consumers','Half-saturation level','location','northeast')


%% Fig. 4 + Fig. S2
type='Area';
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,5,1); binplot(A,X(ind_R),'k',type,3); title('Predators'); ylabel('Density [m^{-3}]')
set(gca,'xlim',[1e6 1e10],'ylim',[1e-9 1e-6],'yscale','log'); hold on;
binplot(A,X_noFlow(ind_R),[0.5 0 0.5],type,3); 
text(2e6,8e-7,sprintf('%.3e',sum(X(ind_R).*V)));
text(2e6,6.25e-7,sprintf('%.2e (%.1f%%)',sum(X_noFlow(ind_R).*V),100*sum(X_noFlow(ind_R).*V)/sum(X(ind_R).*V)),'Color',[0.5 0 0.5]);
subplot(2,5,2); binplot(A,X(ind_S),'k',type,3); title('Shredders'); ylabel('Density [m^{-3}]')
set(gca,'xlim',[1e6 1e10],'ylim',[1e-8 1e-4],'yscale','log'); hold on;
binplot(A,X_noFlow(ind_S),[0.5 0 0.5],type,3); 
text(2e6,8e-5,sprintf('%.3e',sum(X(ind_S).*V)));
text(2e6,5.5e-5,sprintf('%.2e (%.1f%%)',sum(X_noFlow(ind_S).*V),100*sum(X_noFlow(ind_S).*V)/sum(X(ind_S).*V)),'Color',[0.5 0 0.5]);
subplot(2,5,3); binplot(A,X(ind_C),'k',type,3); title('Collectors'); ylabel('Density [m^{-3}]')
set(gca,'xlim',[1e6 1e10],'ylim',[1e-8 1e-4],'yscale','log'); hold on;
binplot(A,X_noFlow(ind_C),[0.5 0 0.5],type,3); 
text(2e6,8e-5,sprintf('%.3e',sum(X(ind_C).*V)));
text(2e6,5.5e-5,sprintf('%.2e (%.1f%%)',sum(X_noFlow(ind_C).*V),100*sum(X_noFlow(ind_C).*V)/sum(X(ind_C).*V)),'Color',[0.5 0 0.5]);
subplot(2,5,4); binplot(A,X(ind_F),'k',type,3); title('Filter feeders'); ylabel('Density [m^{-3}]')
set(gca,'xlim',[1e6 1e10],'ylim',[1e-8 1e-4],'yscale','log'); hold on;
binplot(A,X_noFlow(ind_F),[0.5 0 0.5],type,3); 
text(2e6,8e-5,sprintf('%.3e',sum(X(ind_F).*V)));
text(2e6,5.5e-5,sprintf('%.2e (%.1f%%)',sum(X_noFlow(ind_F).*V),100*sum(X_noFlow(ind_F).*V)/sum(X(ind_F).*V)),'Color',[0.5 0 0.5]);
subplot(2,5,5); l1=binplot(A,X(ind_G),'k',type,3); title('Grazers'); ylabel('Density [m^{-3}]')
set(gca,'xlim',[1e6 1e10],'ylim',[1e-8 1e-4],'yscale','log'); hold on;
l2=binplot(A,X_noFlow(ind_G),[0.5 0 0.5],type,3); 
legend([l1 l2],'w/ Flow','w/o Flow')
text(2e6,8e-5,sprintf('%.3e',sum(X(ind_G).*V)));
text(2e6,5.5e-5,sprintf('%.2e (%.1f%%)',sum(X_noFlow(ind_G).*V),100*sum(X_noFlow(ind_G).*V)/sum(X(ind_G).*V)),'Color',[0.5 0 0.5]);
subplot(2,5,6); binplot(A,X(ind_MC),'k',type,3); title('CPOM'); ylabel('Density [m^{-3}]')
set(gca,'xlim',[1e6 1e10],'ylim',[1e-5 1e-2],'yscale','log'); hold on;
binplot(A,X_noFlow(ind_MC),[0.5 0 0.5],type,3);  
text(2e6,8e-3,sprintf('%.3e',sum(X(ind_MC).*V)));
text(2e6,6.5e-3,sprintf('%.2e (%.1f%%)',sum(X_noFlow(ind_MC).*V),100*sum(X_noFlow(ind_MC).*V)/sum(X(ind_MC).*V)),'Color',[0.5 0 0.5]);
subplot(2,5,7); binplot(A,X(ind_MF),'k',type,3);  title('FPOM'); ylabel('Density [m^{-3}]')
set(gca,'xlim',[1e6 1e10],'ylim',[1e-7 1e-3],'yscale','log'); hold on;
binplot(A,X_noFlow(ind_MF),[0.5 0 0.5],type,3); 
text(2e6,8e-4,sprintf('%.3e',sum(X(ind_MF).*V)));
text(2e6,5.5e-4,sprintf('%.2e (%.1f%%)',sum(X_noFlow(ind_MF).*V),100*sum(X_noFlow(ind_MF).*V)/sum(X(ind_MF).*V)),'Color',[0.5 0 0.5]);
subplot(2,5,8); binplot(A,X(ind_MD),'k',type,3); title('DOM'); ylabel('Density [m^{-3}]')
set(gca,'xlim',[1e6 1e10],'ylim',[9.999e-6 1e-2],'yscale','log'); hold on;
binplot(A,X_noFlow(ind_MD),[0.5 0 0.5],type,3); 
text(1e8,8e-3,sprintf('%.3e',sum(X(ind_MD).*V)));
text(1e8,6.5e-3,sprintf('%.2e (%.1f%%)',sum(X_noFlow(ind_MD).*V),100*sum(X_noFlow(ind_MD).*V)/sum(X(ind_MD).*V)),'Color',[0.5 0 0.5]);
subplot(2,5,9); binplot(A,X(ind_N),'k',type,3); title('Nutrients'); ylabel('Density [m^{-3}]')
set(gca,'xlim',[1e6 1e10],'ylim',[9.999e-6 1e-2],'yscale','log'); hold on;
binplot(A,X_noFlow(ind_N),[0.5 0 0.5],type,3); 
text(1e8,8e-3,sprintf('%.3e',sum(X(ind_N).*V)));
text(1e8,6.5e-3,sprintf('%.2e (%.1f%%)',sum(X_noFlow(ind_N).*V),100*sum(X_noFlow(ind_N).*V)/sum(X(ind_N).*V)),'Color',[0.5 0 0.5]);
subplot(2,5,10); binplot(A,X(ind_P),'k',type,3); title('Producers'); ylabel('Density [m^{-3}]')
set(gca,'xlim',[1e6 1e10],'ylim',[1e-7 1e-4],'yscale','log'); hold on;
binplot(A,X_noFlow(ind_P),[0.5 0 0.5],type,3); 
text(1e8,8e-5,sprintf('%.3e',sum(X(ind_P).*V)));
text(1e8,6.25e-5,sprintf('%.2e (%.1f%%)',sum(X_noFlow(ind_P).*V),100*sum(X_noFlow(ind_P).*V)/sum(X(ind_P).*V)),'Color',[0.5 0 0.5]);


%% Fig. 5
type='Area';
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,5,1); binplot(A,X(ind_R),'k',type,3); title('Predators'); ylabel('Density [m^{-3}]')
set(gca,'xlim',[1e6 1e10],'ylim',[1e-9 1e-6],'yscale','log'); hold on;
binplot(A,X_fast(ind_R),[0.6 0.4 0.1],type,3);
text(2e6,8e-7,sprintf('%.3e',sum(X(ind_R).*V)));
text(2e6,6.25e-7,sprintf('%.3e (%.1f%%)',sum(X_fast(ind_R).*V),100*sum(X_fast(ind_R).*V)/sum(X(ind_R).*V)),'Color',[0.6 0.4 0.1]);
subplot(2,5,2); binplot(A,X(ind_S),'k',type,3); title('Shredders'); ylabel('Density [m^{-3}]')
set(gca,'xlim',[1e6 1e10],'ylim',[1e-8 1e-4],'yscale','log'); hold on;
binplot(A,X_fast(ind_S),[0.6 0.4 0.1],type,3);
text(2e6,8e-5,sprintf('%.3e',sum(X(ind_S).*V)));
text(2e6,5.5e-5,sprintf('%.3e (%.1f%%)',sum(X_fast(ind_S).*V),100*sum(X_fast(ind_S).*V)/sum(X(ind_S).*V)),'Color',[0.6 0.4 0.1]);
subplot(2,5,3); binplot(A,X(ind_C),'k',type,3); title('Collectors'); ylabel('Density [m^{-3}]')
set(gca,'xlim',[1e6 1e10],'ylim',[1e-8 1e-4],'yscale','log'); hold on;
binplot(A,X_fast(ind_C),[0.6 0.4 0.1],type,3);
text(2e6,8e-5,sprintf('%.3e',sum(X(ind_C).*V)));
text(2e6,5.5e-5,sprintf('%.3e (%.1f%%)',sum(X_fast(ind_C).*V),100*sum(X_fast(ind_C).*V)/sum(X(ind_C).*V)),'Color',[0.6 0.4 0.1]);
subplot(2,5,4); binplot(A,X(ind_F),'k',type,3); title('Filter feeders'); ylabel('Density [m^{-3}]')
set(gca,'xlim',[1e6 1e10],'ylim',[1e-8 1e-4],'yscale','log'); hold on;
binplot(A,X_fast(ind_F),[0.6 0.4 0.1],type,3);
text(2e6,8e-5,sprintf('%.3e',sum(X(ind_F).*V)));
text(2e6,5.5e-5,sprintf('%.3e (%.1f%%)',sum(X_fast(ind_F).*V),100*sum(X_fast(ind_F).*V)/sum(X(ind_F).*V)),'Color',[0.6 0.4 0.1]);
subplot(2,5,5); l1=binplot(A,X(ind_G),'k',type,3); title('Grazers'); ylabel('Density [m^{-3}]')
set(gca,'xlim',[1e6 1e10],'ylim',[1e-8 1e-4],'yscale','log'); hold on;
l2=binplot(A,X_fast(ind_G),[0.6 0.4 0.1],type,3);
legend([l1 l2],'slow CPOM, FPOM','fast CPOM, FPOM','location','southeast')
text(2e6,8e-5,sprintf('%.3e',sum(X(ind_G).*V)));
text(2e6,5.5e-5,sprintf('%.3e (%.1f%%)',sum(X_fast(ind_G).*V),100*sum(X_fast(ind_G).*V)/sum(X(ind_G).*V)),'Color',[0.6 0.4 0.1]);
subplot(2,5,6); binplot(A,X(ind_MC),'k',type,3); title('CPOM'); ylabel('Density [m^{-3}]')
set(gca,'xlim',[1e6 1e10],'ylim',[9.99e-6 1e-2],'yscale','log'); hold on;
binplot(A,X_fast(ind_MC),[0.6 0.4 0.1],type,3);   xlabel('Drainage Area [m^2]')
text(2e6,8e-3,sprintf('%.3e',sum(X(ind_MC).*V)));
text(2e6,6.5e-3,sprintf('%.3e (%.1f%%)',sum(X_fast(ind_MC).*V),100*sum(X_fast(ind_MC).*V)/sum(X(ind_MC).*V)),'Color',[0.6 0.4 0.1]);
subplot(2,5,7); binplot(A,X(ind_MF),'k',type,3);  title('FPOM'); ylabel('Density [m^{-3}]')
set(gca,'xlim',[1e6 1e10],'ylim',[1e-7 1e-3],'yscale','log'); hold on;
binplot(A,X_fast(ind_MF),[0.6 0.4 0.1],type,3);   xlabel('Drainage Area [m^2]')
text(2e6,8e-4,sprintf('%.3e',sum(X(ind_MF).*V)));
text(2e6,5.5e-4,sprintf('%.3e (%.1f%%)',sum(X_fast(ind_MF).*V),100*sum(X_fast(ind_MF).*V)/sum(X(ind_MF).*V)),'Color',[0.6 0.4 0.1]);
subplot(2,5,8); binplot(A,X(ind_MD),'k',type,3); title('DOM'); ylabel('Density [m^{-3}]')
set(gca,'xlim',[1e6 1e10],'ylim',[9.999e-6 1e-0],'yscale','log'); hold on;   xlabel('Drainage Area [m^2]')
binplot(A,X_fast(ind_MD),[0.6 0.4 0.1],type,3);
text(1e8,7e-1,sprintf('%.3e',sum(X(ind_MD).*V)));
text(1e8,4e-1,sprintf('%.3e (%.1f%%)',sum(X_fast(ind_MD).*V),100*sum(X_fast(ind_MD).*V)/sum(X(ind_MD).*V)),'Color',[0.6 0.4 0.1]);
subplot(2,5,9); binplot(A,X(ind_N),'k',type,3); title('Nutrients'); ylabel('Density [m^{-3}]')
set(gca,'xlim',[1e6 1e10],'ylim',[9.999e-6 1e-0],'yscale','log'); hold on;
binplot(A,X_fast(ind_N),[0.6 0.4 0.1],type,3);  xlabel('Drainage Area [m^2]')
text(1e8,7e-1,sprintf('%.3e',sum(X(ind_N).*V)));
text(1e8,4e-1,sprintf('%.3e (%.1f%%)',sum(X_noFlow(ind_N).*V),100*sum(X_fast(ind_N).*V)/sum(X(ind_N).*V)),'Color',[0.6 0.4 0.1]);
subplot(2,5,10); binplot(A,X(ind_P),'k',type,3); title('Producers'); ylabel('Density [m^{-3}]')
set(gca,'xlim',[1e6 1e10],'ylim',[1e-7 1e-4],'yscale','log'); hold on;
binplot(A,X_fast(ind_P),[0.6 0.4 0.1],type,3);  xlabel('Drainage Area [m^2]')
text(1e8,8e-5,sprintf('%.3e',sum(X(ind_P).*V)));
text(1e8,6.25e-5,sprintf('%.3e (%.1f%%)',sum(X_fast(ind_P).*V),100*sum(X_fast(ind_P).*V)/sum(X(ind_P).*V)),'Color',[0.6 0.4 0.1]);

%% Fig. S3
type='Area';
figure('units','normalized','outerposition',[0 0 1 0.5])
subplot(1,3,1); binplot(A,X(ind_G),'k',type,3); xlabel('Drainage Area [m^2]'); hold on
set(gca,'xlim',[1e6 1e10],'ylim',[1e-8 1e-4],'yscale','log'); title('Grazers'); ylabel('Density [m^{-3}]')
binplot(A,X_mixed(ind_G),[0 0.5 0.5],type,3);
binplot(A,X_noFlow(ind_G),[0.5 0 0.5],type,3); xlabel('Drainage Area [m^2]'); hold on
binplot(A,X_noFlow_mixed(ind_G),[0.25 0.25 0],type,3);
text(2e6,8e-5,sprintf('%.3e',sum(X(ind_G).*V)));
text(2e6,5.5e-5,sprintf('%.2e (%.1f%%)',sum(X_noFlow(ind_G).*V),100*sum(X_noFlow(ind_G).*V)/sum(X(ind_G).*V)),'Color',[0.5 0 0.5]);
text(2e6,4.2e-5,sprintf('%.2e (%.1f%%)',sum(X_mixed(ind_G).*V),100*sum(X_mixed(ind_G).*V)/sum(X(ind_G).*V)),'Color',[0 0.5 0.5]);
text(2e6,3e-5,sprintf('%.2e (%.1f%%)',sum(X_noFlow_mixed(ind_G).*V),100*sum(X_noFlow_mixed(ind_G).*V)/sum(X(ind_G).*V)),'Color',[0.25 0.25 0]);
subplot(1,3,2); binplot(A,X(ind_P),'k',type,3); xlabel('Drainage Area [m^2]'); hold on
set(gca,'xlim',[1e6 1e10],'ylim',[1e-7 1e-4],'yscale','log'); title('Producers');
binplot(A,X_mixed(ind_P),[0 0.5 0.5],type,3);
binplot(A,X_noFlow(ind_P),[0.5 0 0.5],type,3); xlabel('Drainage Area [m^2]'); hold on
binplot(A,X_noFlow_mixed(ind_P),[0.25 0.25 0],type,3);
text(1e8,8e-5,sprintf('%.3e',sum(X(ind_P).*V)));
text(1e8,5.5e-5,sprintf('%.2e (%.1f%%)',sum(X_noFlow(ind_P).*V),100*sum(X_noFlow(ind_P).*V)/sum(X(ind_P).*V)),'Color',[0.5 0 0.5]);
text(1e8,4.2e-5,sprintf('%.2e (%.1f%%)',sum(X_mixed(ind_P).*V),100*sum(X_mixed(ind_P).*V)/sum(X(ind_P).*V)),'Color',[0 0.5 0.5]);
text(1e8,3e-5,sprintf('%.2e (%.1f%%)',sum(X_noFlow_mixed(ind_P).*V),100*sum(X_noFlow_mixed(ind_P).*V)/sum(X(ind_P).*V)),'Color',[0.25 0.25 0]);
subplot(1,3,3); l1=binplot(A,X(ind_N),'k',type,3); xlabel('Drainage Area [m^2]'); hold on
set(gca,'xlim',[1e6 1e10],'ylim',[9.999e-6 1e-0],'yscale','log'); title('Nutrients');
l2=binplot(A,X_mixed(ind_N),[0 0.5 0.5],type,3);
l3=binplot(A,X_noFlow(ind_N),[0.5 0 0.5],type,3); xlabel('Drainage Area [m^2]'); hold on
l4=binplot(A,X_noFlow_mixed(ind_N),[0.25 0.25 0],type,3);
l5=plot([1e6 1e10],[N0 N0],'--r');
text(1e8,8e-1,sprintf('%.3e',sum(X(ind_N).*V)));
text(1e8,4e-1,sprintf('%.2e (%.1f%%)',sum(X_noFlow(ind_N).*V),100*sum(X_noFlow(ind_N).*V)/sum(X(ind_N).*V)),'Color',[0.5 0 0.5]);
text(1e8,2e-1,sprintf('%.2e (%.1f%%)',sum(X_mixed(ind_N).*V),100*sum(X_mixed(ind_N).*V)/sum(X(ind_N).*V)),'Color',[0 0.5 0.5]);
text(1e8,8e-2,sprintf('%.2e (%.1f%%)',sum(X_noFlow_mixed(ind_N).*V),100*sum(X_noFlow_mixed(ind_N).*V)/sum(X(ind_N).*V)),'Color',[0.25 0.25 0]);
legend([l1 l2 l3 l4 l5],'Default','mixed','no flow','mixed no flow','half-saturation value','location','east')

%% Fig. 2b--e
nodes=[958,1601,1913];
figure('units','centimeters','position',[0 0 8 12]);
subplot(3,1,1); histogram(log10(A),[6:0.25:10]); box off; set(gca,'tickdir','out','xlim',[6 10],'xtick',[6:1:10],'ylim',[0 2000]);
ylabel('Number of nodes'); hold on
for i=1:3; plot(log10(A(nodes(i)) + AS(nodes(i)))*[1 1],[0 2000],'k'); end
subplot(3,1,2); binplot(A,V,[0.5 0.5 0],'Area',3); set(gca,'tickdir','out','yscale','log'); ylabel('Water volume [m^3]')
subplot(3,1,3); binplot(A,light,[0.7 0 0.3],'Area',3);
set(gca,'tickdir','out','ylim',[0 0.6]); xlabel('Drainage Area [m^2]'); ylabel('Light factor [-]')
for i=1:3; plot((A(nodes(i)) + AS(nodes(i)))*[1 1],[0 0.6],'k'); end