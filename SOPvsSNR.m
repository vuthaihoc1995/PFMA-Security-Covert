clc; clear variables; close all;
%% Network - Fading parameters
OmegaSB = 0.5; OmegaSR = 1;
OmegaST = 0.125;
%%   Transmit power &  Noise
SNRdB = -10:2.5:30;
SNRdB1 = -10:0.5:30;

%% Rate requirement
rB = 0.5; rR = 1;
uB = 2^rB -1; uR = 2.^rR -1;

R_B = 0.25;
theta = 2^R_B - 1;
tau = 2^R_B; 
%% Bandwith + Power
% Power allocation
betaB = 0.7; 
betaR = 1-betaB;
% Bandwith allocation
alphaBR = 0.3;
alphaB = (1-alphaBR)/2;
alphaR = (1-alphaBR)/2;
[rhoB,rhoR,deltaB,deltaR,vB,vR,psiB,psiR] = func_para(alphaB,alphaR,alphaBR,betaB,betaR);

%% Channel realization
trial = 5e5;
hSB = sqrt(OmegaSB/2)*(randn(1,trial) + 1i*randn(1,trial));
hSR  = sqrt(OmegaSR/2)*(randn(1,trial) + 1i*randn(1,trial));
hST = sqrt(OmegaST/2)*(randn(1,trial) + 1i*randn(1,trial));

for idx =1:length(SNRdB)
    %% SINR
    snr = db2pow(SNRdB(idx));
    snr_B   =  rhoB*snr*abs(hSB).^2./(   alphaBR*betaR*snr*abs(hSB).^2 + deltaB );           
    snr_R_xR   =  rhoR*snr*abs(hSR).^2./(   alphaBR*betaB*snr*abs(hSR).^2 + deltaR );
    snr_R_xB   =  alphaBR*betaB*snr*abs(hSR).^2/deltaR;
   %% Internal eavesdropper
             % Simulation
             simSOP_int(idx)         =    mean(  snr_R_xR >   uR   &  (1+ snr_B)./(1 + snr_R_xB)  <    tau );

             % NOMA comparison
             snr_B_xB_noma        =  betaB*snr*abs(hSB).^2./(   betaR*snr*abs(hSB).^2 + 1 ); 
             snr_B_xR_noma        =   betaR*snr*abs(hSB).^2; 
             snr_R_xR_noma        =   betaR*snr*abs(hSR).^2;
            simSOP_noma(idx) =  1- mean(  snr_B_xB_noma > uB     & (1+ snr_R_xR_noma)./(1 + snr_B_xR_noma)  <    tau );
    %% External eavesdropper
            % SNR
            snr_sb =  snr*abs(hST).^2;
            snr_sbr =  betaB*snr*abs(hST).^2;
            snr_sbr0 =  betaB*snr*abs(hST).^2./(   betaR*snr*abs(hST).^2 + 1 );
            snr_sbbr=  rhoB*snr*abs(hST).^2/deltaB;
            snr_sbbr0= rhoB*snr*abs(hST).^2./(   alphaBR*betaR*snr*abs(hST).^2 + deltaB );
            % Simulation
            simSOP_ext(1,idx)            =    mean((1  + snr_B)./(1 + snr_sb)  <    tau );
             simSOP_ext(2,idx)        =    mean((1  + snr_B)./(1 + snr_sbbr)  <    tau );
            simSOP_ext(3,idx)          =    mean((1  + snr_B)./(1 + snr_sbr)  <    tau );
            % NOMA comparison
            simSOP_ext_noma(idx)           =    mean((1  + snr_B_xB_noma)./(1 + snr_sbr)  <    tau );
           
          
end

for idx =1:length(SNRdB1)
    snr = db2pow(SNRdB1(idx));
             % Analysis 
             q_func = @(x)  exp( -  psiB/(OmegaSB*snr).*x./(vB - x )  -    psiR/(  OmegaSR*tau*snr).*x  ) ;    
             CDF_SR = @(x)  1 - exp( -  x./OmegaSR);
              
             phi = uR/(vR-uR); 
             sumQx = 0; K=20; 
             for k = 1:K
                 larg_Phi = exp(psiR*theta/(OmegaSR*tau*snr))*psiR/(OmegaSR*tau*snr);
                 eta_k = cos((2*k-1)*pi/2/K);
                 var_k = eta_k*(vB - tau*phi - theta)/2 + (vB + tau*phi + theta)/2;
                 Lamk = larg_Phi*(vB - tau*phi - theta)/2*pi*sqrt(1 - eta_k^2)/K;
                 sumQx = sumQx + Lamk*q_func(var_k);
             end

             anaSOP_int(idx)         =    (1 - CDF_SR(psiR*phi/snr) - sumQx).*(tau*phi + theta < vB)...
                 +(1 - CDF_SR(psiR*phi/snr)).*(tau*phi + theta >= vB);
             % Asymptotic
             asySOP_int(idx)         =    exp(- 1/snr/OmegaSR*max(psiR*uR/(vR-uR),psiR* (vB-theta)/tau  )  );

            % Analysis
      
            K=20; zeta=1;
            [sumPx] = ext_SOP(K,theta,OmegaST,zeta,tau,snr,vB,psiB,OmegaSB);
            anaSOP_ext(1,idx)            =   (1-sumPx).*(vB>theta ) + 1.*(vB<=theta );

            zeta=betaB;
            [sumPx] = ext_SOP(K,theta,OmegaST,zeta,tau,snr,vB,psiB,OmegaSB);
            anaSOP_ext(2,idx)           =   (1-sumPx).*(vB>theta ) + 1.*(vB<=theta );

            zeta=rhoB/deltaB;
            [sumPx] = ext_SOP(K,theta,OmegaST,zeta,tau,snr,vB,psiB,OmegaSB);
            anaSOP_ext(3,idx)            =   (1-sumPx).*(vB>theta ) + 1.*(vB<=theta );

           % Asymptotic
            asySOP_ext(1,idx)             =    exp(-(vB-theta)/(OmegaST*tau*snr));
            asySOP_ext(2,idx)            =    exp(-(vB-theta)/(OmegaST*betaB*tau*snr));
            asySOP_ext(3,idx)          =    exp(-(vB-theta)/(OmegaST*rhoB/deltaB*tau*snr));

end
%% Plot result
figure(1)
% % Plot
subplot(1,2,1)
sim0 = semilogy(SNRdB, simSOP_int,'ro',  'linewidth',1,'MarkerSize',6,'MarkerFaceColor','r'); hold on;
sim1 = semilogy(SNRdB, simSOP_noma,'bs-',  'linewidth',1,'MarkerSize',6,'MarkerFaceColor','b'); hold on;
asy0 = semilogy(SNRdB1, asySOP_int,'k--',  'linewidth',0.5,'MarkerSize',8.5); hold on;
ana0 = semilogy(SNRdB1, anaSOP_int,'k-',  'linewidth',1,'MarkerSize',8.5); hold on;

% % Grid
lgd1= legend([sim0(1),sim1(1),ana0(1),asy0(1)],...
   'PFMA',...
   'NOMA',...
   'Theory',...
   'Bound');
xlim([min(SNRdB)  max(SNRdB)]);
xlabel('SNR $\overline{\gamma}$ [dB]','FontName','Times New Roman','FontSize',15,'Interpreter','latex');
ylabel('Secrecy Outage Probability','FontName','Times New Roman','FontSize',15);

% % Grid
lgd1.NumColumns = 1;
lgd1.FontSize = 13;
set(gca,'fontsize',14);
xlabel('SNR $\overline{\gamma}$ [dB]','FontSize',15,'Interpreter','latex') ;
ylabel('Secrecy Outage Probability','FontSize',15);
axis([min(SNRdB) max(SNRdB) 1e-3 1]);

% % Plot
subplot(1,2,2)

sim5 = plot(SNRdB, simSOP_ext_noma,'b-s',  'linewidth',1,'MarkerSize',6); hold on; 
asy2 = semilogy(SNRdB1, asySOP_ext,'k--',  'linewidth',0.5,'MarkerSize',8); hold on;
ana2 = semilogy(SNRdB1, anaSOP_ext,'k-',  'linewidth',1,'MarkerSize',8); hold on;
sim2 = plot(SNRdB, simSOP_ext(1,:),'ro',  'linewidth',1,'MarkerSize',6,'MarkerFaceColor','r'); hold on; 
sim3 = plot(SNRdB, simSOP_ext(2,:),'g^',  'linewidth',1,'MarkerSize',6,'MarkerFaceColor','g'); hold on; 
sim4 = plot(SNRdB, simSOP_ext(3,:),'md',  'linewidth',1,'MarkerSize',6,'MarkerFaceColor','m'); hold on;
lgd2 = legend([sim2(1),sim3(1),sim4(1),sim5(1),ana2(1),asy2(1)],...
   'SCE-1',...
   'SCE-2',...
   'SCE-3',...
   'NOMA',...
   'Theory',...
   'Bound');

% % Grid
lgd2.NumColumns = 1;
lgd2.FontSize = 13;
set(gca,'fontsize',14);
xlabel('SNR $\overline{\gamma}$ [dB]','FontSize',15,'Interpreter','latex') 
axis([min(SNRdB) max(SNRdB) 0.25 1]);

%%% Define function
function [sumPx] = ext_SOP(K,theta,OmegaST,zeta,tau,snr,vB,psiB,OmegaSB)
    sumPx =0;
    p_func = @(x,zeta)  exp( -  psiB/(OmegaSB*snr).*x./(vB - x )  -    x/(  OmegaST*zeta*tau*snr)  ) ;
    for k = 1:K
         Denta = exp(theta/(OmegaST*zeta*tau*snr))/(OmegaST*zeta*tau*snr);
         eta_k = cos((2*k-1)*pi/2/K);
         var_k = eta_k*(vB- theta)/2 + (vB + theta)/2;
         sumPx = sumPx + Denta*(vB -theta)/2*pi*sqrt(1 - eta_k^2)/K*p_func(var_k,zeta);
     end
end