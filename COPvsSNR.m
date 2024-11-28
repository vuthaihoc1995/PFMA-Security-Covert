clc; clear variables; close all;
%% Network - Fading parameters
OmegaSB = 0.5; OmegaSR = 1;

%%   Transmit power &  Noise
SNRdB = -5:2.5:30;
snr = db2pow(SNRdB);

%% Rate requirement
rB = 0.5; rR = 1;
uB = 2^rB -1; uR = 2.^rR -1;
%% Bandwith + Power
% power allocation
betaB = 0.7; 
betaR = 1-betaB;


%% Channel realization
trial = 5e5;
hSB = sqrt(OmegaSB/2)*(randn(1,trial) + 1i*randn(1,trial));
hSR  = sqrt(OmegaSR/2)*(randn(1,trial) + 1i*randn(1,trial));

%% COP Evaluation
% Plot (2,1,1)
for idx =1:length(snr)
    alphaBR = 0.1:0.1:0.3;
    for ss = 1:length(alphaBR)
        % bandwidth allocation
        alphaB = (1-alphaBR(ss))/2;
        alphaR = (1-alphaBR(ss))/2;
        [rhoB,rhoR,deltaB,deltaR,vB,vR,psiB,psiR] = func_para(alphaB,alphaR,alphaBR(ss),betaB,betaR);
    % Simulation            
        snr_B   =  rhoB*snr(idx)*abs(hSB).^2./(   alphaBR(ss)*betaR*snr(idx)*abs(hSB).^2 + deltaB );           
        snr_R_xR   =  rhoR*snr(idx)*abs(hSR).^2./(   alphaBR(ss)*betaB*snr(idx)*abs(hSR).^2 + deltaR );
        simCOP(ss,idx) =  1- mean(  snr_B > uB     &  snr_R_xR >   uR);

        
    % Analysis
        F_Bx = 1 - exp(- psiB*uB/(OmegaSB*(vB - uB)*snr(idx) ) );
        F_Rx = 1 - exp(- psiR*uR/(OmegaSR*(vR - uR)*snr(idx) ) );
        ProdX = 1 - (1-F_Bx)*(1-F_Rx);
        anaCOP(ss,idx) =  ProdX*( vB > uB &  vR > uR)  + 1*( vB <= uB ||  vR <= uR);
      % Asymtotic
        asyCOP(ss,idx) =  1/snr(idx) * (psiB*uB/(OmegaSB*(vB - uB)) +  psiR*uR/(OmegaSR*(vR - uR))  )*( vB > uB &  vR > uR)  + 1*( vB <= uB ||  vR <= uR);
    end
    % NOMA comparison
        snr_B_noma   =  betaB*snr(idx)*abs(hSB).^2./(   betaR*snr(idx)*abs(hSB).^2 + 1 );           
        snr_R_xB_noma        =  betaB*snr(idx)*abs(hSR).^2./(   betaR*snr(idx)*abs(hSR).^2 + 1 );
        snr_R_xR_noma        =  betaR*snr(idx)*abs(hSR).^2;

        simCOP_noma(idx) =  1- mean(  snr_B_noma > uB     & snr_R_xB_noma > uB     &  snr_R_xR_noma >   uR);
end
 
%% Plot result
figure(1)
% % Plot
sim1 = semilogy(SNRdB, simCOP,'ro',  'linewidth',1,'MarkerSize',8.5); hold on; 
sim2 = semilogy(SNRdB, simCOP_noma,'b-s',  'linewidth',1,'MarkerSize',9.5); hold on; 
ana1 = semilogy(SNRdB, anaCOP,'k-',  'linewidth',1,'MarkerSize',8.5); hold on; 
asy1 = semilogy(SNRdB, asyCOP,'k--',  'linewidth',1,'MarkerSize',8.5); hold on; 
% % Grid
% % Legend
lgd2=legend([sim1(1),sim2(1),ana1(1),asy1(1)],...
   'PFMA (sim.)',...
   'NOMA  (sim.)',...
   'Theory',...
   'Bound');
% xlabel('SNR $\overline{\gamma}$ [dB]','FontName','Times New Roman','FontSize',15,'Interpreter','latex');
% ylabel('Connection Outage Probability','FontName','Times New Roman','FontSize',15);

lgd2.NumColumns = 1;
lgd2.FontSize = 13;
set(gca,'fontsize',14);
xlabel('SNR $\overline{\gamma}$ [dB]','FontSize',15,'Interpreter','latex') 
ylabel('Connection Outage Probability','FontSize',15) 
axis([min(SNRdB) max(SNRdB) 1e-2 1]);

