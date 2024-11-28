clc; clear variables; close all;
%% Network - Fading parameters
OmegaST = 0.125;
%%   Transmit power &  Noise
SNRdB = 0:5:40;
SNRdB1 = 0:1:40;
N0=1;
%% Threshold
lambda = db2pow(1);
%% Bandwith + Power
% Power allocation
betaB = 0.7;  betaR = 1-betaB;  
% Bandwith allocation
alphaBR = 0.3;
alphaR = (1-alphaBR)/2; 
deltaR = alphaR + alphaBR;

%% Channel realization
trial = 1e4;
hST = sqrt(OmegaST/2)*(randn(1,trial) + 1i*randn(1,trial));
nSTsR  = sqrt(alphaR*N0/2)*(randn(1,trial) + 1i*randn(1,trial));
nSTsBR  = sqrt(alphaBR*N0/2)*(randn(1,trial) + 1i*randn(1,trial));
nSTsRpBR  = sqrt(deltaR*N0/2)*(randn(1,trial) + 1i*randn(1,trial));

nSTnoma  = sqrt(N0/2)*(randn(1,trial) + 1i*randn(1,trial));

%% DEP Evaluation

%% Simulation - Binary Hypothesis
for idx =1:length(SNRdB)
    snr = db2pow(SNRdB(idx));
    % scenario 1
    ProFA1 = mean( abs(nSTsR).^2 > lambda);
    ProMD1 = mean( alphaR*snr*abs(hST).^2 + abs(nSTsR).^2 < lambda);
    simDEPsR(idx) = ProFA1 + ProMD1;
    % scenario 2
    ProFA2 = mean( alphaBR*betaB*snr*abs(hST).^2 + abs(nSTsBR).^2 > lambda);
    ProMD2 = mean( alphaBR*snr*abs(hST).^2 + abs(nSTsBR).^2 < lambda);
    simDEPsBR(idx) = ProFA2 + ProMD2;

    % scenario 3
    ProFA3  = mean( alphaBR*betaB*snr*abs(hST).^2 + abs(nSTsRpBR).^2 > lambda);
    ProMD3 = mean( deltaR*snr*abs(hST).^2 + abs(nSTsRpBR).^2 < lambda);
    simDEPsRpBR(idx) = ProFA3 + ProMD3;

    % NOMA comparision
    ProFAnoma = mean( betaB*snr*abs(hST).^2 + abs(nSTnoma).^2 > lambda);
    ProMDnoma = mean( snr*abs(hST).^2 + abs(nSTnoma).^2 < lambda);
    simDEPnoma(idx) = ProFAnoma + ProMDnoma;
end
 %% Analytical - Binary Hypothesis
for idx =1:length(SNRdB1)
    snr = db2pow(SNRdB1(idx));   


   func_Xi = @(z,a,OmegaX,OmegaY) OmegaX*a/(OmegaX*a - OmegaY)*( exp(- z/( OmegaX*a))  - exp(- z/OmegaY));

   % SCENARIO 1 
   anaDEPsR(idx) =  1 - func_Xi(lambda,alphaR*N0*snr,OmegaST,alphaR*N0);

    % SCENARIO 2
    anaDEPsBR(idx) =  1+ func_Xi(lambda,alphaBR*betaB*N0*snr,OmegaST,alphaBR*N0) - func_Xi(lambda,alphaBR*N0*snr,OmegaST,alphaBR*N0);
     
    % SCENARIO 3 
    anaDEPsRpBR(idx) = 1+ func_Xi(lambda,alphaBR*betaB*N0*snr,OmegaST,(alphaBR+alphaR)*N0) - func_Xi(lambda,(alphaBR+alphaR)*snr,OmegaST,(alphaBR+alphaR)*N0);

end

%% Optimal solution

snr_starsR =max(snr);
aptDEPsR =  1 - func_Xi(lambda,alphaR*N0*snr_starsR,OmegaST,alphaR*N0);

snr_starsBR = lambda/(OmegaST*alphaBR*betaB*N0)*(1-betaB)/log(1/betaB);
aptDEPsBR =  1+ func_Xi(lambda,alphaBR*betaB*N0*snr_starsBR,OmegaST,alphaBR*N0) - func_Xi(lambda,alphaBR*N0*snr_starsBR,OmegaST,alphaBR*N0);

snr_starsRpBR = lambda*(1/(alphaBR*betaB)-1/deltaR)/(OmegaST*N0*log(deltaR/(alphaBR*betaB)));
aptDEPsRpBR =  1+ func_Xi(lambda,alphaBR*betaB*N0*snr_starsRpBR,OmegaST,(alphaBR+alphaR)*N0) - func_Xi(lambda,(alphaBR+alphaR)*snr_starsRpBR,OmegaST,(alphaBR+alphaR)*N0);

%% Save result
simDEPvsSNR = [simDEPsR;simDEPsBR;simDEPsRpBR ]; 
snrdBsim = SNRdB;
anaDEPvsSNR = [anaDEPsR;anaDEPsBR;anaDEPsRpBR ]; 
snrdBana = SNRdB1;
simDEPnomavsSNR = simDEPnoma;

optSNRdB = 10*log10([ snr_starsR, snr_starsBR,  snr_starsRpBR]);
optDEPvsSNR= [aptDEPsR, aptDEPsBR,aptDEPsRpBR];

%% Next-Plot
snr = db2pow(20);  lambdadB = -20:5:20;  lambdadB1 = -20:0.5:20; 
simDEPsR=[];simDEPsBR=[];simDEPsRpBR=[];simDEPnoma=[];
anaDEPsR=[];anaDEPsBR=[];anaDEPsRpBR=[];
%% DEP Evaluation
for idx =1:length(lambdadB)
    lambda = db2pow(lambdadB(idx));
    %% Simulation - Binary Hypothesis
    % scenario 1
    ProFA1 = mean( abs(nSTsR).^2 > lambda);
    ProMD1 = mean( alphaR*snr*abs(hST).^2 + abs(nSTsR).^2 < lambda);
    simDEPsR(idx) = ProFA1 + ProMD1;
    % scenario 2
    ProFA2 = mean( alphaBR*betaB*snr*abs(hST).^2 + abs(nSTsBR).^2 > lambda);
    ProMD2 = mean( alphaBR*snr*abs(hST).^2 + abs(nSTsBR).^2 < lambda);
    simDEPsBR(idx) = ProFA2 + ProMD2;

    % scenario 3
    ProFA3  = mean( alphaBR*betaB*snr*abs(hST).^2 + abs(nSTsRpBR).^2 > lambda);
    ProMD3 = mean( deltaR*snr*abs(hST).^2 + abs(nSTsRpBR).^2 < lambda);
    simDEPsRpBR(idx) = ProFA3 + ProMD3;

    % NOMA comparision
    ProFAnoma = mean( betaB*snr*abs(hST).^2 + abs(nSTnoma).^2 > lambda);
    ProMDnoma = mean( snr*abs(hST).^2 + abs(nSTnoma).^2 < lambda);
    simDEPnoma(idx) = ProFAnoma + ProMDnoma;
end
for idx =1:length(lambdadB1)
    lambda = db2pow(lambdadB1(idx));    
   %% Analytical - Binary Hypothesis
   func_Xi = @(z,a,OmegaX,OmegaY) OmegaX*a/(OmegaX*a - OmegaY)*( exp(- z/( OmegaX*a))  - exp(- z/OmegaY));

   % SCENARIO 1 
   anaDEPsR(idx) =  1 - func_Xi(lambda,alphaR*N0*snr,OmegaST,alphaR*N0);

    % SCENARIO 2
    anaDEPsBR(idx) =  1+ func_Xi(lambda,alphaBR*betaB*N0*snr,OmegaST,alphaBR*N0) - func_Xi(lambda,alphaBR*N0*snr,OmegaST,alphaBR*N0);
     
    % SCENARIO 3 
    anaDEPsRpBR(idx) = 1+ func_Xi(lambda,alphaBR*betaB*N0*snr,OmegaST,(alphaBR+alphaR)*N0) - func_Xi(lambda,(alphaBR+alphaR)*snr,OmegaST,(alphaBR+alphaR)*N0);

end
%% Optimal solution
lambda_starsR = OmegaST*alphaR*N0*snr*log(OmegaST*snr)/(OmegaST*snr - 1);
optDEPvslambdasR = 1-OmegaST*snr/(OmegaST*snr - 1)*( exp(- lambda_starsR/( OmegaST*alphaR*N0*snr))  - exp(- lambda_starsR/(alphaR*N0)));
% SCENARIO 2
lambda_starsBR = OmegaST*alphaBR*betaB*N0*snr*log(1/betaB)/(1-betaB);
optDEPvslambdasBR =  1+ func_Xi(lambda_starsBR,alphaBR*betaB*N0*snr,OmegaST,alphaBR*N0) - func_Xi(lambda_starsBR,alphaBR*N0*snr,OmegaST,alphaBR*N0);
% SCENARIO 3
lambda_starsRpBR = log(deltaR/(alphaBR*betaB))/(1/(alphaBR*betaB)-1/deltaR)*OmegaST*N0*snr;
optDEPvslambdasRpBR = 1+ func_Xi(lambda_starsRpBR,alphaBR*betaB*N0*snr,OmegaST,(alphaBR+alphaR)*N0) - func_Xi(lambda_starsRpBR,(alphaBR+alphaR)*snr,OmegaST,(alphaBR+alphaR)*N0);


%% Save result
simDEPvslambda = [simDEPsR;simDEPsBR;simDEPsRpBR ]; 
lambdadBsim =lambdadB;
anaDEPvslambda = [anaDEPsR;anaDEPsBR;anaDEPsRpBR ]; 
lambdadBana = lambdadB1;
simDEPnomavslambda = simDEPnoma;

optlambdadB = 10*log10([lambda_starsR, lambda_starsBR,lambda_starsRpBR]);
optDEPvslambda = [optDEPvslambdasR,optDEPvslambdasBR,optDEPvslambdasRpBR];

%% Plot result
figure(1)

% % Plot
subplot(1,2,1)
plot(snrdBsim, simDEPvsSNR(1,:),'o',  'linewidth',2,'MarkerSize',8.5,'Color',[0.93,0.69,0.13]); hold on;
plot(snrdBsim, simDEPvsSNR(2,:),'bs',  'linewidth',2,'MarkerSize',9.5); hold on;
plot(snrdBsim, simDEPvsSNR(3,:),'g^',  'linewidth',2,'MarkerSize',9.5); hold on;
plot(snrdBsim, simDEPnomavsSNR,'m-d',  'linewidth',2,'MarkerSize',9.5); hold on;

plot(optSNRdB, optDEPvsSNR,'rp',  'linewidth',2,'MarkerSize',9.5); hold on;

plot(snrdBana, anaDEPvsSNR,'k-',  'linewidth',2,'MarkerSize',8.5); hold on;

% 
% xlabel('SNR $\overline{\gamma}$ [dB]','FontName','Times New Roman','FontSize',15,'Interpreter','latex');
% ylabel('Detection Error Probability','FontName','Times New Roman','FontSize',15);
% 
% lgd1.NumColumns = 1;
% lgd1.FontSize = 13;
set(gca,'fontsize',14);
xlabel('SNR $\overline{\gamma}$ [dB]','FontSize',15,'Interpreter','latex') ;
ylabel('Detection Error Probability','FontSize',15);
axis([min(snrdBsim) max(snrdBsim) 0 1]);

% % Plot
subplot(1,2,2)
sim1 = plot(lambdadBsim, simDEPvslambda(1,:),'o',  'linewidth',2,'MarkerSize',8.5,'Color',[0.93,0.69,0.13]); hold on;
sim2 = plot(lambdadBsim, simDEPvslambda(2,:),'bs',  'linewidth',2,'MarkerSize',9.5); hold on;
sim3 = plot(lambdadBsim, simDEPvslambda(3,:),'g^',  'linewidth',2,'MarkerSize',9.5); hold on;
sim4 = plot(lambdadBsim, simDEPnomavslambda,'m-d',  'linewidth',2,'MarkerSize',9.5); hold on;

ana1 = plot(lambdadBana, anaDEPvslambda,'k-',  'linewidth',2,'MarkerSize',8.5); hold on;
opt1 = plot(optlambdadB, optDEPvslambda,'rp',  'linewidth',2,'MarkerSize',9.5); hold on;

% % Grid

lgd2=legend([sim1(1),sim2(1),sim3(1),sim4(1),ana1(1),opt1(1)],...
   'Scenario 1',...
   'Scenario 2',...
   'Scenario 3',...
   'NOMA',...
   'Theory',...
   'Optimal');

% ylim([0 1]);
% xlim([min(lambdadBsim)  max(lambdadBsim)]);
% xlabel('Threshold $\lambda$ [dB]','FontName','Times New Roman','FontSize',15,'Interpreter','latex');

lgd2.NumColumns = 1;
lgd2.FontSize = 13;
set(gca,'fontsize',14);
xlabel('Threshold $\lambda$ [dB]','FontSize',15,'Interpreter','latex') 
axis([min(lambdadBsim) max(lambdadBsim) 0 1]);


% function [out_snr] = find_SNR(lambda,alphaR,N0,OmegaST,snr)
% a = min(snr); b = max(snr)
% 
% func_x= @(z,a,OmegaX,OmegaY) OmegaX*a/(OmegaX*a - OmegaY)*( exp(- z/( OmegaX*a))  - exp(- z/OmegaY));
% % initial value
% xn = db2pow(0);  
% fxn  =  1 - func_x(lambda,alphaR*N0*xn,OmegaST,alphaR*N0)
% gxn = - exp(- lambda/( OmegaST*alphaR*N0*xn))*lambda/( OmegaST*alphaR*N0*xn^2);
% xn1 = xn - fxn/gxn;
% fxn1  =  1 - func_x(lambda,alphaR*N0*xn1,OmegaST,alphaR*N0)
% 
%     while abs(fxn1-fxn)^2/fxn1 > threshold
%         xn = xn1
%          fxn  =  1 - func_x(lambda,alphaR*N0*xn,OmegaST,alphaR*N0)
%          gxn = - exp(- lambda/( OmegaST*alphaR*N0*xn))*lambda/( OmegaST*alphaR*N0*xn^2);
%          xn1 = xn - fxn/gxn
%          fxn1  =  1 - func_x(lambda,alphaR*N0*xn1,OmegaST,alphaR*N0);
%     end
% 
% out_snr = xn;
% 
% end