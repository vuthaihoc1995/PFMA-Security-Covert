function [rhoB,rhoR,deltaB,deltaR,vB,vR,psiB,psiR] = func_para(alphaB,alphaR,alphaBR,betaB,betaR)

rhoB    = alphaB + alphaBR*betaB;
rhoR    = alphaR + alphaBR*betaR;
deltaB = alphaB + alphaBR;
deltaR = alphaR + alphaBR;

vB = rhoB/(alphaBR*betaR);
vR = rhoR/(alphaBR*betaB);
psiB = deltaB/(alphaBR*betaR);
psiR = deltaR/(alphaBR*betaB);
end