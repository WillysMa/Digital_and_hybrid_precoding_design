function EE=EE_SCHBF(BW,bit,SE_set,Nr,Nt,Nrfset)
factor2=1e-3;
P_LNA=25*factor2;
P_PA=60*factor2;
P_SP=19.5*factor2;
P_C=P_SP;
% P_SP=0*factor2;
P_DAC=130*factor2;
% P_M=19*factor2;
% P_LO=5*factor2;
% P_LPF=14*factor2;
% P_BBamp=5*factor2;
P_RF=43*factor2;
P_PS=23*factor2;
% kappa=99;
kappa=494;
P_ADC=kappa*BW*2.^bit*1e-15;

Ptx=Nrfset*(P_RF+2*P_DAC+P_SP)+Nt*(P_C+P_PA+P_PS);
Prx=Nr*(P_LNA+P_RF+2*P_ADC);

P_HBF=Ptx'+Prx;
EE=SE_set./P_HBF;
end