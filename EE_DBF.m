function EE=EE_DBF(BW,bit,SE_set,Nr,Nt)
factor2=1e-3;
P_LNA=25*factor2;
P_PA=60*factor2;
% P_SP=19.5*factor2;
P_DAC=130*factor2;
% P_M=19*factor2;
% P_LO=5*factor2;
% P_LPF=14*factor2;
% P_BBamp=5*factor2;
P_RF=43*factor2;
% kappa=99;
kappa=494;
P_ADC=kappa*BW*2.^bit*1e-15;


Ptx=Nt*(P_PA+P_RF+2*P_DAC);


Prx=Nr*(P_LNA+P_RF+2*P_ADC);



P_DBF=Ptx+Prx;

EE=SE_set./P_DBF;
end