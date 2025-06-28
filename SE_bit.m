function TC=SE_bit
clc;clear all;close all;
tic
importmanopt()
rng(0);
DistortionFactor=load('Distr_factor.mat');
DistortionFactor_list=DistortionFactor.Distr_factor;
% Quantizer_Type='Uniform';
Quantizer_Type='Non-uniform';
if strcmp(Quantizer_Type,'Uniform')
    distortion_list=DistortionFactor_list.Unf; 
   indx=2;
else
    distortion_list=DistortionFactor_list.Nnf;
   indx=1;
end
Nt=64;
[Nt_h, Nt_v] = array_dimension(Nt);

Nr=64;
[Nr_h, Nr_v] = array_dimension(Nr);

Ncl=5;
Nray=10;


data_opt=load('TL_VS_B.mat');


Pb_dBm=30;
Pb=10^(Pb_dBm/10);%mW

SNRx_dB=20;
% SNRx_dB=-10:5:30;
SNRx=10.^(SNRx_dB/10);
noise_power=Pb./SNRx;% noise power mW 0.3634, 


Imax=1e+5;tol_cov=1e-3;
iter_max=5000;In_max=1;pga_iter=1;tol_SE=1e-6;
params_wmmse.tol=tol_SE;
params_wmmse.Imax=iter_max;
params_wmmse.TxPowerBudget=Pb;


Nrf=8;
Ns=Nrf;

Ball=1:8;lenB=length(Ball);
MC=2;


SE_opt_all=zeros(MC,1);

Ini_Mtrx=zeros(MC,lenB);

WF_DBF=Ini_Mtrx;
JPC_DBF=Ini_Mtrx;
MMheur_DBF=Ini_Mtrx;

MM_FCHBF=Ini_Mtrx;
Mjh_FCHBF=Ini_Mtrx;
AO_FCHBF=Ini_Mtrx;

MM_SCHBF=Ini_Mtrx;
Mjh_SCHBF=Ini_Mtrx;
Sic_SCHBF=Ini_Mtrx;

WF_DBFapx=Ini_Mtrx;
JPC_DBFapx=Ini_Mtrx;
MMheur_DBFapx=Ini_Mtrx;

MM_FCHBFapx=Ini_Mtrx;
Mjh_FCHBFapx=Ini_Mtrx;
AO_FCHBFapx=Ini_Mtrx;

MM_SCHBFapx=Ini_Mtrx;
Mjh_SCHBFapx=Ini_Mtrx;
Sic_SCHBFapx=Ini_Mtrx;
for avr=1:MC 
       
    disp([' avr=', num2str(avr)])
%     H=1/sqrt(2)*(randn(Nr,Nt)+1i*randn(Nr,Nt));
    H=Channel_Gen_UPA(Nt_h,Nt_v,Nr_h,Nr_v,Ncl,Nray);
    

    [SE_opt_all(avr),F_wf,pv_wf]=Narrowband_Waterfilling(H,Ns,Pb,noise_power);


        for id=1:lenB
            b=Ball(id);Tb=data_opt.T_0{b,indx};Lb=data_opt.L_0{b,indx};           
            gamma=distortion_list(b);
            Bsg=(1-distortion_list(b))*eye(Nr);

            [CoVar_eta_wf,gamma_hat,ii]=CoV_eta_evaluated(H, F_wf,Ns,noise_power,gamma,Tb,Lb,Imax,tol_cov);% compute numerical CoVar            
            WF_DBF(avr,id)=rate_cal(CoVar_eta_wf,Bsg, noise_power,H,F_wf);
            F=F_wf;
            CoVar_eta_apx=gamma*(1-gamma)*diag(diag(H*F*(H*F)'+noise_power*eye(Nr)));
            WF_DBFapx(avr,id)=rate_cal(CoVar_eta_apx,Bsg, noise_power,H,F);

            [MMheur_DBFapx(avr,id),Fpv,~]=MMheuristic_DBF(pv_wf,H,Ns,Pb,noise_power,Bsg,tol_SE,iter_max);        
            [CoVar_eta_heur,gamma_hat,ii]=CoV_eta_evaluated(H, Fpv,Ns,noise_power,gamma,Tb,Lb,Imax,tol_cov);% compute numerical CoVar  
            MMheur_DBF(avr,id)=rate_cal(CoVar_eta_heur,Bsg, noise_power,H,Fpv);


            [JPC_DBFapx(avr,id), Fjpc, U,Obj]=WMMSE_DBFdesign(params_wmmse,F_wf,H,Ns,noise_power,Bsg);
            [CoVar_eta_jpc,gamma_hat,ii]=CoV_eta_evaluated(H, Fjpc,Ns,noise_power,gamma,Tb,Lb,Imax,tol_cov);% compute numerical CoVar  
            JPC_DBF(avr,id)=rate_cal(CoVar_eta_jpc,Bsg, noise_power,H,Fjpc);

            
            type='FC-HBF';
            [MM_FCHBFapx(avr,id),FmmFC,~]=MM_AltMin_HBF(F_wf,H,Bsg,Pb,noise_power,Nrf,tol_SE,iter_max,In_max,pga_iter,type);
            [CoVar_eta,gamma_hat,ii]=CoV_eta_evaluated(H, FmmFC,Ns,noise_power,gamma,Tb,Lb,Imax,tol_cov);% compute numerical CoVar  
            MM_FCHBF(avr,id)=rate_cal(CoVar_eta,Bsg, noise_power,H,FmmFC);

            [Mjh_FCHBFapx(avr,id),FmjhFC]=HBF_MJH(H,Bsg,Nrf,Ns,Pb,noise_power,type);
            [CoVar_eta,gamma_hat,ii]=CoV_eta_evaluated(H, FmjhFC,Ns,noise_power,gamma,Tb,Lb,Imax,tol_cov);% compute numerical CoVar  
            Mjh_FCHBF(avr,id)=rate_cal(CoVar_eta,Bsg, noise_power,H,FmjhFC);

            [AO_FCHBFapx(avr,id),F_AO,~]=MO_AltOpt_HBF(H,F_wf,Bsg,Pb,noise_power,Nrf,iter_max);
            [CoVar_eta_evl,gamma_hat,ii]=CoV_eta_evaluated(H, F_AO,Ns,noise_power,gamma,Tb,Lb,Imax,tol_cov);% compute numerical CoVar  
            AO_FCHBF(avr,id)=rate_cal(CoVar_eta_evl,Bsg, noise_power,H,F_AO);
            

            type='SC-HBF';
            [MM_SCHBFapx(avr,id),FmmSC,~]=MM_AltMin_HBF(F_wf,H,Bsg,Pb,noise_power,Nrf,tol_SE,iter_max,In_max,pga_iter,type);
            [CoVar_eta,gamma_hat,ii]=CoV_eta_evaluated(H, FmmSC,Ns,noise_power,gamma,Tb,Lb,Imax,tol_cov);% compute numerical CoVar  
            MM_SCHBF(avr,id)=rate_cal(CoVar_eta,Bsg, noise_power,H,FmmSC);


            [Mjh_SCHBFapx(avr,id),FmjhSC]=HBF_MJH(H,Bsg,Nrf,Ns,Pb,noise_power,type);
            [CoVar_eta,gamma_hat,ii]=CoV_eta_evaluated(H, FmjhSC,Ns,noise_power,gamma,Tb,Lb,Imax,tol_cov);% compute numerical CoVar  
            Mjh_SCHBF(avr,id)=rate_cal(CoVar_eta,Bsg, noise_power,H,FmjhSC);

            [Sic_SCHBFapx(avr,id),Fsic]=SIC_precoding(H,Bsg,Nrf,Pb,noise_power);
            [CoVar_eta_evl,gamma_hat,ii]=CoV_eta_evaluated(H, Fsic,Ns,noise_power,gamma,Tb,Lb,Imax,tol_cov);% compute numerical CoVar  
            Sic_SCHBF(avr,id)=rate_cal(CoVar_eta_evl,Bsg, noise_power,H,Fsic);   
          
            ccc=1;
        end
 
% time_mark=datestr(now,'mmmm-dd');
% file_name=['SE_Ns_',time_mark];
% save(file_name,'-v7.3')
end
TC=toc;



SE_opt_avr=mean(SE_opt_all);
WF_DBF_avr=mean(WF_DBF); 
JPC_DBF_avr=mean(JPC_DBF);
MMheur_DBF_avr=mean(MMheur_DBF);

MM_FCHBF_avr=mean(MM_FCHBF);
Mjh_FCHBF_avr=mean(Mjh_FCHBF);
AO_FCHBF_avr=mean(AO_FCHBF);

MM_SCHBF_avr=mean(MM_SCHBF);  
Mjh_SCHBF_avr=mean(Mjh_SCHBF);
Sic_SCHBF_avr=mean(Sic_SCHBF);



time_mark=datestr(now,'mmmm-dd');
file_name=['SE_Ns_',time_mark];
% save(file_name,'-v7.3')
%% draw figure, Evl



lineSpec1={'r-o';'k-s';'b-+';'m-^';'g-d';'c-x';{[0.6,0.2,0],'-*'};{[0.4940 0.1840 0.5560],'-p'};{[0.9290 0.6940 0.1250],'-v'}};
lineSpec2={'r--o';'k--s';'b--+';'m--^';'g--d';'c--x';{[0.6,0.2,0],'--*'};{[0.4940 0.1840 0.5560],'--p'};{[0.9290 0.6940 0.1250],'--v'}};
lineSpec3={'r:o';'k:s';'b:+';'m:^';'g:d';'c:x';{[0.6,0.2,0],':*'};{[0.4940 0.1840 0.5560],':p'};{[0.9290 0.6940 0.1250],':v'}};
lineSpec4={'r-.o';'k-.s';'b-.+';'m-.^';'g-.d';'c-.x';{[0.6,0.2,0],'-.*'};{[0.4940 0.1840 0.5560],'-.p'};{[0.9290 0.6940 0.1250],'-.v'}};
lineSpec0={'ro';'ks';'b+';'m^';'gd';'cx';{[0.6,0.2,0],'*'};{[0.4940 0.1840 0.5560],'p'};{[0.9290 0.6940 0.1250],'v'}};
% % set(h30,'Color',[0.6,0.2,0])%棕色
% % set(h31,'Color',[0.4940 0.1840 0.5560])
% % set(h32,'Color',[0.9290 0.6940 0.1250])
% % set(h0,'Color',[0 0.4470 0.7410])
% % % set(h4,'Color',[0.3010 0.7450 0.9330])

xx=Ball;
one_vector=ones(1,lenB);
figure
lineSpec=lineSpec1;
h0=plot(xx,SE_opt_avr*one_vector,'-h','color',[0 0.4470 0.7410],'LineWidth',1.5);hold on;
h1=plot(xx,WF_DBF_avr,lineSpec{1},'LineWidth',1.5);hold on;
h2=plot(xx,JPC_DBF_avr,lineSpec{2},'LineWidth',1.5);hold on;
h3=plot(xx,MMheur_DBF_avr,lineSpec{3},'LineWidth',1.5);hold on;

lineSpec=lineSpec2;
h4=plot(xx, MM_FCHBF_avr,lineSpec{4},'LineWidth',1.5);hold on;
h5=plot(xx,Mjh_FCHBF_avr,lineSpec{5},'LineWidth',1.5);hold on;
h6=plot(xx, AO_FCHBF_avr,lineSpec{6},'LineWidth',1.5);hold on;

lineSpec=lineSpec3;
h7=plot(xx, MM_SCHBF_avr,lineSpec{7}{2},'color',lineSpec{7}{1},'LineWidth',1.5);hold on;
h8=plot(xx, Mjh_SCHBF_avr,lineSpec{8}{2},'color',lineSpec{8}{1},'LineWidth',1.5);hold on;
h9=plot(xx, Sic_SCHBF_avr,lineSpec{9}{2},'color',lineSpec{9}{1},'LineWidth',1.5);hold on;


legend([h0,h1,h2,h3],'UqOpt','DBF:WF', 'DBF:JPC','DBF:proposed','FontSize',10,'FontWeight','bold')
xlabel('ADC resolution [bit]')
ylabel('SE [bits/s/Hz]')
xticks(xx)
xlim([xx(1) xx(end)])
grid on
box on
axesNew = axes('position',get(gca,'position'),'visible','off');
legend(axesNew,[h4,h5,h6],'FC-HBF:proposed','FC-HBF:SVD-based','FC-HBF:AO','FontSize',10,'FontWeight','bold')

axesNew = axes('position',get(gca,'position'),'visible','off');
legend(axesNew,[h7,h8,h9],'SC-HBF:proposed','SC-HBF:SVD-based','SC-HBF:SIC','FontSize',10,'FontWeight','bold')



%% EE versus bit

BW=1e+9;
xx=Ball;


EE_DBF_avr=EE_DBF(BW,Ball,MMheur_DBF_avr,Nr,Nt);
EE_FCHBF_avr=EE_FCHBF(BW,Ball,MM_FCHBF_avr,Nr,Nt,Nrf);
EE_SCHBF_avr=EE_SCHBF(BW,Ball,MM_SCHBF_avr,Nr,Nt,Nrf);
figure
h=plot(xx,EE_DBF_avr,'r-o',xx,EE_FCHBF_avr,'b--s',xx,EE_SCHBF_avr,'k:x','LineWidth',1.5);
legend([h(1), h(2), h(3)],'DBF:proposed','FC-HBF:proposed','SC-HBF:proposed')
xlabel('ADC resolution [bit]')
ylabel('EE (Gbits/J)')
xticks(xx)
xlim([xx(1) xx(end)])
grid on
box on


cccc=1;

end
