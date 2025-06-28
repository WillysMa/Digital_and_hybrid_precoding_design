function [SE, precoder,power_allocatoin]=Narrowband_Waterfilling(H,Ns,P,noise_power)
episilon=0;
% Heff= 1/sqrt(noise_power)*H;
[Nr,~]=size(H);
[U,S,V]=svd(H);
sigular_vec=diag(S);
eigen_vec=sigular_vec.^2;
if Ns>1
rank_apx=min(Ns,sum(eigen_vec>episilon)); % speed up the water-filling alg
else
   rank_apx=Ns; 
end
eigen_eff=eigen_vec(1:rank_apx);

power_allocatoin=water_filling(P,rank_apx,eigen_eff,noise_power);
power_matrix=diag(power_allocatoin(1:rank_apx));
V_effect=V(:,1:rank_apx);
precoder=V_effect*sqrt(power_matrix);

SE=real(log2( det(eye(Nr) + 1/noise_power*H*precoder*(H*precoder)') ));

end

