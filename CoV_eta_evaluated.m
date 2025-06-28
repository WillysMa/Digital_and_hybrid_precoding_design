function [CoVar_eta,gamma_hat,ii]=CoV_eta_evaluated(H, precoder,Ns,noise_power,distortion_factor,Tb,Lb,Imax,tol)

[Nr,Nt]=size(H);
gamma=distortion_factor;

CoVar_y = H*precoder*(H*precoder)'+noise_power*eye(Nr);
Var_y=real(diag(CoVar_y));

% samples=zeros(Nr,Imax);
% samples_quan=zeros(Nr,Imax,lenB);
CoVar_eta=zeros(Nr,Nr);

gamma_hat=0;
Flag=1;ii=0;
for ii=1:Imax
% while Flag && ii<=Imax
%         ii=ii+1; 
    xs = 1/sqrt(2)*(randn(Ns,1)+1i*randn(Ns,1)); % modulated symbols            
    xp=precoder*xs;% precoded signal
   
    noise_vec = sqrt(noise_power/2)*(randn(Nr,1) + 1i*randn(Nr,1)); % white gaussian noise
    y = H*xp + noise_vec;

%     samples(:,ii)=y;    
    sp_z=zeros(Nr,1);
    for jj=1:Nr    
        sp_z(jj)=scalar_quantization(y(jj),Var_y(jj),Tb,Lb,'complex');  
        tmp_gamma=abs(sp_z(jj)-y(jj))^2/Var_y(jj);
        ccc=1;
    end
    err_lmmse=sp_z-(1-gamma)*y;
    if ii==1
       gamma_hat=tmp_gamma;
       CoVar_eta=err_lmmse*err_lmmse';
    else
       gamma_hat=gamma_hat+1/ii*(tmp_gamma-gamma_hat);
       CoVar_eta = CoVar_eta + 1/ii*(err_lmmse*err_lmmse'-CoVar_eta);
    end
    err=abs(gamma_hat-gamma)./gamma;
    if sum(err)<=tol 
        Flag=0;    
    end      
end%->ii=1:Imax

ccc=1;
end


