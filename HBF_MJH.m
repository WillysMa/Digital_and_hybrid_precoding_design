function [SE,F]=HBF_MJH(H,Bsg,Nrf,Ns,Pb,noise_power,type)
W=SC_HBF_antenna_mapping_matrix(H,Nrf,'Fix');
if type=='SC-HBF'
    Proj=@(X) W.*exp(1j*angle(X));
else
    Proj=@(X) exp(1j*angle(X));
end

[Nr,~]=size(H);
[~,~,V]=svd(H);
Frf=V(:,1:Nrf);
Frf=Proj(Frf);
dis_cache=[];iter_max=100;
for ii=1:iter_max
    Frf_old=Frf;
    

%     temp=(Frf*Frf)^(-0.5);
    Frf=Frf*(Frf'*Frf)^(-0.5);
    Frf=Proj(Frf);
    dis=norm(Frf_old-Frf,'fro')^2;

    dis_cache=[dis_cache dis];
    err=dis/norm(Frf_old,'fro')^2;
    if err<=1e-6
        break
    end
end
% figure
% plot(dis_cache)
% xlabel('Iterations')
% ylabel('Distance')
% grid on
% box on

G=H*Frf;
[SE_opt,Fbb_hat]=Narrowband_Waterfilling(G,Ns,Pb,noise_power);
Fbb=(Frf'*Frf)^(-0.5)*Fbb_hat;

F=Frf*Fbb;
Ce=Bsg*(eye(Nr)-Bsg)*diag(diag(H*F*F'*H'))+noise_power*Bsg;
SE=real(log2(det(eye(Nr)+inv(Ce)*Bsg*H*F*F'*H'*Bsg )));

end