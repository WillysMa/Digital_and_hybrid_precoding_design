function [SE,F]=SIC_precoding(H,Bsg,Nt_RF,Pb,noise_power)
[Nr,Nt]=size(H);
M=Nt/Nt_RF; % number of antennas connected to one RF chains
F=zeros(Nt,Nt_RF);
coef=1/noise_power;
for i=1:Nt_RF
    G=H'*inv(eye(Nr)+coef*H*F(:,1:(i-1))*F(:,1:(i-1))'*H')*H;
    f=zeros(Nt,1);
    temp=G(M*(i-1)+1:M*(i-1)+M,M*(i-1)+1:M*(i-1)+M);
    [~,~,V]=svd(temp);
    op=V(:,1);
    phase=exp(1i*angle(op))/sqrt(M);
    a=(op'*phase+phase'*op)/(2*phase'*phase);
    f(M*(i-1)+1:M*(i-1)+M)=a*phase;
    F(:,i)=f;
end
F=sqrt(Pb/Nt_RF)*F;

Ce=Bsg*(eye(Nr)-Bsg)*diag(diag(H*F*F'*H'))+noise_power*Bsg;
SE=real(log2(det(eye(Nr)+inv(Ce)*Bsg*H*F*F'*H'*Bsg )));
end