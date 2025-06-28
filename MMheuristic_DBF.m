function [SE_mm,Fmm,Obj]=MMheuristic_DBF(pv_wf,H,Ns,Pb,noise_power,Bsg,tol,iter_max)
[~,~,V]=svd(H);
Vd=V(:,1:Ns);
qv=sqrt(pv_wf);
[Nr,Nt]=size(H);
fX=@(q) Bsg*H*Vd*diag(q);
fT=@(q) Bsg*(eye(Nr)-Bsg)*diag(diag( H*Vd*diag(q.^2)*Vd'*H'   )) + noise_power*Bsg;


% iter_max=1000;
% tol=1e-6;
Obj=[];
acc=1e-4;SE_iter=10;
for ii=1:iter_max
    SE_old=SE_iter;
    qv_old=qv;

    X_hat=fX(qv_old);
    T_hat=fT(qv_old);
    tmp=inv(T_hat);
    C= X_hat'*tmp*Bsg*H*Vd;
    Y=tmp-inv(T_hat+X_hat*X_hat');
    J=Vd'*H'*( Bsg*Y*Bsg+ diag(diag(Y))*Bsg*(eye(Nr)-Bsg) )*H*Vd;

    qv=diag(real(C))./diag(real(J));
%     for i=1:Ns
%         qv(i)=real(C(i,i))/real(J(i,i));
%     end
    if norm(qv)^2 <= Pb
        ;
    else
        mu_ub=norm(diag(real(C)),'fro')/sqrt(Pb);
        mu_lb=0;
        while abs(mu_ub-mu_lb)/abs(mu_ub)>=acc
            mu=(mu_ub+mu_lb)/2;
%             for i=1:Ns
%                 qv(i)=real(C(i,i))/(real(J(i,i))+mu);
%             end

            qv=diag(real(C))./diag(real(J)+mu);

            if norm(qv)^2<= Pb
                mu_ub=mu;
            else
                mu_lb=mu;
            end
    
        end
    end
    SE_iter=rate_pv(H,Bsg,Vd,qv,noise_power);

     Obj=[Obj SE_iter];
     err=abs(SE_old-SE_iter)/SE_old;
     if err<=tol
        break
     end
end
SE_mm=SE_iter;
Fmm=Vd*diag(qv);
% x=1:length(Obj);
% figure
% plot(x,Obj,'r-' ,'Linewidth',1.2)
% legend('Proposed MM-based power allocation')
% xlabel('Iteration')
% ylabel('SE (bits/s/Hz)')
% xlim([x(1),x(end)])
% grid on
% box on
cccc=1;
end
%% defined functions
function SE=rate_pv(H,Bsg,Vd,p_sqrt,noise_power)
[Nr,Nt]=size(H);
T=Bsg*(eye(Nr)-Bsg)*diag(diag( H*Vd*diag(p_sqrt.^2)*Vd'*H'   )) + noise_power*Bsg;
X=Bsg*H*Vd*diag(p_sqrt);
SE=real(log2(det( eye(Nr)+ inv(T)*X*X'   )));
end
