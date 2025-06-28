function [SE_end, F, U,Obj]=WMMSE_DBFdesign(params_wmmse,PrecoderIni,H,Ns,noise_power,Bsg)
[Nr,Nt]=size(H);
tol=params_wmmse.tol;
iter_max=params_wmmse.Imax;
Pb=params_wmmse.TxPowerBudget;
% PrecoderIni=params_wmmse.PrecoderIni;

F_func=@(tmp_C,J,F) real(trace(J*F*F'))-2*real(trace(tmp_C*F' ));% define function 
grad_F=@(tmp_C,J,F) J*F-tmp_C;% define function
Proj=@(X) sqrt(Pb)*X/norm(X,'fro');

F=PrecoderIni;
tic

% Bsg=(1-distortion_factor)*eye(Nr);
Obj=[];Cache=[];
acc=1e-4;shrink_factor=1;   
SE_iter=1e-4;tar_func=1e-4;
W=eye(Ns);
 for i=1:iter_max
    tar_func_old=tar_func;
    SE_old=SE_iter;
    W_old=W;

    CoV_e=Bsg*(eye(Nr)-Bsg)*diag(diag(H*F*F'*H'))+noise_power*Bsg;
    U=inv(Bsg*H*F*F'*H'*Bsg+CoV_e)*Bsg*H*F;

%     W=inv(eye(Ns)+ F'*Heff'*inv(CoV_e)*Heff*F);
    W=eye(Ns)+ F'*H'*Bsg*inv(CoV_e)*Bsg*H*F;


    J=H'*(Bsg* U*W*U' + diag(diag(U*W*U'))*(eye(Nr)-Bsg) )*Bsg*H; 
    tmp=H'*Bsg*U*W;


    % func_obj=@(X) F_func(tmp,J,X);% specify default value      
    % func_grad=@(X) grad_F(tmp,J,X);
    % [Obj_Opt,F,iter_data]=PGD(F,Proj,func_obj,func_grad);

    
    F=pinv(J)*tmp;
    if  norm(F,'fro')^2<= Pb
          ;
    else

        mu_ub=norm(tmp,'fro')/sqrt(Pb);
        mu_lb=0;
        while abs(mu_ub-mu_lb)/abs(mu_ub)>=acc*shrink_factor
            mu=(mu_ub+mu_lb)/2;
            F=inv(J+mu*eye(Nt))*tmp;
            if norm(F,'fro')^2<= Pb
                mu_ub=mu;
            else
                mu_lb=mu;
            end

        end
    end
    
    
Ce=Bsg*(eye(Nr)-Bsg)*diag(diag(H*F*F'*H'))+noise_power*Bsg;
SE_iter=real(log2(det(eye(Nr)+inv(Ce)*Bsg*H*F*F'*H'*Bsg )));

MSE=U'*(Bsg*H*F*F'*H'*Bsg+CoV_e)*U+eye(Ns)-U'*Bsg*H*F-(U'*Bsg*H*F)';
tar_func=real(trace(W*MSE))-real(log2(det(W)));

 Obj=[Obj SE_iter];
Cache=[Cache tar_func];


 err=abs(SE_old-SE_iter)/SE_old;
 if err<=tol || SE_iter+0.1<SE_old
    break
 end


ccc=1;
 end
 t_WMMSE=toc;
 Ce=Bsg*(eye(Nr)-Bsg)*diag(diag(H*F*F'*H'))+noise_power*Bsg;
SE_end=real(log2(det(eye(Ns)+pinv(U'*Ce*U)*U'*Bsg*H*F*F'*H'*Bsg*U )));
 
% x=1:length(Obj);
% figure
% plot(x,Obj,'r-' ,'Linewidth',1.2)
% legend('Proposed WMMSE-based DBF')
% xlabel('Iteration')
% ylabel('SE (bits/s/Hz)')
% % xlim([x(1),x(end)])
% grid on
% box on
% 
% figure
% plot(Cache,'b-' ,'Linewidth',1.2)
% xlabel('Iteration')
% ylabel('Objective value')
% grid on
% box on
ccc=1;
end