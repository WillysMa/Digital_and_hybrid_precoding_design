function [SE_iter,F,Cache_out]=MM_AltMin_HBF(F_wf,H,Bsg,Pb,noise_power,Nrf,tol_SE,Out_max,In_max,pga_iter,type)
[Nr,Nt]=size(H);
W=SC_HBF_antenna_mapping_matrix(H,Nrf,'Fix');
if type=='SC-HBF'
    Proj=@(X) W.*exp(1j*angle(X));
else
    Proj=@(X) exp(1j*angle(X));
end

fX=@(F) Bsg*H*F;
fT=@(F) Bsg*(eye(Nr)-Bsg)*diag(diag( H*F*F'*H'   )) + noise_power*Bsg;

[~,~,V]=svd(H);
Frf_ini=V(:,1:Nrf);

Frf=Proj(Frf_ini);
F=F_wf;
Fbb=pinv(Frf)*F;

% Out_max=iter_max;
Cache_out=[];SE_iter=1;

% In_max=100;
acc=1e-4;tol_in=1e-4;

for jj=1:Out_max
    SE_old=SE_iter;
    F_old=F;

    X_hat=fX(F_old);
    T_hat=fT(F_old);
    tmpT=inv(T_hat);
    Y=tmpT-inv(T_hat+X_hat*X_hat');
    K=H'*( Bsg*Y*Bsg+ diag(diag(Y))*Bsg*(eye(Nr)-Bsg) )*H;

    
    S=X_hat'*tmpT*Bsg*H;

    Obj_all=[];obj_iter=1;
    for ii=1:In_max
         %% update Frf
        obj_old=obj_iter;
%         [Frf, cost] = MO_Frf(Frf, Fbb,J,S);
%         [SE_new,Frf,iter_data]=PGD_FrfDesign(J,Fbb,S,Frf);
%         [SE_new,Frf,iter_data]=MM_PGD_FrfSCHBP(J,Fbb,S,Frf,W,);
        [SE_new,Frf,iter_data]=MM_PGD_Frf(K,Fbb,S,Frf,W,pga_iter,type);
        
%         BigMatrix = kron((Fbb*Fbb').',K);
%         Frf_vec = pinv(BigMatrix)*vec(S'*Fbb');
%         Frf_opt = reshape(Frf_vec, Nt,Nrf);
%         Frf=Proj(Frf_opt);
        %% update Fbb
        

%         [SE_new,Fbb,iter_data]=PGD_FbbDesign(params_pga,J,Fbb,S,Frf,Pb);

        J_hat=Frf'*K*Frf;
        tmp=Frf'*S';
        Fbb=pinv(J_hat)*tmp;
        if  norm(Frf*Fbb,'fro')^2<= Pb
              ;
        else
        
            mu_ub=norm((Frf'*Frf)^(-0.5)*tmp,'fro')/sqrt(Pb);
            mu_lb=0;
            while abs(mu_ub-mu_lb)/abs(mu_ub)>=acc
                mu=(mu_ub+mu_lb)/2;
                Fbb=inv(J_hat+mu*Frf'*Frf)*tmp;
                if norm(Frf*Fbb,'fro')^2<= Pb
                    mu_ub=mu;
                else
                    mu_lb=mu;
                end
        
            end
        end
    F=Frf*Fbb;  
    obj_iter=real(trace(K*F*F'))-2*real(trace(S*F));
    Obj_all=[Obj_all obj_iter];
    err=abs(obj_old-obj_iter)/abs(obj_old);
    if ii>=2 && err<= tol_in
        break   
    end
        ccc=1;
    end
%     figure
%     plot(Obj_all)
%     xlabel('Iteration')
%     ylabel('Objective')
%     % xlim([x(1),x(end)])
%     grid on
%     box on
    
    F=Frf*Fbb;
    Ce=Bsg*(eye(Nr)-Bsg)*diag(diag(H*F*F'*H'))+noise_power*Bsg;
    SE_iter=real(log2(det(eye(Nr)+inv(Ce)*Bsg*H*F*F'*H'*Bsg )));
    
    Cache_out=[Cache_out SE_iter];
    
    err_out=abs(SE_iter-SE_old)/SE_old;%
    if err_out<=tol_SE
        break
    end
    ccc=1;
end
% figure
% plot(Cache_out)
% xlabel('Iteration')
% ylabel('SE (bits/s/Hz)')
% % xlim([x(1),x(end)])
% grid on
% box on
% title(type)

cccc=1;
end