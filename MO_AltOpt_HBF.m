function [SE,F,Obj_all]=MO_AltOpt_HBF(H,Fopt,Bsg,Pb,noise_power,N_rf,iter_max)
% only works for FC-HBF
[Nr,Nt]=size(H);
Imax=iter_max;
Obj_all=[];Obj_new=1;
Frf=analog_bf_compute(H,Nt,N_rf);
tol_obj=1e-4;
for ii=1:Imax
    Obj_old=Obj_new;
    Fbb=pinv(Frf)*Fopt;
    [Frf, cost] = MO_Frf_LS(Frf, Fbb,Fopt);

    Obj_new=norm(Fopt-Frf*Fbb,'fro');
    Obj_all=[Obj_all Obj_new];
    err=abs(Obj_new-Obj_old)/Obj_old;
    if err<=tol_obj
        break
    end
end
% figure
% plot(Obj_all)
% xlabel('Iteration')
% ylabel('Objective')
% % xlim([x(1),x(end)])
% grid on
% box on
Fbb=sqrt(Pb)*Fbb/norm(Frf*Fbb,'fro');%scaling power
F=Frf*Fbb;

Ce=Bsg*(eye(Nr)-Bsg)*diag(diag(H*F*F'*H'))+noise_power*Bsg;
SE=real(log2(det(eye(Nr)+inv(Ce)*Bsg*H*F*F'*H'*Bsg )));

ccc=1;
end