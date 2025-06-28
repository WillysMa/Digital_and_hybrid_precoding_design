function [SE_new,Frf_new,iter_data]=MM_PGD_Frf(J,Fbb,S,Frf,W,pga_iter,type)
if type=='SC-HBF'
    Proj=@(X) W.*exp(1j*angle(X));
else
    Proj=@(X) exp(1j*angle(X));
end

tol=1e-4;
iter_max=pga_iter;
alpha=0.1;beta=0.3;

Obj_cal=@(J,Fbb,S,X) real(trace(J*X*Fbb*Fbb'*X'))-2*real(trace( Fbb*S*X ));
% Proj=@(X) W.*exp(1j*angle(X));

fw=Obj_cal(J,Fbb,S,Frf);grad=1;
fw_all=[];
for iter=1:iter_max
%     step_size=1/sqrt(iter+1);

    Frf_old=Frf;
    fw_old=fw;
    grad_old=grad;

    grad=J*Frf*Fbb*Fbb'-S'*Fbb';
    

    grad_unit=grad/norm(grad,'fro');

    step_size=1;counter=0;
    while 1 && counter<=30
        counter=counter+1;

        % x_new=Frf_old-step_size*grad_unit;
%         x_proj=exp(1j*angle(x_new));
        x_proj=Proj(Frf_old-step_size*grad_unit);
        obj_new=Obj_cal(J,Fbb,S,x_proj);
        obj_benchmark=fw-alpha*step_size*norm(grad_unit,"fro")^2;

        if obj_new<obj_benchmark
            break
        else
            step_size=beta*step_size;
        end
    end%backtracking line search


    Frf=Proj(Frf-step_size*grad_unit);
    fw=Obj_cal(J,Fbb,S,Frf);

    % grad_new=J*Frf*Fbb*Fbb'-S'*Fbb';

    fw_all=[fw_all,fw];
%     err=fw-fw_old;
    err=norm(grad-grad_old,'fro')^2/norm(grad_old,'fro')^2;
    if err<=tol
        break
    end
end
Frf_new=Frf;
SE_new=fw;
iter_data=fw_all;

% figure
% plot(fw_all)
% ylabel('Objective')
% xlabel('iteration')
% grid on
% box on
% % xlim([1,length(fw_all)])
% title('PGD Algorithm')
cccc=1;
end




