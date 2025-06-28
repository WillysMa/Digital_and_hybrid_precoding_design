function analog_bf=analog_bf_compute(H_avr,N,N_rf)
% C=randn(N_rf,5*N_rf) + 1j*randn(N_rf,5*N_rf);
% B=C*C';% B can be any invertible matrix
[U,~,V]=svd(H_avr);


Wb_ini=V(:,1:N_rf);
W_rf=zeros(N,N_rf);
for row_id=1:N
    for col_id=1:N_rf
        W_rf(row_id,col_id)=exp(1j*angle(Wb_ini(row_id,col_id)));%1/sqrt(N)*
    end
end
analog_bf=W_rf;
end

%% subfunctions

