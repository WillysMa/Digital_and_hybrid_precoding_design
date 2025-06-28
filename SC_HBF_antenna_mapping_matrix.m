function W=SC_HBF_antenna_mapping_matrix(H,N_rf,type)
[Nr,Nt]=size(H);
M=Nt/N_rf;
W=zeros(Nt,N_rf);

if type=='Dyn'

    Hstrength=zeros(1,Nr);
    for rol=1:Nr
    Hstrength(rol)=norm(H(rol,:));
    end
    [v,ind]=sort(Hstrength,'descend');
    Htarg=H(ind(1:N_rf),:);
    Hamp=abs(Htarg)';
    
    for m=1:M
        for n=1:N_rf
            [v,id]=max(Hamp(:,n));
            W(id,n)=1;
            Hamp(id,:)=0;
            cccc=1;
        end
    end

else

    for i=1:N_rf
        index=(i-1)*M+1:i*M;
        W(index,i)=1;
    end
end


end