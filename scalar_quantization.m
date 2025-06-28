function x_qt=scalar_quantization(x_in,Var,Tb,Lb,type)
%'''
%x_in:  signal
%Tb: quantizer threshold list
%Lb: quantizer codebook
%type: either 'complex' or 'real'
%x_qt: output a quantized complex signal
%'''
if strcmp(type,'complex')
     for t=1:length(Tb)
         if real(x_in)/sqrt(Var/2)<Tb(t)
             out_n_real=Lb(t);
             break
         else
             out_n_real=Lb(end);
         end
     end
     for t=1:length(Tb)
         if imag(x_in)/sqrt(Var/2)<Tb(t)
             out_n_imag=Lb(t);
             break
         else
             out_n_imag=Lb(end);
         end
     end 
    
    x_qt=sqrt(Var/2)*(out_n_real+1i*out_n_imag); 
else
 for t=1:length(Tb)
     if x_in/sqrt(Var)<Tb(t)
         out=Lb(t);
         break
     else
         out=Lb(end);
     end
 end
 x_qt=sqrt(Var)*out; 
end

