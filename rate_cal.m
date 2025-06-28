function SE=rate_cal(CoVar_eta,Bsg, noise_power,H,precoder)
[Nr,Nt]=size(H);
 Cov_eff=CoVar_eta+noise_power*Bsg^2;
 SE=real(log2(det( eye(Nr)+inv(Cov_eff)*Bsg*H*precoder*(H*precoder)'*Bsg )));
end