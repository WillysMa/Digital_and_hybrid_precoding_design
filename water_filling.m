function power_allocatoin=water_filling(P,sum_rank,eigen_all_vec,noise_power)

level_upper=P+noise_power/min(eigen_all_vec);
level_lower=min(0,noise_power/max(eigen_all_vec));
power_all=zeros(sum_rank,1);
while abs(level_upper-level_lower)/level_upper > 1e-4
    water_level=(level_upper+level_lower)/2;
    for i=1:sum_rank        
    power_all(i)=max(0,water_level-noise_power/eigen_all_vec(i));
         if power_all(i)>=P
              power_all(i)=P;
          end
    end   
    if sum(power_all)>=P
        level_upper=water_level;
    else
        level_lower=water_level;
    end
end
power_allocatoin=power_all;
end