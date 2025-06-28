function H=Channel_Gen_UPA(Nt_h,Nt_v,Nr_h,Nr_v,Ncl,Nray)
% mmWave channel model
% L=10; % path number

Nt=Nt_h*Nt_v;
Nr=Nr_h*Nr_v;
% Lp=Ncl*Nray;
var=1;
path_gain_all=sqrt(var/2)*(randn(Ncl,Nray)+1i*randn(Ncl,Nray));


AoD_H_mean=(2*rand(1,Ncl)-1)*pi; AoD_V_mean=(2*rand(1,Ncl)-1)*pi/2;
AoA_H_mean=(2*rand(1,Ncl)-1)*pi;AoA_V_mean=(2*rand(1,Ncl)-1)*pi/2;
std_H=10/180*pi; std_V=3/180*pi;

AoD_H_all=zeros(Ncl,Nray);
AoA_H_all=zeros(Ncl,Nray);
AoD_V_all=zeros(Ncl,Nray);
AoA_V_all=zeros(Ncl,Nray);
for i=1:Ncl
    AoD_H_all(i,:)=laprnd(1,Nray,AoD_H_mean(i),std_H);
    AoA_H_all(i,:)=laprnd(1,Nray,AoA_H_mean(i),std_H);
    AoD_V_all(i,:)=laprnd(1,Nray,AoD_V_mean(i),std_V);
    AoA_V_all(i,:)=laprnd(1,Nray,AoA_V_mean(i),std_V);
end

H=0;
for i=1:Ncl
    for l=1:Nray

            ar=array_steering_dictionary(AoA_H_all(i,l),AoA_V_all(i,l),Nr_h,Nr_v);
            at=array_steering_dictionary(AoD_H_all(i,l),AoD_V_all(i,l),Nt_h,Nt_v);
            H=H+path_gain_all(i,l)*ar*at';
    end
end

beta=sqrt(Nt*Nr/(Ncl*Nray));
H=beta*H;

end       

% ccc=1;

%% defined functions

function array_vector=array_steering_dictionary(Angle_H,Agngle_V,Nh,Nv)

spatial_h=sin(Angle_H).*sin(Agngle_V);
spatial_v=cos(Agngle_V);
factor_h=0:Nh-1;
factor_v=0:Nv-1;
steering_vector_h = 1/sqrt(Nh)*exp(-1i*pi*factor_h'*spatial_h);
steering_vector_v = 1/sqrt(Nv)*exp(-1i*pi*factor_v'*spatial_v);
array_vector= kron(steering_vector_h,steering_vector_v);

ccc=1;
end