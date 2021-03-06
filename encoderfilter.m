clear all
%load data
load sig21b
load sig21c
load sig21d
load sig21f

%fix switched
encoderc1=sig21c{16};
encoderc2=sig21c{17};
velocity_FL=sig21c{12};
sig21c{17}=[];
velocity_FR=sig21c{13};
sig21c{12}=[];
sig21c{13}=[];
sig21c{16}=[];
sig21c{17}=[];
sig21c{12}=encoderc1;
sig21c{13}=encoderc2;
sig21c{16}=velocity_FL;
sig21c{17}=velocity_FR;

encoderf1=sig21f{16};
encoderf2=sig21f{17};
velocity_FL=sig21f{12};
sig21f{17}=[];
velocity_FR=sig21f{13};
sig21f{12}=[];
sig21f{13}=[];
sig21f{16}=[];
sig21f{17}=[];
sig21f{12}=encoderf1;
sig21f{13}=encoderf2;
sig21f{16}=velocity_FL;
sig21f{17}=velocity_FR;

%concatenate
for i=1:19
sig21{i}.Data=[sig21b{i}.Data;sig21c{i}.Data;sig21d{i}.Data;sig21f{i}.Data];
sig21{i}.Time=[sig21b{i}.Time;sig21c{i}.Time;sig21d{i}.Time;sig21f{i}.Time];
end 

%encoder time: only first two
sig21b{10}.Time= (sig21b{10}.Time - sig21b{10}.Time(1)).*86400;
sig21b{12}.Time= (sig21b{12}.Time - sig21b{12}.Time(1)).*86400;
sig21b{13}.Time= (sig21b{13}.Time - sig21b{13}.Time(1)).*86400;

sig21c{10}.Time= (sig21c{10}.Time - sig21c{10}.Time(1)).*86400;
sig21c{12}.Time= (sig21c{12}.Time - sig21c{12}.Time(1)).*86400;
% sig21c{12}.Time= sig21c{12}.Time +
% sig21b{12}.Time(length(sig21b{12}.Time)); % useful for concatenation!
sig21c{13}.Time= (sig21c{13}.Time - sig21c{13}.Time(1)).*86400;
% sig21c{13}.Time= sig21c{13}.Time +
% sig21b{13}.Time(length(sig21b{13}.Time)); % useful for concatenation!

sig21c{10}.Time= (sig21c{10}.Time - sig21c{10}.Time(1)).*86400;
sig21d{12}.Time= (sig21d{12}.Time - sig21d{12}.Time(1)).*86400;
sig21d{13}.Time= (sig21d{13}.Time - sig21d{13}.Time(1)).*86400;

sig21c{10}.Time= (sig21c{10}.Time - sig21c{10}.Time(1)).*86400;
sig21f{12}.Time= (sig21f{12}.Time - sig21f{12}.Time(1)).*86400;
sig21f{13}.Time= (sig21f{13}.Time - sig21f{13}.Time(1)).*86400;

%computer velocity from each encoder
l=1;

encoderFL=sig21b{12};
while l+i<length(encoderFL.Data)
        i=1;
        while encoderFL.Data(l)==encoderFL.Data(l+i) && l+i<length(encoderFL.Data)
            i=i+1;
        end
        v_encoderFL.Data(l)=(((encoderFL.Data(l+i)-encoderFL.Data(l))/8)*2*pi*.025)/(encoderFL.Time(l+i)-encoderFL.Time(l));
        v_encoderFL.Time(l)=encoderFL.Time(l);
        l=l+1;
end

l=1;
encoderFR=sig21b{13};
while l+i<length(encoderFR.Data)
        i=1;
        while encoderFR.Data(l)==encoderFR.Data(l+i) && l+i<length(encoderFR.Data)
            i=i+1;
        end
        v_encoderFR.Data(l)=(((encoderFR.Data(l+i)-encoderFR.Data(l))/8)*2*pi*.025)/(encoderFR.Time(l+i)-encoderFR.Time(l));
        v_encoderFR.Time(l)=encoderFR.Time(l);
        l=l+1;
end

%smooth 
sw=40;

v_encoderFR.Data=movmean(v_encoderFR.Data,sw);


v_encoderFL.Data=movmean(v_encoderFL.Data,sw);


%time derivative: acceleration 
% acc_encoderFR=ts_derivative_edit(v_encoderFR);
acc_encoderFR.Data=[];
acc_encoderFR.Data(1)=0;
for j=2:length(v_encoderFR.Data)
  acc_encoderFR.Data(j)=(v_encoderFR.Data(j)-v_encoderFR.Data(j-1))/(v_encoderFR.Time(j)-v_encoderFR.Time(j-1));
end
acc_encoderFR.Time=v_encoderFR.Time;

acc_encoderFL.Data=[];
acc_encoderFL.Data(1)=0;
for j=2:length(v_encoderFL.Data)
  acc_encoderFL.Data(j)=(v_encoderFL.Data(j)-v_encoderFL.Data(j-1))/(v_encoderFL.Time(j)-v_encoderFL.Time(j-1));
end
acc_encoderFL.Time=v_encoderFL.Time;

sw=30;
accFR_smooth=movmean(acc_encoderFR.Data,sw);

% acc_encoderFL=ts_derivative_edit(v_encoderFL);
accFL_smooth=movmean(acc_encoderFL.Data,sw);

%average acceleration for each test
v_ave=zeros(max([length(acc_encoderFL.Data) length(acc_encoderFR.Data)]),1); 
acc_ave=zeros(max([length(acc_encoderFL.Data) length(acc_encoderFR.Data)]),1); 
for i=1:max([length(acc_encoderFL.Data) length(acc_encoderFR.Data)])
    if i<=min([length(acc_encoderFL.Data) length(acc_encoderFR.Data)])
        v_ave(i)=(v_encoderFL.Data(i)+v_encoderFR.Data(i))/2; 
        acc_ave(i)=(acc_encoderFL.Data(i)+acc_encoderFR.Data(i))/2; 
    else
        if length(acc_encoderFR.Data)==max([length(acc_encoderFL.Data) length(acc_encoderFR.Data)])
           acc_ave(i) = acc_encoderFR.Data(i); 
           v_ave(i)=v_encoderFR.Data(i);
        else
            acc_ave(i) = acc_encoderFL.Data(i);
            v_ave(i)=v_encoderFL.Data(i);
        end
    end     
end
acc_ave_smooth=movmean(acc_ave,sw);

if length(v_encoderFR.Time)>length(v_encoderFL.Time)
    vel_structb.Time= v_encoderFR.Time;
    acc_structb.Time= v_encoderFR.Time;
else
    vel_structb.Time= v_encoderFL.Time;
    acc_structb.Time= v_encoderFL.Time;
end
acc_structb.Data=acc_ave_smooth';
vel_structb.Data=v_ave';
motor_pwm=timeseries;

motor_pwm.Data= sig21b{10}.Data-1500;
motor_pwm.Time= sig21b{10}.Time;
% motor_short=timeseries;

resample_motor_pwm=resample(motor_pwm,acc_structb.Time,'linear');


for i=1:length(acc_structb.Time)
    v(i,1)=vel_structb.Data(i);
    v(i,2)=resample_motor_pwm.Data(i);
end

A=v\acc_ave_smooth;

first_order=tf([A(2)],[1 -A(1)]);

v_est.Time=vel_structb.Time;
for i=1:length(v_est.Time)
    v_est.Data(i)=(acc_ave_smooth(i)-A(2)*resample_motor_pwm.Data(i))/A(1);
end
plot(vel_structb.Time,vel_structb.Data)
hold on
plot(v_est.Time,v_est.Data)
% vdot=A(1).*v(:,1)+A(2).*v(:,2);
% plot(v_est.Time,vdot)
% hold on
% plot(acc_structb.Time,acc_structb.Data)


%concatenate all v, a to matrix

%get transfer function

%vdot ode

%compare v ave from encoder to vdot integral
