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
sig21b{12}.Time= (sig21b{12}.Time - sig21b{12}.Time(1)).*86400;
sig21b{13}.Time= (sig21b{13}.Time - sig21b{13}.Time(1)).*86400;

sig21c{12}.Time= (sig21c{12}.Time - sig21c{12}.Time(1)).*86400;
sig21c{13}.Time= (sig21c{13}.Time - sig21c{13}.Time(1)).*86400;

sig21d{12}.Time= (sig21d{12}.Time - sig21d{12}.Time(1)).*86400;
sig21d{13}.Time= (sig21d{13}.Time - sig21d{13}.Time(1)).*86400;

sig21f{12}.Time= (sig21f{12}.Time - sig21f{12}.Time(1)).*86400;
sig21f{13}.Time= (sig21f{13}.Time - sig21f{13}.Time(1)).*86400;

%computer velocity from each encoder
l=1;

encoderFL=sig21b{12};
while l+i<length(encoderFL.Data)
        i=1;
        while encoderFL.Data(l)==encoderFL.Data(l+i) && l+i<length(encoderFL.Data)
%         if encoderFR.Data(l)~=encoderFR.Data(l+1)
%             v_encoderFR(l)=(encoderFR.Data(l+1)-encoderFR.Data(l))/(encoderFR.Time(l+1)-encoderFR.Time(l));
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
v_encoderFR.Name=v_encoderFR;

v_encoderFL.Data=movmean(v_encoderFL.Data,sw);
v_encoderFL.Name=v_encoderFL;

%time derivative: acceleration 
acc_encoderFR=ts_derivative(v_encoderFR);
sw=30; 
accFR=movmean(acc_encoderFR,sw);

acc_encoderFL=ts_derivative(v_encoderFL);
accFL=movmean(acc_encoderFL,sw);

%average acceleration for each test
v_ave=zeros(max([length(acc_encoderFL) length(acc_encoderFR)]),1); 
acc_ave=zeros(max([length(acc_encoderFL) length(acc_encoderFR)]),1); 
for i=1:max([length(acc_encoderFL) length(acc_encoderFR)])
    if i<=min([length(acc_encoderFL) length(acc_encoderFR)])
        v_ave(i)=(v_encoderFL.Data(i)+v_encoderFR.Data(i))/2; 
        acc_ave(i)=(acc_encoderFL(i)+acc_encoderFR(i))/2; 
    else
        if length(acc_encoderFR)==max([length(acc_encoderFL) length(acc_encoderFR)])
           acc_ave(i) = acc_encoderFR(i); 
           v_ave(i)=v_encoderFR.Data(i);
        else
            acc_ave(i) = acc_encoderFL(i);
            v_ave(i)=v_encoderFL.Data(i);
        end
    end     
end
acc_ave_smooth=movmean(acc_ave,sw);
%concatenate all v, a to matrix

%get transfer function

%vdot ode

%compare v ave from encoder to vdot integral
