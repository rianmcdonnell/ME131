%% concatenate 
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

%computer acceleration from encoders
% m12=floor(max(sig21b{12}.Data)/8); %why max? could use length?
% final_full=find(sig21b{12}.Data==m12*8);
% percent=(max(sig21b{12}.Data)/8 - m12);
% i=1;
% v_b12=[];
% count=0;
% for i=1:length(sig21b{12}.Data)
%     if i<=m12*8
%         for k=1:8:final_full
%             ind12=find(sig21b{12}.Data==k);
%             ind12=ind12(1);
%             ind12_2=find(sig21b{12}.Data==k+8);
%             ind12_2=ind12_2(1);
%             v_b12=[v_b12 2*pi*.025/(sig21b{12}.Time(ind12_2)-sig21b{12}.Time(ind12))];
%             i=ind12_2;
%         end
%     else
%         v_b12=[v_b12 percent*2*pi*.025/(sig21b{12}.Time(end)-sig21b{12}.Time(i))];
%     end
%     count=count+1;
% end

% diff encoder counter
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
%     ind13=find(sig21b{13}.Data==i);
%     ind13=ind13(1);
%     ind13_2=find(sig21b{13}.Data==i+8);
%     ind13_2=ind13_2(1);
% 
%     v_b12(i)=2*pi*.025/(sig21b{12}.Time(ind12)-sig21b{12}.Time(ind12_2));
%     v_b13(i)=2*pi*.025/(sig21b{13}.Time(ind13)-sig21b{13}.Time(ind13_2));
%     if length(nonzeros([v_b12(i) v_b13(i)]))>=1
%     v_b(i)=v_b12(i)+v_b12(i)/length(nonzeros([v_b12(i) v_b13(i)])); 
%     else 
%         v_b(i)= 0; 
%     end
% end


% while i<=length(sig21c{12}.Data)-8
% v_c12(i)=2*pi*.025/(sig21c{12}.Time(i+8)-sig21c{12}.Time(i))^2;
% v_c13(i)=2*pi*.025/(sig21c{13}.Time(i+8)-sig21c{13}.Time(i))^2;
% end
% 
% while i<=length(sig21d{12}.Data)-8
% v_d12(i)=2*pi*.025/(sig21d{12}.Time(i+8)-sig21d{12}.Time(i))^2;
% v_d13(i)=2*pi*.025/(sig21d{13}.Time(i+8)-sig21d{13}.Time(i))^2;
% end
% 
% while i<=length(sig21f{12}.Data)-8
% v_f12(i)=2*pi*.025/(sig21f{12}.Time(i+8)-sig21f{12}.Time(i))^2;
% v_f13(i)=2*pi*.025/(sig21f{13}.Time(i+8)-sig21f{13}.Time(i))^2;
% end
% 
% 
% Linear_acceleration_x=sig21{7}.Data;
% motor_pwm=sig21{10}.Data;
% velocity_FL=sig21{16}.Data;
% velocity_FR=sig21{17}.Data;
% %velocity_BL=sig21{18}.Data;
% %velocity_BR=sig21{19}.Data;
% 
% 
% %average
% Vave=zeros(1,length(velocity_FR));
% for j=1:length(velocity_FR)
%     if velocity_FR(j)+velocity_FL(j)==velocity_FL(j)
%         Vave(j)=(velocity_FR(j)+velocity_FL(j))/1;
%     elseif velocity_FR(j)+velocity_FL(j)==velocity_FR(j)
%         Vave(j)=(velocity_FR(j)+velocity_FL(j))/1;
%     elseif velocity_FR(j)+velocity_FL(j)==0
%         Vave(j)=0;
%     else 
%         Vave(j)=(velocity_FR(j)+velocity_FL(j))/2;
%     end
% end
% Vave_clean=Vave;
% ind=find(Vave_clean==0);
% sw=150;
% ax=movmean(sig21{7}.Data,sw);
% %ax=smooth(sig21{7}.Data);
% 
% 
% ts_motor_pwm = timeseries(motor_pwm(:,1), sig21{10}.Time, 'name', 'motor_pwm');
% ts_ax = timeseries(ax(:,1), sig21{7}.Time, 'name', 'ax');
% ts_Vave = timeseries(Vave_clean(1,:)', sig21{16}.Time, 'name', 'Vave_clean');
% ts_Vave_clean = delsample(ts_Vave,'Index',ind);
% 
% resample_ts_motor_pwm=resample(ts_motor_pwm,ts_ax.Time,'linear');
% resample_ts_motor_pwm.Data=resample_ts_motor_pwm.Data - 1500; 
% resample_ts_Vave_clean=resample(ts_Vave_clean,ts_ax.Time,'linear');
% 
% % sw=10; 
% % ts_Vaveclean=movmean(resample_ts_Vave_clean.Data,sw);
% % ts_matrix=[ts_Vaveclean,resample_ts_motor_pwm.Data];
% ts_matrix=[resample_ts_Vave_clean.Data,resample_ts_motor_pwm.Data];
% 
% 
% ts_matrix(isnan(ts_matrix))=0;
% A=ts_matrix\ts_ax.Data;
% %Task 2_2
% plot(ts_motor_pwm.Data)
% ylabel('pwm')
% xlabel('time')
% title('pwm vs. time')
% set(gcf,'color','w');
% 
% %Task 2_1
% first_order=tf([A(2)],[1 -A(1)]);
% 
% %Task 2_3
% vdot=A(1).*ts_matrix(:,1)+A(2).*ts_matrix(:,2);
% plot(vdot)
% hold on
% plot(ts_ax.Data)