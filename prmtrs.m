close all
clear 
clc

[rs,fs]=audioread('d:\Users\ataka\OneDrive\Masaüstü\bildiri impulse responses\S2R4_sweep4000.wav');
%%  sum channels to mono
szr=size(rs); % checking the column number of signal

if szr(2)==2   % if the signal is stereo
r=(rs(:,1)+rs(:,2)); %mix channels
else            %if mono, 
    r=rs;   %keep it the same 
end
%% filter
fm=[63 125 250 500]; 
%     1000 2000 4000 8000]; 
results= [];

for f_index=1:length(fm)
    F0=fm(f_index);
B  = 1;      % Bands per octave
N  = 6;      % Order
    % Center frequency
Fs = 48000;  % Sampling Frequency
    
h = fdesign.octave(B, 'Class 0', 'N,F0', N, F0, Fs);
Hd = design(h, 'butter', ...
        'SOSScaleNorm', 'Linf');  
set(Hd,'PersistentMemory',true);
y = filter(Hd,r);
audioOut=y;
%% hilbert transform
env=abs(hilbert(audioOut)); %envelope
%% dB conversion
dbl=20*log10(env/max(env)); %db conversion
lpl=length(dbl); %length of all
T= 1/fs;
tpl=0:T:(lpl-1)/fs; %time vector of all

rtz=find(dbl==max(dbl)); %cut beginning of ır to direct sound
db=dbl((rtz-250):end);
lp=length(db);
tp=0:T:(lp-1)/fs;
nfft=lp;
env=env((rtz-250):end);
%% moving average filter
a=7001;
num=(1/a)*ones(1,a);
den=(1);
db_smoothed=filter(num,den,db);
%% noise level estimation, index of intersection% 
[abrch]=ischange(db,'linear','MaxNumChanges',1); %find a point signal changed abruptly, probably near to noise level, steadty-state position
nl=find(abrch==1); % common point at reveberation line and noise level will be used schcroeder integral formula as upper integral limit
%% polinom
p = polyfit(tp,db_smoothed,1); %the best fit polinom function of moving average filtered signal
fit = polyval(p,tp); 
%% figure
% figure(1)
% plot(tp,20*log10(r))
%% parameters - rt60,edt,t20,t30 from moving average filter
rt60i=find(abs(db_smoothed+60) < 0.003); 
if rt60i~=0  
rt60i=rt60i(1);
rt60=tp(rt60i);
else
    rt60=0;
end

rt5i=find(abs(db_smoothed+5) < 0.01);
rt5i=rt5i(1);

rt35i=find(abs(db_smoothed+35) < 0.01);

TF = isempty(rt35i);
if TF==1
    rt35i=find(abs(db_smoothed+25) < 0.01);
    rt35i=rt35i(1);
else
    rt35i=rt35i(1);
end

p = polyfit(tp(rt5i:rt35i),db_smoothed(rt5i:rt35i),1); %the best fit line between 5-35 dB to estimate 60 dB LEVEL
fitrt35i = polyval(p,tp(rt5i:rt35i)); 

rt30 = interp1(fitrt35i,tp(rt5i:rt35i),-60, 'linear', 'extrap');

if fit(end)>-20
    rt20='incomputable';
else
rt25i=find(abs(db_smoothed+25) < 0.01);
rt25i=rt25i(1);

p = polyfit(tp(1:rt25i),db_smoothed(1:rt25i),1); 
fitrt25i = polyval(p,tp(rt5i:rt25i)); 

rt20 = interp1(fitrt25i,tp(rt5i:rt25i),-60, 'linear', 'extrap');
end

edti=find(abs(db_smoothed+10) < 0.01);
edti=edti(1);

p = polyfit(tp(1:edti),db_smoothed(1:edti),1); 
fitedt = polyval(p,tp(1:edti)); 

edt = interp1(fitedt,tp(1:edti),-60, 'linear', 'extrap');

t15i=find(abs(db_smoothed+15) < 0.01);
t15i=t15i(1);

p = polyfit(tp(rt5i:t15i),db_smoothed(rt5i:t15i),1); 
fitt10 = polyval(p,tp(rt5i:t15i)); 
t10i = interp1(fitt10,tp(rt5i:t15i),-60, 'linear', 'extrap');

if rt60==0
   rt60=tp(rt35i)*2;
end
%% figure
txt = [num2str(F0),' Hz'];
figure(f_index)

% sgtitle(txt,'fontname','times')

% subplot(3,3,1)
% plot(tp,db)
% hold on
% plot(tp,db_smoothed,'k')
% hold on
% set(gca,'fontname','times','fontsize',10) 
% if rt60==tp(rt35i)*2
%     rt60i=rt35i;
% end
% 
% plot(tp(rt60i),db_smoothed(rt60i),'r*');
% title('RT60 from moving average filter')
% set(gca,'fontname','times','fontsize',10) 
% subplot(3,3,2)
% plot(tp(1:rt35i),db_smoothed(1:rt35i), 'k')
% hold on
% plot(tp(rt5i:rt35i),fitrt35i,'r','Linewidth',1.5,'LineStyle','--')
% title('T30 from moving average filter')
% set(gca,'fontname','times','fontsize',10) 
% subplot(3,3,3)
% plot(tp(1:rt25i),db_smoothed(1:rt25i), 'k')
% hold on
% plot(tp(rt5i:rt25i),fitrt25i,'r','Linewidth',1.5,'LineStyle','--')
% title('T20 from moving average filter')
% set(gca,'fontname','times','fontsize',10) 
% subplot(3,3,4)
% plot(tp(1:edti),db_smoothed(1:edti), 'k')
% hold on
% plot(tp(1:edti),fitedt,'r','Linewidth',1.5,'LineStyle','--')
% title('EDT from moving average filter')

%% 2. METHOD schroeder curve method, input will be envelope signal
sch=10*log10(cumsum(env(nl:-1:1).^2)/sum(env(1:nl).^2)); %db conversion + integration, nl=noise level (upper limit of integral)
schend=10*log10(cumsum(env(end:-1:1).^2)/sum(env(1:end).^2)); %db conversion + integration, nl=noise level (upper limit of integral)
sch=sch(end:-1:1); %fliplr
schend=schend(end:-1:1);
%% rt60 regression from schroeder
p60=polyfit(tp(1:nl),sch(1:nl),1); 
fitsch60 = polyval(p60,tp(1:nl)); 
rt60sch = interp1(fitsch60,tp(1:nl),-60, 'linear', 'extrap');

subplot(2,2,1)
% plot(tp,db,'k')
% hold on
plot(tp(1:nl),sch,'k','Linewidth',1.3)
hold on
plot(tp(1:nl),fitsch60,'r','Linewidth',1.7,'LineStyle','--')
% plot(tp,schend,'b','Linewidth',1.3)
% legend('RIR','∞ = Geçiş değeri','∞ = RIR sonu ')
% 'RT60 Regression')
title('RT60 Regresyon Analizi')
set(gca,'fontname','times','fontsize',14) 
xlabel('Zaman (s)')
ylabel('Genlik (dB)')
%% edt regression from schroeder
edtschi=find(abs(sch+10) < 0.01);
edtschi=edtschi(1);
p = polyfit(tp(1:edtschi),sch(1:edtschi),1); 
fitschedt = polyval(p,tp(1:edtschi)); 
edtsch = interp1(fitschedt,tp(1:edtschi),-60, 'linear', 'extrap');
subplot(2,2,2)
plot(tp(1:edtschi),sch(1:edtschi), 'k','Linewidth',1.3)
hold on
plot(tp(1:edtschi),fitschedt,'r','Linewidth',1.7,'LineStyle','--')
title('EDT Regresyon Analizi')
set(gca,'fontname','times','fontsize',14)
xlabel('Zaman (s)')
ylabel('Genlik (dB)')

%% t20-t30 schroeder calculation
t5schi=find(abs(sch+5) < 0.01);
t25schi=find(abs(sch+25) < 0.01);
t35schi=find(abs(sch+35) < 0.01);

p20=polyfit(tp(t5schi:t25schi),sch(t5schi:t25schi),1); 
p30=polyfit(tp(t5schi:t35schi),sch(t5schi:t35schi),1); 
fitsch20 = polyval(p20,tp(t5schi:t25schi)); 
fitsch30 = polyval(p30,tp(t5schi:t35schi)); 

t20sch = interp1(fitsch20,tp(t5schi:t25schi),-60, 'linear', 'extrap');
t30sch = interp1(fitsch30,tp(t5schi:t35schi),-60, 'linear', 'extrap');
% title('Noise level - Intersection point')
%% noiselevel from ischange function
% subplot(3,3,7)
% figure(53)
% plot(tp,db,'k')
% hold on
% plot(tp(nl),db(nl),'r*','Linewidth',2)
% xline(nl/fs, 'r:','Linewidth',2)
% legend('RIR','Geçiş Değeri','Kesim çizgisi')
% set(gca,'fontname','times','fontsize',16) 
% xlabel('Zaman (s)')
% ylabel('Genlik (dB)')
%% t20-t30 plot from Schroeder
subplot(2,2,3)
plot(tp(1:t25schi),sch(1:t25schi), 'k','Linewidth',1.3)
hold on
plot(tp(t5schi:t25schi),fitsch20,'r','Linewidth',1.7,'LineStyle','--')
title('T20 Regresyon Analizi')
set(gca,'fontname','times','fontsize',14) 
xlabel('Zaman (s)')
ylabel('Genlik (dB)')
subplot(2,2,4)
plot(tp(1:t35schi),sch(1:t35schi), 'k','Linewidth',1.3)
hold on
plot(tp(t5schi:t35schi),fitsch30,'r','Linewidth',1.7,'LineStyle','--')
title('T30 Regresyon Analizi')
set(gca,'fontname','times','fontsize',14) 
xlabel('Zaman (s)')
ylabel('Genlik (dB)')

%% result matrix store
result=[rt60;rt30;rt20;edt;edtsch;t30sch;t20sch;rt60sch];
results=[results result];
end

%% table
T = table({'RT60 moving average';'T30 moving average';'T20 moving average';'EDT Moving Average';'EDT Schroeder Curve';'T30 SCH';'T20 SCH';'RT60 SCH'},[results(:,1)], ...
          [results(:,2)],[results(:,3)],[results(:,4)],[results(:,5)],[results(:,6)],[results(:,7)],[results(:,8)]);
T.Properties.VariableNames = [{'Parameters'},'63','125','250','500','1000','2000','4000','8000']


figure(50)
plot(tpl,audioOut,'k');
hold on
plot(tpl,env,'r');
xlabel('Zaman (s)')
ylabel('Genlik')
set(gca,'fontname','times','fontsize',16)
xlim([0.15 0.8])
legend('RIR','Sinyal Zarfı')


