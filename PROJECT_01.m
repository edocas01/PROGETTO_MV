%% Import of data
clc;
clear all;
close all;


data = importdata('data.txt');
time = data(:,1);
acc1 = data(:,3);
acc2 = data(:,4);
force = data(:,2);

%% Plot of the data as they are
figure
subplot(211)
plot(time, force)
leg = legend("External force",'Location','northeast');
xlabel('Acquisition time [s]')
ylabel('Force [N]')
set(gca,'FontSize',20);
subplot(212)
plot(time,acc1)
hold on
plot(time,acc2)
leg = legend("acc1","acc2",'Location','northeast');
xlabel('Acquisition time [s]')
ylabel('Aceleration [m * s^{-2}]')
set(gca,'FontSize',20);
hold off
%% Filtering the acceleration with a moving mean
acc1 = movmean(acc1,15);
acc2 = movmean(acc2,15);

%% Selecting the desired interval

initialt = 1.2;
finalt = 1.7;

indexinitial = find(time == initialt);
indexfinal = find(time == finalt);

time = time(indexinitial:indexfinal);
acc1 = acc1(indexinitial:indexfinal);
acc2 = acc2(indexinitial:indexfinal);
force = force(indexinitial:indexfinal);

%% Mean acceleration 
meanacc = (acc1 - acc2)/2;

figure_mean_acc = figure('Position', get(0, 'Screensize')); 
plot(time, meanacc,'LineWidth',1.2);
xlabel('Acquisition time [s]')
ylabel('Mean acceleration [m * s^{-2}]')
title {Mean acceleration as function of time}
set(gca,'FontSize',20);
grid on
saveas(figure_mean_acc, 'C:\\EDOARDO_01\\Università---------------------------\CORSI\VIBRATION\PROGETTO\RELAZIONE\IMMAGINI\MEAN_ACC.png','png');
%% Numerical transfer function
[ Tf , Fr ] = tfestimate(force, meanacc, [], [], [], 51200);
Tm = abs(Tf);

% ideal values of resonance frequency
fideal1 = 124.403; 
fideal2 = 497.612; 
fideal3 = 1119.63;

figure_TF = figure('Position', get(0, 'Screensize')); 
plot(Fr,Tm,'LineWidth',1.2);
xlim([0 2000]);
xlabel('Frequency [Hz]');
ylabel('Transfer function [(m * s^{-2})/N]');
title {Transfer function between force and acceleration};
set(gca,'FontSize',20);
grid on
hold on
% Extracting the values of the ideal TF
valuesTF = readmatrix('dataTF.txt');
Fraxis = 0:1:2000;
plot(Fraxis,valuesTF,'LineWidth',1.2);
xline(fideal1, "LineStyle","--");
xline(fideal2, "LineStyle","--");
xline(fideal3, "LineStyle","--");
leg = legend("Experimental TF","Ideal TF",'Location','east');

hold off
saveas(figure_TF, 'C:\\EDOARDO_01\\Università---------------------------\CORSI\VIBRATION\PROGETTO\RELAZIONE\IMMAGINI\COMPARE_TF.png','png');
%% Damping ratio estimate

dampideal1 = 0.05;
dampideal3 = 0.01;

% interpolating function of Tm
TF = @(x) interp1(Fr, abs(Tm), x, 'spline' );
TFplot = TF(Fr);

% first damping ratio
peak_freq_1 = max(TF(0:0.01:2*fideal1))/sqrt(2);
TFtmp = @(x) interp1(Fr, abs(Tm), x, 'spline' ) - peak_freq_1;

inters1 = fsolve(TFtmp,fideal1 - 5);
inters2 = fsolve(TFtmp,fideal1 + 100);

freqreal1 = (inters2 + inters1)/2;
dampreal1 = (inters2 - inters1)/(freqreal1*2);

% second damping ratio
peak_freq_2 = max(TF(fideal2-100:0.01:fideal2+100))/sqrt(2);
TFtmp = @(x) interp1(Fr, abs(Tm), x, 'spline' ) - peak_freq_2;

inters3 = fsolve(TFtmp,fideal2 - 20);
inters4 = fsolve(TFtmp,fideal2 + 5);

freqreal2 = (inters3 + inters4)/2;
dampreal2 = (inters4 - inters3)/(freqreal2*2);

% third damping ratio
peak_freq_3 = max(TF(fideal2:0.01:fideal3+300))/sqrt(2);
TFtmp = @(x) interp1(Fr, abs(Tm), x, 'spline' ) - peak_freq_3;

inters5 = fsolve(TFtmp,fideal3 - 50);
inters6 = fsolve(TFtmp,fideal3 + 5);

freqreal3 = (inters5 + inters6)/2;
dampreal3 = (inters6 - inters5)/(freqreal3*2);

%% Plot of half power method
fontdim = 18;
width = .2;

figure_HP = figure('Position', get(0, 'Screensize'));
plot(Fr,TFplot,'LineWidth',1.4);
hold on
xlim([000 1400]);

% First peak
yline(peak_freq_1, "LineStyle","--", "Color","Red","LineWidth",width);
xline(inters1, "LineStyle","--", "Color","Red","LineWidth",width);
xline(inters2, "LineStyle","--", "Color","Red","LineWidth",width);
plot(freqreal1,peak_freq_1*sqrt(2)-0.4,'.',"MarkerSize",20,"Color","Black");
plot(inters1,peak_freq_1,'.',"MarkerSize",20,"Color","Black");
plot(inters2,peak_freq_1,'.',"MarkerSize",20,"Color","Black");

A = "A";
text(inters1 - 20,peak_freq_1 + 1, A,'Color','red','FontSize',fontdim);
B = "B";
text(freqreal1 +7 ,peak_freq_1*sqrt(2) + 1.2, B,'Color','red','FontSize',fontdim);
C = "C";
text(inters2 + 10,peak_freq_1 + 1, C,'Color','red','FontSize',fontdim);
txtpos = [180 peak_freq_1*sqrt(2)+14];
strb ="A: " + "(" + num2str(inters1) + ", " + num2str(peak_freq_1) + ")";
text(txtpos(1),txtpos(2), strb,'Color','red','FontSize',fontdim);
strb ="B: " + "(" + num2str(freqreal1) + ", " + num2str(peak_freq_1*sqrt(2)) + ")";
text(txtpos(1),txtpos(2)-2, strb,'Color','red','FontSize',fontdim);
strb ="C: " + "(" + num2str(inters2) + ", " + num2str(peak_freq_1) + ")";
text(txtpos(1),txtpos(2)-4, strb,'Color','red','FontSize',fontdim);
strdfreq = "Freq = " + num2str(freqreal1) + " Hz";
text(txtpos(1),txtpos(2)-8, strdfreq,'Color','red','FontSize',fontdim);
strdamp = "\zeta = " + num2str(dampreal1);
text(txtpos(1),txtpos(2)-10, strdamp,'Color','red','FontSize',fontdim);

% Second peak
yline(peak_freq_2, "LineStyle","--", "Color","Black","LineWidth",width);
xline(inters3, "LineStyle","--", "Color","Black","LineWidth",width);
xline(inters4, "LineStyle","--", "Color","Black","LineWidth",width);
plot(freqreal2,peak_freq_2*sqrt(2),'.',"MarkerSize",20,"Color","Black");
plot(inters3,peak_freq_2,'.',"MarkerSize",20,"Color","Black");
plot(inters4,peak_freq_2,'.',"MarkerSize",20,"Color","Black");

A = "A";
text(inters3 - 30,peak_freq_2 + 1, A,'Color',"Black",'FontSize',fontdim);
B = "B";
text(freqreal2 -2 ,peak_freq_2*sqrt(2) + 1.2, B,'Color',"Black",'FontSize',fontdim);
C = "C";
text(inters4 + 20,peak_freq_2 + 1, C,'Color',"Black",'FontSize',fontdim);
txtpos = [550 peak_freq_1*sqrt(2)+4];
strb ="A: " + "(" + num2str(inters3) + ", " + num2str(peak_freq_2) + ")";
text(txtpos(1),txtpos(2), strb,'Color',"Black",'FontSize',fontdim);
strb ="B: " + "(" + num2str(freqreal2) + ", " + num2str(peak_freq_2*sqrt(2)) + ")";
text(txtpos(1),txtpos(2)-2, strb,'Color',"Black",'FontSize',fontdim);
strb ="C: " + "(" + num2str(inters4) + ", " + num2str(peak_freq_2) + ")";
text(txtpos(1),txtpos(2)-4, strb,'Color',"Black",'FontSize',fontdim);
strdfreq = "Freq = " + num2str(freqreal2) + " Hz";
text(txtpos(1),txtpos(2)-8, strdfreq,'Color',"Black",'FontSize',fontdim);
strdamp = "\zeta = " + num2str(dampreal2);
text(txtpos(1),txtpos(2)-10, strdamp,'Color',"Black",'FontSize',fontdim);

% Third peak
yline(peak_freq_3, "LineStyle","--", "Color",'#77AC30',"LineWidth",width);
xline(inters5, "LineStyle","--", "Color",'#77AC30',"LineWidth",width);
xline(inters6, "LineStyle","--", "Color",'#77AC30',"LineWidth",width);
plot(freqreal3,peak_freq_3*sqrt(2),'.',"MarkerSize",20,"Color","Black");
plot(inters5,peak_freq_3,'.',"MarkerSize",20,"Color","Black");
plot(inters6,peak_freq_3,'.',"MarkerSize",20,"Color","Black");



A = "A";
text(inters5 - 20,peak_freq_3 + 2, A,'Color','#77AC30','FontSize',fontdim);
B = "B";
text(freqreal3 -2 ,peak_freq_3*sqrt(2) + 1.2, B,'Color','#77AC30','FontSize',fontdim);
C = "C";
text(inters6 + 10,peak_freq_3 + 2, C,'Color','#77AC30','FontSize',fontdim);
txtpos = [700 peak_freq_3*sqrt(2)];
strb ="A: " + "(" + num2str(inters5) + ", " + num2str(peak_freq_3) + ")";
text(txtpos(1),txtpos(2), strb,'Color','#77AC30','FontSize',fontdim);
strb ="B: " + "(" + num2str(freqreal3) + ", " + num2str(peak_freq_3*sqrt(2)) + ")";
text(txtpos(1),txtpos(2)-2, strb,'Color','#77AC30','FontSize',fontdim);
strb ="C: " + "(" + num2str(inters6) + ", " + num2str(peak_freq_3) + ")";
text(txtpos(1),txtpos(2)-4, strb,'Color','#77AC30','FontSize',fontdim);
strdfreq = "Freq = " + num2str(freqreal3) + " Hz";
text(txtpos(1),txtpos(2)-8, strdfreq,'Color','#77AC30','FontSize',fontdim);
strdamp = "\zeta = " + num2str(dampreal3);
text(txtpos(1),txtpos(2)-10, strdamp,'Color','#77AC30','FontSize',fontdim);

xlabel('Frequency [Hz]');
ylabel('Transfer function [(m * s^{-2})/N]');
title {Half power points method};
set(gca,'FontSize',20);

leg = legend("Experimental TF",'Location','east');
hold off
saveas(figure_HP, 'C:\\EDOARDO_01\\Università---------------------------\CORSI\VIBRATION\PROGETTO\RELAZIONE\IMMAGINI\HALF_POWER.png','png');
%%
force = force - mean(force);
Y = fft(force);
L = length(force);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
Fs = 51200;
f = Fs*(0:(L/2))/L;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of the external force')
xlabel('Frequency (Hz)')
ylabel('|Force|')
set(gca,'FontSize',20);
xlim([0 2000]);