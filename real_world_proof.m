close all, clearvars, clc

load benchmark_signals_time_var.mat fiterr
file = fullfile("data/2022-07-05.json");
fid = fopen(file); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
val = jsondecode(str);
steps = val.steps.data;
hr = val.heart_rate.data;
TT.Time = seconds([hr.offset]);
TT.heart_rate = [hr.value]';
TT = timetable(TT.Time', TT.heart_rate, 'VariableNames',{'heart_rate'});
ST.Time = seconds([steps.offset]) + hours(2);
ST.steps = [steps.value]';

fc = 0.2;
Fs = seconds(TT.Time(2)-TT.Time(1));
[b,a] = butter(6,fc/(Fs/2));
butterFilt = TT.heart_rate;
butterFilt = filtfilt(b,a, butterFilt);
shr = butterFilt;
ref = array2timetable(shr, 'RowTimes', TT.Time, 'VariableNames',{'heart_rate'});

m_avg = smooth(TT.heart_rate, 25);
m_avg = round(m_avg);
rmse_avg = sqrt(nanmean((m_avg-ref.heart_rate).^2));

%%
wlen = 240;
kstd = 10;

[uHatSmooth, sigmaSmooth, stdSmooth, idxWindStarts] = bayesian_smoothing(TT.heart_rate, wlen, kstd, ...
    'noisecorrpar', fiterr.A, 'showplots', false);

uHatSmooth = round(uHatSmooth);
smooth = timetable(TT.Time, uHatSmooth, stdSmooth, VariableNames={'uHatSmooth', 'stdSmooth'});
%% 
TT.Time.Format = 'hh:mm:ss';
ax = figure('WindowState','maximized');
ax(1) = subplot(4,1,1:3);
plot(TT.Time, TT.heart_rate,'k.-','linewidth',1); hold on;
fill([TT.Time; flipud(TT.Time)], [(smooth.uHatSmooth)+stdSmooth;  flipud((smooth.uHatSmooth)-smooth.stdSmooth)],'y','FaceAlpha',0.3, 'EdgeAlpha',.1);
plot(TT.Time,smooth.uHatSmooth,'b-','LineWidth',1.5)
%plot(ref.Time, ref.heart_rate, 'r-', 'LineWidth',1.5)
%plot(TT.Time(idxWindStarts), TT.heart_rate(idxWindStarts), 'bd', 'LineWidth', 2)
%stem(ST.Time, ST.steps./10)
%title(['Bayesian smoothing with wlen: ', num2str(wlen), ', reconciliation: ', num2str(kstd)] )
legend('Raw signal', ...
    'CI (+/-SD)',...
    'Bayesian filter','Butterwort filter'    )%, ...
    %'Window start')
axis tight;

axes('position',[0.17,0.7,0.2,0.2])
box on % put box around new pair of axes
indexOfInterest = timerange('11:10:00','12:10:00'); % range of t near perturbation
plot(TT.Time(indexOfInterest), TT.heart_rate(TT.Time(indexOfInterest)),'k.-','linewidth',1); hold on;
fill([TT.Time(indexOfInterest); flipud(TT.Time(indexOfInterest))], [(smooth.uHatSmooth(TT.Time(indexOfInterest)))+smooth.stdSmooth(TT.Time(indexOfInterest));  
    flipud((smooth.uHatSmooth(TT.Time(indexOfInterest)))-smooth.stdSmooth(TT.Time(indexOfInterest)))],'y','FaceAlpha',0.3, 'EdgeAlpha',.1);
plot(TT.Time(indexOfInterest),smooth.uHatSmooth(indexOfInterest),'b-','LineWidth',1.5)
%plot(ref.Time(indexOfInterest), ref.heart_rate(indexOfInterest), 'r-', 'LineWidth',1.5)

axis tight;

ax(2) = subplot(4,1,4);
plot(TT.Time, sigmaSmooth','r','LineWidth',2)
title('Estimated noise variance')
setFonts(ax(2))
% xlabel('[MM/DD/YY]')
linkaxes(ax,'x')
rmse = sqrt(nanmean((smooth.uHatSmooth-ref.heart_rate).^2));
