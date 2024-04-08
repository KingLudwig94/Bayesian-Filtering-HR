% this file performs the analysis of the BF against the MA on the generated signal
close all, clearvars, clc

load benchmark_signals_time_var.mat
noisy.heart_rate = round(noisy.heart_rate);
ref.heart_rate = round(ref.heart_rate);

%% smoothing moving average

m_avg = smooth(noisy.heart_rate, 25);
m_avg = round(m_avg);
ma = timetable(ref.Time, m_avg, VariableNames="m_avg");
%% bayesian

wlen = 240;
kstd = 10;
[uHatSmooth, sigmaSmooth, stdSmooth, idxWindStarts] = bayesian_smoothing(noisy.heart_rate, wlen, kstd, ...
    'noisecorrpar', fiterr.A, 'showplots', false);
uHatSmooth = round(uHatSmooth);
uHat = timetable(ref.Time, uHatSmooth, stdSmooth, 'VariableNames',{'uHat', 'std'});
save results_benchmark

%% figures
close all
load results_benchmark
%% bayes full
ref.Time.Format = 'hh:mm:ss';
figure('WindowState','maximized')
ax(1) = subplot(4,1,1:3);
plot(ref.Time, ref.heart_rate,'r-.','LineWidth',1.5); hold on;
fill([ref.Time; flipud(ref.Time)],[(uHatSmooth)+stdSmooth;  flipud((uHatSmooth)-stdSmooth)],'y','FaceAlpha',0.3, 'EdgeAlpha',.1);
%plot(ref.Time, noisy.heart_rate, 'g.-', 'linewidth', 1)
plot(ref.Time,uHatSmooth,'b.-','LineWidth',1.5)
plot(ref.Time(idxWindStarts), ref.heart_rate(idxWindStarts), 'bd', 'LineWidth', 2)
%title(['Bayesian smoothing with wlen: ', num2str(wlen), ', reconciliation: ', num2str(kstd)] )
legend('Reference', ...    
    'CI (+/-SD)',...
    'Bayesian')
%'Noisy', ...
axis tight;
ylim([30 160])

axes('position',[.18 .7 .2 .2])
box on % put box around new pair of axes
indexOfInterest = timerange('02:00:00','05:30:00'); % range of t near perturbation
plot(ref.Time(indexOfInterest), ref.heart_rate(ref.Time(indexOfInterest)),'r.-','linewidth',1.5); hold on;
fill([ref.Time(indexOfInterest); flipud(ref.Time(indexOfInterest))],[(uHat.uHat(indexOfInterest))+uHat.std(indexOfInterest); ...
flipud((uHat.uHat(indexOfInterest))-uHat.std(indexOfInterest))],'y','FaceAlpha',0.3, 'EdgeAlpha',.1);
%plot(ref.Time(indexOfInterest),noisy.heart_rate(noisy.Time(indexOfInterest)),'g.-','LineWidth',1)
plot(ref.Time(indexOfInterest),uHatSmooth(ismember(ref.Time,ref.Time(indexOfInterest))),'b.-','LineWidth',1)
axis tight;
ylim([40 65])


ax(2) = subplot(4,1,4);
plot(ref.Time, sigmaSmooth','r','LineWidth',2)
sigmas = sigma.*ones(720,1);
sigmas = sigmas(:);
hold on
plot(ref.Time, sigmas(1:length(ref.Time)),'r--')
%title('Estimated noise variance vs simulated noise')
% xlabel('[MM/DD/YY]')
legend('Estimated variance', 'Variance of generated noise')
linkaxes(ax,'x')


%% residuals
smt = timetable(ref.Time, uHatSmooth, VariableNames="uHatSmooth");

figure('WindowState','maximized')
err_b = sqrt((uHatSmooth-ref.heart_rate).^2);
err_avg = sqrt((m_avg-ref.heart_rate).^2);
rmse_b = sqrt(nanmean((uHatSmooth-ref.heart_rate).^2));
rmse_avg = sqrt(nanmean((m_avg-ref.heart_rate).^2));
mard_b = 100*nanmean(abs((uHatSmooth-ref.heart_rate)./ref.heart_rate));
mard_avg = 100*nanmean(abs((m_avg-ref.heart_rate)./ref.heart_rate));
mae_b = nanmean(abs((uHatSmooth-ref.heart_rate)));
mae_avg = nanmean(abs((m_avg-ref.heart_rate)));

indexOfInterest = timerange('05:30:00','23:59:00'); % range of t near perturbation
rmse_day = sqrt(nanmean((smt.uHatSmooth(indexOfInterest)-ref.heart_rate(indexOfInterest)).^2));
rmse_day_avg = sqrt(nanmean((ma.m_avg(indexOfInterest)-ref.heart_rate(indexOfInterest)).^2));

indexOfInterest = timerange('00:00:00','05:30:00'); % range of t near perturbation
rmse_night = sqrt(nanmean((smt.uHatSmooth(indexOfInterest)-ref.heart_rate(indexOfInterest)).^2));
rmse_night_avg = sqrt(nanmean((ma.m_avg(indexOfInterest)-ref.heart_rate(indexOfInterest)).^2));

ax = axes();
plot(ref.Time, err_b);
hold on
plot(ref.Time, err_avg)
legend({'Bayesian smoothing', 'Moving average'})
%title('Residuals')

%% boxplot
figure('WindowState','maximized')
ax = axes();
boxplot([err_b, err_avg], {'Bayesian' 'Moving_average'})


%% signals
figure('WindowState','maximized')
ax = axes();
plot(ref.Time, ref.heart_rate,'k.-','linewidth',1.5); hold on;
plot(ref.Time,noisy.heart_rate,'g.-','LineWidth',1)
ylabel('Heart rate [bpm]')
axis tight;
legend('Original', 'Noisy')

%% comparison
figure('WindowState','maximized')
ax = axes();
plot(ref.Time, ref.heart_rate,'r.-','linewidth',1.5); hold on;
plot(ref.Time,noisy.heart_rate,'.-','LineWidth',1, 'Color', [0.39,0.83,0.07])
plot(ref.Time, m_avg, 'c.-', 'linewidth', 1)
plot(ref.Time,uHatSmooth,'b.-','LineWidth',1)
ylabel('Heart rate [bpm]')
axis tight;
legend('Ground truth', 'Noisy', 'Moving average', 'Bayesian filtering')

axes('position',[0.07,0.75,0.2,0.2])
box on % put box around new pair of axes
indexOfInterest = timerange('05:30:00','07:30:00'); % range of t near perturbation
plot(ref.Time(indexOfInterest), ref.heart_rate(ref.Time(indexOfInterest)),'r.-','linewidth',1.5); hold on;
plot(ref.Time(indexOfInterest),noisy.heart_rate(noisy.Time(indexOfInterest)),'.-','LineWidth',1, 'Color', [0.39,0.83,0.07])
plot(ref.Time(indexOfInterest), m_avg(ismember(ref.Time,ref.Time(indexOfInterest))), 'c.-', 'linewidth', 1)
plot(ref.Time(indexOfInterest),uHatSmooth(ismember(ref.Time,ref.Time(indexOfInterest))),'b.-','LineWidth',1)
axis tight;

%% ttest

% [H, pval] = ttest2(err_b, err_avg)
% d = meanEffectSize(err_b, err_avg)
% figure('WindowState','maximized')
% gardnerAltmanPlot(err_b,err_avg,Effect="mediandiff");


%% 
% 
% clear mard_avg rmse_avg
% for i = 10:50
%     m_avg = smooth(noisy.heart_rate, i);
%     mard_avg = [mard_avg 100*nanmean(abs((m_avg-ref.heart_rate)./ref.heart_rate))];
%     rmse_avg = [rmse_avg sqrt(nanmean((m_avg-ref.heart_rate).^2))];
% end
% 
% figure
% plot(rmse_avg)
% hold on 
% plot(mard_avg)