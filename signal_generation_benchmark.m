% this code generates the ground truth signal and characterizes the HR noise 

close all, clearvars, clc


file = fullfile("data/2022-05-12.json");
fid = fopen(file); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
val = jsondecode(str);
hr = val.heart_rate.data;
offset_hr = [hr.offset];
value = [hr.value]';
offset = seconds(offset_hr(1):15:offset_hr(end));
timestamp = offset;
[~, idx] = intersect(offset, seconds(offset_hr), 'stable');
heart_rate = nan(length(offset), 1);
heart_rate(idx) = value;
TT = timetable(timestamp', heart_rate);

%% frequency check
L = length(value);
Y = fft(value);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
Fs = seconds(TT.Time(2)-TT.Time(1));
f = Fs*(0:(L/2))/L;
plot(f,P1) 
%xlim([0 0.4])
ylim([0 10])
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%% smoothing

% automatic smoothing with moving average
hr = smooth(value);
hr = hr(:);

% tt.heart_rate = hr(:).';
% tt.timestamp = TT.timestamp;

fc = 0.4;
[b,a] = butter(6,fc/(Fs/2));

%freqz(b,a,[],Fs)

butterFilt = hr;
butterFilt = filtfilt(b,a, butterFilt);

shr = nan(length(offset), 1);
shr(idx) = butterFilt;

ref = array2timetable(shr, 'RowTimes', TT.Time, 'VariableNames',{'heart_rate'});

errsm = butterFilt - hr;


%% noise characterization

% zero mean smoothing error
errsm = errsm - mean(errsm);
% find best AR order
% kmax = 160;
% AIC = zeros(kmax, 1);
% for k = 1:kmax
%     ar estimation
%     fiterr = ar(errsm, k);
%     AIC(k) = aic(fiterr);
% end
% figure
% plot(AIC)
% xlabel('ord')
% title AIC

% PACF
[pacf, lags] = parcorr(errsm, 'NumLags', 200);


%ord = find(AIC==min(AIC));
ord = 20; %152

% % Anderson test to find small enough order whith AIC close to minimum
% h = 0;
% ordtest = ord+1;
% while h == 0
%     ordtest = ordtest -1;
%     fiterr = ar(errsm, ordtest);
%     errfit = filter(fiterr.A, 1, errsm);
%     h = adtest(errfit, "Alpha",0.05);
% end
% ord = ordtest;

fiterr = ar(errsm, ord);

%% add generated noise

time_varying = true;
mu = 0;
rng(1)

if time_varying
    wlen = 720;
    n = ceil(L/wlen);
    noise = zeros(1, L);
    sigma = [4 3.1 2.7 0.4 1.3 1 0.5 1 1.7 2 1.6].*fiterr.NoiseVariance;
    %sigma = [.8 .6 .3 .1 .4 .6 .5 1 .4 1.3 .3];
    %sigma = fliplr(sigma);
    sigma = sigma(1:n);
    for i = 1:n
        l = wlen;
        idxmax = i*wlen;
        if idxmax > L
            idxmax = L;
            l = L-(i-1)*wlen;
        end
        noisein = sqrt(sigma(i)) * randn(1, l)+mu;
        noise((i-1)*wlen+1:idxmax) = filter(1, fiterr.A,noisein);
    end
else
    sigma = fiterr.NoiseVariance;
    noisein = sqrt(sigma) * randn(1, L)+mu;
    noise = filter(1, fiterr.A,noisein);
end

noise = noise(:);
%disp(['coefficienti AR (',num2str(ord),'): ', num2str(fiterr.A)])

errfit = filter(fiterr.A, 1, errsm);


noisyhr = butterFilt;
noisyhr = noisyhr + noise;
nhr = nan(length(offset), 1);
nhr(idx) = noisyhr;
noisy = array2timetable(nhr, 'RowTimes', TT.Time, 'VariableNames',{'heart_rate'});


ax=figure;
ref.Time.Format = 'hh:mm:ss';
TT.Time.Format = 'hh:mm:ss';
noisy.Time.Format = 'hh:mm:ss';
subplot(3,1,1)
plot(TT.Time, TT.heart_rate, 'k', 'Linewidth',1.5)
hold on
plot(ref.Time,ref.heart_rate, 'r', 'Linewidth',1.5)
legend("Raw", "Ground truth")
subplot(3,1,2)
plot(noisy.Time, noisy.heart_rate, 'Color', [0.39,0.83,0.07])
hold on
plot(ref.Time, ref.heart_rate, 'r', 'Linewidth',1.5)
legend("Noisy", "Ground truth")
%title 'Signal with added noise'
subplot(3,1,3)
sigmaplot = reshape(repmat(sigma, wlen,1),1,[]);
sigmaplot = sigmaplot(1:length(TT.Time));
plot(TT.Time, sigmaplot, 'r--')
legend("Sigma of generated noise")

%%
% if time_varying
%     save benchmark_signals_time_var.mat ref noisy fiterr sigma
% else
%     save benchmark_signals.mat ref noisy fiterr sigma
% end
