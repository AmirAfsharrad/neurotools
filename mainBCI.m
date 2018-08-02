clear
clc
close

cd('./dataset')
load('subject_1.mat')
cd('..')
fs = 2400;

%% Loading Preprocssed Data
% clear
% load Data.mat

%% Loading Subject's Data

Data.exe.Arm.signal     = data.train{1};
Data.exe.Thumb.signal   = data.train{2};
Data.exe.Leg.signal     = data.train{3};
Data.exe.Idle.signal    = data.train{4};

Data.exe.test.signal = data.test;

% EEG Signals
FData.exe.Arm.signal     = Data.exe.Arm.signal(1:63,:,:);
FData.exe.Thumb.signal   = Data.exe.Thumb.signal(1:63,:,:);
FData.exe.Leg.signal     = Data.exe.Leg.signal(1:63,:,:);
FData.exe.Idle.signal    = Data.exe.Idle.signal(1:63,:,:);

% Start of Movement's Index
[FData.exe.Arm.index, ~]   = find(squeeze(Data.exe.Arm.signal(64,:,:)));
[FData.exe.Thumb.index, ~] = find(squeeze(Data.exe.Thumb.signal(64,:,:)));
[FData.exe.Leg.index, ~]   = find(squeeze(Data.exe.Leg.signal(64,:,:)));
[FData.exe.Idle.index, ~]  = find(squeeze(Data.exe.Idle.signal(64,:,:)));


FData.exe.test.signal    = Data.exe.test.signal(:,:,:);

clear Data

%% Removing DC Component

s_Arm   = size(FData.exe.Arm.signal);
s_Leg   = size(FData.exe.Leg.signal);
s_Thumb = size(FData.exe.Thumb.signal);
s_Idle  = size(FData.exe.Idle.signal);
s_Test  = size(FData.exe.test.signal);


mean_arm     = reshape(repmat(squeeze(mean(FData.exe.Arm.signal, 2)), [1, 1, s_Arm(2)]), [s_Arm(1), s_Arm(2), s_Arm(3)]);
mean_leg     = reshape(repmat(squeeze(mean(FData.exe.Leg.signal, 2)), [1, 1, s_Leg(2)]), [s_Leg(1), s_Leg(2), s_Leg(3)]);
mean_thumb   = reshape(repmat(squeeze(mean(FData.exe.Thumb.signal, 2)), [1, 1, s_Thumb(2)]), [s_Thumb(1), s_Thumb(2), s_Thumb(3)]);
mean_idle    = reshape(repmat(squeeze(mean(FData.exe.Idle.signal, 2)), [1, 1, s_Idle(2)]), [s_Idle(1), s_Idle(2), s_Idle(3)]);

mean_test = reshape(repmat(squeeze(mean(FData.exe.test.signal, 2)), [1, 1, s_Test(2)]), [s_Test(1), s_Test(2), s_Test(3)]);

FData.exe.Arm.signal   = FData.exe.Arm.signal   - mean_arm;
FData.exe.Leg.signal   = FData.exe.Leg.signal   - mean_leg;
FData.exe.Thumb.signal = FData.exe.Thumb.signal - mean_thumb;
FData.exe.Idle.signal  = FData.exe.Idle.signal  - mean_idle;      
FData.exe.test.signal  = FData.exe.test.signal  - mean_test;      

clear mean_arm mean_idle mean_leg mean_test mean_thumb data

%% Variance

FData.exe.Arm.var    =  squeeze(var(FData.exe.Arm.signal,0, 2));
FData.exe.Leg.var    =  squeeze(var(FData.exe.Leg.signal,0, 2));
FData.exe.Idle.var   =  squeeze(var(FData.exe.Idle.signal,0, 2));
FData.exe.Thumb.var  =  squeeze(var(FData.exe.Thumb.signal,0, 2));        
FData.exe.test.var   =  squeeze(var(FData.exe.test.signal,0, 2));


%% Skewness

FData.exe.Arm.skew    =  squeeze(skewness(FData.exe.Arm.signal,0, 2));
FData.exe.Leg.skew    =  squeeze(skewness(FData.exe.Leg.signal,0, 2));
FData.exe.Idle.skew   =  squeeze(skewness(FData.exe.Idle.signal,0, 2));
FData.exe.Thumb.skew  =  squeeze(skewness(FData.exe.Thumb.signal,0, 2));        
FData.exe.test.skew   =  squeeze(skewness(FData.exe.test.signal,0, 2));

%% Root Mean Square

FData.exe.Arm.RMS    =  squeeze(rms(FData.exe.Arm.signal,2));
FData.exe.Leg.RMS    =  squeeze(rms(FData.exe.Leg.signal,2));
FData.exe.Idle.RMS   =  squeeze(rms(FData.exe.Idle.signal,2));
FData.exe.Thumb.RMS  =  squeeze(rms(FData.exe.Thumb.signal,2));
FData.exe.test.RMS   =  squeeze(rms(FData.exe.test.signal,2));


%% Mean Freq

for channel = 1 : 63
    FData.exe.Arm.meanFreq(channel,:)    =  meanfreq(squeeze(FData.exe.Arm.signal(channel,:,:)), fs);
    FData.exe.Leg.meanFreq(channel,:)    =  meanfreq(squeeze(FData.exe.Leg.signal(channel,:,:)), fs);
    FData.exe.Idle.meanFreq(channel,:)   =  meanfreq(squeeze(FData.exe.Idle.signal(channel,:,:)), fs);
    FData.exe.Thumb.meanFreq(channel,:)  =  meanfreq(squeeze(FData.exe.Thumb.signal(channel,:,:)), fs);
    FData.exe.test.meanFreq(channel, :)  =  meanfreq(squeeze(FData.exe.test.signal(channel, :, :)), fs);
end

%% Med Freq

for channel = 1 : 63
    FData.exe.Arm.medFreq(channel,:)    =  medfreq(squeeze(FData.exe.Arm.signal(channel,:,:)), fs);
    FData.exe.Leg.medFreq(channel,:)    =  medfreq(squeeze(FData.exe.Leg.signal(channel,:,:)), fs);
    FData.exe.Idle.medFreq(channel,:)   =  medfreq(squeeze(FData.exe.Idle.signal(channel,:,:)), fs);
    FData.exe.Thumb.medFreq(channel,:)  =  medfreq(squeeze(FData.exe.Thumb.signal(channel,:,:)), fs);
    FData.exe.test.medFreq(channel, :)  =  medfreq(squeeze(FData.exe.test.signal(channel, :, :)), fs);
end

%% Mode Freq

for trials = 1 : 18
    for channels = 1 : 64
        Y = fft(FData.exe.Arm.signal(channels, :, trials));
        L = length(Y);
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        [~, I] = max(P1);
        f = 120*(0:(L/2))/L;
        FData.exe.Arm.modeFreq(channels, trials) = f(I);
        
        
        Y = fft(FData.exe.Leg.signal(channels, :, trials));
        L = length(Y);
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        [~, I] = max(P1);
        f = 120*(0:(L/2))/L;
        FData.exe.Leg.modeFreq(channels, trials) = f(I);
        
        
        Y = fft(FData.exe.Thumb.signal(channels, :, trials));
        L = length(Y);
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        [~, I] = max(P1);
        f = 120*(0:(L/2))/L;
        FData.exe.Thumb.modeFreq(channels, trials) = f(I);
        
        
        Y = fft(FData.exe.Idle.signal(channels, :, trials));
        L = length(Y);
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        [~, I] = max(P1);
        f = 120*(0:(L/2))/L;
        FData.exe.Idle.modeFreq(channels, trials) = f(I);
    
       
    end
end
for trials = 1 : Nexe
    for channels = 1 : 63
        Y = fft(FData.exe.test.signal(channels, :, trials));
        L = length(Y);
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        [~, I] = max(P1);
        f = 120*(0:(L/2))/L;
        FData.exe.test.modeFreq(channels, trials) = f(I);
    end
end


clear Y f L P1 P2 I


%% Sine Transform
for trials = 1 : 20
    for channels = 1 : 64
        FData.exe.Arm.DST  (channels, :, trials)  =  dst(FData.exe.Arm.signal(channels, :, trials));
        FData.exe.Leg.DST  (channels, :, trials)  =  dst(FData.exe.Leg.signal(channels, :, trials));
        FData.exe.Idle.DST (channels, :, trials)  =  dst(FData.exe.Idle.signal(channels, :, trials));
        FData.exe.Thumb.DST(channels, :, trials)  =  dst(FData.exe.Thumb.signal(channels, :, trials));
        
        FData.img.Arm.DST  (channels, :, trials)  =  dst(FData.img.Arm.signal(channels, :, trials));
        FData.img.Leg.DST  (channels, :, trials)  =  dst(FData.img.Leg.signal(channels, :, trials));
        FData.img.Thumb.DST(channels, :, trials)  =  dst(FData.img.Thumb.signal(channels, :, trials));
    end
end
for trials = 1 : Nexe
    for channels = 1 : 64
        FData.exe.test.DST(channels, :, trials)  =  dst(FData.exe.test.signal(channels, :, trials));
    end
end


for trials = 1 : Nimg
    for channels = 1 : 64
        FData.img.test.DST(channels, :, trials)  =  dst(FData.img.test.signal(channels, :, trials));
    end
end
%% Cosine Transform
for trials = 1 : 20
    for channels = 1 : 64
        FData.exe.Arm.DCT  (channels, :, trials)  =  dct(FData.exe.Arm.signal(channels, :, trials));
        FData.exe.Leg.DCT  (channels, :, trials)  =  dct(FData.exe.Leg.signal(channels, :, trials));
        FData.exe.Idle.DCT (channels, :, trials)  =  dct(FData.exe.Idle.signal(channels, :, trials));
        FData.exe.Thumb.DCT(channels, :, trials)  =  dct(FData.exe.Thumb.signal(channels, :, trials));
        
        FData.img.Arm.DCT  (channels, :, trials)  =  dct(FData.img.Arm.signal(channels, :, trials));
        FData.img.Leg.DCT  (channels, :, trials)  =  dct(FData.img.Leg.signal(channels, :, trials));
        FData.img.Thumb.DCT(channels, :, trials)  =  dct(FData.img.Thumb.signal(channels, :, trials));
    end
end
for trials = 1 : Nexe
    for channels = 1 : 64
        FData.exe.test.DCT(channels, :, trials)  =  dct(FData.exe.test.signal(channels, :, trials));
    end
end


for trials = 1 : Nimg
    for channels = 1 : 64
        FData.img.test.DCT(channels, :, trials)  =  dct(FData.img.test.signal(channels, :, trials));
    end
end
%% Band Pass Filters
% Frequency Bands According to the HW :))
h_alpha = BPF(360, 7.5  , 13.5, 120);
h_beta =  BPF(360, 13.5 , 20, 120);
h_theta = BPF(360, 3.5  ,  7.5, 120);
h_delta = BPF(360, eps,  3.5, 120);
close all 
%% Filtering Frequency Bands
for trials = 1 : 20
    for channels = 1 : 64
        FData.exe.Arm.alpha_band(channels,:,trials)  = doFilt(h_alpha, FData.exe.Arm.signal(channels, :, trials));
        FData.exe.Arm.beta_band(channels,:, trials)  = doFilt(h_beta,  FData.exe.Arm.signal(channels, :, trials));
        FData.exe.Arm.theta_band(channels,:, trials) = doFilt(h_theta, FData.exe.Arm.signal(channels, :, trials));
        FData.exe.Arm.delta_band(channels,:, trials) = doFilt(h_delta, FData.exe.Arm.signal(channels, :, trials));
        
        FData.exe.Leg.alpha_band(channels,:,trials)  = doFilt(h_alpha, FData.exe.Leg.signal(channels, :, trials));
        FData.exe.Leg.beta_band(channels,:, trials)  = doFilt(h_beta,  FData.exe.Leg.signal(channels, :, trials));
        FData.exe.Leg.theta_band(channels,:, trials) = doFilt(h_theta, FData.exe.Leg.signal(channels, :, trials));
        FData.exe.Leg.delta_band(channels,:, trials) = doFilt(h_delta, FData.exe.Leg.signal(channels, :, trials));
        
        FData.exe.Thumb.alpha_band(channels,:,trials)  = doFilt(h_alpha, FData.exe.Thumb.signal(channels, :, trials));
        FData.exe.Thumb.beta_band(channels,:, trials)  = doFilt(h_beta,  FData.exe.Thumb.signal(channels, :, trials));
        FData.exe.Thumb.theta_band(channels,:, trials) = doFilt(h_theta, FData.exe.Thumb.signal(channels, :, trials));
        FData.exe.Thumb.delta_band(channels,:, trials) = doFilt(h_delta, FData.exe.Thumb.signal(channels, :, trials));
        
        FData.exe.Idle.alpha_band(channels,:,trials)  = doFilt(h_alpha, FData.exe.Idle.signal(channels, :, trials));
        FData.exe.Idle.beta_band(channels,:, trials)  = doFilt(h_beta,  FData.exe.Idle.signal(channels, :, trials));
        FData.exe.Idle.theta_band(channels,:, trials) = doFilt(h_theta, FData.exe.Idle.signal(channels, :, trials));
        FData.exe.Idle.delta_band(channels,:, trials) = doFilt(h_delta, FData.exe.Idle.signal(channels, :, trials));
        
        
        FData.img.Arm.alpha_band(channels,:,trials)  = doFilt(h_alpha, FData.img.Arm.signal(channels, :, trials));
        FData.img.Arm.beta_band(channels,:, trials)  = doFilt(h_beta,  FData.img.Arm.signal(channels, :, trials));
        FData.img.Arm.theta_band(channels,:, trials) = doFilt(h_theta, FData.img.Arm.signal(channels, :, trials));
        FData.img.Arm.delta_band(channels,:, trials) = doFilt(h_delta, FData.img.Arm.signal(channels, :, trials));
        
        FData.img.Leg.alpha_band(channels,:,trials)  = doFilt(h_alpha, FData.img.Leg.signal(channels, :, trials));
        FData.img.Leg.beta_band(channels,:, trials)  = doFilt(h_beta,  FData.img.Leg.signal(channels, :, trials));
        FData.img.Leg.theta_band(channels,:, trials) = doFilt(h_theta, FData.img.Leg.signal(channels, :, trials));
        FData.img.Leg.delta_band(channels,:, trials) = doFilt(h_delta, FData.img.Leg.signal(channels, :, trials));
        
        FData.img.Thumb.alpha_band(channels,:,trials)  = doFilt(h_alpha, FData.img.Thumb.signal(channels, :, trials));
        FData.img.Thumb.beta_band(channels,:, trials)  = doFilt(h_beta,  FData.img.Thumb.signal(channels, :, trials));
        FData.img.Thumb.theta_band(channels,:, trials) = doFilt(h_theta, FData.img.Thumb.signal(channels, :, trials));
        FData.img.Thumb.delta_band(channels,:, trials) = doFilt(h_delta, FData.img.Thumb.signal(channels, :, trials));
        
    end
end

for trials = 1 : Nexe
    for channels = 1 : 64
        FData.exe.test.alpha_band(channels,:, trials) = doFilt(h_alpha, FData.exe.test.signal(channels, :, trials));
        FData.exe.test.beta_band(channels,:, trials) = doFilt(h_beta, FData.exe.test.signal(channels, :, trials));
        FData.exe.test.theta_band(channels,:, trials) = doFilt(h_theta, FData.exe.test.signal(channels, :, trials));
        FData.exe.test.delta_band(channels,:, trials) = doFilt(h_delta, FData.exe.test.signal(channels, :, trials));
    end
end


for trials = 1 : Nimg
    for channels = 1 : 64
        FData.img.test.alpha_band(channels,:, trials) = doFilt(h_alpha, FData.img.test.signal(channels, :, trials));
        FData.img.test.beta_band(channels,:, trials) = doFilt(h_beta, FData.img.test.signal(channels, :, trials));
        FData.img.test.theta_band(channels,:, trials) = doFilt(h_theta, FData.img.test.signal(channels, :, trials));
        FData.img.test.delta_band(channels,:, trials) = doFilt(h_delta, FData.img.test.signal(channels, :, trials));
    end
end

clear h_alpha h_beta h_delta h_theta

%% Checking The Procedure
figure
plotFFT(Data(1).img.leg(3,:,16), 120, 0, 60, 'FFT','|fft|',10);
hold on
plotFFT(FData.img.Leg.alpha_band(3,:,16), 120, 0, 60, 'FFT of img.Idle Channel 3 Trial 16 ','',10);
plotFFT(FData.img.Leg.beta_band(3,:,16), 120, 0, 60, 'FFT of img.Idle Channel 3 Trial 16 ','',10);
plotFFT(FData.img.Leg.theta_band(3,:,16), 120, 0, 60, 'FFT of img.Idle Channel 3 Trial 16 ','',10);
plotFFT(FData.img.Leg.delta_band(3,:,16), 120, 0, 60, 'FFT of img.Idle Channel 3 Trial 16 ','',10);
legend('full','alpha','beta','theta','delta');

figure
plotFFT(Data(1).exe.test(2,:,16), 120, 0, 60, 'FFT','|fft|',10);
hold on
plotFFT(FData.exe.test.alpha_band(2,:,16), 120, 0, 60, 'FFT of exe.test Channel 2 Trial 16 ','',10);
plotFFT(FData.exe.test.beta_band(2,:,16), 120, 0, 60, 'FFT of exe.test Channel 2 Trial 16 ','',10);
plotFFT(FData.exe.test.theta_band(2,:,16), 120, 0, 60, 'FFT of exe.test Channel 2 Trial 16 ','',10);
plotFFT(FData.exe.test.delta_band(2,:,16), 120, 0, 60, 'FFT of exe.test Channel 2 Trial 16 ','',10);
legend('full','alpha','beta','theta','delta');

%% Freq Band Energy
for trials = 1 : 20
    for channels = 1 : 64
        % EXE Data
        FData.exe.Arm.alphaEnergy(channels, trials) = norm(FData.exe.Arm.alpha_band(channels, :, trials)).^2;
        FData.exe.Arm.betaEnergy(channels, trials) = norm(FData.exe.Arm.beta_band(channels, :, trials)).^2;
        FData.exe.Arm.thetaEnergy(channels, trials) = norm(FData.exe.Arm.theta_band(channels, :, trials)).^2;
        FData.exe.Arm.deltaEnergy(channels, trials) = norm(FData.exe.Arm.delta_band(channels, :, trials)).^2;
        
        FData.exe.Leg.alphaEnergy(channels, trials) = norm(FData.exe.Leg.alpha_band(channels, :, trials)).^2;
        FData.exe.Leg.betaEnergy(channels, trials) = norm(FData.exe.Leg.beta_band(channels, :, trials)).^2;
        FData.exe.Leg.thetaEnergy(channels, trials) = norm(FData.exe.Leg.theta_band(channels, :, trials)).^2;
        FData.exe.Leg.deltaEnergy(channels, trials) = norm(FData.exe.Leg.delta_band(channels, :, trials)).^2;
        
        FData.exe.Thumb.alphaEnergy(channels, trials) = norm(FData.exe.Thumb.alpha_band(channels, :, trials)).^2;
        FData.exe.Thumb.betaEnergy(channels, trials) = norm(FData.exe.Thumb.beta_band(channels, :, trials)).^2;
        FData.exe.Thumb.thetaEnergy(channels, trials) = norm(FData.exe.Thumb.theta_band(channels, :, trials)).^2;
        FData.exe.Thumb.deltaEnergy(channels, trials) = norm(FData.exe.Thumb.delta_band(channels, :, trials)).^2;
        
        FData.exe.Idle.alphaEnergy(channels, trials) = norm(FData.exe.Idle.alpha_band(channels, :, trials)).^2;
        FData.exe.Idle.betaEnergy(channels,  trials) = norm(FData.exe.Idle.beta_band(channels, :, trials)).^2;
        FData.exe.Idle.thetaEnergy(channels, trials) = norm(FData.exe.Idle.theta_band(channels, :, trials)).^2;
        FData.exe.Idle.deltaEnergy(channels, trials) = norm(FData.exe.Idle.delta_band(channels, :, trials)).^2;
        
        % IMG Data
        
        FData.img.Arm.alphaEnergy(channels, trials) = norm(FData.img.Arm.alpha_band(channels, :, trials)).^2;
        FData.img.Arm.betaEnergy(channels, trials)  = norm(FData.img.Arm.beta_band(channels, :, trials)).^2;
        FData.img.Arm.thetaEnergy(channels, trials) = norm(FData.img.Arm.theta_band(channels, :, trials)).^2;
        FData.img.Arm.deltaEnergy(channels, trials) = norm(FData.img.Arm.delta_band(channels, :, trials)).^2;
        
        FData.img.Leg.alphaEnergy(channels, trials) = norm(FData.img.Leg.alpha_band(channels, :, trials)).^2;
        FData.img.Leg.betaEnergy(channels, trials)  = norm(FData.img.Leg.beta_band(channels, :, trials)).^2;
        FData.img.Leg.thetaEnergy(channels, trials) = norm(FData.img.Leg.theta_band(channels, :, trials)).^2;
        FData.img.Leg.deltaEnergy(channels, trials) = norm(FData.img.Leg.delta_band(channels, :, trials)).^2;
        
        FData.img.Thumb.alphaEnergy(channels, trials) = norm(FData.img.Thumb.alpha_band(channels, :, trials)).^2;
        FData.img.Thumb.betaEnergy(channels, trials)  = norm(FData.img.Thumb.beta_band(channels, :, trials)).^2;
        FData.img.Thumb.thetaEnergy(channels, trials) = norm(FData.img.Thumb.theta_band(channels, :, trials)).^2;
        FData.img.Thumb.deltaEnergy(channels, trials) = norm(FData.img.Thumb.delta_band(channels, :, trials)).^2;
        
    end
end

for trials = 1 : Nexe
    for channels = 1 : 64
        FData.exe.test.alphaEnergy(channels, trials) = norm(FData.exe.test.alpha_band(channels, :, trials)).^2;
        FData.exe.test.betaEnergy(channels, trials)  = norm(FData.exe.test.beta_band(channels, :, trials)).^2;
        FData.exe.test.thetaEnergy(channels, trials) = norm(FData.exe.test.theta_band(channels, :, trials)).^2;
        FData.exe.test.deltaEnergy(channels, trials) = norm(FData.exe.test.delta_band(channels, :, trials)).^2;
    end
end


for trials = 1 : Nimg
    for channels = 1 : 64
        FData.img.test.alphaEnergy(channels, trials) = norm(FData.img.test.alpha_band(channels, :, trials)).^2;
        FData.img.test.betaEnergy(channels, trials)  = norm(FData.img.test.beta_band(channels, :, trials)).^2;
        FData.img.test.thetaEnergy(channels, trials) = norm(FData.img.test.theta_band(channels, :, trials)).^2;
        FData.img.test.deltaEnergy(channels, trials) = norm(FData.img.test.delta_band(channels, :, trials)).^2;
    end
end

%% STFT
% clear FData
% load('FData.mat')
for trials = 1 : 20
    for channels = 1 : 64
        FData.exe.Arm.STFT  (channels, :, :, trials)  =  abs(spectrogram(FData.exe.Arm.signal(channels, :, trials),10,0,15,120));
        FData.exe.Leg.STFT  (channels, :, :, trials)  =  abs(spectrogram(FData.exe.Leg.signal(channels, :, trials),10,0,15,120));
        FData.exe.Idle.STFT (channels, :, :, trials)  =  abs(spectrogram(FData.exe.Idle.signal(channels, :, trials),10,0,15,120));
        FData.exe.Thumb.STFT(channels, :, :, trials)  =  abs(spectrogram(FData.exe.Thumb.signal(channels, :, trials),10,0,15,120));
        
        FData.img.Arm.STFT (channels, :, :, trials)  =  abs(spectrogram(FData.img.Arm.signal(channels, :, trials),10,0,15,120));
        FData.img.Leg.STFT  (channels, :, :, trials)  =  abs(spectrogram(FData.img.Leg.signal(channels, :, trials),10,0,15,120));
        FData.img.Thumb.STFT(channels, :, :, trials)  =  abs(spectrogram(FData.img.Thumb.signal(channels, :, trials),10,0,15,120));
    end
end

for trials = 1 : Nexe
    for channels = 1 : 64
        FData.exe.test.STFT(channels, :, :, trials)  =  abs(spectrogram(FData.exe.test.signal(channels, :, trials),10,0,15,120));
    end
end


for trials = 1 : Nimg
    for channels = 1 : 64
        FData.img.test.STFT(channels, :, :, trials)  =  abs(spectrogram(FData.img.test.signal(channels, :, trials),10,0,15,120));
    end
end

%% Alpha and Beta Band Filtering
h = BPF(360, 7  , 30, 120);

for trials = 1 : 20
    for channels = 1 : 64
        FData.exe.Arm.mainbands(channels,:,trials)    = doFilt(h, FData.exe.Arm.signal(channels, :, trials));
        FData.exe.Leg.mainbands(channels,:,trials)    = doFilt(h, FData.exe.Leg.signal(channels, :, trials));
        FData.exe.Thumb.mainbands(channels,:,trials)  = doFilt(h, FData.exe.Thumb.signal(channels, :, trials));
        FData.exe.Idle.mainbands(channels,:,trials)   = doFilt(h, FData.exe.Idle.signal(channels, :, trials));
        
        FData.img.Arm.mainbands(channels,:,trials)    = doFilt(h, FData.img.Arm.signal(channels, :, trials));
        FData.img.Leg.mainbands(channels,:,trials)    = doFilt(h, FData.img.Leg.signal(channels, :, trials));
        FData.img.Thumb.mainbands(channels,:,trials)  = doFilt(h, FData.img.Thumb.signal(channels, :, trials));
    end
end

for trials = 1 : Nexe
    for channels = 1 : 64
        FData.exe.test.mainbands(channels,:,trials)    = doFilt(h, FData.exe.test.signal(channels, :, trials));
    end
end


for trials = 1 : Nimg
    for channels = 1 : 64
        FData.img.test.mainbands(channels,:,trials)    = doFilt(h, FData.img.test.signal(channels, :, trials));
    end
end


clear h
%%
figure
hold on
plotFFT(FData.exe.Leg.mainbands(11,:,12), 120, 0, 60, '', '', 10);
plotFFT(FData.exe.Leg.signal(11,:,12), 120, 0, 60, '', '', 10);
legend('full','\alpha and \beta');
%% Common Spatial Patterns
FData.exe.covMat(1, :, :) = cov(mean(FData.exe.Arm.mainbands, 3)');
FData.exe.covMat(2, :, :) = cov(mean(FData.exe.Leg.mainbands, 3)');
FData.exe.covMat(3, :, :) = cov(mean(FData.exe.Thumb.mainbands, 3)');
FData.exe.covMat(4, :, :) = cov(mean(FData.exe.Idle.mainbands, 3)');



FData.img.covMat(1, :, :) = cov(mean(FData.img.Arm.mainbands, 3)');
FData.img.covMat(2, :, :) = cov(mean(FData.img.Leg.mainbands, 3)');
FData.img.covMat(3, :, :) = cov(mean(FData.img.Thumb.mainbands, 3)');


spatialFilter_EXE = MulticlassCSP(squeeze(FData.exe.covMat), 2);
spatialFilter_IMG = MulticlassCSP(squeeze(FData.img.covMat), 2);


for trials = 1 : 20
    FData.exe.Arm.CSP1(:, trials)    = squeeze(spatialFilter_EXE(1,:))*squeeze(FData.exe.Arm.signal(:,:,trials));
    FData.exe.Arm.CSP2(:, trials)    = squeeze(spatialFilter_EXE(2,:))*squeeze(FData.exe.Arm.signal(:,:,trials));
    FData.exe.Leg.CSP1(:, trials)    = squeeze(spatialFilter_EXE(1,:))*squeeze(FData.exe.Leg.signal(:,:,trials));
    FData.exe.Leg.CSP2(:, trials)    = squeeze(spatialFilter_EXE(2,:))*squeeze(FData.exe.Leg.signal(:,:,trials));
    FData.exe.Idle.CSP1(:, trials)   = squeeze(spatialFilter_EXE(1,:))*squeeze(FData.exe.Idle.signal(:,:,trials));
    FData.exe.Idle.CSP2(:, trials)   = squeeze(spatialFilter_EXE(2,:))*squeeze(FData.exe.Idle.signal(:,:,trials));
    FData.exe.Thumb.CSP1(:, trials)  = squeeze(spatialFilter_EXE(1,:))*squeeze(FData.exe.Thumb.signal(:,:,trials));
    FData.exe.Thumb.CSP2(:, trials)  = squeeze(spatialFilter_EXE(2,:))*squeeze(FData.exe.Thumb.signal(:,:,trials));
    
    
    FData.img.Arm.CSP1(:, trials)    = squeeze(spatialFilter_IMG(1,:))*squeeze(FData.img.Arm.signal(:,:,trials));
    FData.img.Arm.CSP2(:, trials)    = squeeze(spatialFilter_IMG(2,:))*squeeze(FData.img.Arm.signal(:,:,trials));
    FData.img.Leg.CSP1(:, trials)    = squeeze(spatialFilter_IMG(1,:))*squeeze(FData.img.Leg.signal(:,:,trials));
    FData.img.Leg.CSP2(:, trials)    = squeeze(spatialFilter_IMG(2,:))*squeeze(FData.img.Leg.signal(:,:,trials));
    FData.img.Thumb.CSP1(:, trials)  = squeeze(spatialFilter_IMG(1,:))*squeeze(FData.img.Thumb.signal(:,:,trials));
    FData.img.Thumb.CSP2(:, trials)  = squeeze(spatialFilter_IMG(2,:))*squeeze(FData.img.Thumb.signal(:,:,trials));
end

for trials = 1 : Nexe
    for channels = 1 : 64
        FData.exe.test.CSP1(:, trials)  = squeeze(spatialFilter_EXE(1,:))*squeeze(FData.exe.test.signal(:,:,trials));
        FData.exe.test.CSP2(:, trials)  = squeeze(spatialFilter_EXE(2,:))*squeeze(FData.exe.test.signal(:,:,trials));
    end
end


for trials = 1 : Nimg
    for channels = 1 : 64
        FData.img.test.CSP1(:, trials)  = squeeze(spatialFilter_EXE(1,:))*squeeze(FData.img.test.signal(:,:,trials));
        FData.img.test.CSP2(:, trials)  = squeeze(spatialFilter_EXE(2,:))*squeeze(FData.img.test.signal(:,:,trials));
    end
end
%% Checking CSP

figure
hold on
plot(FData.exe.Arm.CSP1(:,1))
plot(FData.exe.Arm.signal(1:3,:,1)')
xlabel('t');
legend('CSP-Filtered Signal','Original Signal Channel 01',...
    'Original Signal Channel 02','Original Signal Channel 03')

figure
hold on
plotFFT(FData.exe.Arm.CSP1(:,10),120, 0, 60, '', '|fft|', 10)
plotFFT(spatialFilter_EXE(1,:)*FData.exe.Arm.mainbands(:,:,10), 120, 0, 60,'','|fft|',10)
plotFFT(FData.exe.Arm.signal(1,:,10), 120, 0, 60, '','|fft|',10)
legend('a Sample CSP','a Sample CSP(\alpha and \beta)','Original Signal')

%% CSP Band Energy

for trials = 1 : 20
    FData.exe.Arm.CSP1_alpha  (trials) = bandpower(FData.exe.Arm.CSP1(:,trials),120,[7.5,13.5]);
    FData.exe.Leg.CSP1_alpha (trials)  = bandpower(FData.exe.Leg.CSP1(:,trials),120,[7.5,13.5]);
    FData.exe.Thumb.CSP1_alpha(trials) = bandpower(FData.exe.Thumb.CSP1(:,trials),120,[7.5,13.5]);
    FData.exe.Idle.CSP1_alpha (trials) = bandpower(FData.exe.Idle.CSP1(:,trials),120,[7.5,13.5]);
    
    FData.exe.Arm.CSP2_alpha(trials)   = bandpower(FData.exe.Arm.CSP2(:,trials),120,[7.5,13.5]);
    FData.exe.Leg.CSP2_alpha(trials)   = bandpower(FData.exe.Leg.CSP2(:,trials),120,[7.5,13.5]);
    FData.exe.Thumb.CSP2_alpha(trials) = bandpower(FData.exe.Thumb.CSP2(:,trials),120,[7.5,13.5]);
    FData.exe.Idle.CSP2_alpha(trials)  = bandpower(FData.exe.Idle.CSP2(:,trials),120,[7.5,13.5]);
    
    
    FData.exe.Arm.CSP1_beta  (trials) = bandpower(FData.exe.Arm.CSP1(:,trials),120,[13.5,20]);
    FData.exe.Leg.CSP1_beta (trials)  = bandpower(FData.exe.Leg.CSP1(:,trials),120,[13.5,20]);
    FData.exe.Thumb.CSP1_beta(trials) = bandpower(FData.exe.Thumb.CSP1(:,trials),120,[13.5,20]);
    FData.exe.Idle.CSP1_beta (trials) = bandpower(FData.exe.Idle.CSP1(:,trials),120,[13.5,20]);
    
    FData.exe.Arm.CSP2_beta(trials)   = bandpower(FData.exe.Arm.CSP2(:,trials),120,[13.5,20]);
    FData.exe.Leg.CSP2_beta(trials)   = bandpower(FData.exe.Leg.CSP2(:,trials),120,[13.5,20]);
    FData.exe.Thumb.CSP2_beta(trials) = bandpower(FData.exe.Thumb.CSP2(:,trials),120,[13.5,20]);
    FData.exe.Idle.CSP2_beta(trials)  = bandpower(FData.exe.Idle.CSP2(:,trials),120,[13.5,20]);



    FData.img.Arm.CSP1_beta  (trials) = bandpower(FData.img.Arm.CSP1(:,trials),120,[13.5,20]);
    FData.img.Leg.CSP1_beta (trials)  = bandpower(FData.img.Leg.CSP1(:,trials),120,[13.5,20]);
    FData.img.Thumb.CSP1_beta(trials) = bandpower(FData.img.Thumb.CSP1(:,trials),120,[13.5,20]);

    FData.img.Arm.CSP2_beta(trials)   = bandpower(FData.img.Arm.CSP2(:,trials),120,[13.5,20]);
    FData.img.Leg.CSP2_beta(trials)   = bandpower(FData.img.Leg.CSP2(:,trials),120,[13.5,20]);
    FData.img.Thumb.CSP2_beta(trials) = bandpower(FData.img.Thumb.CSP2(:,trials),120,[13.5,20]);

    
    FData.img.Arm.CSP1_alpha  (trials) = bandpower(FData.img.Arm.CSP1(:,trials),120,[7.5,13.5]);
    FData.img.Leg.CSP1_alpha (trials)  = bandpower(FData.img.Leg.CSP1(:,trials),120,[7.5,13.5]);
    FData.img.Thumb.CSP1_alpha(trials) = bandpower(FData.img.Thumb.CSP1(:,trials),120,[7.5,13.5]);

    FData.img.Arm.CSP2_alpha(trials)   = bandpower(FData.img.Arm.CSP2(:,trials),120,[7.5,13.5]);
    FData.img.Leg.CSP2_alpha(trials)   = bandpower(FData.img.Leg.CSP2(:,trials),120,[7.5,13.5]);
    FData.img.Thumb.CSP2_alpha(trials) = bandpower(FData.img.Thumb.CSP2(:,trials),120,[7.5,13.5]);    
    
end

for trials = 1 : Nexe
	FData.exe.test.CSP1_alpha  (trials) = bandpower(FData.exe.test.CSP1(:,trials),120,[7.5,13.5]);
    FData.exe.test.CSP2_alpha (trials)  = bandpower(FData.exe.test.CSP2(:,trials),120,[7.5,13.5]);
	FData.exe.test.CSP1_beta  (trials) = bandpower(FData.exe.test.CSP1(:,trials),120,[13.5,20]);
    FData.exe.test.CSP2_beta (trials)  = bandpower(FData.exe.test.CSP2(:,trials),120,[13.5,20]);
end


for trials = 1 : Nimg
	FData.img.test.CSP1_alpha  (trials) = bandpower(FData.img.test.CSP1(:,trials),120,[7.5,13.5]);
    FData.img.test.CSP2_alpha (trials)  = bandpower(FData.img.test.CSP2(:,trials),120,[7.5,13.5]);
	FData.img.test.CSP1_beta  (trials) = bandpower(FData.img.test.CSP1(:,trials),120,[13.5,20]);
    FData.img.test.CSP2_beta (trials)  = bandpower(FData.img.test.CSP2(:,trials),120,[13.5,20]);
end


%% CSP Mean Freq

for trials = 1 : 20
    FData.exe.Arm.CSP1_meanFreq  (trials)  =  meanfreq(FData.exe.Arm.CSP1(:, trials), 120);
    FData.exe.Leg.CSP1_meanFreq  (trials)  =  meanfreq(FData.exe.Leg.CSP1(:, trials), 120);
    FData.exe.Idle.CSP1_meanFreq (trials)  =  meanfreq(FData.exe.Idle.CSP1(:, trials), 120);
    FData.exe.Thumb.CSP1_meanFreq(trials)  =  meanfreq(FData.exe.Thumb.CSP1(:, trials), 120);

    FData.img.Arm.CSP1_meanFreq  (trials)  =  meanfreq(FData.img.Arm.CSP1(:, trials), 120);
    FData.img.Leg.CSP1_meanFreq  (trials)  =  meanfreq(FData.img.Leg.CSP1(:, trials), 120);
    FData.img.Thumb.CSP1_meanFreq(trials)  =  meanfreq(FData.img.Thumb.CSP1(:, trials), 120);
 
    
    FData.exe.Arm.CSP2_meanFreq  (trials)  =  meanfreq(FData.exe.Arm.CSP2(:, trials), 120);
    FData.exe.Leg.CSP2_meanFreq  (trials)  =  meanfreq(FData.exe.Leg.CSP2(:, trials), 120);
    FData.exe.Idle.CSP2_meanFreq (trials)  =  meanfreq(FData.exe.Idle.CSP2(:, trials), 120);
    FData.exe.Thumb.CSP2_meanFreq(trials)  =  meanfreq(FData.exe.Thumb.CSP2(:, trials), 120);

    FData.img.Arm.CSP2_meanFreq  (trials)  =  meanfreq(FData.img.Arm.CSP2(:, trials), 120);
    FData.img.Leg.CSP2_meanFreq  (trials)  =  meanfreq(FData.img.Leg.CSP2(:, trials), 120);
    FData.img.Thumb.CSP2_meanFreq(trials)  =  meanfreq(FData.img.Thumb.CSP2(:, trials), 120);
end

for trials = 1 : Nexe
    FData.exe.test.CSP1_meanFreq  (trials)  =  meanfreq(FData.exe.test.CSP1(:, trials), 120);
	FData.exe.test.CSP2_meanFreq  (trials)  =  meanfreq(FData.exe.test.CSP2(:, trials), 120);
end


for trials = 1 : Nimg
    FData.img.test.CSP1_meanFreq  (trials)  =  meanfreq(FData.img.test.CSP1(:, trials), 120);
	FData.img.test.CSP2_meanFreq  (trials)  =  meanfreq(FData.img.test.CSP2(:, trials), 120);
end

%% CSP Median Freq

for trials = 1 : 20
    FData.exe.Arm.CSP1_medFreq  (trials)  =  medfreq(FData.exe.Arm.CSP1(:, trials), 120);
    FData.exe.Leg.CSP1_medFreq  (trials)  =  medfreq(FData.exe.Leg.CSP1(:, trials), 120);
    FData.exe.Idle.CSP1_medFreq (trials)  =  medfreq(FData.exe.Idle.CSP1(:, trials), 120);
    FData.exe.Thumb.CSP1_medFreq(trials)  =  medfreq(FData.exe.Thumb.CSP1(:, trials), 120);

    FData.img.Arm.CSP1_medFreq  (trials)  =  medfreq(FData.img.Arm.CSP1(:, trials), 120);
    FData.img.Leg.CSP1_medFreq  (trials)  =  medfreq(FData.img.Leg.CSP1(:, trials), 120);
    FData.img.Thumb.CSP1_medFreq(trials)  =  medfreq(FData.img.Thumb.CSP1(:, trials), 120);
 
    
    FData.exe.Arm.CSP2_medFreq  (trials)  =  medfreq(FData.exe.Arm.CSP2(:, trials), 120);
    FData.exe.Leg.CSP2_medFreq  (trials)  =  medfreq(FData.exe.Leg.CSP2(:, trials), 120);
    FData.exe.Idle.CSP2_medFreq (trials)  =  medfreq(FData.exe.Idle.CSP2(:, trials), 120);
    FData.exe.Thumb.CSP2_medFreq(trials)  =  medfreq(FData.exe.Thumb.CSP2(:, trials), 120);

    FData.img.Arm.CSP2_medFreq  (trials)  =  medfreq(FData.img.Arm.CSP2(:, trials), 120);
    FData.img.Leg.CSP2_medFreq  (trials)  =  medfreq(FData.img.Leg.CSP2(:, trials), 120);
    FData.img.Thumb.CSP2_medFreq(trials)  =  medfreq(FData.img.Thumb.CSP2(:, trials), 120);
end

for trials = 1 : Nexe
    FData.exe.test.CSP1_medFreq  (trials)  =  medfreq(FData.exe.test.CSP1(:, trials), 120);
	FData.exe.test.CSP2_medFreq  (trials)  =  medfreq(FData.exe.test.CSP2(:, trials), 120);
end


for trials = 1 : Nimg
    FData.img.test.CSP1_medFreq  (trials)  =  medfreq(FData.img.test.CSP1(:, trials), 120);
	FData.img.test.CSP2_medFreq  (trials)  =  medfreq(FData.img.test.CSP2(:, trials), 120);
end
%% Wavelet
% for trials = 1 : 20
%     for channels = 1 : 64
%         
%         FData.exe.Arm.DWT(channels, trials) = dwt(FData.exe.Arm.alpha_band(channels, :, trials));
%         FData.exe.Leg.DWT(channels, trials) = dwt(FData.exe.Arm.alpha_band(channels, :, trials));
%         FData.exe.Thumb.DWT(channels, trials) = dwt(FData.exe.Arm.alpha_band(channels, :, trials));
%         FData.exe.Idle.DWT(channels, trials) = dwt(FData.exe.Arm.alpha_band(channels, :, trials));
%         
%     end
% end


%% JValue Matrix 

% J_Index s are not correct :(

clear J_EXE J_IMG
[J_EXE(:,1),~] = Jvalue(FData.exe.Arm, FData.exe.Leg);
[J_EXE(:,2),~] = Jvalue(FData.exe.Arm, FData.exe.Thumb);
[J_EXE(:,3),~] =  Jvalue(FData.exe.Arm, FData.exe.Idle);
[J_EXE(:,4),~] = Jvalue(FData.exe.Leg, FData.exe.Thumb);
[J_EXE(:,5),~] = Jvalue(FData.exe.Leg, FData.exe.Idle);
[J_EXE(:,6),J_EXE_Index,Feat] = Jvalue(FData.exe.Thumb, FData.exe.Idle);


[J_IMG(:,1),~] = Jvalue(FData.img.Arm, FData.img.Leg);
[J_IMG(:,2),~] = Jvalue(FData.img.Arm, FData.img.Thumb);
[J_IMG(:,3), J_IMG_Index] = Jvalue(FData.img.Leg, FData.img.Thumb);

%% 5 Fold Cross Validation
score = [];
num = [];
j = 0;
for i = 1 : 10 : 1000
    j = j + 1;
    [TrainSet, labels, feature_list2] = Feature_Selector(FData,'exe',i);
    score(j) = 100-Final_CrossVal(TrainSet, labels, 5, 10);
    num(j) = i;
end
for i = 1001 : 100 : 10000
    j = j + 1;
    [TrainSet, labels, feature_list2] = Feature_Selector(FData,'exe',i);
    score(j) = 100-Final_CrossVal(TrainSet, labels, 5, 10);
    num(j) = i;
end
for i = 10001 : 2000 : 20000
    j = j + 1;
    [TrainSet, labels, feature_list2] = Feature_Selector(FData,'exe',100);
    score(j) = 100-Final_CrossVal(TrainSet, labels, 5, 10);
    num(j) = i;
end

figure
semilogx(num, score);
title('5-Fold Cross Validation Score - Motor Execution');
ylabel('% of Correct Classification');
xlabel('Number of Features');

%%
score = [];
num = [];
j = 0;
for i = 1 : 10 : 1000
    j = j + 1;
    [TrainSet, labels, feature_list2] = Feature_Selector(FData,'exe',i);
    score(j) = 100-Final_CrossVal(TrainSet, labels, 5, 10);
    num(j) = i;
end
for i = 1001 : 100 : 10000
    j = j + 1;
    [TrainSet, labels, feature_list2] = Feature_Selector(FData,'exe',i);
    score(j) = 100-Final_CrossVal(TrainSet, labels, 5, 10);
    num(j) = i;
end
for i = 10001 : 2000 : 20000
    j = j + 1;
    [TrainSet, labels, feature_list2] = Feature_Selector(FData,'exe',100);
    score(j) = 100-Final_CrossVal(TrainSet, labels, 5, 10);
    num(j) = i;
end
%%
figure
semilogx(num, score);
title('5-Fold Cross Validation Score - Motor Imagery');
ylabel('% of Correct Classification');
xlabel('Number of Features');


%% A K-Fold like algorithm
DataSet = [];
labels = [];
[DataSet, labels, ~] = Feature_Selector(FData,'exe',100);
 
for i = 1 : 400
    I = randperm(80);
    TestIndex    = I(1:16);
    TrainIndex   = I(17:80);        
    result = MultiSVM(DataSet(TrainIndex,:), labels(TrainIndex), DataSet(TestIndex,:));
    correct(i) = sum(labels(TestIndex) == result')/length(result)*100;
end
correct = mean(correct)

%%

[TrainSet_exe, labels_exe, feature_list_exe, TestSet_exe] = Feature_Selector(FData, 'exe', 100);
sub11_exe = MultiSVM(TrainSet_exe,labels_exe, TestSet_exe);

[TrainSet_img, labels_img, feature_list_img, TestSet_img] = Feature_Selector(FData, 'img', 100);
sub11_img = MultiSVM(TrainSet_img,labels_img, TestSet_img);

%% Performing PCA on Features - EXE

[TrainSet, labels, feature_list2] = Feature_Selector(FData,'exe',1000);
[coeff,score,latent] = pca(TrainSet);
figure
stem(cumsum(latent)./sum(latent))
ylabel('Cumm. Sum of Cov Eigenvalue;');
xlabel('Eigenvector number');

Y = TrainSet* coeff;
figure
plot3(Y(1:20,1), Y(1:20,2), Y(1:20,3),'.');
hold on
plot3(Y(21:40,1), Y(21:40,2), Y(21:40,3),'.');
plot3(Y(41:60,1), Y(41:60,2), Y(41:60,3),'.');
plot3(Y(61:80,1), Y(61:80,2), Y(61:80,3),'.');
legend('Arm', 'Leg', 'Thumb', 'Idle');
xlabel('PCA 1')
ylabel('PCA 2')
zlabel('PCA 3')
title('PCA on exe Data');

labels(1:20)  = 1; labels(21:40) = 2;
labels(41:60) = 3; labels(61:80) = 4;
correct_mean = [];
correct = [];
correct_var = [];

for N = 1 : 79
    for i = 1 : 10
        I = randperm(80);
        TestIndex    = I(1:16);
        TrainIndex   = I(17:80);        
        result = MultiSVM(Y(TrainIndex,1:N), labels(TrainIndex), Y(TestIndex,1:N));
        correct(i) = sum(labels(TestIndex) == result')/length(result)*100;
    end
    correct_mean(N) = mean(correct);
    correct_var(N) = var(correct);
end

%%
figure
hold on
plot(correct_mean + sqrt(correct_var))
plot(correct_mean)
plot(correct_mean - sqrt(correct_var))

ylabel('Correct Percentage');
xlabel('Number of eigenvalues');
legend('Mean - STD','Mean', 'Mean + STD');
title('PCA Based Classification Cross Validation Results - EXE');
%% Performing PCA on Features - Img

[TrainSet, labels, feature_list2] = Feature_Selector(FData,'img',1000);
[coeff,score,latent] = pca(TrainSet);
figure
stem(cumsum(latent)./sum(latent))
ylabel('Cumm. Sum of Cov Eigenvalue;');
xlabel('Eigenvector number');

Y = TrainSet* coeff;
figure
plot3(Y(1:20,1), Y(1:20,2), Y(1:20,3),'.');
hold on
plot3(Y(21:40,1), Y(21:40,2), Y(21:40,3),'.');
plot3(Y(41:60,1), Y(41:60,2), Y(41:60,3),'.');
legend('Arm', 'Leg', 'Thumb');
xlabel('PCA 1')
ylabel('PCA 2')
zlabel('PCA 3')
title('PCA on img Data');

labels(1:20)  = 1; labels(21:40) = 2;
labels(41:60) = 3; 
correct_mean = [];
correct = [];
correct_var = [];

for N = 1 : 59
    for i = 1 : 10
        I = randperm(60);
        TestIndex    = I(1:12);
        TrainIndex   = I(13:60);        
        result = MultiSVM(Y(TrainIndex,1:N), labels(TrainIndex), Y(TestIndex,1:N));
        correct(i) = sum(labels(TestIndex) == result')/length(result)*100;
    end
    correct_mean(N) = mean(correct);
    correct_var(N) = var(correct);
end
%%
figure
hold on
plot(correct_mean + sqrt(correct_var))
plot(correct_mean)
plot(correct_mean - sqrt(correct_var))

ylabel('Correct Percentage');
xlabel('Number of eigenvalues');
legend('Mean - STD','Mean', 'Mean + STD');
title('PCA Based Classification Cross Validation Results - IMG');

%% Classify Test
Ytest = [];
Ytrain = [];

[TrainSet_img, labels_img, feature_list_img, TestSet_img] = Feature_Selector(FData, 'img', 100);
[coeff,score,latent] = pca(TrainSet_img);
Ytest  = TestSet_img* coeff;
Ytrain = TrainSet_img* coeff;
sub_img = MultiSVM(Ytrain(:,1:20), labels_img, Ytest(:,1:20));
save('PCA_sub03_img', 'sub_img');
%% Classify Test
Ytest = [];
Ytrain = [];

[TrainSet_exe, labels_exe, feature_list_exe, TestSet_exe] = Feature_Selector(FData, 'exe', 100);
[coeff,score,latent] = pca(TrainSet_exe);
Ytest  = TestSet_exe* coeff;
Ytrain = TrainSet_exe* coeff;
sub_exe = MultiSVM(Ytrain(:,1:20), labels_exe, Ytest(:,1:20));
save('PCA_sub03_exe', 'sub_exe');