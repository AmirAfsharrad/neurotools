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

% EEG Signals
FData.exe.Arm.signal     = Data.exe.Arm.signal(1:63,1:20:end,:);
FData.exe.Thumb.signal   = Data.exe.Thumb.signal(1:63,1:20:end,:);
FData.exe.Leg.signal     = Data.exe.Leg.signal(1:63,1:20:end,:);
FData.exe.Idle.signal    = Data.exe.Idle.signal(1:63,1:20:end,:);
fs = fs/20;

% Start of Movement's Index
[FData.exe.Arm.index, ~]   = find(squeeze(Data.exe.Arm.signal(64,:,:)));
[FData.exe.Thumb.index, ~] = find(squeeze(Data.exe.Thumb.signal(64,:,:)));
[FData.exe.Leg.index, ~]   = find(squeeze(Data.exe.Leg.signal(64,:,:)));
[FData.exe.Idle.index, ~]  = find(squeeze(Data.exe.Idle.signal(64,:,:)));


clear Data

%% Removing DC Component

s_Arm   = size(FData.exe.Arm.signal);
s_Leg   = size(FData.exe.Leg.signal);
s_Thumb = size(FData.exe.Thumb.signal);
s_Idle  = size(FData.exe.Idle.signal);



mean_arm     = reshape(repmat(squeeze(mean(FData.exe.Arm.signal, 2)), [1, 1, s_Arm(2)]), [s_Arm(1), s_Arm(2), s_Arm(3)]);
mean_leg     = reshape(repmat(squeeze(mean(FData.exe.Leg.signal, 2)), [1, 1, s_Leg(2)]), [s_Leg(1), s_Leg(2), s_Leg(3)]);
mean_thumb   = reshape(repmat(squeeze(mean(FData.exe.Thumb.signal, 2)), [1, 1, s_Thumb(2)]), [s_Thumb(1), s_Thumb(2), s_Thumb(3)]);
mean_idle    = reshape(repmat(squeeze(mean(FData.exe.Idle.signal, 2)), [1, 1, s_Idle(2)]), [s_Idle(1), s_Idle(2), s_Idle(3)]);


FData.exe.Arm.signal   = FData.exe.Arm.signal   - mean_arm;
FData.exe.Leg.signal   = FData.exe.Leg.signal   - mean_leg;
FData.exe.Thumb.signal = FData.exe.Thumb.signal - mean_thumb;
FData.exe.Idle.signal  = FData.exe.Idle.signal  - mean_idle;      


clear mean_arm mean_idle mean_leg mean_test mean_thumb data

%% Variance

FData.exe.Arm.var    =  squeeze(var(FData.exe.Arm.signal,0, 2));
FData.exe.Leg.var    =  squeeze(var(FData.exe.Leg.signal,0, 2));
FData.exe.Idle.var   =  squeeze(var(FData.exe.Idle.signal,0, 2));
FData.exe.Thumb.var  =  squeeze(var(FData.exe.Thumb.signal,0, 2));        


%% Skewness

FData.exe.Arm.skew    =  squeeze(skewness(FData.exe.Arm.signal,0, 2));
FData.exe.Leg.skew    =  squeeze(skewness(FData.exe.Leg.signal,0, 2));
FData.exe.Idle.skew   =  squeeze(skewness(FData.exe.Idle.signal,0, 2));
FData.exe.Thumb.skew  =  squeeze(skewness(FData.exe.Thumb.signal,0, 2));        

%% Root Mean Square

FData.exe.Arm.RMS    =  squeeze(rms(FData.exe.Arm.signal,2));
FData.exe.Leg.RMS    =  squeeze(rms(FData.exe.Leg.signal,2));
FData.exe.Idle.RMS   =  squeeze(rms(FData.exe.Idle.signal,2));
FData.exe.Thumb.RMS  =  squeeze(rms(FData.exe.Thumb.signal,2));

%% Mean Freq

for channel = 1 : 63
    FData.exe.Arm.meanFreq(channel,:)    =  meanfreq(squeeze(FData.exe.Arm.signal(channel,:,:)), fs);
    FData.exe.Leg.meanFreq(channel,:)    =  meanfreq(squeeze(FData.exe.Leg.signal(channel,:,:)), fs);
    FData.exe.Idle.meanFreq(channel,:)   =  meanfreq(squeeze(FData.exe.Idle.signal(channel,:,:)), fs);
    FData.exe.Thumb.meanFreq(channel,:)  =  meanfreq(squeeze(FData.exe.Thumb.signal(channel,:,:)), fs);
end

%% Med Freq

for channel = 1 : 63
    FData.exe.Arm.medFreq(channel,:)    =  medfreq(squeeze(FData.exe.Arm.signal(channel,:,:)), fs);
    FData.exe.Leg.medFreq(channel,:)    =  medfreq(squeeze(FData.exe.Leg.signal(channel,:,:)), fs);
    FData.exe.Idle.medFreq(channel,:)   =  medfreq(squeeze(FData.exe.Idle.signal(channel,:,:)), fs);
    FData.exe.Thumb.medFreq(channel,:)  =  medfreq(squeeze(FData.exe.Thumb.signal(channel,:,:)), fs);    
end

%% Sine Transform

for channel = 1 : 63
    FData.exe.Arm.DST  (channel, :, :)  =  dst(squeeze(FData.exe.Arm.signal(channel, :, :)));
    FData.exe.Leg.DST  (channel, :, :)  =  dst(squeeze(FData.exe.Leg.signal(channel, :, :)));
    FData.exe.Idle.DST (channel, :, :)  =  dst(squeeze(FData.exe.Idle.signal(channel, :, :)));
    FData.exe.Thumb.DST(channel, :, :)  =  dst(squeeze(FData.exe.Thumb.signal(channel, :, :)));
end

%% Cosine Transform

for channel = 1 : 63
    FData.exe.Arm.DCT  (channel, :, :)  =  dct(squeeze(FData.exe.Arm.signal(channel, :, :)));
    FData.exe.Leg.DCT  (channel, :, :)  =  dct(squeeze(FData.exe.Leg.signal(channel, :, :)));
    FData.exe.Idle.DCT (channel, :, :)  =  dct(squeeze(FData.exe.Idle.signal(channel, :, :)));
    FData.exe.Thumb.DCT(channel, :, :)  =  dct(squeeze(FData.exe.Thumb.signal(channel, :, :)));
end

%% Band Pass Filters

h_alpha = BPF(7200, 7.5  , 13.5, fs);
h_beta =  BPF(7200, 13.5 , 20, fs);
h_theta = BPF(7200, 3.5  ,  7.5, fs);
h_delta = BPF(7200, eps,  3.5, fs);

%% Filtering Frequency Bands

for channels = 1 : 63
    FData.exe.Arm.alpha_band(channels,:,:)  = doFilt(h_alpha, squeeze(FData.exe.Arm.signal(channels, :, :)));
    FData.exe.Arm.beta_band(channels,:, :)  = doFilt(h_beta,  squeeze(FData.exe.Arm.signal(channels, :, :)));
    FData.exe.Arm.theta_band(channels,:, :) = doFilt(h_theta, squeeze(FData.exe.Arm.signal(channels, :, :)));
    FData.exe.Arm.delta_band(channels,:, :) = doFilt(h_delta, squeeze(FData.exe.Arm.signal(channels, :, :)));

    FData.exe.Leg.alpha_band(channels,:,:)  = doFilt(h_alpha, squeeze(FData.exe.Leg.signal(channels, :, :)));
    FData.exe.Leg.beta_band(channels,:, :)  = doFilt(h_beta,  squeeze(FData.exe.Leg.signal(channels, :, :)));
    FData.exe.Leg.theta_band(channels,:, :) = doFilt(h_theta, squeeze(FData.exe.Leg.signal(channels, :, :)));
    FData.exe.Leg.delta_band(channels,:, :) = doFilt(h_delta, squeeze(FData.exe.Leg.signal(channels, :, :)));

    FData.exe.Thumb.alpha_band(channels,:,:)  = doFilt(h_alpha, squeeze(FData.exe.Thumb.signal(channels, :, :)));
    FData.exe.Thumb.beta_band(channels,:, :)  = doFilt(h_beta,  squeeze(FData.exe.Thumb.signal(channels, :, :)));
    FData.exe.Thumb.theta_band(channels,:, :) = doFilt(h_theta, squeeze(FData.exe.Thumb.signal(channels, :, :)));
    FData.exe.Thumb.delta_band(channels,:, :) = doFilt(h_delta, squeeze(FData.exe.Thumb.signal(channels, :, :)));

    FData.exe.Idle.alpha_band(channels,:,:)  = doFilt(h_alpha, squeeze(FData.exe.Idle.signal(channels, :, :)));
    FData.exe.Idle.beta_band(channels,:, :)  = doFilt(h_beta,  squeeze(FData.exe.Idle.signal(channels, :, :)));
    FData.exe.Idle.theta_band(channels,:, :) = doFilt(h_theta, squeeze(FData.exe.Idle.signal(channels, :, :)));
    FData.exe.Idle.delta_band(channels,:, :) = doFilt(h_delta, squeeze(FData.exe.Idle.signal(channels, :, :)));        
end

clear h_alpha h_beta h_delta h_theta

%% Checking The Procedure

figure
plotFFT(FData.exe.Arm.signal(1,:,16), fs, 0, 60, 'FFT','|fft|',10);
hold on
plotFFT(FData.exe.Arm.alpha_band(1,:,16), fs, 0, 60, 'FFT of exe.test Channel 2 Trial 16 ','',10);
plotFFT(FData.exe.Arm.beta_band(1,:,16), fs, 0, 60, 'FFT of exe.test Channel 2 Trial 16 ','',10);
plotFFT(FData.exe.Arm.theta_band(1,:,16), fs, 0, 60, 'FFT of exe.test Channel 2 Trial 16 ','',10);
plotFFT(FData.exe.Arm.delta_band(1,:,16), fs, 0, 60, 'FFT of exe.test Channel 2 Trial 16 ','',10);
legend('full','alpha','beta','theta','delta');

%% Freq Band Energy

FData.exe.Arm.alphaEnergy = squeeze(sum(FData.exe.Arm.alpha_band.^2, 2));
FData.exe.Arm.betaEnergy  = squeeze(sum(FData.exe.Arm.beta_band.^2, 2));
FData.exe.Arm.thetaEnergy = squeeze(sum(FData.exe.Arm.theta_band.^2, 2));
FData.exe.Arm.deltaEnergy = squeeze(sum(FData.exe.Arm.delta_band.^2, 2));

FData.exe.Leg.alphaEnergy = squeeze(sum(FData.exe.Leg.alpha_band.^2, 2));
FData.exe.Leg.betaEnergy  = squeeze(sum(FData.exe.Leg.beta_band.^2, 2));
FData.exe.Leg.thetaEnergy = squeeze(sum(FData.exe.Leg.theta_band.^2, 2));
FData.exe.Leg.deltaEnergy = squeeze(sum(FData.exe.Leg.delta_band.^2, 2));

FData.exe.Thumb.alphaEnergy = squeeze(sum(FData.exe.Thumb.alpha_band.^2, 2));
FData.exe.Thumb.betaEnergy  = squeeze(sum(FData.exe.Thumb.beta_band.^2, 2));
FData.exe.Thumb.thetaEnergy = squeeze(sum(FData.exe.Thumb.theta_band.^2, 2));
FData.exe.Thumb.deltaEnergy = squeeze(sum(FData.exe.Thumb.delta_band.^2, 2));

FData.exe.Idle.alphaEnergy = squeeze(sum(FData.exe.Idle.alpha_band.^2, 2));
FData.exe.Idle.betaEnergy  = squeeze(sum(FData.exe.Idle.beta_band.^2, 2));
FData.exe.Idle.thetaEnergy = squeeze(sum(FData.exe.Idle.theta_band.^2, 2));
FData.exe.Idle.deltaEnergy = squeeze(sum(FData.exe.Idle.delta_band.^2, 2));

%% STFT

for channels = 1 : 63
    FData.exe.Arm.STFT  (channels, :, :, :)  =  abs(STFT(squeeze(FData.exe.Arm.signal(channels, :, :)),10,0,15,fs));
    FData.exe.Leg.STFT  (channels, :, :, :)  =  abs(STFT(squeeze(FData.exe.Leg.signal(channels, :, :)),10,0,15,fs));
    FData.exe.Idle.STFT (channels, :, :, :)  =  abs(STFT(squeeze(FData.exe.Idle.signal(channels, :, :)),10,0,15,fs));
    FData.exe.Thumb.STFT(channels, :, :, :)  =  abs(STFT(squeeze(FData.exe.Thumb.signal(channels, :, :)),10,0,15,fs));
end

%% Common Spatial Patterns

FData.exe.covMat(1, :, :) = cov(mean(FData.exe.Arm.alpha_band + FData.exe.Arm.beta_band, 3)');
FData.exe.covMat(2, :, :) = cov(mean(FData.exe.Leg.alpha_band + FData.exe.Leg.beta_band, 3)');
FData.exe.covMat(3, :, :) = cov(mean(FData.exe.Idle.alpha_band + FData.exe.Idle.beta_band, 3)');
FData.exe.covMat(4, :, :) = cov(mean(FData.exe.Thumb.alpha_band + FData.exe.Thumb.beta_band, 3)');

spatialFilter_EXE = MulticlassCSP(squeeze(FData.exe.covMat), 2);

for trials = 1 : s_Arm(3)
    FData.exe.Arm.CSP1(:, trials)    = squeeze(spatialFilter_EXE(1,:))*squeeze(FData.exe.Arm.signal(:,:,trials));
    FData.exe.Arm.CSP2(:, trials)    = squeeze(spatialFilter_EXE(2,:))*squeeze(FData.exe.Arm.signal(:,:,trials));
end    
for trials = 1 : s_Leg(3)
    FData.exe.Leg.CSP1(:, trials)    = squeeze(spatialFilter_EXE(1,:))*squeeze(FData.exe.Leg.signal(:,:,trials));
    FData.exe.Leg.CSP2(:, trials)    = squeeze(spatialFilter_EXE(2,:))*squeeze(FData.exe.Leg.signal(:,:,trials));
end
for trials = 1 : s_Idle(3)
    FData.exe.Idle.CSP1(:, trials)   = squeeze(spatialFilter_EXE(1,:))*squeeze(FData.exe.Idle.signal(:,:,trials));
    FData.exe.Idle.CSP2(:, trials)   = squeeze(spatialFilter_EXE(2,:))*squeeze(FData.exe.Idle.signal(:,:,trials));
end    
for trials = 1 : s_Thumb(3)
    FData.exe.Thumb.CSP1(:, trials)  = squeeze(spatialFilter_EXE(1,:))*squeeze(FData.exe.Thumb.signal(:,:,trials));
    FData.exe.Thumb.CSP2(:, trials)  = squeeze(spatialFilter_EXE(2,:))*squeeze(FData.exe.Thumb.signal(:,:,trials));
end

%% Checking CSP

figure
plot(FData.exe.Arm.CSP1(:,1)')
plot(FData.exe.Arm.signal(1:3,:,1)')
xlabel('t');
legend('CSP-Filtered Signal','Original Signal Channel 01',...
    'Original Signal Channel 02','Original Signal Channel 03')

figure
hold on
plotFFT(FData.exe.Arm.CSP1(:,10),fs, 0, 600, '', '|fft|', 10)
plotFFT(FData.exe.Arm.signal(1,:,10), fs, 0, 600, '','|fft|',10)
legend('a Sample CSP','a Sample CSP(\alpha and \beta)','Original Signal')

%% Feature Matrix as a Whole

Feature = [];
label = [];
for i = 1 : s_Arm(3)
    Feature(i,:) = [reshape(FData.exe.Arm.signal(:,:,i),1,[]) ...
                    FData.exe.Arm.var(:,i)' FData.exe.Arm.skew(:,i)' ...
                    FData.exe.Arm.RMS(:,i)' FData.exe.Arm.meanFreq(:,i)' ...
                    FData.exe.Arm.medFreq(:,i)'...
                    reshape(FData.exe.Arm.DST(:,:,i),1,[]) ...
                    reshape(FData.exe.Arm.DCT(:,:,i),1,[]) ...
                    reshape(FData.exe.Arm.alpha_band(:,:,i),1,[]) ...
                    reshape(FData.exe.Arm.beta_band(:,:,i),1,[]) ...
                    reshape(FData.exe.Arm.theta_band(:,:,i),1,[]) ...
                    reshape(FData.exe.Arm.delta_band(:,:,i),1,[]) ...
                    FData.exe.Arm.alphaEnergy(:,i)'...
                    FData.exe.Arm.betaEnergy(:,i)'...
                    FData.exe.Arm.thetaEnergy(:,i)'...
                    FData.exe.Arm.deltaEnergy(:,i)'...
                    reshape(FData.exe.Arm.STFT(:,:,:,i),1,[]) ...
                    FData.exe.Arm.CSP1(:,i)' ...
                    FData.exe.Arm.CSP2(:,i)' ...
                    ];
    label(i) = 1;    
end

for i = 1 : s_Thumb(3)
    Feature(i+s_Arm(3),:) = ...
                    [reshape(FData.exe.Thumb.signal(:,:,i),1,[]) ...
                    FData.exe.Thumb.var(:,i)' FData.exe.Thumb.skew(:,i)' ...
                    FData.exe.Thumb.RMS(:,i)' FData.exe.Thumb.meanFreq(:,i)' ...
                    FData.exe.Thumb.medFreq(:,i)'...
                    reshape(FData.exe.Thumb.DST(:,:,i),1,[]) ...
                    reshape(FData.exe.Thumb.DCT(:,:,i),1,[]) ...
                    reshape(FData.exe.Thumb.alpha_band(:,:,i),1,[]) ...
                    reshape(FData.exe.Thumb.beta_band(:,:,i),1,[]) ...
                    reshape(FData.exe.Thumb.theta_band(:,:,i),1,[]) ...
                    reshape(FData.exe.Thumb.delta_band(:,:,i),1,[]) ...
                    FData.exe.Thumb.alphaEnergy(:,i)'...
                    FData.exe.Thumb.betaEnergy(:,i)'...
                    FData.exe.Thumb.thetaEnergy(:,i)'...
                    FData.exe.Thumb.deltaEnergy(:,i)'...
                    reshape(FData.exe.Thumb.STFT(:,:,:,i),1,[]) ...
                    FData.exe.Thumb.CSP1(:,i)' ...
                    FData.exe.Thumb.CSP2(:,i)' ...
                    ];
    label(i+s_Arm(3)) = 2;
end

for i = 1 : s_Leg(3)
    Feature(i+s_Arm(3)+s_Thumb(3),:) = ...
                    [reshape(FData.exe.Leg.signal(:,:,i),1,[]) ...
                    FData.exe.Leg.var(:,i)' FData.exe.Leg.skew(:,i)' ...
                    FData.exe.Leg.RMS(:,i)' FData.exe.Leg.meanFreq(:,i)' ...
                    FData.exe.Leg.medFreq(:,i)'...
                    reshape(FData.exe.Leg.DST(:,:,i),1,[]) ...
                    reshape(FData.exe.Leg.DCT(:,:,i),1,[]) ...
                    reshape(FData.exe.Leg.alpha_band(:,:,i),1,[]) ...
                    reshape(FData.exe.Leg.beta_band(:,:,i),1,[]) ...
                    reshape(FData.exe.Leg.theta_band(:,:,i),1,[]) ...
                    reshape(FData.exe.Leg.delta_band(:,:,i),1,[]) ...
                    FData.exe.Leg.alphaEnergy(:,i)'...
                    FData.exe.Leg.betaEnergy(:,i)'...
                    FData.exe.Leg.thetaEnergy(:,i)'...
                    FData.exe.Leg.deltaEnergy(:,i)'...
                    reshape(FData.exe.Leg.STFT(:,:,:,i),1,[]) ...
                    FData.exe.Leg.CSP1(:,i)' ...
                    FData.exe.Leg.CSP2(:,i)' ...
                    ];
    label(i+s_Arm(3)+s_Thumb(3)) = 3;
end


for i = 1 : s_Idle(3)
    Feature(i+s_Arm(3)+s_Thumb(3)+s_Leg(3),:) =...
                    [reshape(FData.exe.Idle.signal(:,:,i),1,[]) ...
                    FData.exe.Idle.var(:,i)' FData.exe.Idle.skew(:,i)' ...
                    FData.exe.Idle.RMS(:,i)' FData.exe.Idle.meanFreq(:,i)' ...
                    FData.exe.Idle.medFreq(:,i)'...
                    reshape(FData.exe.Idle.DST(:,:,i),1,[]) ...
                    reshape(FData.exe.Idle.DCT(:,:,i),1,[]) ...
                    reshape(FData.exe.Idle.alpha_band(:,:,i),1,[]) ...
                    reshape(FData.exe.Idle.beta_band(:,:,i),1,[]) ...
                    reshape(FData.exe.Idle.theta_band(:,:,i),1,[]) ...
                    reshape(FData.exe.Idle.delta_band(:,:,i),1,[]) ...
                    FData.exe.Idle.alphaEnergy(:,i)'...
                    FData.exe.Idle.betaEnergy(:,i)'...
                    FData.exe.Idle.thetaEnergy(:,i)'...
                    FData.exe.Idle.deltaEnergy(:,i)'...
                    reshape(FData.exe.Idle.STFT(:,:,:,i),1,[]) ...
                    FData.exe.Idle.CSP1(:,i)' ...
                    FData.exe.Idle.CSP2(:,i)' ...
                    ];
    label(i+s_Arm(3)+s_Thumb(3)+s_Leg(3)) = 4;
end

label = label';
Feature = Feature';

%% Feature Names

FeatureName = [repmat({'Signal'}, 1, length(reshape(FData.exe.Arm.signal(:,:,1),1,[]))) ...
               repmat({'Var'}, 1, length( FData.exe.Arm.var(:,1)' )) ...
               repmat({'Skew'}, 1, length( FData.exe.Arm.skew(:,1)' )) ...
               repmat({'RMS'}, 1, length( FData.exe.Arm.RMS(:,1)' )) ...
               repmat({'meanFreq'}, 1, length( FData.exe.Arm.meanFreq(:,1)' )) ...
               repmat({'medFreq'}, 1, length( FData.exe.Arm.medFreq(:,1)'  )) ...
               repmat({'DST'}, 1, length( reshape(FData.exe.Arm.DST(:,:,1),1,[]))) ...
               repmat({'DCT'}, 1, length( reshape(FData.exe.Arm.DCT(:,:,1),1,[])  )) ...
               repmat({'alphaSignal'}, 1, length( reshape(FData.exe.Arm.alpha_band(:,:,1),1,[]) )) ...
               repmat({'betaSignal'}, 1, length( reshape(FData.exe.Arm.beta_band(:,:,1),1,[]) )) ...
               repmat({'thetaSignal'}, 1, length( reshape(FData.exe.Arm.theta_band(:,:,1),1,[]) )) ...
               repmat({'deltaSignal'}, 1, length( reshape(FData.exe.Arm.delta_band(:,:,1),1,[]) )) ...
               repmat({'alphaEnergy'}, 1, length( FData.exe.Arm.alphaEnergy(:,1)')) ...
               repmat({'betaEnergy'}, 1, length( FData.exe.Arm.betaEnergy(:,1)')) ...
               repmat({'thetaEnergy'}, 1, length( FData.exe.Arm.thetaEnergy(:,1)' )) ...
               repmat({'deltaEnergy'}, 1, length( FData.exe.Arm.deltaEnergy(:,1)' )) ...
               repmat({'STFT'}, 1, length( reshape(FData.exe.Arm.STFT(:,:,:,1),1,[]) )) ...
               repmat({'CSP1'}, 1, length( FData.exe.Arm.CSP1(:,1)' )) ...
               repmat({'CSP2'}, 1, length( FData.exe.Arm.CSP2(:,1)' )) ...
];
                  
%% ANOVA for Hold-Out Validation

I = randperm(77);
test_I  = I(1:20); 
train_I = I(21:end);

train_Set = Feature(:, train_I);
test_Set = Feature(:, test_I);

train_Label = label(train_I);
test_Label = label(test_I);

for i = 1 : size(train_Set, 1)
   pVal(i) = anova1(train_Set(i,:)', train_Label, 'off');
   if(mod(i, 10000) == 0)
       i
   end
end

%% Hold-Out Validation

Obj = fitcdiscr(train_Set(pVal<0.0001,:)', train_Label)
predict(Obj, test_Set(pVal<0.0001,:)')'
test_Label'

