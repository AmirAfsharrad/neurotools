clear
clc
close

cd('./dataset')
load('subject_1.mat')
cd('..')
fs = 2400;

%% Loading Subject's Data

signal = cat(3, data.train{1}, data.train{2}, data.train{3}, data.train{4});
signal = signal(1:63, :, :);
signal = cat(3, signal, data.test);
signal = signal(:, 1:20:end, :);
fs = fs/20;
FData.signal = signal;
clear signal

%% Removing DC Component

matSize   = size(FData.signal);
mean     = reshape(repmat(squeeze(mean(FData.signal, 2)), [1, 1, matSize(2)]), [matSize(1), matSize(2), matSize(3)]);
FData.signal  = FData.signal - mean;      
clear mean

%% Feature Extraction

FData.var    =  squeeze(var(FData.signal,0, 2));       % Variance
FData.skew   =  squeeze(skewness(FData.signal,0, 2));  % Skewness
FData.RMS    =  squeeze(rms(FData.signal,2));          % RMS Value

% Mean Frequency
for channel = 1 : 63
    FData.meanFreq(channel,:)    =  meanfreq(squeeze(FData.signal(channel,:,:)), fs);
end

% Median Frequency
for channel = 1 : 63
    FData.medFreq(channel,:)    =  medfreq(squeeze(FData.signal(channel,:,:)), fs);
end

% Sine Transform
for channel = 1 : 63
    FData.DST  (channel, :, :)  =  dst(squeeze(FData.signal(channel, :, :)));
end

% Cosine Transform
for channel = 1 : 63
    FData.DCT  (channel, :, :)  =  dct(squeeze(FData.signal(channel, :, :)));
end

% Frequncy Band Filters
h_alpha = BPF(7200, 7.5  , 13.5, fs);
h_beta =  BPF(7200, 13.5 , 20, fs);
h_theta = BPF(7200, 3.5  ,  7.5, fs);
h_delta = BPF(7200, eps,  3.5, fs);


% Filtering
for channel = 1 : 63
    FData.alpha_band(channel,:,:)  = doFilt(h_alpha, squeeze(FData.signal(channel, :, :)));
    FData.beta_band (channel,:, :) = doFilt(h_beta,  squeeze(FData.signal(channel, :, :)));
    FData.theta_band(channel,:, :) = doFilt(h_theta, squeeze(FData.signal(channel, :, :)));
    FData.delta_band(channel,:, :) = doFilt(h_delta, squeeze(FData.signal(channel, :, :)));
end

clear h_alpha h_beta h_delta h_theta


FData.alphaEnergy = squeeze(sum(FData.alpha_band.^2, 2));
FData.betaEnergy  = squeeze(sum(FData.beta_band.^2, 2));
FData.thetaEnergy = squeeze(sum(FData.theta_band.^2, 2));
FData.deltaEnergy = squeeze(sum(FData.delta_band.^2, 2));


for channel = 1 : 63
    FData.STFT  (channel, :, :, :)  =  abs(STFT(squeeze(FData.signal(channel, :, :)),10,0,15,fs));
end

clear channel

%% Train Feature Matrix

Feature = [];
label = [ones(1, size(data.train{1}, 3)), 2*ones(1, size(data.train{2}, 3)), ...
    3*ones(1, size(data.train{3}, 3)), 4*ones(1, size(data.train{4}, 3)), zeros(1, size(data.test, 3))];

for i = 1 : matSize(3)
    Feature(i,:) = [reshape(FData.signal(:,:,i),1,[]) ...
                    FData.var(:,i)' FData.skew(:,i)' ...
                    FData.RMS(:,i)' FData.meanFreq(:,i)' ...
                    FData.medFreq(:,i)'...
                    reshape(FData.DST(:,:,i),1,[]) ...
                    reshape(FData.DCT(:,:,i),1,[]) ...
                    reshape(FData.alpha_band(:,:,i),1,[]) ...
                    reshape(FData.beta_band(:,:,i),1,[]) ...
                    reshape(FData.theta_band(:,:,i),1,[]) ...
                    reshape(FData.delta_band(:,:,i),1,[]) ...
                    FData.alphaEnergy(:,i)'...
                    FData.betaEnergy(:,i)'...
                    FData.thetaEnergy(:,i)'...
                    FData.deltaEnergy(:,i)'...
                    reshape(FData.STFT(:,:,:,i),1,[]) ...
                    ];  
end

clear i

%% Feature Names

FeatureName = [repmat({'Signal'}, 1, length(reshape(FData.signal(:,:,1),1,[]))) ...
               repmat({'Var'}, 1, length( FData.var(:,1)' )) ...
               repmat({'Skew'}, 1, length( FData.skew(:,1)' )) ...
               repmat({'RMS'}, 1, length( FData.RMS(:,1)' )) ...
               repmat({'meanFreq'}, 1, length( FData.meanFreq(:,1)' )) ...
               repmat({'medFreq'}, 1, length( FData.medFreq(:,1)'  )) ...
               repmat({'DST'}, 1, length( reshape(FData.DST(:,:,1),1,[]))) ...
               repmat({'DCT'}, 1, length( reshape(FData.DCT(:,:,1),1,[])  )) ...
               repmat({'alphaSignal'}, 1, length( reshape(FData.alpha_band(:,:,1),1,[]) )) ...
               repmat({'betaSignal'}, 1, length( reshape(FData.beta_band(:,:,1),1,[]) )) ...
               repmat({'thetaSignal'}, 1, length( reshape(FData.theta_band(:,:,1),1,[]) )) ...
               repmat({'deltaSignal'}, 1, length( reshape(FData.delta_band(:,:,1),1,[]) )) ...
               repmat({'alphaEnergy'}, 1, length( FData.alphaEnergy(:,1)')) ...
               repmat({'betaEnergy'}, 1, length( FData.betaEnergy(:,1)')) ...
               repmat({'thetaEnergy'}, 1, length( FData.thetaEnergy(:,1)' )) ...
               repmat({'deltaEnergy'}, 1, length( FData.deltaEnergy(:,1)' )) ...
               repmat({'STFT'}, 1, length( reshape(FData.STFT(:,:,:,1),1,[]) )) ...
];

%% Seperating Train and Test

Train_Feature = Feature(label ~= 0, :);
Train_Label = label(label ~= 0);

Test_Feature = Feature(label == 0, :);

%% ANOVA

for i = 1 : size(Train_Feature, 2)
    pVal(i) = anova1(Train_Feature(:,i), Train_Label, 'off');
    if(mod(i, 10000) == 0)
        i
    end
end

%% Choosen Features

figure
pie(categorical(FeatureName(pVal<0.0001)))
figure
histogram(categorical(FeatureName(pVal<0.0001)))

%% Train LDA

Obj = fitcdiscr(Train_Feature(:, pVal<0.0001), Train_Label)

%% Predict Test Data

predict(Obj, Test_Feature(:,pVal<0.0001))

