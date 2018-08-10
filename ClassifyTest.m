clc
clear all
load('complete_train.mat');

cd('./dataset')
load('subject_1.mat')
cd('..')

Test.exe.signal = data.test(:,1:20:end, :);

s_Test  = size(Test.exe.signal);
mean_test    = reshape(repmat(squeeze(mean(Test.exe.signal, 2)), [1, 1, s_Test(2)]), [s_Test(1), s_Test(2), s_Test(3)]);


Test.exe.signal  = Test.exe.signal  - mean_test;      
Test.exe.var   =  squeeze(var(Test.exe.signal,0, 2));
Test.exe.skew   =  squeeze(skewness(Test.exe.signal,0, 2));
Test.exe.RMS   =  squeeze(rms(Test.exe.signal,2));

for channel = 1 : 63
    Test.exe.meanFreq(channel, :)  =  meanfreq(squeeze(Test.exe.signal(channel, :, :)), fs);
end

for channel = 1 : 63
    Test.exe.medFreq(channel, :)  =  medfreq(squeeze(Test.exe.signal(channel, :, :)), fs);
end

for channel = 1 : 63
    Test.exe.DST (channel, :, :)  =  dst(squeeze(Test.exe.signal(channel, :, :)));
end

for channel = 1 : 63
    Test.exe.DCT (channel, :, :)  =  dct(squeeze(Test.exe.signal(channel, :, :)));
end

h_alpha = BPF(7200, 7.5  , 13.5, fs);
h_beta =  BPF(7200, 13.5 , 20, fs);
h_theta = BPF(7200, 3.5  ,  7.5, fs);
h_delta = BPF(7200, eps,  3.5, fs);

for channels = 1 : 63    
    Test.exe.alpha_band(channels,:,:)  = doFilt(h_alpha, squeeze(Test.exe.signal(channels, :, :)));
    Test.exe.beta_band(channels,:, :)  = doFilt(h_beta,  squeeze(Test.exe.signal(channels, :, :)));
    Test.exe.theta_band(channels,:, :) = doFilt(h_theta, squeeze(Test.exe.signal(channels, :, :)));
    Test.exe.delta_band(channels,:, :) = doFilt(h_delta, squeeze(Test.exe.signal(channels, :, :)));        
end

clear h_alpha h_beta h_delta h_theta


Test.exe.alphaEnergy = squeeze(sum(Test.exe.alpha_band.^2, 2));
Test.exe.betaEnergy  = squeeze(sum(Test.exe.beta_band.^2, 2));
Test.exe.thetaEnergy = squeeze(sum(Test.exe.theta_band.^2, 2));
Test.exe.deltaEnergy = squeeze(sum(Test.exe.delta_band.^2, 2));

for channels = 1 : 63
    Test.exe.STFT (channels, :, :, :)  =  abs(STFT(squeeze(Test.exe.signal(channels, :, :)),10,0,15,fs));
end


for trials = 1 : s_Test(3)
    Test.exe.CSP1(:, trials)  = squeeze(spatialFilter_EXE(1,:))*squeeze(Test.exe.signal(:,:,trials));
    Test.exe.CSP2(:, trials)  = squeeze(spatialFilter_EXE(2,:))*squeeze(Test.exe.signal(:,:,trials));
end

%%
for i = 1 : s_Test(3)
    Test_Feature(i,:) =...
                    [reshape(Test.exe.signal(:,:,i),1,[]) ...
                    Test.exe.var(:,i)' Test.exe.skew(:,i)' ...
                    Test.exe.RMS(:,i)' Test.exe.meanFreq(:,i)' ...
                    Test.exe.medFreq(:,i)'...
                    reshape(Test.exe.DST(:,:,i),1,[]) ...
                    reshape(Test.exe.DCT(:,:,i),1,[]) ...
                    reshape(Test.exe.alpha_band(:,:,i),1,[]) ...
                    reshape(Test.exe.beta_band(:,:,i),1,[]) ...
                    reshape(Test.exe.theta_band(:,:,i),1,[]) ...
                    reshape(Test.exe.delta_band(:,:,i),1,[]) ...
                    Test.exe.alphaEnergy(:,i)'...
                    Test.exe.betaEnergy(:,i)'...
                    Test.exe.thetaEnergy(:,i)'...
                    Test.exe.deltaEnergy(:,i)'...
                    reshape(Test.exe.STFT(:,:,:,i),1,[]) ...
                    Test.exe.CSP1(:,i)' ...
                    Test.exe.CSP2(:,i)' ...
                    ];
end

%% Classify Test Signals

predicted_Label = predict(Obj, Test_Feature(:, pVal<0.0001))
