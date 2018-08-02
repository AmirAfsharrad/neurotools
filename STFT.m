function [output] = STFT(signals, a, b, c, fs)
s = size(signals);
for trials = 1 : s(2)
    output(:,:,trials)  =  abs(spectrogram(squeeze(signals(:, trials)),a, b, c, fs));
end
end

