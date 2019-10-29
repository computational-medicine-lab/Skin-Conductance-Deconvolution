function y_filtered = LowPassFilter(y, Fs, wc, filter_order)

% y_filtered = LowPassFilter(y, Fs, wc, filter_order)
% y : Input signal 
% Fs: sampling Frequency
% wc: Cutoff frequency
% filter_order: Filter Order


y = y(:); %ensure colum vector
Fs = Fs/2;
rng default;
% Design a 64th order lowpass FIR filter with cutoff frequency of 75 Hz.

Fnorm = wc/(Fs/2);           % Normalized frequency
df = designfilt('lowpassfir','FilterOrder',filter_order,'CutoffFrequency',Fnorm);
%grpdelay(df,2048,Fs);   % plot group delay
D = mean(grpdelay(df)); %mean group delay
temp = filter(df,[y(1) * ones(D,1); y ; y(end) * ones(D,1)]); % Append D zeros to the input data
temp = temp(2*D+1:end);                  % Shift data to compensate for delay
y_filtered = temp(:);
end