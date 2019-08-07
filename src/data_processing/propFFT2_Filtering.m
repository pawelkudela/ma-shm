function [Filtered_FFTData] = propFFT2_Filtering(Data,FilterMask)
[rows cols samples] = size(Data);
Filtered_FFTData = zeros(rows,cols,samples);
for k = 1:samples
     Filtered_FFTData(:,:,k) = Data(:,:,k).*FilterMask;
end