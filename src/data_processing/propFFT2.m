function [FFTData] = propFFT2(Data)
[rows cols samples] = size(Data);
FFTData = zeros(rows,cols,samples);
for k = 1:samples
     FFTData(:,:,k) = fftshift(fftn(Data(:,:,k)));
end

