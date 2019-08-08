function [IDATA] = propIFFT2(Data)
[rows cols samples] = size(Data);
IDATA = zeros(rows,cols,samples);
for k = 1:samples
     IDATA(:,:,k) = ifftn(ifftshift(Data(:,:,k)));
end
