function [RMS] = propRMS(Data)
RMS = abs(sqrt(sum(Data(:,:,1:end).^2,3)));