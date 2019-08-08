function [FData] = propMedian(Data,Nmed)
[rows cols samples] = size(Data);
FData = zeros(rows,cols,samples);
for k = 1:samples
     FData(:,:,k) = medfilt2(Data(:,:,k),[Nmed Nmed],'symmetric'); 
end

