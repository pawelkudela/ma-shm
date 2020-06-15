function [FilterMask,CDF_x,CDF,h] = propMask_cum(Data,NFrames,thresh)

h1 = figure('Name','CFD');
[rows cols ~] = size(Data);
AvgFrames = zeros(rows,cols);

for i = 1:length(NFrames)
    AvgFrames = AvgFrames + (Data(:,:,NFrames(i)));
end

FilterMask = ones(rows,cols); 
%%%
MinX = min(min(abs(AvgFrames)));
MaxX = max(max(abs(AvgFrames)));
x = MinX:(MaxX-MinX)/(1000-1):MaxX;
A = reshape(AvgFrames,cols*rows,1);
%   figure(1)
%    hist(A,100)
n_elements = histc(abs(A),x);
c_elements = cumsum(n_elements);
CDF = c_elements/cols/rows;
CDF_x = x;
    %figure(2)
    %bar(x,c_elements,'BarWidth',0.33)
    plot(CDF_x,CDF,'black','LineWidth',2);
    axis tight
    xlabel('Magnitude')
    ylabel('Cumulative distribution function')
    ylim([0.5 1.05])
    xlim([-max(x)*0.02 max(x)]);
    
thresho = (100-thresh)/100*cols*rows;
n = 1;
while c_elements(n) < thresho
    n = n + 1;
end
threshold = x(n);

    hold on
    line([threshold threshold], [0 (100-thresh)/100],'LineWidth',1,'LineStyle','--','Color',[.1 .1 1]);
    line([-max(x)*0.01 threshold], [(100-thresh)/100 (100-thresh)/100],'LineWidth',1,'LineStyle','--','Color',[.1 .1 1] )
     
for x = 1:cols
    for y = 1:rows
        if abs((AvgFrames(y,x))) > threshold;
            FilterMask(y,x) = 0;
        end
    end
end
