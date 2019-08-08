function [FilterMask] = propMask(Data,NFrames,thresh)
% thersold = thresh * mean
[rows cols ~] = size(Data);
AvgFrames = zeros(rows,cols);

for i = 1:length(NFrames)
    AvgFrames = AvgFrames + Data(:,:,NFrames(i));
end

FilterMask = ones(rows,cols);         
threshold = mean(mean(abs(real(AvgFrames))))*thresh;
   
for x = 1:cols
    for y = 1:rows
        if abs((AvgFrames(y,x))) > threshold;
            FilterMask(y,x) = 0;
        end
    end
end

%figure
%[Cmap] = colormapwlasny(0,0.1);
%colormap(Cmap)
%shading interp; view(2); axis off; axis equal; set(gcf,'Renderer','zbuffer');