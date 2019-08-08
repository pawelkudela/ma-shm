function [H2D] = Hann2D(cols,rows,mrg)
%% Space window
H = hanning(floor(rows*mrg));
W = hanning(floor(rows*mrg)); %colss dla proporcjonalnych
A = [H(1:floor(size(H,1)/2))' ones(rows-size(H,1),1)' H(floor(size(H,1)/2)+1:end)']';
B = [W(1:floor(size(W,1)/2))' ones(cols-size(W,1),1)' W(floor(size(W,1)/2)+1:end)'];
H2D = A*B;