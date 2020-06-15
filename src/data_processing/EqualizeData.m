function [EData,E,P,E2] = EqualizeData(Data,time,Vv,Vh,trs)

[rows cols samples] = size(Data);

E = zeros(samples,1);
P = zeros(samples,1);
EData = zeros(rows,cols,samples);

% Energy
for n = 1:samples
    E(n) = sqrt(sum(sum(abs(Data(:,:,n).^2))));
    P(n) = pi*time(n)*(3/2*(Vh+Vv)-sqrt(Vv*Vh));
    EData(:,:,n) = Data(:,:,n)/E(n)*sqrt(P(n));
end

[maxx ~] = max(E);

strt = 1;

while strt < samples
    if E(strt) > trs*maxx  %0.98 normalnie
        break
    else
     strt = strt + 1;
    end
end


E2 = E;
n = 1;

while n < strt
    E2(n) = maxx;
    n = n+1;
end