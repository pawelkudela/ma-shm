function [H] = Gauss(hsize,sigma)
% sigma - standard deviation

n1=hsize;
n2=hsize;
H = zeros(n1,n2);
a = 1;
b = n1/n2;
    for x = 1:n2
       for y = 1:n1        
           d = sqrt( (((n2/2)-x)^2)*a + (((n1/2)-y)^2)*b);           
          % H(y,x) = exp(-(x^2+y^2)/(2*sigma^2));
           H(y,x) = exp(-(d^2)/(2*sigma^2));
        end
    end
end