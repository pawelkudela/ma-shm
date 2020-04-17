function [Am,A] = amplification_factor(Data_exp,Data_num,N)
% AMPLIFICATION_FACTOR   Estimate amplification factor for numerical signals
%    such that Data_exp ~= A*Data_num 
%    time samples between experimental and numerical data must match
%    experimental data should be cleansed in respect to outliers
% 
% Syntax: [A] = amplification_factor(Data_exp,Data_num) 
% 
% Inputs: 
%    Data_exp - Experimental wavefield data, double matrix, dimensions [nx,ny,nft] 
%    Data_num - Numerical wavefield data, double matrix, dimensions [nx,ny,nft]
%    nx - number of points in x direction
%    ny - number of points in y direction
%    nft - number of time steps
%    N - number of test points (optional argument, default N=100)
%     
% 
% Outputs: 
%    A - Amplification factor
% 
% Example: [A] = amplification_factor(Data_exp,Data_num)
%              [A] = amplification_factor(Data_exp,Data_num,150)
% 
% Other m-files required: none 
% Subfunctions: none 
% MAT-files required: none 
% See also:  
% 

% Author: Pawel Kudela, D.Sc., Ph.D., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pk@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 

%---------------------- BEGIN CODE---------------------- 

% check arguments
if nargin == 2
    N=100;
end
Amax_num = max(max(max(abs(Data_num))));
Amax_exp = max(max(max(abs(Data_exp))));
Am = Amax_exp/Amax_num;
Arange= linspace(Am/4,4*Am,N);
F=zeros(length(Arange),1);
G=zeros(length(Arange),1);
for k=1:length(Arange)
    [k,N]
    F(k)=sum(sum(sum((max(Data_exp,[],3) - Arange(k)*max(Data_num,[],3)).^2)));
    %G(k)=sum(sum(sum(((Data_exp) + Arange(k)*(Data_num)).^2)));
end
[Fmin,I] = min(F);
% [Gmin,J] = min(G);
% if(Fmin < Gmin)
%     A= Arange(I);
% else
%     A= -1*Arange(J);
% end
A= Arange(I);

%---------------------- END OF CODE---------------------- 

% ================ [amplification_factor.m] ================  
