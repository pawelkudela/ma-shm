function [IG1,IG2,IG3,IG4,IG5,IG6,IG7,IG8,IG9,IG10,IG11,IG12,IL1,IL2,IL3,IL4,IL5,IL6,IL7,IL8,IL9,IL10,IL11,IL12] = global_local_12baskets_5dofs(IG1,IG2,IG3,IG4,IG5,IG6,IG7,IG8,IG9,IG10,IG11,IG12,IL1,IL2,IL3,IL4,IL5,IL6,IL7,IL8,IL9,IL10,IL11,IL12,mode)
% GLOBAL_LOCAL_12BASKETS_5DOFS   convert nodal numbers into degrees of freedom numbers 
% 
% Syntax: [IG1,IG2,IG3,IG4,IG5,IG6,IG7,IG8,IG9,IG10,IG11,IG12,IL1,IL2,IL3,IL4,IL5,IL6,IL7,IL8,IL9,IL10,IL11,IL12] = global_local_12baskets_5dofs(IG1,IG2,IG3,IG4,IG5,IG6,IG7,IG8,IG9,IG10,IG11,IG12,IL1,IL2,IL3,IL4,IL5,IL6,IL7,IL8,IL9,IL10,IL11,IL12,mode)
% 
% Inputs: 
%    IG - global node numbers, integer, dimensions [m, 1]
%    IL - local node numbers, interger, dimensions [m, 1]
%    mode - 'cpu' or 'gpu', string
% 
% Outputs: 
%    IG - global degree of freedom numbers, integer, dimensions [5*m, 1]
%    IL - local degree of freedom numbers, interger, dimensions [5*m, 1]
%         m = fen*NofElNodes/12
% Example: 
%    [IG1,IG2,IG3,IG4,IG5,IG6,IG7,IG8,IG9,IG10,IG11,IG12,IL1,IL2,IL3,IL4,IL5,IL6,IL7,IL8,IL9,IL10,IL11,IL12] = global_local_12baskets_5dofs(IG1,IG2,IG3,IG4,IG5,IG6,IG7,IG8,IG9,IG10,IG11,IG12,IL1,IL2,IL3,IL4,IL5,IL6,IL7,IL8,IL9,IL10,IL11,IL12,mode)
% 
% Other m-files required: none 
% Subfunctions: none 
% MAT-files required: none 
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2 
% 

% Author: Pawel Kudela, D.Sc., Ph.D., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pk@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 

%---------------------- BEGIN CODE---------------------- 

Itemp=zeros(5*length(IG1),1);
Itemp(1:5:end)=5*IG1-4;
Itemp(2:5:end)=5*IG1-3;
Itemp(3:5:end)=5*IG1-2;
Itemp(4:5:end)=5*IG1-1;
Itemp(5:5:end)=5*IG1;
IG1=Itemp;

Itemp=zeros(5*length(IG2),1);
Itemp(1:5:end)=5*IG2-4;
Itemp(2:5:end)=5*IG2-3;
Itemp(3:5:end)=5*IG2-2;
Itemp(4:5:end)=5*IG2-1;
Itemp(5:5:end)=5*IG2;
IG2=Itemp;

Itemp=zeros(5*length(IG3),1);
Itemp(1:5:end)=5*IG3-4;
Itemp(2:5:end)=5*IG3-3;
Itemp(3:5:end)=5*IG3-2;
Itemp(4:5:end)=5*IG3-1;
Itemp(5:5:end)=5*IG3;
IG3=Itemp;

Itemp=zeros(5*length(IG4),1);
Itemp(1:5:end)=5*IG4-4;
Itemp(2:5:end)=5*IG4-3;
Itemp(3:5:end)=5*IG4-2;
Itemp(4:5:end)=5*IG4-1;
Itemp(5:5:end)=5*IG4;
IG4=Itemp;

Itemp=zeros(5*length(IG5),1);
Itemp(1:5:end)=5*IG5-4;
Itemp(2:5:end)=5*IG5-3;
Itemp(3:5:end)=5*IG5-2;
Itemp(4:5:end)=5*IG5-1;
Itemp(5:5:end)=5*IG5;
IG5=Itemp;

Itemp=zeros(5*length(IG6),1);
Itemp(1:5:end)=5*IG6-4;
Itemp(2:5:end)=5*IG6-3;
Itemp(3:5:end)=5*IG6-2;
Itemp(4:5:end)=5*IG6-1;
Itemp(5:5:end)=5*IG6;
IG6=Itemp;

Itemp=zeros(5*length(IG7),1);
Itemp(1:5:end)=5*IG7-4;
Itemp(2:5:end)=5*IG7-3;
Itemp(3:5:end)=5*IG7-2;
Itemp(4:5:end)=5*IG7-1;
Itemp(5:5:end)=5*IG7;
IG7=Itemp;

Itemp=zeros(5*length(IG8),1);
Itemp(1:5:end)=5*IG8-4;
Itemp(2:5:end)=5*IG8-3;
Itemp(3:5:end)=5*IG8-2;
Itemp(4:5:end)=5*IG8-1;
Itemp(5:5:end)=5*IG8;
IG8=Itemp;

Itemp=zeros(5*length(IG9),1);
Itemp(1:5:end)=5*IG9-4;
Itemp(2:5:end)=5*IG9-3;
Itemp(3:5:end)=5*IG9-2;
Itemp(4:5:end)=5*IG9-1;
Itemp(5:5:end)=5*IG9;
IG9=Itemp;

Itemp=zeros(5*length(IG10),1);
Itemp(1:5:end)=5*IG10-4;
Itemp(2:5:end)=5*IG10-3;
Itemp(3:5:end)=5*IG10-2;
Itemp(4:5:end)=5*IG10-1;
Itemp(5:5:end)=5*IG10;
IG10=Itemp;

Itemp=zeros(5*length(IG11),1);
Itemp(1:5:end)=5*IG11-4;
Itemp(2:5:end)=5*IG11-3;
Itemp(3:5:end)=5*IG11-2;
Itemp(4:5:end)=5*IG11-1;
Itemp(5:5:end)=5*IG11;
IG11=Itemp;

Itemp=zeros(5*length(IG12),1);
Itemp(1:5:end)=5*IG12-4;
Itemp(2:5:end)=5*IG12-3;
Itemp(3:5:end)=5*IG12-2;
Itemp(4:5:end)=5*IG12-1;
Itemp(5:5:end)=5*IG12;
IG12=Itemp;

Itemp=zeros(5*length(IL1),1);
Itemp(1:5:end)=5*IL1-4;
Itemp(2:5:end)=5*IL1-3;
Itemp(3:5:end)=5*IL1-2;
Itemp(4:5:end)=5*IL1-1;
Itemp(5:5:end)=5*IL1;
IL1=Itemp;

Itemp=zeros(5*length(IL2),1);
Itemp(1:5:end)=5*IL2-4;
Itemp(2:5:end)=5*IL2-3;
Itemp(3:5:end)=5*IL2-2;
Itemp(4:5:end)=5*IL2-1;
Itemp(5:5:end)=5*IL2;
IL2=Itemp;

Itemp=zeros(5*length(IL3),1);
Itemp(1:5:end)=5*IL3-4;
Itemp(2:5:end)=5*IL3-3;
Itemp(3:5:end)=5*IL3-2;
Itemp(4:5:end)=5*IL3-1;
Itemp(5:5:end)=5*IL3;
IL3=Itemp;

Itemp=zeros(5*length(IL4),1);
Itemp(1:5:end)=5*IL4-4;
Itemp(2:5:end)=5*IL4-3;
Itemp(3:5:end)=5*IL4-2;
Itemp(4:5:end)=5*IL4-1;
Itemp(5:5:end)=5*IL4;
IL4=Itemp;

Itemp=zeros(5*length(IL5),1);
Itemp(1:5:end)=5*IL5-4;
Itemp(2:5:end)=5*IL5-3;
Itemp(3:5:end)=5*IL5-2;
Itemp(4:5:end)=5*IL5-1;
Itemp(5:5:end)=5*IL5;
IL5=Itemp;

Itemp=zeros(5*length(IL6),1);
Itemp(1:5:end)=5*IL6-4;
Itemp(2:5:end)=5*IL6-3;
Itemp(3:5:end)=5*IL6-2;
Itemp(4:5:end)=5*IL6-1;
Itemp(5:5:end)=5*IL6;
IL6=Itemp;

Itemp=zeros(5*length(IL7),1);
Itemp(1:5:end)=5*IL7-4;
Itemp(2:5:end)=5*IL7-3;
Itemp(3:5:end)=5*IL7-2;
Itemp(4:5:end)=5*IL7-1;
Itemp(5:5:end)=5*IL7;
IL7=Itemp;

Itemp=zeros(5*length(IL8),1);
Itemp(1:5:end)=5*IL8-4;
Itemp(2:5:end)=5*IL8-3;
Itemp(3:5:end)=5*IL8-2;
Itemp(4:5:end)=5*IL8-1;
Itemp(5:5:end)=5*IL8;
IL8=Itemp;

Itemp=zeros(5*length(IL9),1);
Itemp(1:5:end)=5*IL9-4;
Itemp(2:5:end)=5*IL9-3;
Itemp(3:5:end)=5*IL9-2;
Itemp(4:5:end)=5*IL9-1;
Itemp(5:5:end)=5*IL9;
IL9=Itemp;

Itemp=zeros(5*length(IL10),1);
Itemp(1:5:end)=5*IL10-4;
Itemp(2:5:end)=5*IL10-3;
Itemp(3:5:end)=5*IL10-2;
Itemp(4:5:end)=5*IL10-1;
Itemp(5:5:end)=5*IL10;
IL10=Itemp;

Itemp=zeros(5*length(IL11),1);
Itemp(1:5:end)=5*IL11-4;
Itemp(2:5:end)=5*IL11-3;
Itemp(3:5:end)=5*IL11-2;
Itemp(4:5:end)=5*IL11-1;
Itemp(5:5:end)=5*IL11;
IL11=Itemp;

Itemp=zeros(5*length(IL12),1);
Itemp(1:5:end)=5*IL12-4;
Itemp(2:5:end)=5*IL12-3;
Itemp(3:5:end)=5*IL12-2;
Itemp(4:5:end)=5*IL12-1;
Itemp(5:5:end)=5*IL12;
IL12=Itemp;

clear Itemp;

switch mode
    case 'gpu'
    IG1=gpuArray(IG1);
    IG2=gpuArray(IG2);
    IG3=gpuArray(IG3);
    IG4=gpuArray(IG4);
    IG5=gpuArray(IG5);
    IG6=gpuArray(IG6);
    IG7=gpuArray(IG7);
    IG8=gpuArray(IG8);
    IG9=gpuArray(IG9);
    IG10=gpuArray(IG10);
    IG11=gpuArray(IG11);
    IG12=gpuArray(IG12);

    IL1=gpuArray(IL1);
    IL2=gpuArray(IL2);
    IL3=gpuArray(IL3);
    IL4=gpuArray(IL4);
    IL5=gpuArray(IL5);
    IL6=gpuArray(IL6);
    IL7=gpuArray(IL7);
    IL8=gpuArray(IL8);
    IL9=gpuArray(IL9);
    IL10=gpuArray(IL10);
    IL11=gpuArray(IL11);
    IL12=gpuArray(IL12);

    case 'cpu'
    
end

%---------------------- END OF CODE---------------------- 

% ================ [global_local_12baskets_5dofs.m] ================  
