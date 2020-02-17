function [spec_element_nodes,spec_coords] = quad2spec(q,quad_coords,N)
% QUAD2SPEC   convert quad nodes mesh into spectral element mesh 
%    only for linear quad elements 
% 
% Syntax: [spec_element_nodes,spec_coords] = quad2spec(q,quad_coords,N) 
% 
% Inputs: 
%    q - Quad nodes topology (element nodes), integer, dimensions [nQuadElements, 4]
%    quad_coords - coordinates of quad element nodes, double, dimensions [nQuadNodes, 3], Units: m 
%    N - element approximation order, integer (N=3,4,5,6,7,8,9)
% 
% Outputs: 
%    spec_element_nodes - spectral elements topology (element nodes), integer, dimensions [nQuadElements,(N+1)^2]
%    spec_coords - coordinates of spectral element nodes, integer, dimensions [nSpecNodes, 3], Units: m 
% 
% Example: 
%    [spec_element_nodes,spec_coords] = quad2spec(q,quad_coords,N)
% 
% Other m-files required: gll.m, spec_nodes.m 
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

n=N+1; % number of nodes in 1D element
% coordinates of the nodes and weights
[ksi,~]=gll(n);

nSpecNodes = n^2;

nQuadElements=length(q);
nQuadNodes=size(quad_coords,1);

% connectivity matrix
A=sparse(4,4);
A(1,2)=1;A(1,4)=1;A(2,3)=1;A(3,4)=1;A(4,1)=1;
A(2,1)=1;A(3,2)=1;A(4,3)=1;
K=sparse(nQuadNodes,nQuadNodes);
%disp('calculate connectivity matrix...')
for k=1:nQuadElements
    iout=sparse(4,nQuadNodes);
    iout(1,q(k,1))=1;
    iout(2,q(k,2))=1;
    iout(3,q(k,3))=1;
    iout(4,q(k,4))=1;
    K=K+iout'*A*iout;
end
spec_element_nodes=zeros(nQuadElements,nSpecNodes);
spec_coords=zeros(nQuadElements*nSpecNodes,3);

division=n-1;
nodeTabOld=zeros(division+1, division+1);
cc=0;
for i=1:division+1
    for j=1:division+1
        cc=cc+1;
        nodeTabOld(i,j)=cc;
    end
end
% corner node numbers
mc(1,1)=nodeTabOld(1,1);
mc(1,2)=nodeTabOld(1,division+1);
mc(1,3)=nodeTabOld(division+1,division+1);
mc(1,4)=nodeTabOld(division+1,1);

% outer node numbers without corner node numbers
mo=[nodeTabOld(1,2:division),nodeTabOld(2:division,division+1)',nodeTabOld(division+1,division:-1:2),nodeTabOld(division:-1:2,1)'];
% inner node numbers
cc=0;
mi=zeros(1,(division-1)*(division-1));

for i=2:division
    for j=2:division
        cc=cc+1;
        mi(1,cc)=nodeTabOld(i,j);
    end
end

mol=length(mo)/4;
mil=length(mi);
% corner nodes
spec_coords(1:nQuadNodes,:)=quad_coords;
spec_element_nodes(1:nQuadElements,mc)=q;
cc=0;ck=length(quad_coords);
m=[1,2,3,4,1];
Kp(nQuadNodes,nQuadNodes).q=[1:mol];
%disp('mesh conversion...');
for k=1:nQuadElements
   % coordinates of the nodes of the element 
   [X,Y]=spec_nodes(ksi,quad_coords(q(k,:),1),quad_coords(q(k,:),2));
        cc=cc+1;
        % edge nodes and coordinates
     for j=1:4
          c=(j-1)*mol;
         switch K(q(k,m(j)),q(k,m(j+1)))     
             case 1
             for i=1:mol
                 c=c+1;
                 ck=ck+1;
                 spec_element_nodes(cc,mo(c))=ck;
                 spec_coords(ck,1)=X(mo(c));spec_coords(ck,2)=Y(mo(c));
             end
             case 2
             K(q(k,m(j+1)),q(k,m(j)))=3;
             for i=1:mol
                 c=c+1;
                 ck=ck+1;
                 spec_element_nodes(cc,mo(c))=ck;
                 spec_coords(ck,1)=X(mo(c));spec_coords(ck,2)=Y(mo(c));
                 Kp(q(k,m(j+1)),q(k,m(j))).q(i)=ck;
             end
             case 3
             c1=mol;
             for i=1:mol
                 c=c+1;
                 spec_element_nodes(cc,mo(c))=Kp(q(k,m(j)),q(k,m(j+1))).q(c1);
                 c1=c1-1;
             end
         end
     end
     % inner nodes
     cm=4*mol+4;
     for i=1:mil
         ck=ck+1;
         cm=cm+1;
         spec_element_nodes(cc,mi(i))=ck;
         spec_coords(ck,1)=X(mi(i)); spec_coords(ck,2)=Y(mi(i));
     end
end
spec_coords=(spec_coords(1:ck,:));


%---------------------- END OF CODE---------------------- 

% ================ [quad2spec.m] ================  
