function [nodes3D_new,coords3D_new]=connect_3Dsolid_bottom2(NofElNodesx,NofElNodesy,NofElNodesz,solidEl,nodes3D,coords3D,solid_t)
% CONNECT_3DSOLID_TOP   One line description of what the function or script performs (H1 line) 
%    optional: more details about the function than in the H1 line 
%    optional: more details about the function than in the H1 line 
%    optional: more details about the function than in the H1 line 
% 
% Syntax: [output1,output2] = connect_3Dsolid_top(input1,input2,input3) 
% 
% Inputs: 
%    input1 - Description, string, dimensions [m, n], Units: ms 
%    input2 - Description, logical, dimensions [m, n], Units: m 
%    input3 - Description, double, dimensions [m, n], Units: N 
% 
% Outputs: 
%    output1 - Description, integer, dimensions [m, n], Units: - 
%    output2 - Description, double, dimensions [m, n], Units: m/s^2 
% 
% Example: 
%    [output1,output2] = connect_3Dsolid_top(input1,input2,input3) 
%    [output1,output2] = connect_3Dsolid_top(input1,input2) 
%    [output1] = connect_3Dsolid_top(input1,input2,input3) 
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


[ksi,wi]=gll(NofElNodesz);
nsolid=length(solidEl);
[fen,NofElNodes]=size(nodes3D);
Nodemax=max(max(nodes3D));
cc=0;
 if(isempty(solidEl))
     
 else
     nodes3Dsolidlay=[];c=0;
     nodes3Dsolid=zeros(nsolid,NofElNodesx*NofElNodesy*NofElNodesz);
     for m=1:nsolid
         ne=solidEl(m);
         c=c+1;
         nodes3Dsolidlay=[nodes3Dsolidlay,nodes3D(ne,1:NofElNodesx*NofElNodesy)];
         %nodes3Dsolid(c,1:NofElNodesx*NofElNodesy)=nodes3D(ne,1:NofElNodesx*NofElNodesy);
         nodes3Dsolid(c,NofElNodesx*NofElNodesy*NofElNodesz-NofElNodesx*NofElNodesy+1:NofElNodesx*NofElNodesy*NofElNodesz)=nodes3D(ne,1:NofElNodesx*NofElNodesy);
     end
     [Bnodes3D,I,J]=unique(nodes3Dsolidlay);
     doflay=length(Bnodes3D);
     for k=NofElNodesz-1:-1:1
        
        Bnodes3D=[1:doflay]+Nodemax;
        nodes3Dt=Bnodes3D(J);
        c=0;
        for m=1:nsolid
            ne=solidEl(m);
         c=c+1;
         nodes3Dsolid(c,NofElNodesx*NofElNodesy*(k-1)+1:NofElNodesx*NofElNodesy*k)=nodes3Dt((c-1)*NofElNodesx*NofElNodesy+1:c*NofElNodesx*NofElNodesy);
        end
        Nodemax=Nodemax+doflay;
        
     end
 end
Nodemax=Nodemax-doflay;
coords3D_new=zeros(Nodemax,3);
coords3D_new(1:length(coords3D),:)=coords3D;
nodes3D_new=[nodes3D;nodes3Dsolid];
c=0;
for m=1:nsolid
    ne=solidEl(m);c=c+1;
    coords3D_new(nodes3D_new(fen+c,:),1)=coords3D_new(nodes3D_new(ne,:),1);
    coords3D_new(nodes3D_new(fen+c,:),2)=coords3D_new(nodes3D_new(ne,:),2);
   cc=0;
    for k=length(ksi):-1:1
        cc=cc+1;
        %coords3D_new(nodes3D_new(fen+c,(NofElNodesx*NofElNodesy)*(k-1)+1:k*NofElNodesx*NofElNodesy),3)=0-(ksi(k)+1)/2*solid_t;
        
        coords3D_new(nodes3D_new(fen+c,(NofElNodesx*NofElNodesy)*(k-1)+1:k*NofElNodesx*NofElNodesy),3)=coords3D(nodes3D(ne,1:NofElNodesx*NofElNodesy),3)-(ksi(cc)+1)/2*solid_t;
    end
    
end


%---------------------- END OF CODE---------------------- 

% ================ [connect_3Dsolid_top.m] ================  
