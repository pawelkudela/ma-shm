function [nodes3D,coords3D]=extrude_mesh(nodes2D,coords2D,NofElNodesx,NofElNodesy,NofElNodesz,h,Nz)
% EXTRUDE_MESH   One line description of what the function or script performs (H1 line) 
%    optional: more details about the function than in the H1 line 
%    optional: more details about the function than in the H1 line 
%    optional: more details about the function than in the H1 line 
% 
% Syntax: [output1,output2] = extrude_mesh(input1,input2,input3) 
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
%    [output1,output2] = extrude_mesh(input1,input2,input3) 
%    [output1,output2] = extrude_mesh(input1,input2) 
%    [output1] = extrude_mesh(input1,input2,input3) 
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

% function extrude_mesh(meshfile,h,Nz)
% meshfile.mat should contain nodes and coords2D (topology and point coordinates)
% h - thickness of extruded solid
% Nz - number of layers in z direction



[fen,NofElNodes]=size(nodes2D);

[ksi,wi]=gll(NofElNodesz);
dofl=max(max(nodes2D));
nodes3D=zeros(Nz*fen,NofElNodesx*NofElNodesy*NofElNodesz); % allocate table for element node numbers (grid topology)
nodes3D(1:fen,1:NofElNodesx*NofElNodesy)=nodes2D;
dof_lay=0;
for i=1:Nz   
    dof_lay=dof_lay-dofl;
    for k=1:NofElNodesz
    dof_lay=dof_lay+dofl;
    nodes3D((i-1)*fen+1:(i-1)*fen+fen,(k-1)*NofElNodesx*NofElNodesy+1:k*NofElNodesx*NofElNodesy)=nodes3D(1:fen,1:NofElNodesx*NofElNodesy)+dof_lay;
    end
end
dof=max(max(nodes3D));
coords3D=zeros(dof,3);


hlay=0;
for iz=1:Nz
    dz=h(iz);
for ine=1:fen
    ne=(iz-1)*fen+ine;
    cc=0;
  
    zz=dz*(ksi+1)/2+hlay;
    
    for jj=1:NofElNodesz   
        c=0;
    for k=1:NofElNodesx
        for n=1:NofElNodesy  
            c=c+1;
            cc=cc+1;
            coords3D(nodes3D(ne,cc),1)=coords2D(nodes2D(ine,c),1);
            coords3D(nodes3D(ne,cc),2)=coords2D(nodes2D(ine,c),2);
            coords3D(nodes3D(ne,cc),3)=zz(jj);
        end
    end
    end
end
hlay=dz;
end

%---------------------- END OF CODE---------------------- 

% ================ [extrude_mesh.m] ================  
