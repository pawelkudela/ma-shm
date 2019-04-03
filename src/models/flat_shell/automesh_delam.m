function [nodes,coords,den_under,den_above,IG1,IG2,IG3,IG4,IG5,IG6,IG7,IG8,IG9,IG10,IG11,IG12,IL1,IL2,IL3,IL4,IL5,IL6,IL7,IL8,IL9,IL10,IL11,IL12] = automesh_delam...
    (L,W,a,b,xCenter,yCenter,rotAngle,r,xpzt,ypzt,shape_order,CharacteristicLengthFactor,CharacteristicLengthMin,CharacteristicLengthMax,Smoothing,mesh_filename,modelfolder)
% AUTOMESH_DELAM   automatic mesh generation including 1 delamination region 
%    it uses gmsh external software 
% 
% Syntax: [nodes,coords,den_under,den_above,IG1,IG2,IG3,IG4,IG5,IG6,IG7,IG8,IG9,IG10,IG11,IG12,IL1,IL2,IL3,IL4,IL5,IL6,IL7,IL8,IL9,IL10,IL11,IL12] = automesh_delam...
%    (L,W,a,b,xCenter,yCenter,rotAngle,r,xpzt,ypzt,shape_order,CharacteristicLengthFactor,CharacteristicLengthMin,CharacteristicLengthMax,Smoothing,mesh_filename,modelfolder)
% 
% Inputs: 
%    L - plate length, double, Units: m 
%    W - plate width, double, Units: m 
%    a - delamination semi-major axis, double, Units: m
%    b - delamination semi-minor axis, double, Units: m
%    xCenter -  delamination x coordinate, double, Units: m
%    yCenter -  delamination y coordinate, double, Units: m
%    rotAngle - delamination rotation angle [0:180), Units: deg
%    r - pzt radius, double, Units: m
%    xpzt - pzt x coordinate, double, Units: m
%    ypzt - pzt y coordinate, double, Units: m
%    shape_order - element shape function order, integer, Number of nodes in one direction is shape_order+1
%    CharacteristicLengthFactor - see gmsh manual 
%    CharacteristicLengthMin - see gmsh manual 
%    CharacteristicLengthMax  - see gmsh manual 
%    Smoothing - see gmsh manual
%    mesh_filename - file name for mesh, string
%    modelfolder - corresponding model folder name
% 
% Outputs: 
%    nodes - nodes of spectral elements (topology), integer, dimensions [NofElements, NofNodes] 
%    coords - coordinates of spectral element nodes, double, dimensions [NofNodes, 3], Units: m 
%    den_under - delaminated element numbers under split interface
%    den_above - delaminated element numbers above split interface
%    IGi,ILi - global and corresponding local node number for parallel computation
% 
% Example: 
%    [nodes,coords,den_under,den_above,IG1,IG2,IG3,IG4,IG5,IG6,IG7,IG8,IG9,IG10,IG11,IG12,IL1,IL2,IL3,IL4,IL5,IL6,IL7,IL8,IL9,IL10,IL11,IL12] = automesh_delam...
%    (L,W,a,b,xCenter,yCenter,rotAngle,r,xpzt,ypzt,shape_order,CharacteristicLengthFactor,CharacteristicLengthMin,CharacteristicLengthMax,Smoothing,mesh_filename,modelfolder)
% 
% Other m-files required: quad2spec_boundary.m, split_delam_nodes_flat_shell.m, parallel_LG_nodes_Modified_Zb.m 
% Subfunctions: none 
% MAT-files required: none 
% geo-files required: all delam1*.geo in geo folder
% See also:
% 

% Author: Pawel Kudela, D.Sc., Ph.D., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pk@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 

%---------------------- BEGIN CODE---------------------- 

% load projectroot path
load project_paths projectroot src_path;

gmsh_path = fullfile(projectroot,'bin','external','gmsh','gmsh ');
mesh_geometry_path = fullfile(projectroot,'src','models',modelfolder,'geo',filesep);
mesh_output_path = fullfile(projectroot,'src','models',modelfolder,'gmsh_out',filesep);
spec_mesh_output_path = fullfile(projectroot,'src','models',modelfolder,'mesh',filesep);
gmsh_options = ' -2 -format auto -v 1 -o '; % non-verbose
%gmsh_options = ' -2 -format auto -o '; % verbose
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ExpectedMaxElementSize = CharacteristicLengthFactor*CharacteristicLengthMax*L;
fprintf('Expected max element size: %f  [mm]\n', ExpectedMaxElementSize*1e3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check delam position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta = rotAngle/180*pi;
d1 = (-yCenter)*sin(theta)*cos(theta);
d2 = (L-xCenter)*sin(theta)*cos(theta);
d3 = (W-yCenter)*sin(theta)*cos(theta);
d4 = (-xCenter)*sin(theta)*cos(theta);
A2 = (sin(theta)^2)/a^2 + (cos(theta)^2)/b^2;
A3 = (cos(theta)^2)/a^2 + (sin(theta)^2)/b^2;
B1 = -2*xCenter*A3+2*d1/a^2-2*d1/b^2;
B2 = -2*yCenter*A2+2*d2/a^2-2*d2/b^2;
B3 = -2*xCenter*A3+2*d3/a^2-2*d3/b^2;
B4 = -2*yCenter*A2+2*d4/a^2-2*d4/b^2;
C1 = (-yCenter)^2*A2 -2*xCenter*d1/a^2 +2*xCenter*d1/b^2 +xCenter^2*A3-1;
C2 = (L-xCenter)^2*A3 - 2*yCenter*d2/a^2 +2*yCenter*d2/b^2 +yCenter^2*A2-1;
C3 = (W-yCenter)^2*A2 -2*xCenter*d3/a^2 +2*xCenter*d3/b^2 +xCenter^2*A3-1;
C4 = (-xCenter)^2*A3 - 2*yCenter*d4/a^2 +2*yCenter*d4/b^2 +yCenter^2*A2-1;
delta1 = B1^2-4*A3*C1;
delta2 = B2^2-4*A2*C2;
delta3 = B3^2-4*A3*C3;
delta4 = B4^2-4*A2*C4;
%[delta1;delta2;delta3;delta4]

delam_position='';
if(delta1>0 && delta2<0 && delta3<0 && delta4<0)
    delam_position = 'edge1';
end
if(delta1<0 && delta2>0 && delta3<0 && delta4<0)
    delam_position = 'edge2';
end
if(delta1<0 && delta2<0 && delta3>0 && delta4<0)
    delam_position = 'edge3';
end
if(delta1<0 && delta2<0 && delta3<0 && delta4>0)
    delam_position = 'edge4';
end
if(delta1>0 && delta2<0 && delta3<0 && delta4>0)
    delam_position = 'corner1';
end
if(delta1>0 && delta2>0 && delta3<0 && delta4<0)
    delam_position = 'corner2';
end
if(delta1<0 && delta2>0 && delta3>0 && delta4<0)
    delam_position = 'corner3';
end
if(delta1<0 && delta2<0 && delta3>0 && delta4>0)
    delam_position = 'corner4';
end
if(delta1<0 && delta2<0 && delta3<0 && delta4<0 && xCenter<L && xCenter>0 && yCenter<W && yCenter>0)
    delam_position = 'inner';
end
%delam_position

switch delam_position
    case 'inner'
        mesh_template = 'delam_1_inner_0';
    case 'edge1'
        if (rotAngle==0)
            mesh_template = 'delam_1_edge1_0'; % mesh template for delamination angle 0 deg (edge1)
        else
            mesh_template = 'delam_1_edge1_1_179'; % mesh template for delamination angles 1-179 deg (edge 1)
        end
    case 'edge2'
        if(rotAngle<90)
            mesh_template = 'delam_1_edge2_0_89'; % mesh template for delamination angles 1-89 deg (edge 2)
        end
        if(rotAngle==90)
            mesh_template = 'delam_1_edge2_90'; % mesh template for delamination angles 90 deg (edge 2)
        end
        if(rotAngle>90)
            mesh_template = 'delam_1_edge2_91_179'; % mesh template for delamination angles 91-179 deg (edge 2)
        end
    case 'edge3'
        if(rotAngle==0)
            mesh_template = 'delam_1_edge3_0'; % mesh template for delamination angle 0 deg (edge3)
        end
        if(rotAngle>0)
            mesh_template = 'delam_1_edge3_1_179';% mesh template for delamination angles 1-179 deg (edge 3)
        end
    case 'edge4'
        if(rotAngle<90)
            mesh_template = 'delam_1_edge4_0_89';% mesh template for delamination angles 1-89 deg (edge 4)
        end
        if(rotAngle==90)
            mesh_template = 'delam_1_edge4_90';% mesh template for delamination angles 90 deg (edge 4)
        end
        if(rotAngle>90)
            mesh_template = 'delam_1_edge4_91_179';% mesh template for delamination angles 91-179 deg (edge 4)
        end
    case 'corner1'
        if(rotAngle==0)
            mesh_template = 'delam_1_corner1_0'; % mesh template for delamination angle 0 deg (corner1)
        end
        if(rotAngle>0 && rotAngle<90)
            mesh_template = 'delam_1_corner1_1_89'; % mesh template for delamination angle 1-89 deg (corner1)
        end
        if(rotAngle==90)
            mesh_template = 'delam_1_corner1_90'; % mesh template for delamination angle 90 deg (corner1)
        end
        if(rotAngle>90)
            mesh_template = 'delam_1_corner1_91_179'; % mesh template for delamination angle 91-179 deg (corner1)
        end
    case 'corner2'
        if(rotAngle>=0 && rotAngle<90)
            mesh_template = 'delam_1_corner2_0_89'; % mesh template for delamination angle 0-89 deg (corner2)
        end
        if(rotAngle==90)
            mesh_template = 'delam_1_corner2_90'; % mesh template for delamination angle 90 deg (corner2)
        end
        if(rotAngle>90)
            mesh_template = 'delam_1_corner2_91_179'; % mesh template for delamination angle 0-89 deg (corner2)
        end
    case 'corner3'
        mesh_template = 'delam_1_corner3_0_179'; % mesh template for delamination angle 0-179 deg (corner3)
    case 'corner4'
        if(rotAngle==0)
            mesh_template = 'delam_1_corner4_0'; % mesh template for delamination angle 0 deg (corner4)
        end
        if(rotAngle>0)
            mesh_template = 'delam_1_corner4_1_179'; % mesh template for delamination angle 0 deg (corner4)
        end
end


gmsh_command = [gmsh_path, mesh_geometry_path, mesh_filename,'.geo', gmsh_options, mesh_output_path, mesh_filename,'.m'];


% read template mesh file
text = fileread([mesh_geometry_path,mesh_template,'.geo']);
fid = fopen([mesh_geometry_path,mesh_filename,'.geo'],'w');
fprintf(fid,'SetFactory("OpenCASCADE");\n');
fprintf(fid,'L=%f; \n',L);
fprintf(fid,'W=%f; \n',W);
fprintf(fid,'a=%f; \n',a);
fprintf(fid,'b=%f; \n',b);
fprintf(fid,'x=%f; \n',xCenter);
fprintf(fid,'y=%f; \n',yCenter);
fprintf(fid,'alpha=%f *Pi; \n',rotAngle/180);
fprintf(fid,'r=%f; \n',r);
fprintf(fid,'xpzt=%f; \n',xpzt);
fprintf(fid,'ypzt=%f; \n',ypzt);
% append  template mesh file
fprintf(fid,text);
fprintf(fid,'\n\n');
fprintf(fid,'Mesh.Algorithm = 6; // Frontal\n');
fprintf(fid,'Mesh.CharacteristicLengthFactor = %f;\n',CharacteristicLengthFactor);
fprintf(fid,'Mesh.CharacteristicLengthMin = %f;\n',CharacteristicLengthMin);
fprintf(fid,'Mesh.CharacteristicLengthMax = %f;\n',CharacteristicLengthMax);
fprintf(fid,'Mesh.RecombineAll = 1; // Apply recombination algorithm\n');
fprintf(fid,'Mesh.SubdivisionAlgorithm = 1; // 1: all quadrangles\n');
fprintf(fid,'Mesh.Smoothing = %d;',Smoothing);
fclose(fid);

% run external gmsh.exe - make mesh file from geometry file
status = system(gmsh_command);
if(status)
    disp('Mesh generation in gmsh failed');
end

% load mesh into matlab
run([mesh_output_path, mesh_filename,'.m']);
plot_mesh(msh);
disp('Quad to spectral mesh conversion...');
%[nodes,coords] = quad2spec(msh.QUADS(:,1:4),msh.POS,shape_order); % my implementation
%[nodes,coords,boundary_nodes] = quad2spec_boundary(msh.QUADS(:,1:4),msh.POS,shape_order);

%  % Piotr Fiborek implementation
%[nodes_pl,coords_pl]= quad2spectral_Fiborek(msh.QUADS(:,1:4),msh.POS,shape_order);
[nodes,coords,boundary_nodes]= quad2spectral_Fiborek(msh.QUADS(:,1:4),msh.POS,shape_order);
% plot(coords(boundary_nodes,1),coords(boundary_nodes,2),'m.');

% delete gmsh out m file
delete([mesh_output_path, mesh_filename,'.m']);
nRegions = max(msh.QUADS(:,5));
nbElements = size(msh.QUADS,1);
Regions = cell(nRegions,1);
c=1;
for k=1:nRegions
    c1=c;
    while(c <= nbElements && msh.QUADS(c,5)==k)
        c=c+1;
    end 
    Regions{k} = c1:c-1;
end
pztEl = Regions{1};
delamEl1 = Regions{2};
mesh_min=msh.MIN;
mesh_max=msh.MAX;
% split delamination nodes
NofElNodesx = shape_order +1;
NofElNodesy = shape_order +1;

%plot(coords(boundary_nodes,1),coords(boundary_nodes,2),'m.');
[nodes,coords,den_under1,den_above1] = split_delam_nodes_flat_shell(nodes,coords,delamEl1,NofElNodesx,NofElNodesy,boundary_nodes);
den_under{1} = den_under1;
den_above{1} = den_above1;
delamEl = [den_under1, den_above1];
disp('12 baskets: calculating local and global node numbers...');
[IG1,IG2,IG3,IG4,IG5,IG6,IG7,IG8,IG9,IG10,IG11,IG12,IL1,IL2,IL3,IL4,IL5,IL6,IL7,IL8,IL9,IL10,IL11,IL12]=parallel_LG_nodes_Modified_Zb(nodes);

save([spec_mesh_output_path,mesh_filename,'.mat'],'nodes','coords','pztEl','delamEl','IG1','IG2','IG3','IG4','IG5','IG6','IG7','IG8','IG9','IG10','IG11','IG12','IL1','IL2','IL3','IL4','IL5','IL6','IL7','IL8','IL9','IL10','IL11','IL12','mesh_min','mesh_max','shape_order','den_under','den_above');

%---------------------- END OF CODE---------------------- 

% ================ [automesh_delam.m] ================  
