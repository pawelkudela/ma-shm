% testing automeshing with gmsh

% load projectroot path
load project_paths projectroot src_path;

% mesh_template = 'delam_1_inner_0';
% mesh_filename = 'delam_1_inner_1';

%mesh_template = 'delam_1_edge1_1_179'; % mesh template for delamination angles 1-179 deg (edge 1)
% mesh_template = 'delam_1_edge1_0'; % mesh template for delamination angle 0 deg (edge1)
% mesh_filename = 'delam_1_edge1_1';
%mesh_template = 'delam_1_edge2_0_89'; % mesh template for delamination angles 1-89 deg (edge 2)
%mesh_template = 'delam_1_edge2_90'; % mesh template for delamination angles 90 deg (edge 2)
%mesh_template = 'delam_1_edge2_91_179'; % mesh template for delamination angles 91-179 deg (edge 2)
%mesh_filename = 'delam_1_edge2_1';
%mesh_template = 'delam_1_edge3_0'; % mesh template for delamination angle 0 deg (edge3)
%mesh_template = 'delam_1_edge3_1_179';% mesh template for delamination angles 1-179 deg (edge 3)
%mesh_filename = 'delam_1_edge3_1';
%mesh_template = 'delam_1_edge4_0_89';% mesh template for delamination angles 1-89 deg (edge 4)
%mesh_template = 'delam_1_edge4_90';% mesh template for delamination angles 90 deg (edge 4)
% mesh_template = 'delam_1_edge4_91_179';% mesh template for delamination angles 91-179 deg (edge 4)
% mesh_filename = 'delam_1_edge4_1';


%mesh_template = 'delam_1_corner1_0'; % mesh template for delamination angle 0 deg (corner1)
%mesh_template = 'delam_1_corner1_1_89'; % mesh template for delamination angle 1-89 deg (corner1)
%mesh_template = 'delam_1_corner1_90'; % mesh template for delamination angle 90 deg (corner1)
%mesh_template = 'delam_1_corner1_91_179'; % mesh template for delamination angle 91-179 deg (corner1)
% mesh_filename = 'delam_1_corner1_1';

%mesh_template = 'delam_1_corner2_0_89'; % mesh template for delamination angle 0-89 deg (corner2)
%mesh_template = 'delam_1_corner2_90'; % mesh template for delamination angle 0-89 deg (corner2)
% mesh_template = 'delam_1_corner2_91_179'; % mesh template for delamination angle 0-89 deg (corner2)

%mesh_template = 'delam_1_corner3_0_179'; % mesh template for delamination angle 0-179 deg (corner3)

%mesh_template = 'delam_1_corner4_0'; % mesh template for delamination angle 0 deg (corner4)
mesh_template = 'delam_1_corner4_1_179'; % mesh template for delamination angle 0 deg (corner4)

mesh_filename = 'delam_1_corner4_1';

modelfolder = 'flat_shell';
gmsh_path = fullfile(projectroot,'bin','external','gmsh','gmsh ');
mesh_geometry_path = fullfile(projectroot,'src','models',modelfolder,'geo',filesep);
mesh_output_path = fullfile(projectroot,'src','models',modelfolder,'gmsh_out',filesep);
%gmsh_options = ' -2 -format auto -v 1 -o '; % non-verbose
gmsh_options = ' -2 -format auto -o '; % verbose
gmsh_command = [gmsh_path, mesh_geometry_path, mesh_filename,'.geo', gmsh_options, mesh_output_path, mesh_filename,'.m'];

%%%%%%%%%%%%%%%%%%%%%%%%%%
a=0.02; % delamination semi-major axis
b=0.01; % delamination semi-minor axis
x=0.;   % delamination x coordinate
y=0.5;   % delamination y coordinate
alpha = 179; % delamination rotation
r=0.005; % pzt radius
xpzt=0.25; % pzt x coordinate
ypzt=0.25; % pzt y coordinate

% mesh parameters
CharacteristicLengthFactor = 0.08; 
CharacteristicLengthMin = 0.001; 
CharacteristicLengthMax = 0.2;
Smoothing = 1;

% read template mesh file
text = fileread([mesh_geometry_path,mesh_template,'.geo']);
fid = fopen([mesh_geometry_path,mesh_filename,'.geo'],'w');
fprintf(fid,'SetFactory("OpenCASCADE");\n');
fprintf(fid,'a=%f; \n',a);
fprintf(fid,'b=%f; \n',b);
fprintf(fid,'x=%f; \n',x);
fprintf(fid,'y=%f; \n',y);
fprintf(fid,'alpha=%f *Pi; \n',alpha/180);
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
pause;
close all;
% figfilename = 'test_mesh';
% print([figfilename],'-dpng', '-r600'); 