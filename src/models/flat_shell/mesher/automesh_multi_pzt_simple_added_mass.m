function [nodes,coords,addedMassEl,I_G,I_L] = automesh_multi_pzt_simple_added_mass(nPZT,shape_order,mesh_filename,modelfolder,figfilename,isAddedMassOn)
% AUTOMESH_MULTI_PZT_SIMPLE   automatic mesh generation including 1 delamination region based on geo file
%    for the case of multi pzts; it uses gmsh external software 
% 
% Syntax: [nodes,coords,den_under,den_above,I_G,I_L] = automesh_multi_pzt_simple_added_mass(mesh_filename,modelfolder,figfilename)
% 
% Inputs: 
%    shape_order - element shape function order, integer, Number of nodes in one direction is shape_order+1
%    mesh_filename - file name for mesh, string
%    modelfolder - corresponding model folder name
%    figfilename - file name for png figure of mesh
%    isAddedMassOn - boolean, if true - save element no in which added mass should be placed
% 
% Outputs: 
%    nodes - nodes of spectral elements (topology), integer, dimensions [NofElements, NofNodes] 
%    coords - coordinates of spectral element nodes, double, dimensions [NofNodes, 3], Units: m 
%    I_G,I_L - global and corresponding local node number for parallel computation
%    addedMassEl - added mass element number, integer
% 
% Example: 
%    [nodes,coords,den_under,den_above,I_G,I_L] = automesh_multi_pzt_simple_added_mass(mesh_filename,modelfolder,figfilename)
% 
% Other m-files required: quad2spectral_Fiborek.m, split_delam_nodes_flat_shell.m, nodesMaps_Fiborek.m 
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

gmsh_command = [gmsh_path, mesh_geometry_path, mesh_filename,'.geo', gmsh_options, mesh_output_path, mesh_filename,'.m'];

% run external gmsh.exe - make mesh file from geometry file
status = system(gmsh_command);
if(status)
    disp('Mesh generation in gmsh failed');
end

% load mesh into matlab
run([mesh_output_path, mesh_filename,'.m']);
[nodes1,coords1]=change_turn_quad(msh.QUADS(:,1:4),msh.POS(:,1:2));
msh.QUADS(:,1:4) = nodes1;
plot_mesh(msh,figfilename);

close all;
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
for k=1:nPZT
    pztEl{k} = Regions{k};
end
if(isAddedMassOn)
    addedMassEl = Regions{k+1};
else
    addedMassEl = [];
end
mesh_min=msh.MIN;
mesh_max=msh.MAX;
% split delamination nodes
NofElNodesx = shape_order +1;
NofElNodesy = shape_order +1;

%plot(coords(boundary_nodes,1),coords(boundary_nodes,2),'m.');

    delamEl=[];
    den_under=[];
    den_above=[];

disp('12 baskets: calculating local and global node numbers...');
%[IG1,IG2,IG3,IG4,IG5,IG6,IG7,IG8,IG9,IG10,IG11,IG12,IL1,IL2,IL3,IL4,IL5,IL6,IL7,IL8,IL9,IL10,IL11,IL12]=parallel_LG_nodes_Modified_Zb(nodes); % Zhibo implementation
n_z=1; % number of nodes in thickness direction
[I_G,I_L]=nodesMaps_Fiborek(nodes,size(coords,1),shape_order,n_z); % Piotr Fiborek implementation

%save([spec_mesh_output_path,mesh_filename,'.mat'],'nodes','coords','pztEl','delamEl','IG1','IG2','IG3','IG4','IG5','IG6','IG7','IG8','IG9','IG10','IG11','IG12','IL1','IL2','IL3','IL4','IL5','IL6','IL7','IL8','IL9','IL10','IL11','IL12','mesh_min','mesh_max','shape_order','den_under','den_above');
save([spec_mesh_output_path,mesh_filename,'added_mass.mat'],'nodes','coords','pztEl','delamEl','I_G','I_L','mesh_min','mesh_max','shape_order','den_under','den_above','addedMassEl','msh');

%---------------------- END OF CODE---------------------- 

% ================ [automesh_multi_pzt_simple_added_mass.m] ================  
