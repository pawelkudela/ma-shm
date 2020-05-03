clear all; close all;
% generate mesh for single geometry file (1 delamination, square pzt)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare output directories
% allow overwriting existing results if true
%overwrite=false;
overwrite=true;
% retrieve model name based on running file and folder
currentFile = mfilename('fullpath');
[pathstr,name,ext] = fileparts( currentFile );
idx = strfind( pathstr,filesep );
modelfolder = pathstr(idx(end)+1:end); % name of folder
modelname = name; 
% prepare model output path
%image_label_path = prepare_model_paths('raw','num',modelfolder,modelname);
figure_output_path = prepare_figure_paths(modelfolder,modelname);
%%  INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mesh parameters
shape_order = 5; % element shape function order, Number of nodes in one direction is shape_order+1
NofElNodesz = 3; % number of nodes through thickness per layer
nPZT = 18; % number of piezoelectric transducers
Nz = 2; % number of layers
delamination_layer=1;
h = [1.5/4,1.5*3/4]/1000; % thickness of layers [m]
isDelamOn = true; % baseline
stiffener_thickness = 1.5/1000; % [m]
pzt_thickness=0.000254; %[m]
rivet_head_thickness = 3/1000; % [m]
node_separation_dist=0.0001; % 0.1 mm
% Physical Surface in gmsh geo file
Physical_Surface_pzt=[1:18];
Physical_Surface_delam=19;
Physical_Surface_rivet_pin = [20:2:36];
Physical_Surface_rivet_head = [21:2:37];
Physical_Surface_stiffener = [38];
Physical_Surface_plate_left = [39];
Physical_Surface_plate_right = [40];
%%

mesh_filename = 'sensor_nonopt_pzt_stiffener_low_rivets_delam'; 

figfilename = [figure_output_path,mesh_filename];

if(overwrite||(~overwrite && ~exist([mesh,filesep,mesh_filename,'.mat'], 'file')))
         %% RUN AUTOMESH
     try
        disp(mesh_filename);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
 % move later to automesh_multi_pzt_complex_3d    
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
        %  % Piotr Fiborek implementation
        [nodes,coords,boundary_nodes]= quad2spectral_Fiborek(msh.QUADS(:,1:4),msh.POS,shape_order);
        % plot(coords(boundary_nodes,1),coords(boundary_nodes,2),'m.');
        % delete gmsh out m file
        delete([mesh_output_path, mesh_filename,'.m']);
        nRegions = max(msh.QUADS(:,5));
        nbElements = size(msh.QUADS,1);
        
        pztEl_all=[];
        for k=1:length(Physical_Surface_pzt)
            IndexReg = find(msh.QUADS(:,5)==Physical_Surface_pzt(k));
            pztEl{k} = IndexReg;
            pztEl_all=[pztEl_all;pztEl{k}];
        end
        delamEll_all =[];
        if(isDelamOn)
            for k=1:length(Physical_Surface_delam)
                delamEl{k} = find(msh.QUADS(:,5)==Physical_Surface_delam(k));
                delamEll_all=[delamEll_all;delamEl{k}];
            end
        end
        rivet_pinEl_all=[];
        for k=1:length(Physical_Surface_rivet_pin)
            IndexReg = find(msh.QUADS(:,5)==Physical_Surface_rivet_pin(k));
            rivet_pinEl{k} = IndexReg;
            rivet_pinEl_all=[rivet_pinEl_all;rivet_pinEl{k}];
        end  
        rivet_headEl_all=[];
        for k=1:length(Physical_Surface_rivet_head)
            IndexReg = find(msh.QUADS(:,5)==Physical_Surface_rivet_head(k));
            rivet_headEl{k} = IndexReg;
            rivet_headEl_all=[ rivet_headEl_all;rivet_headEl{k}];
        end  
        stiffenerEl_all=[];
        for k=1:length(Physical_Surface_stiffener)
            IndexReg = find(msh.QUADS(:,5)==Physical_Surface_stiffener(k));
            stiffenerEl{k} = IndexReg;
            stiffenerEl_all=[stiffenerEl_all;stiffenerEl{k}];
        end 
        plate_leftEl_all=[];
        for k=1:length(Physical_Surface_plate_left)
            IndexReg = find(msh.QUADS(:,5)==Physical_Surface_plate_left(k));
            plate_leftEl{k} = IndexReg;
            plate_leftEl_all=[plate_leftEl_all; plate_leftEl{k}];
        end 
        plate_rightEl_all=[];
        for k=1:length(Physical_Surface_plate_right)
            IndexReg = find(msh.QUADS(:,5)==Physical_Surface_plate_right(k));
            plate_rightEl{k} = IndexReg;
            plate_rightEl_right=[plate_rightEl_all;plate_rightEl{k} ];
        end 
        % extrude mesh
        NofElNodesx = shape_order+1;
        NofElNodesy = shape_order+1;
       
        [nodes3D,coords3D]=extrude_mesh(nodes,coords,NofElNodesx,NofElNodesy,NofElNodesz,h,Nz);
        fen = size(nodes3D,1);
        fen_per_lay = fen/Nz;
        % inser delamination
        if(isDelamOn)
            dln=zeros(1,length(delamEll_all))+delamination_layer;% delamination layer
            [nodes3D,coords3D,den_under,den_above]=split_delam_nodes2(nodes3D,coords3D,node_separation_dist,delamEll_all,dln,fen_per_lay,NofElNodesx,NofElNodesy);
        end 
        disp('connect stiffener, pzt and rivets ');
        % connect stiffener
        solidEl= [stiffenerEl_all;rivet_pinEl_all;rivet_headEl_all]+(Nz-1)*fen_per_lay;
        [nodes3D,coords3D]=connect_3Dsolid_top(NofElNodesx,NofElNodesy,NofElNodesz,solidEl,nodes3D,coords3D,stiffener_thickness);
        stiffenerEllength=length(stiffenerEl_all);
        solidEllength=length(solidEl);
        [m,n]=size(nodes3D);
        solidElnew=m-solidEllength+1:m;
        stiffenerEl_new = solidElnew(1:length(stiffenerEl_all));
        %plot3Dmesh_pzt(nodes3D,coords3D,NofElNodesx,NofElNodesy,NofElNodesz,stiffenerEl_new);
        
        % connect pzt
        pztEl_all = pztEl_all + (Nz-1)*fen_per_lay;
        [nodes3D,coords3D]=connect_3Dsolid_top(NofElNodesx,NofElNodesy,NofElNodesz,pztEl_all,nodes3D,coords3D, pzt_thickness);
        [m,n]=size(nodes3D);
        pztEllength=length(pztEl_all);
        pztElnew=m-pztEllength+1:m;
        %plot3Dmesh_pzt(nodes3D,coords3D,NofElNodesx,NofElNodesy,NofElNodesz,pztElnew);
     
        % rivet pin
        rivet_pinEl_new = [];
        for k=1:Nz
             rivet_pinEl_new = [rivet_pinEl_new;(k-1)*fen_per_lay+rivet_pinEl_all]; % through host plate
        end
        % through stiffener
        rivet_pinEl_stiffener = [solidElnew(length(stiffenerEl_all)+1:length(stiffenerEl_all)+length(rivet_pinEl_all))]';
        % all rivet pins
        rivet_pinEl_new = [rivet_pinEl_new;rivet_pinEl_stiffener];
        % rivet head elements for extrusion
        rivet_headEl_stiffener = [solidElnew(length(stiffenerEl_all)+length(rivet_pinEl_all)+1:length(solidElnew))]';
        rivet_headEl_top = [rivet_pinEl_stiffener;rivet_headEl_stiffener];
        %plot3Dmesh_pzt(nodes3D,coords3D,NofElNodesx,NofElNodesy,NofElNodesz, rivet_headEl_top);
        % connect rivet top
        [nodes3D,coords3D]=connect_3Dsolid_top(NofElNodesx,NofElNodesy,NofElNodesz,rivet_headEl_top,nodes3D,coords3D, rivet_head_thickness);
        [m,n]=size(nodes3D);
        rivet_headEl_top_new=m-length(rivet_headEl_top)+1:m;
        %plot3Dmesh_pzt(nodes3D,coords3D,NofElNodesx,NofElNodesy,NofElNodesz, rivet_headEl_top_new);

        % connect rivet bottom
        rivet_headEl_bottom = [rivet_pinEl_all;rivet_headEl_all];
        [nodes3D,coords3D]=connect_3Dsolid_bottom2(NofElNodesx,NofElNodesy,NofElNodesz,rivet_headEl_bottom,nodes3D,coords3D,rivet_head_thickness);
        [m,n]=size(nodes3D);
        rivet_headEl_bottom_new=m-length(rivet_headEl_bottom)+1:m;
        rivetEl=[rivet_headEl_bottom_new,rivet_headEl_top_new,rivet_pinEl_new'];
        plot3Dmesh_pzt(nodes3D,coords3D,NofElNodesx,NofElNodesy,NofElNodesz, rivetEl);
        %plot3Dmesh_pzt(nodes3D,coords3D,NofElNodesx,NofElNodesy,NofElNodesz, den_above);
        %plot3Dmesh_pzt(nodes3D,coords3D,NofElNodesx,NofElNodesy,NofElNodesz, den_under);
        
        disp('12 baskets: calculating local and global node numbers...');
        %[IG1,IG2,IG3,IG4,IG5,IG6,IG7,IG8,IG9,IG10,IG11,IG12,IL1,IL2,IL3,IL4,IL5,IL6,IL7,IL8,IL9,IL10,IL11,IL12]=parallel_LG_nodes_Modified_Zb(nodes); % Zhibo implementation
        n_z=NofElNodesz; % number of nodes in thickness direction
        [I_G,I_L]=nodesMaps_Fiborek(nodes3D,size(coords3D,1),shape_order,n_z); % Piotr Fiborek implementation
        save([spec_mesh_output_path,mesh_filename,'.mat'],'nodes3D','coords3D','pztEl','delamEl','rivetEl','I_G','I_L','shape_order','den_under','den_above','msh');

     catch
        fprintf(['Meshing failed:', mesh_filename,' \n']);
     end
else
    fprintf(['Mesh:', mesh_filename,' already exist\n']);
end
                
 