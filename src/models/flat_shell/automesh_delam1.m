clear all; close all;
% generate meshes for various delamination placements, sizes and orientations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare output directories
% allow overwriting existing results if true
overwrite=false;
%overwrite=true;
% retrieve model name based on running file and folder
currentFile = mfilename('fullpath');
[pathstr,name,ext] = fileparts( currentFile );
idx = strfind( pathstr,filesep );
modelfolder = pathstr(idx(end)+1:end); % name of folder
modelname = name; 
% prepare model output path
image_label_path = prepare_model_paths('raw','num',modelfolder,modelname);
figure_output_path = prepare_figure_paths(modelfolder,modelname);
%%  INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=500; % image size in pixels
n=10; % mesh size of delamination positions
% input for automesh_delam
L=0.5; % plate length % for now fixed
W=0.5; % plate width % for now fixed
a_set=[0.01,0.015]; % set of delamination semi-major axis
b_set=[0.005,0.01]; % set of delamination semi-minor axis
rotAngle_set = [0,30,60,90,120,150]; % set of delamination rotation angles [0:180)
r=0.005; % pzt radius
xpzt=0.25; % pzt x coordinate
ypzt=0.25; % pzt y coordinate

% mesh parameters
shape_order = 5; % element shape function order, Number of nodes in one direction is shape_order+1
CharacteristicLengthFactor = 0.08; 
CharacteristicLengthMin = 0.001; 
CharacteristicLengthMax = 0.2;
Smoothing = 1;
%%
% delamination scenario grid
[xn,yn]=meshgrid(0:n/(n-1):n,0:n/(n-1):n);
% x=xn*N/n; % coordinates in pixels
% y=yn*N/n; % coordinates in pixels
xc=xn*L/n; % coordinates in meters
yc=yn*W/n; % coordinates in meters
% plot(x,y,'ro');
% figure;
% plot(xc,yc,'ro');
mesh_parameters = struct('position_no',{},'xCenter',{},'yCenter',{},'a',{},'b',{},'rotAngle',{},'meshfile',{});
counter = 0;
% delamination position: only one quarter
for j=n/2+1:n 
    for i=n/2+1:n
        position_no=(j-1)*n+i; % delamination position number
        % coordinates of delamination position
        xCenter = xc(j,i);
        yCenter = yc(j,i);    
        % mionor axis
        for ib=1:length(b_set)
            b = b_set(ib);
            % major axis
            for ia=1:length(a_set)
                a = a_set(ia);
                if(a == b)
                    counter = counter + 1;
                    rotAngle = 0;
                    mesh_filename = ['m_',num2str(counter),'_delam1_position_no_',num2str(position_no),'_a_',num2str(a*1e3),'mm_b_',num2str(b*1e3),'mm_angle_',num2str(rotAngle)]; 
                    mesh_parameters(counter).position_no = position_no;
                    mesh_parameters(counter).xCenter = xCenter;
                    mesh_parameters(counter).yCenter = yCenter;
                    mesh_parameters(counter).a = a;
                    mesh_parameters(counter).b = b;
                    mesh_parameters(counter).rotAngle = rotAngle;
                    mesh_parameters(counter).meshfile = mesh_filename;
                    figfilename = [figure_output_path,mesh_filename];
                    image_label_filename = [image_label_path,filesep,mesh_filename];
                    delam_image_label(N,xCenter*N/L,yCenter*N/W,a*N/L,b*N/W,rotAngle,image_label_filename);
                    if(overwrite||(~overwrite && ~exist([figfilename,'.png'], 'file')))
                             %% RUN AUTOMESH
                         try
                            disp(mesh_filename);
                            [nodes,coords,den_under,den_above,IG1,IG2,IG3,IG4,IG5,IG6,IG7,IG8,IG9,IG10,IG11,IG12,IL1,IL2,IL3,IL4,IL5,IL6,IL7,IL8,IL9,IL10,IL11,IL12] = automesh_delam...
                            (L,W,a,b,xCenter,yCenter,rotAngle,r,xpzt,ypzt,shape_order,CharacteristicLengthFactor,CharacteristicLengthMin,CharacteristicLengthMax,Smoothing,mesh_filename,modelfolder);
                             %                  
                            print(figfilename,'-dpng', '-r300'); 
                            close all;
                         catch
                            fprintf(['Meshing failed:', mesh_filename,' \n']);
                         end
                    else
                        fprintf(['Mesh:', mesh_filename,' already exist\n']);
                    end
                else
                    for irot = 1:length(rotAngle_set) % rotation angle
                        counter = counter + 1;
                        rotAngle = rotAngle_set(irot);
                        mesh_filename = ['m_',num2str(counter),'_delam1_position_no_',num2str(position_no),'_a_',num2str(a*1e3),'mm_b_',num2str(b*1e3),'mm_angle_',num2str(rotAngle)]; % a, b in mm                       
                        mesh_parameters(counter).position_no = position_no;
                        mesh_parameters(counter).xCenter = xCenter;
                        mesh_parameters(counter).yCenter = yCenter;
                        mesh_parameters(counter).a = a;
                        mesh_parameters(counter).b = b;
                        mesh_parameters(counter).rotAngle = rotAngle;
                        mesh_parameters(counter).meshfile = mesh_filename;
                        figfilename = [figure_output_path,mesh_filename];
                        image_label_filename = [image_label_path,filesep,mesh_filename];
                        delam_image_label(N,xCenter*N/L,yCenter*N/W,a*N/L,b*N/W,rotAngle,image_label_filename);
                        if(overwrite||(~overwrite && ~exist([figfilename,'.png'], 'file')))
                                  %% RUN AUTOMESH
                            try
                                disp(mesh_filename);
                                [nodes,coords,den_under,den_above,IG1,IG2,IG3,IG4,IG5,IG6,IG7,IG8,IG9,IG10,IG11,IG12,IL1,IL2,IL3,IL4,IL5,IL6,IL7,IL8,IL9,IL10,IL11,IL12] = automesh_delam...
                                (L,W,a,b,xCenter,yCenter,rotAngle,r,xpzt,ypzt,shape_order,CharacteristicLengthFactor,CharacteristicLengthMin,CharacteristicLengthMax,Smoothing,mesh_filename,modelfolder);
                                 %   
                                print(figfilename,'-dpng', '-r300'); 
                                close all;
                            catch
                                fprintf(['Meshing failed:', mesh_filename,' \n']);
                            end
                        else
                            fprintf(['Mesh:', mesh_filename,' already exist\n']);
                        end
                    end
                end
            end
        end
    end
end
save([image_label_path,filesep,'mesh_parameters'],'mesh_parameters');