
% modelfolder='flat_shell';
% spec_mesh_output_path = fullfile(projectroot,'src','models',modelfolder,'mesh',filesep);

input_no = 56;
run(fullfile('inputs',['input',num2str(input_no)]));

mesh_filename = 'rect_70x70'; 

meshfile=fullfile('mesh',mesh_filename);
load(meshfile);
[fen,NofNodes] = size(nodes);
% load pressure

%load('/pkudela_odroid_laser/acoustic_excitation/Pa_1ACT_25mm_30stp'); %
%Data, dx [m], dy [m], dt [1e3*s]
load('/pkudela_odroid_laser/acoustic_excitation/Pa_1ACT_25mm_0stp'); % Data, dx [m], dy [m], dt [1e3*s]
dt_comsol = dt/1e3; % [s] from comsol
%
dt_sem=tt/(nft-1);   % calculation time step [s]
t_sem=([1:nft]-1)*dt_sem;

% spatial interpolation
disp('spatial interpolation');
x = linspace(0,0.5,size(Data,1));
y = linspace(0,0.5,size(Data,2));
[X,Y] = meshgrid(x,y);
p=zeros(size(coords,1),size(Data,3));
for frame_no=1:size(Data,3)
    v1=reshape(Data(:,:,frame_no),[size(Data,1)*size(Data,2),1]);
    x1=reshape(X,[size(Data,1)*size(Data,2),1]);
    y1=reshape(Y,[size(Data,1)*size(Data,2),1]);
    F = scatteredInterpolant(x1,y1,v1);
    % interpolated pressure
    p(:,frame_no) = F(coords(:,1),coords(:,2));
end

% time interpolation
t1=linspace(0,(size(Data,3)-1)*dt_comsol,size(Data,3)); % from comsol
%pt=zeros(size(coords,1),nft); 
pt=zeros(size(coords,1),1);
t_sem=5000*dt_sem;
[A,I]=min(abs(tc-t1));
if(t_sem<t1(end))
    for kp=1:size(coords,1)

        pt(kp,1) = interp1(t1,p(kp,:),t_sem,'spline');

    end
end
return;
% plot
% po=6;
% Xe=zeros(po,po);
% Ye=zeros(po,po);
% Ze=zeros(po,po);
% for ne=1:fen
%    c=0;
%        
%         cc=0;
%     for j=1:po
%         for i=1:po
%             cc=cc+1;c=c+1;
%             Xe(i,j)=coords(nodes(ne,cc),1);
%             Ye(i,j)=coords(nodes(ne,cc),2);
%             Ze(i,j)=p(nodes(ne,cc),200);
%         end
%     end
%     %plot3(X,Y,Z,'r.');hold on;
%     surf(X,Y,Z); hold on;
%     %pause(0.001);
%     shading interp;
% end

