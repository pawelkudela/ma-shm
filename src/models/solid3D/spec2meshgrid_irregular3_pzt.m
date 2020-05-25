close all;clear all;

for k_test=[853];

%Nx=400; %Nx=152;% grid size
%Ny=400; %Ny=300;
% 29
Nx=500;%
Ny=500;

surface='bottom'; % options surface='top';surface='bottom';surface='middle';

os='win';
switch os
    case 'win'
        eval(['run inputs\input',num2str(k_test)]); 
        data_filename=['outputs\output',num2str(k_test),'\plate_w_',num2str(k_test),'_',num2str(Nx),'x',num2str(Ny),surface,'.mat'];
    case 'linux'
        eval(['run inputs/input',num2str(k_test)]); 
	 data_filename=['outputs/output',num2str(k_test),'/plate_w_',num2str(k_test),'_',num2str(Nx),'x',num2str(Ny),surface,'.mat'];

    otherwise
        disp('Unknown system');
end
%nft=nft/2; % for testing (half of the time)
%frm_int=2*frm_int;
Data=zeros(Ny,Nx,nft/frm_int);
load(meshfile);
[fen,no]=size(nodes);
%fen=fen-2*length(pztnum);
fen=fen-1*length(pztnum);% no glue
L=max(coords(:,1))-min(coords(:,1));
xi=min(coords(:,1)):L/(Nx-1):max(coords(:,1));
B=max(coords(:,2))-min(coords(:,2));
yi=min(coords(:,2)):B/(Ny-1):max(coords(:,2));
[XI,YI]=meshgrid(xi,yi);
ZA=zeros(max(max(nodes(1:fen,1:nx*ny))),1);
XA=zeros(max(max(nodes(1:fen,1:nx*ny))),1);
YA=zeros(max(max(nodes(1:fen,1:nx*ny))),1);
switch surface
   case 'bottom'
   for ne=1:fen/lay
      XA(nodes(ne,1:nx*ny))=coords(nodes(ne,1:nx*ny),1);
      YA(nodes(ne,1:nx*ny))=coords(nodes(ne,1:nx*ny),2);
   end
   case 'top'
   for ne=(fen/lay)*(lay-1)+1:fen
      XA(nodes(ne,1:nx*ny))=coords(nodes(ne,nx*ny*nz-nx*ny+1:nx*ny*nz),1);
      YA(nodes(ne,1:nx*ny))=coords(nodes(ne,nx*ny*nz-nx*ny+1:nx*ny*nz),2);
   end
   case 'middle'
   for ne=(fen/lay)*(lay-3)+1:fen
      XA(nodes(ne,1:nx*ny))=coords(nodes(ne,nx*ny*nz-nx*ny+1:nx*ny*nz),1);
      YA(nodes(ne,1:nx*ny))=coords(nodes(ne,nx*ny*nz-nx*ny+1:nx*ny*nz),2);
   end
end

sample=0;
for n=frm_int:frm_int:nft; 
    sample=sample+1;
    [sample,nft/frm_int]
switch os
    case 'win'
        filename=['outputs\output',num2str(k_test),'\frame',num2str(n),'.txt'];
        figfilename=['outputs\output',num2str(k_test),'\frame',num2str(n),'.png'];
        figfilename2=['outputs\output',num2str(k_test),'\frame',num2str(n),'.tiff'];
    case 'linux'
        filename=['outputs/output',num2str(k_test),'/frame',num2str(n),'.txt'];
        figfilename=['outputs/output',num2str(k_test),'/frame',num2str(n),'.png'];
        figfilename2=['outputs/output',num2str(k_test),'/frame',num2str(n),'.tiff'];
    otherwise
        disp('Unknown system');
end
Z=load(filename,'-ascii');
%A=sqrt(Z(1:3:end).^2+Z(2:3:end).^2+Z(3:3:end).^2);
%A=sqrt(Z(1:3:end).^2+Z(2:3:end).^2);
%A=Z(1:3:end);%u
%A=Z(2:3:end);%v
A=Z(3:3:end);%w
%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch surface
   case 'bottom'
   for ne=1:fen/lay
      ZA(nodes(ne,1:nx*ny))=A(nodes(ne,1:nx*ny)); % bottom surface
   end
   case 'top'
   for ne=(fen/lay)*(lay-1)+1:fen
      ZA(nodes(ne,1:nx*ny))=A(nodes(ne,nx*ny*nz-nx*ny+1:nx*ny*nz)); % top surface
   end
   case 'middle'
   for ne=(fen/lay)*(lay-3)+1:fen
      ZA(nodes(ne,1:nx*ny))=A(nodes(ne,nx*ny*nz-nx*ny+1:nx*ny*nz)); % middle surface
   end
end

ZI = griddata(XA,YA,ZA,XI,YI,'cubic');    
Data(:,:,sample)=ZI;
end

save(data_filename,'Data');
%save(data_filename,'Data','-v7.3'); % for data files >2GB (compression is used)
pause(0.1);
end