function [nodes,ksi,fen,xe,ye]=cracked_mesh(nx,ny,L,B,NofElNodes,Ncrack,ncpx,ncpy)
% rectangular grid generator

disp('RECTANGULAR GRID GENERATOR...');
colordef white
%***************************************************************************************
%INPUT
%L=tl; % plate length
%B=tb; % plate width
% number of elements in x and y direction
Nx=nx; 
dx=L/Nx;
%Ny=round(B/dx);if(Ny<=0) Ny=1; end;
Ny=ny;
dy=B/Ny;
%NofElNodes=36;
po=sqrt(NofElNodes);
fen=Nx*Ny; % number of elements
% number of "cracked"elements Ncrack
% number of double vertices
if(Ncrack==0) Ndouble=0; else Ndouble=Ncrack*(po-2)+(Ncrack-1); end;
% crack vertice placement
%ncpx=24;ncpy=16;
crackVertStart = [ncpx*dx, ncpy*dy];
crackVertEnd = crackVertStart + [Ncrack*dx,0];
crackX=[crackVertStart(1,1);crackVertEnd(1,1)];
crackY=[crackVertStart(1,2);crackVertEnd(1,2)];
clear crackVertStart crackVertEnd
% element interior points
if(NofElNodes==100) ksi=[-1,-0.9195339081664589,-0.738773865105505,-0.4779249498104445,-0.16527895766638692,0.16527895766638692,0.4779249498104445,0.738773865105505,0.9195339081664589, 1]; end;
if(NofElNodes==81) ksi=[-1,-0.899758,-0.677186, -0.363117, 0.0, 0.363117, 0.677186, 0.899758, 1]; end;
if(NofElNodes==64) ksi=[-1,-0.87174,-0.5917, -0.209299, 0.209299, 0.5917, 0.87174, 1]; end;
if(NofElNodes==49) ksi=[-1,-0.830224,-0.468849, 0.0, 0.468849, 0.830224, 1]; end;
if(NofElNodes==36) ksi=[-1,-0.765055,-0.285232, 0.285232, 0.765055, 1]; end;
if(NofElNodes==25) ksi=[-1,-0.654654, 0.0, 0.654654, 1]; end;
if(NofElNodes==16) ksi=[-1,-0.447214, 0.447214, 1]; end;
if(NofElNodes==9)  ksi=[-1, 0.0, 1]; end;
disp('allocating tables...');
nodes=zeros(fen,NofElNodes); % allocate table for element nodes numbers (grid topology)
NofNodes=(Nx*(NofElNodes-po)+po)*Ny-(Ny-1)*(Nx*(po-1)+1)+Ndouble;
Coord=zeros(NofNodes,2); % allocate table for nodal coordinate
%***************************************************************************************
%                       GRID GENERATOR
%***************************************************************************************
disp('generating mesh...');
delta=po-1;
countNodes=0;
%***************************************************************************************
% TOPOLOGY
%***************************************************************************************
% one element nodes numbers
for n=1:po
    for m=1:po
        nr=m+(n-1)*po;
        nodes(1,nr)= m+(n-1)*(po-1)*Nx+n-1;
    end;  
end;
% one row element nodes numbers
count=1;
    for ne=2:Nx
        for n=1:NofElNodes
            nodes(ne,n) = nodes(ne-1,n) + delta;
        end;
    count=count+1;
    end;
% remaining row element nodes numbers
delta=(po-1)*((po-1)*Nx+1);    
for m=2:Ny
    for ne=1:Nx
        if((m-1)*Nx+ne>ncpy*Nx) crackT=Ndouble; else crackT=0; end;
        if(m>ncpy+1)
            for n=1:po
                nodes((m-1)*Nx+ne,n) = nodes(ne,n) + (m-1)*delta+crackT;
            end;
        else
            for n=1:po
                nodes((m-1)*Nx+ne,n) = nodes(ne,n) + (m-1)*delta;
            end;
        end;
        for n=po+1:NofElNodes
            nodes((m-1)*Nx+ne,n) = nodes(ne,n) + (m-1)*delta + crackT;
        end;
    count=count+1;
    end;
end;
% modify cracked element numbers - double nodes
eleref=ncpy*Nx+ncpx;
cc=0;
cDelta=(Nx-Ncrack-ncpx)*(po-1)+Ndouble+1;
for n=1:Ncrack
    for m = 1:po
        cc=cc+1;
        if(n==1& m>1) nodes(eleref+n,m) = nodes(eleref+n,m)+cDelta;end;
        if(cc<2*po & n~=1 & m<po) nodes(eleref+n,m) = nodes(eleref+n,m)+cDelta;end;
    end;  
end;
% %***************************************************************************************
% % CORNER ELEMENTS TOPOLOGY - only for control drawing 
% %***************************************************************************************
% CornerNodes(1,1)= 1;
% CornerNodes(1,2)= 2;
% CornerNodes(1,3)= 3+Nx;
% CornerNodes(1,4)= 2+Nx;
% % remaining element nodes numbers
% count=1;
% for ne=2:Nx*Ny
%     if(count==Nx) delta=(2-1)*((2-1)*Nx+1) - (Nx-1)*(2-1);count=0; else delta=2-1; end;
%     for n=1:4
%         CornerNodes(ne,n) = CornerNodes(ne-1,n) + delta;
%     end;
%     count=count+1;
% end;
% %***************************************************************************************
% % COORDINATES
% %***************************************************************************************
% % corners element coordinates
for iy=1:Ny;
  for ix=1:Nx;
    ne=Nx*(iy-1)+ix;
    xe(ne)=(ix-1)*dx;
    ye(ne)=(iy-1)*dy;
  end;
end;
% % nodal corner element coordinates - only for control drawing 
% countNodes=0;
% for m=1:Ny+1
%     for n=1:Nx+1
%         countNodes=countNodes+1;
%         CornElx(countNodes,1)=(n-1)*dx;
%         CornEly(countNodes,1)=(m-1)*dy;      
%     end;
% end
% % nodes global coordinate for each element
% %x=zeros(fen,po);
% %y=zeros(fen,po);
% %for ne=1:fen
% %    x(ne,:)=dx*(ksi+1)/2+xe(ne);
% %    y(ne,:)=dy*(ksi+1)/2+ye(ne);
% %end;
% % convert coordinates matrices into vectors
% for ne=1:fen
%     cc=0;
%     xx=dx*(ksi+1)/2+xe(ne);
%     yy=dy*(ksi+1)/2+ye(ne);
%     for k=1:po
%         for n=1:po
%             cc=cc+1;
%             Coord(nodes(ne,cc),1)=xx(n);
%             Coord(nodes(ne,cc),2)=yy(k);
%         end;
%     end;
% end;
% 
% % plot nodes mesh
% patch('Vertices',[CornElx,CornEly],'Faces',CornerNodes,'FaceColor','y','EdgeColor','r');axis('equal');
% title('rectangular mesh');
% hold on;
% line(crackX,crackY,'Color','g','LineWidth',3); 
% hold on;
% plot(Coord(:,1),Coord(:,2),'.'); 
% %***************************************************************************************
% %   END OF CORNER ELEMENTS TOPOLOGY - only for control drawing 
% %***************************************************************************************

clear B L cc cDelta Nx Ny eleref Coord Ndouble Ncrack count countNodes crackT delta
clear ncpx ncpy nr po

disp('END OF RECTANGULAR GRID GENERATOR...');
return;