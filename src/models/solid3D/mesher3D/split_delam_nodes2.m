function [nodes,coords,den_under,den_above]=split_delam_nodes2(nodes,coords,node_separation_dist,den,dln,fen,NofElNodesx,NofElNodesy)

% function split nodes in elements den in such a way that delamination is
% created between 3D spectral elements (lens shape) 
% node_separation_dist is the meaximum distance between separated nodes [m]
% dln - delaminated layer number 0<dln(i)<Nz
% where Nz  - number of layers in z direction
% fen - number of elements in one layer of the mesh
% nodes - node numbers of each element
% coords - coordinates of nodes
% function is developed only for the case of 3 nodes across the thickness 
% of spectral elements

% corner nodes
mc=[1,NofElNodesx,NofElNodesx*NofElNodesy-NofElNodesx+1,NofElNodesx*NofElNodesy];
% edge nodes
Edge1=2:NofElNodesx-1;
for k=1:NofElNodesy-2
    Edge2(k)=(k+1)*NofElNodesx;
end
Edge3=NofElNodesx*NofElNodesy-NofElNodesx+2:NofElNodesx*NofElNodesy-1;
for k=1:NofElNodesy-2
    Edge4(k)=k*NofElNodesx+1;
end

den_under=zeros(1,length(den));
den_above=zeros(1,length(den));
for k=1:length(den)
    den_under(k)=(dln(k)-1)*fen+den(k);
    den_above(k)=dln(k)*fen+den(k);
end
delamnodes=zeros(length(den)*NofElNodesx*NofElNodesy,1);
for k=1:length(den)
    delamnodes((k-1)*NofElNodesx*NofElNodesy+1:k*NofElNodesx*NofElNodesy)=nodes(den_above(k),1:NofElNodesx*NofElNodesy);
end
delamnodes=unique(delamnodes);
%plot3Dmesh_pzt_new_all(nodes,coords,6,6,3,den_above);
c=0;cc=0;
for k=1:length(den)
    common_nodes=nodes(den_above(k),1:NofElNodesx*NofElNodesy);
    % corner nodes
    for i=1:4
        cn=common_nodes(mc(i));
        % check neighbours
        neighbours=0;
        n=setdiff(1:length(den),k );
        for j=n
            for jc=1:4
                if(nodes(den_above(j),mc(jc))==cn) 
                    neighbours=neighbours+1;
                end
            end
        end
        if(neighbours<=1) 
            cc=cc+1;
            outernodes(cc)=cn;
        end  
    end
    % edge nodes
    % edge 1
    for i=1:length(Edge1) 
        en=common_nodes(Edge1(i));
        % check neighbours
        neighbours=0;
        n=setdiff(1:length(den),k );
        for j=n
            % edge1
            for jc=1:length(Edge1) 
                if(nodes(den_above(j),Edge1(jc))==en) 
                    neighbours=neighbours+1;
                end
            end
            % edge2
            for jc=1:length(Edge2) 
                if(nodes(den_above(j),Edge2(jc))==en) 
                    neighbours=neighbours+1;
                end
            end
             % edge3
            for jc=1:length(Edge3) 
                if(nodes(den_above(j),Edge3(jc))==en) 
                    neighbours=neighbours+1;
                end
            end
             % edge4
            for jc=1:length(Edge4) 
                if(nodes(den_above(j),Edge4(jc))==en) 
                    neighbours=neighbours+1;
                end
            end
        end
        if(neighbours==0) 
            cc=cc+1;
            outernodes(cc)=en;
        end  
    end
    % edge 2
    for i=1:length(Edge2) 
        en=common_nodes(Edge2(i));
        % check neighbours
        neighbours=0;
        n=setdiff(1:length(den),k );
        for j=n
            % edge1
            for jc=1:length(Edge1) 
                if(nodes(den_above(j),Edge1(jc))==en) 
                    neighbours=neighbours+1;
                end
            end
            % edge2
            for jc=1:length(Edge2) 
                if(nodes(den_above(j),Edge2(jc))==en) 
                    neighbours=neighbours+1;
                end
            end
             % edge3
            for jc=1:length(Edge3) 
                if(nodes(den_above(j),Edge3(jc))==en) 
                    neighbours=neighbours+1;
                end
            end
             % edge4
            for jc=1:length(Edge4) 
                if(nodes(den_above(j),Edge4(jc))==en) 
                    neighbours=neighbours+1;
                end
            end
        end
        if(neighbours==0) 
            cc=cc+1;
            outernodes(cc)=en;
        end  
    end
    % edge 3
    for i=1:length(Edge3) 
        en=common_nodes(Edge3(i));
        % check neighbours
        neighbours=0;
        n=setdiff(1:length(den),k );
        for j=n
            % edge1
            for jc=1:length(Edge1) 
                if(nodes(den_above(j),Edge1(jc))==en) 
                    neighbours=neighbours+1;
                end
            end
            % edge2
            for jc=1:length(Edge2) 
                if(nodes(den_above(j),Edge2(jc))==en) 
                    neighbours=neighbours+1;
                end
            end
             % edge3
            for jc=1:length(Edge3) 
                if(nodes(den_above(j),Edge3(jc))==en) 
                    neighbours=neighbours+1;
                end
            end
             % edge4
            for jc=1:length(Edge4) 
                if(nodes(den_above(j),Edge4(jc))==en) 
                    neighbours=neighbours+1;
                end
            end
        end
        if(neighbours==0) 
            cc=cc+1;
            outernodes(cc)=en;
        end  
    end
    % edge 4
    for i=1:length(Edge4) 
        en=common_nodes(Edge4(i));
        % check neighbours
        neighbours=0;
        n=setdiff(1:length(den),k );
        for j=n
            % edge1
            for jc=1:length(Edge1) 
                if(nodes(den_above(j),Edge1(jc))==en) 
                    neighbours=neighbours+1;
                end
            end
            % edge2
            for jc=1:length(Edge2) 
                if(nodes(den_above(j),Edge2(jc))==en) 
                    neighbours=neighbours+1;
                end
            end
             % edge3
            for jc=1:length(Edge3) 
                if(nodes(den_above(j),Edge3(jc))==en) 
                    neighbours=neighbours+1;
                end
            end
             % edge4
            for jc=1:length(Edge4) 
                if(nodes(den_above(j),Edge4(jc))==en) 
                    neighbours=neighbours+1;
                end
            end
        end
        if(neighbours==0) 
            cc=cc+1;
            outernodes(cc)=en;
        end  
    end
end
outernodes=unique(outernodes);
innernodes=setdiff(delamnodes,outernodes);
% split nodes
nodemax=max(max(nodes));
nm=nodemax;
c=0;
for k=1:length(innernodes)
    nodemax=nodemax+1;c=c+1;
    for i=1:length(den)
        for j=1:NofElNodesx*NofElNodesy
            if( nodes(den_above(i),j)==innernodes(k) )
                nodes(den_above(i),j)=nodemax;
                %nodepair(c)=nodes(den_above(i),j);
            end
        end
    end
end
% update coordinates
coords_new=zeros(nodemax,3);
coords_new(1:nm,:)=coords;
coords=coords_new;
minx=min(coords(outernodes,1));maxx=max(coords(outernodes,1));
miny=min(coords(outernodes,2));maxy=max(coords(outernodes,2));
a=(maxx-minx)/2;
b=(maxy-miny)/2;
x0=(maxx+minx)/2;
y0=(maxy+miny)/2;
x=coords(innernodes,1);y=coords(innernodes,2);
z=(((x-x0).^2)./(a^2)+((y-y0).^2)./(b^2) -1)*node_separation_dist/2;
    for k1=1:length(z)
        if(z(k1)>0) z(k1)=0; end;
    end
coords(nm+1:nodemax,3)=coords(innernodes,3)-z;
coords(nm+1:nodemax,1)=coords(innernodes,1);
coords(nm+1:nodemax,2)=coords(innernodes,2);
coords(innernodes,3)=coords(innernodes,3)+z;
for ne=den_under
    coords(nodes(ne,NofElNodesx*NofElNodesy+1:2*NofElNodesx*NofElNodesy),3)=(coords(nodes(ne,1:NofElNodesx*NofElNodesy),3)+coords(nodes(ne,2*NofElNodesx*NofElNodesy+1:3*NofElNodesx*NofElNodesy),3))/2;
end
for ne=den_above
    coords(nodes(ne,NofElNodesx*NofElNodesy+1:2*NofElNodesx*NofElNodesy),3)=(coords(nodes(ne,1:NofElNodesx*NofElNodesy),3)+coords(nodes(ne,2*NofElNodesx*NofElNodesy+1:3*NofElNodesx*NofElNodesy),3))/2;
end
% figure;
% for ne=den_under
%     plot3(coords(nodes(ne,:),1),coords(nodes(ne,:),2),coords(nodes(ne,:),3),'ro');hold on;
% end
% for ne=den_above
%     plot3(coords(nodes(ne,:),1),coords(nodes(ne,:),2),coords(nodes(ne,:),3),'bo');hold on;
% end