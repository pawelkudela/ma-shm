function [nodes,coords,den_under,den_above] = split_delam_nodes_flat_shell(nodes,coords,den,NofElNodesx,NofElNodesy,boundary_nodes)
% SPLIT_DELAM_NODES_FLAT_SHELL   split nodes in elements den
%    zero-distance double nodes are created 
%    outer nodes are connected with undelaminated elements  
% 
% Syntax: [nodes,coords,den_under,den_above] = split_delam_nodes_flat_shell(nodes,coords,den,NofElNodesx,NofElNodesy)
% 
% Inputs: 
%    nodes - node numbers of elements (topology), integer, dimensions [fen, NofElNodesx*NofElNodesy]
%    coords - coordinates of nodes, double, dimensions [NofNodes, 3], Units: m 
%    den - delamination element numbers, integer list
%    NofElNodesx,NofElNodesy - number of element nodes in x and y direction, respectively
%    boundary_nodes - list of node numbers lying on boundary, integer
% 
% Outputs: 
%    nodes - node numbers of elements (topology) after adding delamination, integer, dimensions [fen, NofElNodesx*NofElNodesy]
%    coords - coordinates of nodes after adding delamination, double, dimensions [NofNodes, 3], Units: m 
%    den_under - delamination element numbers under the split interface
%    den_above - delamination element numbers above the split interface
% 
% Example: 
%    [nodes,coords,den_under,den_above] = split_delam_nodes_flat_shell(nodes,coords,den,NofElNodesx,NofElNodesy)
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

fen = size(nodes,1);
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

den_under = den;
den_above = fen+[1:length(den)];

delamnodes = unique(nodes(den,:));

% plot(coords(reshape(nodes(:,mc),[],1),1),coords(reshape(nodes(:,mc),[],1),2),'b.');
% hold on;
% plot(coords(delamnodes,1),coords(delamnodes,2),'r.');
nodemax=max(max(nodes));
nm=nodemax;
nodes = [nodes; nodes(den,:)];
outernodes = zeros(1,length(den)*NofElNodesx);
free_corner =  zeros(1,length(den));
c=0;cc=0;
ifc=0; % index for free corner nodes 
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
        if(neighbours==1 || neighbours==0) 
            cc=cc+1;
            outernodes(cc)=cn;
        end  
        if(neighbours==0) 
            ifc=ifc+1;
            free_corner(ifc)=cn;
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
free_corner = free_corner(1:ifc);
outernodes=outernodes(1:cc);
special_corner = intersect(free_corner,boundary_nodes);
if(length(special_corner)==3)
    special_corner=special_corner(2:3);
end
outernodes=unique(outernodes);
outernodes=setdiff(outernodes,boundary_nodes);
outernodes=[outernodes,special_corner'];
innernodes=setdiff(delamnodes,outernodes);
% split nodes
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
coords(nodes(den_above,:),:) = coords(nodes(den_under,:),:);
% plot(coords(outernodes,1),coords(outernodes,2),'ko');
% plot(coords(innernodes,1),coords(innernodes,2),'r.');
% plot(coords(special_corner,1),coords(special_corner,2),'mo');
% plot(coords(special_corner(1),1),coords(special_corner(1),2),'mo');
% plot(coords(special_corner(2),1),coords(special_corner(2),2),'co');
%---------------------- END OF CODE---------------------- 

% ================ [split_delam_nodes_flat_shell.m] ================  
