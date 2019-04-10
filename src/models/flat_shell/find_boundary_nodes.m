function [boundary_nodes] = find_boundary_nodes(nodes,NofElNodesx,NofElNodesy)
% FIND_BOUNDARY_NODES   find nodes lying on boundary
% 
% Syntax: [boundary_nodes] = find_boundary_nodes(nodes,coords,NofElNodesx,NofElNodesy)
% 
% Inputs: 
%    nodes - node numbers of elements (topology), integer, dimensions [fen, NofElNodesx*NofElNodesy]
%    NofElNodesx,NofElNodesy - number of element nodes in x and y direction, respectively
% 
% Outputs: 
%    boundary_nodes - list of node numbers lying on boundary, integer
% 
% Example: 
%    [boundary_nodes] = find_boundary_nodes(nodes,coords,NofElNodesx,NofElNodesy) 
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

all_elements = 1:fen;

all_nodes = unique(nodes(all_elements,:));

% plot(coords(reshape(nodes(:,mc),[],1),1),coords(reshape(nodes(:,mc),[],1),2),'b.');
% hold on;
% plot(coords(all_nodes,1),coords(all_nodes,2),'r.');

c=0;cc=0;
for k=1:fen
    %[k,fen]
    common_nodes=nodes(all_elements(k),1:NofElNodesx*NofElNodesy);
    % corner nodes
    for i=1:4
        cn=common_nodes(mc(i));
        % check neighbours
        neighbours=0;
        n=setdiff(1:length(all_elements),k );
        for j=n
            for jc=1:4
                if(nodes(all_elements(j),mc(jc))==cn) 
                    neighbours=neighbours+1;
                end
            end
        end
        if(neighbours==1 || neighbours==0) 
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
        n=setdiff(1:length(all_elements),k );
        for j=n
            % edge1
            for jc=1:length(Edge1) 
                if(nodes(all_elements(j),Edge1(jc))==en) 
                    neighbours=neighbours+1;
                end
            end
            % edge2
            for jc=1:length(Edge2) 
                if(nodes(all_elements(j),Edge2(jc))==en) 
                    neighbours=neighbours+1;
                end
            end
             % edge3
            for jc=1:length(Edge3) 
                if(nodes(all_elements(j),Edge3(jc))==en) 
                    neighbours=neighbours+1;
                end
            end
             % edge4
            for jc=1:length(Edge4) 
                if(nodes(all_elements(j),Edge4(jc))==en) 
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
        n=setdiff(1:length(all_elements),k );
        for j=n
            % edge1
            for jc=1:length(Edge1) 
                if(nodes(all_elements(j),Edge1(jc))==en) 
                    neighbours=neighbours+1;
                end
            end
            % edge2
            for jc=1:length(Edge2) 
                if(nodes(all_elements(j),Edge2(jc))==en) 
                    neighbours=neighbours+1;
                end
            end
             % edge3
            for jc=1:length(Edge3) 
                if(nodes(all_elements(j),Edge3(jc))==en) 
                    neighbours=neighbours+1;
                end
            end
             % edge4
            for jc=1:length(Edge4) 
                if(nodes(all_elements(j),Edge4(jc))==en) 
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
        n=setdiff(1:length(all_elements),k );
        for j=n
            % edge1
            for jc=1:length(Edge1) 
                if(nodes(all_elements(j),Edge1(jc))==en) 
                    neighbours=neighbours+1;
                end
            end
            % edge2
            for jc=1:length(Edge2) 
                if(nodes(all_elements(j),Edge2(jc))==en) 
                    neighbours=neighbours+1;
                end
            end
             % edge3
            for jc=1:length(Edge3) 
                if(nodes(all_elements(j),Edge3(jc))==en) 
                    neighbours=neighbours+1;
                end
            end
             % edge4
            for jc=1:length(Edge4) 
                if(nodes(all_elements(j),Edge4(jc))==en) 
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
        n=setdiff(1:length(all_elements),k );
        for j=n
            % edge1
            for jc=1:length(Edge1) 
                if(nodes(all_elements(j),Edge1(jc))==en) 
                    neighbours=neighbours+1;
                end
            end
            % edge2
            for jc=1:length(Edge2) 
                if(nodes(all_elements(j),Edge2(jc))==en) 
                    neighbours=neighbours+1;
                end
            end
             % edge3
            for jc=1:length(Edge3) 
                if(nodes(all_elements(j),Edge3(jc))==en) 
                    neighbours=neighbours+1;
                end
            end
             % edge4
            for jc=1:length(Edge4) 
                if(nodes(all_elements(j),Edge4(jc))==en) 
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
boundary_nodes=unique(outernodes); 

%---------------------- END OF CODE---------------------- 

% ================ [find_boundary_nodes.m] ================  
