function [coords,nodes]=remove_free_nodes(coords,nodes)
%  This function looks for free nodes (nodal coordinates) and removes then renumbers the nodal coordinates
%  in coords so that the numbering is continous.  Also the connectivity
%  matricies in nodes are adjusted for the new node numbering
%    
  
  [row,col]=size(nodes);
  usedNodes=[];

  for n=1:col
    newNodes=unique(nodes(:,n));
    usedNodes=[usedNodes;newNodes];
  end
  
  usedNodes=unique(usedNodes);
  numnode=size(usedNodes,1);
  
  nodeMap=zeros(size(coords,1),1);
  nodeMap(usedNodes)=(1:numnode)';
   
  % remove unused nodes from node 
  coords=coords(usedNodes,:);  
  
  % renumber connectivity lists
  for n=1:col
    for e=1:size(nodes,1)
      conn(e,n)=nodeMap(nodes(e,n))';
    end
  end
  nodes=conn;
  % Done !!!
  
    