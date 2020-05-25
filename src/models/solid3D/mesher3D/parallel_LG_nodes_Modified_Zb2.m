function [I_G,I_L]=parallel_LG_nodes_Modified_Zb2(nodes)

% calculate local and global node numbers necessary for GPU parallel computation
% I_G - vector of global numbers
% I_L - vector of local numbers

NofNodes=max(max(nodes));
[fen,NofElNodes]=size(nodes);
NodesMap(NofNodes).index=0;
NodesMap(NofNodes).element_number=[];
NodesMap(NofNodes).local_number=[];
for k=1:NofNodes-1
    NodesMap(k).index=0;
    NodesMap(k).element_number=[];
    NodesMap(k).local_number=[];
end
for ne=1:fen
    for i=1:NofElNodes
        nod=nodes(ne,i);
        NodesMap(nod).index=NodesMap(nod).index+1;
        idx=NodesMap(nod).index;
        NodesMap(nod).element_number(idx)=ne;
        NodesMap(nod).local_number(idx)=i;
    end
end
%% division of nodes into 12 baskets
Basket_global_list=zeros(12,fen*NofElNodes/12);
Basket_local_list=zeros(12,fen*NofElNodes/12);
Basket_id=zeros(1,12);
%disp('12 baskets: calculating local and global node numbers...')
% loop over all nodes
for nod=1:NofNodes 
    [Y,I]=sort(Basket_id);
    
    for idx=1:NodesMap(nod).index
        ne=NodesMap(nod).element_number(idx);
        index=Basket_id(I(idx))+1;
        Basket_id(I(idx))=index;
        Basket_global_list(I(idx),index)=nod;
        node_num=NodesMap(nod).local_number(idx);
        Basket_local_list(I(idx),index)=(ne-1)*NofElNodes+node_num;
    end
end
I_G=Basket_global_list';
I_L=Basket_local_list';
