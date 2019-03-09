function [IG1,IG2,IG3,IG4,IG5,IG6,IG7,IG8,IG9,IG10,IG11,IG12,IL1,IL2,IL3,IL4,IL5,IL6,IL7,IL8,IL9,IL10,IL11,IL12]=parallel_LG_nodes_Modified_Zb(nodes)

% calculate local and global node numbers necessary for GPU parallel computation
% IG1:IG12 - vector of global numbers
% IL1:IL12 - vector of local numbers

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

IG1=Basket_global_list(1,:);   IG2=Basket_global_list(2,:);   IG3=Basket_global_list(3,:);
IG4=Basket_global_list(4,:);   IG5=Basket_global_list(5,:);   IG6=Basket_global_list(6,:);
IG7=Basket_global_list(7,:);   IG8=Basket_global_list(8,:);   IG9=Basket_global_list(9,:);
IG10=Basket_global_list(10,:); IG11=Basket_global_list(11,:); IG12=Basket_global_list(12,:);

IL1=Basket_local_list(1,:);   IL2=Basket_local_list(2,:);   IL3=Basket_local_list(3,:);
IL4=Basket_local_list(4,:);   IL5=Basket_local_list(5,:);   IL6=Basket_local_list(6,:);
IL7=Basket_local_list(7,:);   IL8=Basket_local_list(8,:);   IL9=Basket_local_list(9,:);
IL10=Basket_local_list(10,:); IL11=Basket_local_list(11,:); IL12=Basket_local_list(12,:);