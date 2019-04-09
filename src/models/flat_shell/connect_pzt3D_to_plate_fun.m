function [nodes2,coords2]=connect_pzt3D_to_plate_fun(NofElNodesx,NofElNodesy,NofElNodesz,pztEl,nodes,pqs,pzt_t)

[ksi,wi]=gll(NofElNodesz);
npzt=length(pztEl);
[fen,NofElNodes]=size(nodes);
Nodemax=max(max(nodes));
cc=0;
 if(isempty(pztEl))
     
 else
     nodespztlay=[];c=0;
     nodespzt=zeros(npzt,NofElNodesx*NofElNodesy*NofElNodesz);
     for m=1:npzt
         ne=pztEl(m);
         c=c+1;
         nodespztlay=[nodespztlay,nodes(ne,1:NofElNodesx*NofElNodesy)];
         nodespzt(c,1:NofElNodesx*NofElNodesy)=nodes(ne,1:NofElNodesx*NofElNodesy);
     end
     [Bnodes,I,J]=unique(nodespztlay);
     doflay=length(Bnodes);
     for k=2:NofElNodesz
        
        Bnodes=[1:doflay]+Nodemax;
        nodest=Bnodes(J);
        c=0;
        for m=1:npzt
            ne=pztEl(m);
         c=c+1;
         nodespzt(c,NofElNodesx*NofElNodesy*(k-1)+1:NofElNodesx*NofElNodesy*k)=nodest((c-1)*NofElNodesx*NofElNodesy+1:c*NofElNodesx*NofElNodesy);
        end
        Nodemax=Nodemax+doflay;
        
     end
 end
Nodemax=Nodemax-doflay;
coords=zeros(Nodemax,3);
coords(1:length(pqs),:)=pqs;
nodes2=[nodes;nodespzt];
c=0;
for m=1:npzt
    ne=pztEl(m);c=c+1;
    coords(nodes(fen+c,:),1)=coords(nodes(ne,:),1);
    coords(nodes(fen+c,:),2)=coords(nodes(ne,:),2);
   
    for k=1:length(ksi)
        coords(nodes(fen+c,(NofElNodesx*NofElNodesy)*(k-1)+1:k*NofElNodesx*NofElNodesy),3)=coords(nodes(ne,NofElNodesx*NofElNodesy*NofElNodesz-NofElNodesx*NofElNodesy+1:NofElNodesx*NofElNodesy*NofElNodesz),3)+(ksi(k)+1)/2*pzt_t;
    end
    
end
