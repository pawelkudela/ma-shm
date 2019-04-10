function [nodespzt,coordspzt]=connect_pzt3D_to_plate_fun(NofElNodesx,NofElNodesy,NofElNodesz,pztEl,nodes,coords,pzt_t,h)

[ksi,wi]=gll(NofElNodesz);
npzt=length(pztEl);
[fen,NofElNodes]=size(nodes);
Nodemax=max(max(nodes));

 if(isempty(pztEl))
     
 else
     c=0;
     nodespzt=zeros(npzt,NofElNodesx*NofElNodesy*NofElNodesz);
     nodespzt(1:length(pztEl),1:NofElNodesx*NofElNodesy)=nodes(pztEl,1:NofElNodesx*NofElNodesy);
     [Bnodes,I,J]=unique(nodespzt(1:length(pztEl),1:NofElNodesx*NofElNodesy)');
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

coordspzt=zeros(Nodemax,3);

c=0;
for m=1:npzt
    ne=pztEl(m);c=c+1;
    for k=1:length(ksi)
        coordspzt(nodespzt(c,NofElNodesx*NofElNodesy*(k-1)+1:k*NofElNodesx*NofElNodesy),1)=coords(nodes(ne,:),1);
        coordspzt(nodespzt(c,NofElNodesx*NofElNodesy*(k-1)+1:k*NofElNodesx*NofElNodesy),2)=coords(nodes(ne,:),2);
        coordspzt(nodespzt(c,(NofElNodesx*NofElNodesy)*(k-1)+1:k*NofElNodesx*NofElNodesy),3)=coords(nodes(ne,1:NofElNodesx*NofElNodesy),3)+(ksi(k)+1)/2*pzt_t+h/2;
    end
    
end

[coordspzt,nodespzt]=remove_free_spec_nodes(coordspzt,nodespzt);