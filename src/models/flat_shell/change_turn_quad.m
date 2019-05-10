function [nodes,pqs]=change_turn_quad(nodes,pqs)
% change order of node numbers
% only for 4-node elements!!

% additional m3 vectors must be defined!


[fen,NofElNodes]=size(nodes);
n=2; 
m=[1,2,4,3];
m3=[1,4,3,2];
ksi(1)=-1;ksi(2)=1;
wi(1)=1;wi(2)=1;


% Vandermonde matrix
for i=1:n
    for j=1:n
        V(i,j)=ksi(j)^(i-1);
    end
end
Q=inv(V);

% calculate values of the shape function at points ksi
% calculate values of the first derivative of the shape function at points ksi
% calculate values of the second derivative of the shape function at points ksi
V=zeros(n,n);
Vx=zeros(n,n);
Vxx=zeros(n,n);
for k=1:n
    kt=ksi(k);
    for i=1:n
%         for j=1:n
%             V(k,i)=V(k,i)+Q(i,j)*kt^(j-1);
%         end;
        for j=2:n
            Vx(k,i)=Vx(k,i)+Q(i,j)*(j-1)*ksi(k)^(j-2);
        end
         for j=3:n
            Vxx(k,i)=Vxx(k,i)+Q(i,j)*(j-2)*(j-1)*ksi(k)^(j-3);
         end
    end
end
V=zeros(n,n);
for k=1:n
        V(k,k)=1;
end
 for k=2:n-1
     Vx(k,k)=0;
 end

% grid generator 


for ne=1:fen
        %[ne]
    J=zeros(2,2,n,n);
    % global coordinates of the element
    cc=0;
    for i=1:n
        for j=1:n
            cc=cc+1;
            x(cc)=pqs(nodes(ne,m(cc)),1);
            y(cc)=pqs(nodes(ne,m(cc)),2);
        end
    end
    % Jacobians J at Gauss points (2D case)
    inJ=zeros(2,2,n,n);
    for j=1:n % eta
        for i=1:n % ksi       
            k=0;
             for kx=1:n % shape function number
                  k=(j-1)*n+kx;
                  J(1,1,i,j)=J(1,1,i,j)+Vx(i,kx)*x(k); 
                  J(1,2,i,j)=J(1,2,i,j)+Vx(i,kx)*y(k);
             end
             for ky=1:n % shape function number
                  k=(ky-1)*n+i;
                 J(2,2,i,j)=J(2,2,i,j)+Vx(j,ky)*y(k);
                 J(2,1,i,j)=J(2,1,i,j)+Vx(j,ky)*x(k);
             end
        end
    end

    cm=0;det_id=0;
    for ky=1:n % eta
         for kx=1:n % ksi     
             cc=0;cm=cm+1;       
              detJ=J(1,1,kx,ky)*J(2,2,kx,ky)-J(1,2,kx,ky)*J(2,1,kx,ky);
              if(detJ<=0)
                  det_id=1; break;
              end     
         end
         if(det_id==1) break; end
    end
    if(det_id==1) disp('error: negative determinant of Jacobi matrix'); 
    disp(['Changing turn in element no ',num2str(ne)]);
                  nodest=nodes(ne,:);
                  for k=1:n*n
                     nodes(ne,k)=nodest(m3(k));
                  end
    end
end
