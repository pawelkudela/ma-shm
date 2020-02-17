function [shapeFunction_P,ownerElement,nodes_order] = ...
    spectral2meshgrid(n_x,n_y,nodeCoordinates_main,elementNodes_main,nodeCoordinates_interface)
Eps = 1e-7;
x_p = nodeCoordinates_interface(:,1);
y_p = nodeCoordinates_interface(:,2);
ownerElement = findOwnerElement(x_p,y_p,nodeCoordinates_main,elementNodes_main,n_x,n_y);


[ksi,~]=gll(n_x);
[eta,~]=gll(n_y);

[Ksi, Eta] = meshgrid(ksi,eta);
Ksi = repmat(reshape(Ksi',1,[]),length(x_p),1);
Eta = repmat(reshape(Eta',1,[]),length(x_p),1);

elementNodes_owner = elementNodes_main(ownerElement,:);
x_0 = nodeCoordinates_main(elementNodes_owner,1);
x_0 = reshape(x_0,[],n_x*n_y);
y_0 = nodeCoordinates_main(elementNodes_owner,2);
y_0 = reshape(y_0,[],n_x*n_y);

[~,point_no] = min(sqrt(bsxfun(@minus, x_0,x_p).^2+bsxfun(@minus, y_0,y_p).^2),[],2);
x_0 = x_0(sub2ind(size(x_0),(1:length(point_no))',point_no));
y_0 = y_0(sub2ind(size(y_0),(1:length(point_no))',point_no));

ksi_0 = Ksi(sub2ind(size(Ksi),ones(length(point_no),1),point_no));
eta_0 = Eta(sub2ind(size(Eta),ones(length(point_no),1),point_no));

n1 = (ownerElement-1)*n_x*n_y+1;
[~,Jacob_P11inv,Jacob_P12inv,Jacob_P21inv,Jacob_P22inv,~,~,~,~,~,~,...
    ~,~,~,~] = Jacob_NbN(5,n_x,n_y,1,elementNodes_main,nodeCoordinates_main,[]);
ksi_p = ksi_0 + Jacob_P11inv(n1).*(x_p - x_0) + Jacob_P21inv(n1).*(y_p - y_0);
eta_p = eta_0 + Jacob_P12inv(n1).*(x_p - x_0) + Jacob_P22inv(n1).*(y_p - y_0);
[ksi_i,~] = gll(n_x);
[eta_i,~] = gll(n_y);
[Q_ksi] = Vandermonde_v2(ksi_i,n_x);
[Q_eta] = Vandermonde_v2(eta_i,n_y);

[shapeFunction_P,~,~] = shapeFunction_2D(Q_ksi,Q_eta,ksi_p,eta_p,'s');
% coordinates of element nodes
nodes_order = reshape(elementNodes_owner',[],1);

x_e = nodeCoordinates_main(nodes_order,1);
y_e = nodeCoordinates_main(nodes_order,2);
    
x_k = shapeFunction_P*x_e;% interpolated values
y_k = shapeFunction_P*y_e;% interpolated values

[A,~] = max(sqrt((x_p-x_k).^2+(y_p-y_k).^2));
k=0;
while A > Eps
    k = k +1;
    disp(k)
    ksi_0 = ksi_p;
    eta_0 = eta_p;
    x_0 = x_k;
    y_0 = y_k;
    % shape function derivatives
    [~,naturalDerivativesX_P,naturalDerivativesY_P] = ...
    shapeFunction_2D(Q_ksi,Q_eta,ksi_0,eta_0,'d');

    % jacobians
    J11 = naturalDerivativesX_P*x_e;
    J12 = naturalDerivativesX_P*y_e;
    J21 = naturalDerivativesY_P*x_e;
    J22 = naturalDerivativesY_P*y_e;
    [Jacob_P11inv,Jacob_P12inv,Jacob_P21inv,Jacob_P22inv] =...
        Jacob_inv(J11,J12,J21,J22);
    % next iteration
    ksi_p = ksi_0 + Jacob_P11inv.*(x_p-x_0) + Jacob_P21inv.*(y_p-y_0);
    eta_p = eta_0 + Jacob_P12inv.*(x_p-x_0) + Jacob_P22inv.*(y_p-y_0);
    
    [shapeFunction_P,~,~] = shapeFunction_2D(Q_ksi,Q_eta,ksi_p,eta_p,'s');
    x_k = shapeFunction_P*x_e;% interpolated values
    y_k = shapeFunction_P*y_e;% interpolated values
    
    [A,~] = max(sqrt((x_p-x_k).^2+(y_p-y_k).^2));
end
shapeFunction_P(abs(shapeFunction_P)<1e-8) = 0;
shapeFunction_P = round(shapeFunction_P*1e8)*1e-8;

