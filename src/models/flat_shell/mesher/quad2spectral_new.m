function [spec_element_nodes,spec_coords,boundary_nodes] = ...
    quad2spectral_new(elementNodes,nodeCoordinates,N_x,N_y)
% convert quad nodes mesh into spectral element mesh 
%    only for linear quad elements 
% 
% Syntax: [spec_element_nodes,spec_coords,boundary_nodes] = ...
%                                   quad2spectral(elementNodes_fem,quad_coords,N) 
% 
% Inputs: 
%    elementNodes_fem - Quad nodes topology (element nodes), integer, dimensions [nQuadElements, 4]
%    nodeCoordinates_fem - coordinates of quad element nodes, double, dimensions [nQuadNodes, 3], Units: m 
%    N - element approximation order, integer (N=3,4,5,6,7,8,9)
% 
% Outputs: 
%    spec_element_nodes - spectral elements topology (element nodes), integer, dimensions [nQuadElements,(N+1)^2]
%    spec_coords - coordinates of spectral element nodes, integer, dimensions [nSpecNodes, 3], Units: m 
%    boundary_nodes - all nodes lying on the boundary of the structure
%  13    14   15   16
%   O----O----O----O   
%   |    |    |    |
%  9O--10O----O11--O12
%   |    |    |    | 
%  5O---6O----O7---O8
%   |    |    |    |
%   O----O----O----O
%   1    2    3    4
%
% Example: 
%    [spec_element_nodes,spec_coords] = quad2spec(q,quad_coords,N)
% 
% Other m-files required: gllnaw.m 
% Subfunctions: edgegeometry 
% MAT-files required: none 
% See also:
% 

% Author: Piotr Fiborek, M.Sc., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pfiborek@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 

%---------------------- BEGIN CODE---------------------- 
    n_x = N_x+1;
    n_y = N_y+1;
    if n_x<3 || n_y<3
        spec_element_nodes = elementNodes(:, [1 2 4 3]);
        spec_coords = nodeCoordinates;
        nodesOnedgeElement = reshape(1:n_x*n_y,n_x,n_y);
        nodesOnedgeElement = unique([nodesOnedgeElement(:,1)' nodesOnedgeElement(:,end)' ...
        nodesOnedgeElement(1,:) nodesOnedgeElement(end,:)]);
        counts = histcounts(spec_element_nodes(:,nodesOnedgeElement), nodesInStr);
        boundary_nodes = nodesInStr(counts==1);
        return
    end
    format long
    ksi = gll(n_x);
    eta = gll(n_y);
    [ksi,eta] = meshgrid(ksi,eta);
    ksi = reshape(ksi',1,[]);
    eta = reshape(eta',1,[]);
    cornersX = reshape(nodeCoordinates(elementNodes',1),4,[])';
    cornersY = reshape(nodeCoordinates(elementNodes',2),4,[])';
    Eps = 1e-10;
    psi = [(1-ksi).*(1-eta);(ksi+1).*(1-eta);(ksi+1).*(1+eta);(1-ksi).*(eta+1)];
    x = 1/4*cornersX*psi; y = 1/4*cornersY*psi;
    xy = [reshape(x',[],1),reshape(y',[],1)];
    [~,ia,ic] = unique(round(xy/Eps)*Eps,'rows');
    spec_coords = xy(ia,:);
    spec_coords(:,3) = 0;
    spec_element_nodes = reshape(ic,n_x*n_y,[])';
    nodesInStr = unique(spec_element_nodes);
    nodesInElemen = reshape(1:n_x^2,n_x,n_x);
    nodesOnedgeElement = unique([nodesInElemen(2:end-1,1)' nodesInElemen(2:end-1,end)' nodesInElemen(1,2:end-1) nodesInElemen(end,2:end-1)]);
    nodesInCorners = unique([nodesInElemen(1,1) nodesInElemen(1,end)' nodesInElemen(end,1) nodesInElemen(end,end)]);
    countsE = histcounts(spec_element_nodes(:,nodesOnedgeElement), nodesInStr);
    countsC = histcounts(spec_element_nodes(:,nodesInCorners), nodesInStr);
    boundary_nodes = [nodesInStr(countsE==1); nodesInStr(countsC<=2&countsC~=0)];
end
%---------------------- END OF CODE--------------------------

% ================ [quad2spectral_Fiborek.m] ================