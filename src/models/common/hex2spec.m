function [nodeCoordinates3D, elementNodes3D,surfStruct] = hex2spec(mesh_filename,parentFolder,modelfolder, N, plotOpt)
% HEX2SPEC   One line description of what the function or script performs (H1 line) 
%    optional: more details about the function than in the H1 line 
%    optional: more details about the function than in the H1 line 
%    optional: more details about the function than in the H1 line 
% 
% Syntax: [output1,output2] = hex2spec(input1,input2,input3) 
% 
% Inputs: 
%    msh - 3D mesh from GMSH structure 
%    N - element approximation order, integer (N=3,4,5,6,7,8,9)
%    plotOpt - option to plot a mesh, logical
% 
% Outputs: 
%    nodeCoordinates3D - coordinates of the nodes, double, dimensions [x, y, z], Units: [m] 
%    elementNodes3D - spectral elements topology (element nodes), integer, dimensions [nHexElements,n_x^] 
% 
% Example: 
%    [output1,output2] = hex2spec(input1,input2,input3) 
%    [output1,output2] = hex2spec(input1,input2) 
%    [output1] = hex2spec(input1,input2,input3) 
% 
% Other m-files required: gll 
% Subfunctions: none 
% MAT-files required: none 
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2 
% 

% Author: Piotr Fiborek, M.Sc., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pfiborek@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 

%---------------------- BEGIN CODE---------------------- 

% corner nodes 
if nargin < 4, disp('Not enough input arguments.'); help(mfilename); return; end
if nargin == 4, plotOpt = false; end


figfilename = fullfile(parentFolder,'src','models',modelfolder,'input','mesh',[mesh_filename,'.png']);
gmsh_path = fullfile(parentFolder,'bin','external','gmsh','gmsh ');
mesh_geometry_path = fullfile(parentFolder,'src','models',modelfolder,'input','mesh','geo',filesep);
mesh_output_path = fullfile(parentFolder,'src','models',modelfolder,'input','mesh','gmsh_out',filesep);
spec_mesh_output_path = fullfile(parentFolder,'src','models',modelfolder,'input','mesh',filesep);
gmsh_options = ' -3 -format auto -v 1 -o '; % non-verbose
%gmsh_options = ' -2 -format auto -o '; % verbose
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gmsh_command = [gmsh_path, mesh_geometry_path, mesh_filename,'.geo', gmsh_options, mesh_output_path, mesh_filename,'.m'];

% run external gmsh.exe - make mesh file from geometry file
status = system(gmsh_command);
if(status)
    disp('Mesh generation in gmsh failed');
end

% load mesh into matlab
run([mesh_output_path, mesh_filename,'.m']);

n_x = N + 1; % number of nodes in 1D element
ksi = gllnaw(n_x);
n_int = n_x-2;
nodeCoordinates = round(msh.POS*1e8)*1e-8;
elementNodes = msh.HEXAS(:,1:8);
subSect = msh.HEXAS(:,9);
nodeID = 1 : length(nodeCoordinates);
emptyNodes = nodeID(~ismember(nodeID',unique(elementNodes)));
while ~isempty(emptyNodes)
   elementNodes(elementNodes>emptyNodes(1)) = elementNodes(elementNodes>emptyNodes(1) ) - 1;
   nodeCoordinates(emptyNodes(1),:) = [];
   nodeID = 1 : length(nodeCoordinates);
   emptyNodes = nodeID(~ismember(nodeID,unique(elementNodes)));
end

%noNodes = size(nodeCoordinates,1);
[noElements,~] = size(elementNodes);
%nodes = struct('x', num2cell(nodeCoordinates(:,1)),'y', num2cell(nodeCoordinates(:,2)),'z', num2cell(nodeCoordinates(:,3)));

%nodeID = 1 : noNodes;
elementID = 1 : noElements;

% nodesInElement = repmat(reshape(elementNodes',1,[]),noNodes,1)==repmat(nodeID',1,noElements*NodsinElement);
% elementsPerNode = sum(nodesInElement,2);

% for k = unique(elementsPerNode)'
%     clc
%     disp('nodes')
%     disp([num2str(k), ' from ', num2str(max(unique(elementsPerNode)))])
%     row = nodeID(elementsPerNode==k);
%     cornerElement = repmat(repmat((1:8),1,noElements),length(row),1)';
%     element = repmat(reshape(repmat(elementID,8,1),1,[]),length(row),1)';
%     
%     corner = num2cell(reshape(cornerElement(nodesInElement(elementsPerNode==k,:)'),k,[])',2);
%     element = num2cell(reshape(element(nodesInElement(elementsPerNode==k,:)'),k,[])',2);
%     [nodes(row).element] = element{:};
%     [nodes(row).corner] = corner{:};
%     
% end
%% 
% edge internal nodes ksi(2:n-1)
edgeNodes = [1 2 2 3 3 4 4 1 5 6 6 7 7 8 8 5 1 5 2 6 3 7 4 8];
elementEdgeNodes = reshape(elementNodes(:,edgeNodes)',2,[])';
elementEdgeNodeSort = sort(elementEdgeNodes,2);

[~, m, n] = unique(elementEdgeNodeSort,'rows','first');
edge = elementEdgeNodes(m,:);
edges = struct('nodes',num2cell(edge,2));

% edgeZeta = false(length(edges),1);
edgeNodes = cell(length(edges),1);
[edgeNodes{:}] = deal(edges.nodes);
edgeNodes = cell2mat(edgeNodes);

ksi_int=repmat(ksi(2:end-1),length(edges),1);

ksi_X = 0.5*(repmat(nodeCoordinates(edgeNodes(:,2),1) - nodeCoordinates(edgeNodes(:,1),1),1,n_int).*ksi_int + ...
    repmat(nodeCoordinates(edgeNodes(:,2),1) + nodeCoordinates(edgeNodes(:,1),1),1,n_int));
ksi_Y = 0.5*(repmat(nodeCoordinates(edgeNodes(:,2),2) - nodeCoordinates(edgeNodes(:,1),2),1,n_int).*ksi_int + ...
    repmat(nodeCoordinates(edgeNodes(:,2),2) + nodeCoordinates(edgeNodes(:,1),2),1,n_int));
ksi_Z = 0.5*(repmat(nodeCoordinates(edgeNodes(:,2),3) - nodeCoordinates(edgeNodes(:,1),3),1,n_int).*ksi_int + ...
    repmat(nodeCoordinates(edgeNodes(:,2),3) + nodeCoordinates(edgeNodes(:,1),3),1,n_int));

X_edge = reshape(ksi_X',[],1); Y_edge = reshape(ksi_Y',[],1); Z_edge = reshape(ksi_Z',[],1);

nC_edge = [nodeCoordinates; X_edge, Y_edge,Z_edge];
nodesEdgeInt = reshape((1:numel(ksi_X))'+ size(nodeCoordinates,1),n_int,[])';

edgeID = num2cell(1 : length(edges))';
edgeInElements = cellfun(@(x) ceil(find(x==n)/12),edgeID, 'UniformOutput',false);
[edges.type] = deal([]);
[edges.condition] = deal([]);
[edges.element] = deal(edgeInElements{:});

edgeNo = {[1 2], [2 3], [4 3], [1 4], [5 6], [6 7], [8 7], [5 8], [1 5], [2 6], [3 7], [4 8]};
edgeComb = {[1 2], [2 1];
            1:n_int, n_int:-1:1};
elementNodesEdgeInt = cell(1,12);
[elementNodesEdgeInt{:}] = deal(zeros(noElements,n_int));

for iEdges = 1 : length(edges)
    if ~mod(iEdges,1000)
        clc
        disp(length(edges)-iEdges)
    end
    for iK = 1 : length(edges(iEdges).element)
        element = edges(iEdges).element(iK);
        edgeNodes = edges(iEdges).nodes;
        col = [find(elementNodes(element,:)==edgeNodes(1)),...
            find(elementNodes(element,:)==edgeNodes(2))];
        for iSide = 1 : length(edgeNo)
            for iComb = 1 : size(edgeComb,2)
                if all(col==edgeNo{iSide}(edgeComb{1,iComb}))
                    elementNodesEdgeInt{iSide}(element,:) = ...
                        nodesEdgeInt(iEdges,edgeComb{2,iComb});
                end
            end
        end
    end
end

%%
% surface internal nodes [ksi(2:n-1),eta(2:n-1)]
surfNodes = [1 2 3 4 5 6 7 8 1 2 6 5 2 3 7 6 4 3 7 8 1 4 8 5];
elementSurfNodes = reshape(elementNodes(:,surfNodes)',4,[])';
elementSurfNodeSort = sort(elementSurfNodes,2);

[~, m, n,] = unique(elementSurfNodeSort,'rows','first');
surfS = elementSurfNodes(m,:);

surfStruct = struct('nodes',num2cell(surfS,2));
surfNodes = cell(length(surfStruct),1);
[surfNodes{:}] = deal(surfStruct.nodes);
surfNodes = cell2mat(surfNodes);


ksi_int = repmat(repmat(ksi(2:end-1)',n_int,1),length(surfStruct),1);

x1 = reshape(repmat(nodeCoordinates(surfNodes(:,1),1)',n_int^2,1),[],1);
x2 = reshape(repmat(nodeCoordinates(surfNodes(:,2),1)',n_int^2,1),[],1);
x3 = reshape(repmat(nodeCoordinates(surfNodes(:,3),1)',n_int^2,1),[],1);
x4 = reshape(repmat(nodeCoordinates(surfNodes(:,4),1)',n_int^2,1),[],1);

y1 = reshape(repmat(nodeCoordinates(surfNodes(:,1),2)',n_int^2,1),[],1);
y2 = reshape(repmat(nodeCoordinates(surfNodes(:,2),2)',n_int^2,1),[],1);
y3 = reshape(repmat(nodeCoordinates(surfNodes(:,3),2)',n_int^2,1),[],1);
y4 = reshape(repmat(nodeCoordinates(surfNodes(:,4),2)',n_int^2,1),[],1);

z1 = reshape(repmat(nodeCoordinates(surfNodes(:,1),3)',n_int^2,1),[],1);
z2 = reshape(repmat(nodeCoordinates(surfNodes(:,2),3)',n_int^2,1),[],1);
z3 = reshape(repmat(nodeCoordinates(surfNodes(:,3),3)',n_int^2,1),[],1);
z4 = reshape(repmat(nodeCoordinates(surfNodes(:,4),3)',n_int^2,1),[],1);


X11 = 0.5*((x2-x1).*ksi_int + (x2+x1));
Y11 = 0.5*((y2-y1).*ksi_int + (y2+y1));
Z11 = 0.5*((z2-z1).*ksi_int + (z2+z1));

X12 = 0.5*((x3-x4).*ksi_int + (x3+x4));
Y12 = 0.5*((y3-y4).*ksi_int + (y3+y4));
Z12 = 0.5*((z3-z4).*ksi_int + (z3+z4));

eta_int = repmat(reshape(repmat(ksi(2:end-1),n_int,1),[],1),length(surfStruct),1);

X21 = 0.5*((x4-x1).*eta_int + (x4+x1));
Y21 = 0.5*((y4-y1).*eta_int + (y4+y1));
Z21 = 0.5*((z4-z1).*eta_int + (z4+z1));

X22 = 0.5*((x3-x2).*eta_int + (x3+x2));
Y22 = 0.5*((y3-y2).*eta_int + (y3+y2));
Z22 = 0.5*((z3-z2).*eta_int + (z3+z2));


ABC_0 = cross([x2-x1,y2-y1,z2-z1],[x3-x1,y3-y1,z3-z1]);
A0 = ABC_0(:,1);
B0 = ABC_0(:,2);
C0 = ABC_0(:,3);
D0 = A0.*x1 + B0.*y1 + C0.*z1;

ABC_1 = cross([A0 B0 C0],[X12 - X11,Y12 - Y11,Z12 - Z11]);
A1 = ABC_1(:,1);
B1 = ABC_1(:,2);
C1 = ABC_1(:,3);
D1 = A1.*X11 + B1.*Y11 + C1.*Z11;

ABC_2 = cross([A0 B0 C0],[X22 - X21,Y22 - Y21,Z22 - Z21]);
A2 = ABC_2(:,1);
B2 = ABC_2(:,2);
C2 = ABC_2(:,3);
D2 = A2.*X21 + B2.*Y21 + C2.*Z21;

W = A0.*B1.*C2 + A1.*B2.*C0 + A2.*B0.*C1 - A2.*B1.*C0 - A0.*B2.*C1 - A1.*B0.*C2;
Wx = D0.*B1.*C2 + D1.*B2.*C0 + D2.*B0.*C1 - D2.*B1.*C0 - D0.*B2.*C1 - D1.*B0.*C2;
Wy = A0.*D1.*C2 + A1.*D2.*C0 + A2.*D0.*C1 - A2.*D1.*C0 - A0.*D2.*C1 - A1.*D0.*C2;
Wz = A0.*B1.*D2 + A1.*B2.*D0 + A2.*B0.*D1 - A2.*B1.*D0 - A0.*B2.*D1 - A1.*B0.*D2;

X_surf = Wx./W; Y_surf = Wy./W; Z_surf = Wz./W;
nC_surf = [nC_edge; X_surf, Y_surf, Z_surf];
nodesSurfInt = reshape((1:numel(X_surf))'+ size(nC_edge,1),n_int^2,[])';


surfID = num2cell(1 : length(surfStruct))';
surfInElements = cellfun(@(x) ceil(find(x==n)/6),surfID, 'UniformOutput',false);
[surfStruct.element] = deal(surfInElements{:});

[surfStruct(:).type] = deal([]);
[surfStruct(:).condition] = deal([]);
[surfStruct(:).side] = deal([]);
[surfStruct(:).boundary] = deal([]);
surfH = cell(6,1);
surfH{1} = 1:n_x^2;
surfH{2} = n_x^2*(n_x-1)+1 : n_x^3;
surfH{3} = reshape((repmat(1:n_x,n_x,1) + repmat(n_x^2*(0:n_x-1)',1,n_x))',1,[]);
surfH{4} = reshape((repmat(n_x:n_x:n_x^2,n_x,1) + repmat(n_x^2*(0:n_x-1)',1,n_x))',1,[]);
surfH{5} = reshape((repmat(1+n_x*(n_x-1):n_x^2,n_x,1) + repmat(n_x^2*(0:n_x-1)',1,n_x))',1,[]);
surfH{6} = reshape((repmat(1:n_x:n_x^2,n_x,1) + repmat(n_x^2*(0:n_x-1)',1,n_x))',1,[]);


elementNodesSurfInt = cell(1,6);
[elementNodesSurfInt{:}] = deal(zeros(noElements,(n_int).^2));
surfNo = {[1 2 3 4], [5 6 7 8], [1 2 6 5], [2 3 7 6], [4 3 7 8], [1 4 8 5]};
surfComb = {1:4, [2 1 4 3], 4:-1:1, [3 4 1 2]
    1:n_int^2, reshape(flipud(reshape(1:n_int^2,n_int,n_int)),1,n_int^2), n_int^2:-1:1, ...
    reshape(flipud(reshape(n_int^2:-1:1,n_int,n_int)),1,n_int^2)};


for iSurfes= 1 : length(surfStruct)
    if ~mod(iSurfes,1000)
        clc
        disp(length(surfStruct)-iSurfes)
    end
    for iK = 1 : length(surfStruct(iSurfes).element)
        
        element = surfStruct(iSurfes).element(iK);
        surfNodes = surfStruct(iSurfes).nodes;
        col = [find(elementNodes(element,:)==surfNodes(1)),...
            find(elementNodes(element,:)==surfNodes(2)),...
            find(elementNodes(element,:)==surfNodes(3)),...
            find(elementNodes(element,:)==surfNodes(4))];
        for iSide = 1 : length(surfNo)
            for iComb = 1 : size(surfComb,2)
                if all(col==surfNo{iSide}(surfComb{1,iComb}))
                    elementNodesSurfInt{iSide}(element,:) = ...
                        nodesSurfInt(iSurfes,surfComb{2,iComb});
                    if length(surfStruct(iSurfes).element)==1
                        surfStruct(iSurfes).boundary = surfH{iSide};
                    end
                end
            end
        end
    end
end
%%
% volume internal nodes [ksi(2:n-1),eta(2:n-1),zeta(2:n-1)]



x11 = reshape(nC_surf(elementNodesSurfInt{1}(:,1:n_int:n_int^2)',1),n_int,[])';
y11 = reshape(nC_surf(elementNodesSurfInt{1}(:,1:n_int:n_int^2)',2),n_int,[])';
z11 = reshape(nC_surf(elementNodesSurfInt{1}(:,1:n_int:n_int^2)',3),n_int,[])';

x11 = reshape(repmat(reshape(repmat(x11',n_int,1),1,[]),n_int,1),[],1);
y11 = reshape(repmat(reshape(repmat(y11',n_int,1),1,[]),n_int,1),[],1);
z11 = reshape(repmat(reshape(repmat(z11',n_int,1),1,[]),n_int,1),[],1);

x12 = reshape(nC_surf(elementNodesSurfInt{1}(:,n_int:n_int:n_int^2)',1),n_int,[])';
y12 = reshape(nC_surf(elementNodesSurfInt{1}(:,n_int:n_int:n_int^2)',2),n_int,[])';
z12 = reshape(nC_surf(elementNodesSurfInt{1}(:,n_int:n_int:n_int^2)',3),n_int,[])';

x12 = reshape(repmat(reshape(repmat(x12',n_int,1),1,[]),n_int,1),[],1);
y12 = reshape(repmat(reshape(repmat(y12',n_int,1),1,[]),n_int,1),[],1);
z12 = reshape(repmat(reshape(repmat(z12',n_int,1),1,[]),n_int,1),[],1);

x13 = reshape(nC_surf(elementNodesSurfInt{2}(:,1:n_int:n_int^2)',1),n_int,[])';
y13 = reshape(nC_surf(elementNodesSurfInt{2}(:,1:n_int:n_int^2)',2),n_int,[])';
z13 = reshape(nC_surf(elementNodesSurfInt{2}(:,1:n_int:n_int^2)',3),n_int,[])';

x13 = reshape(repmat(reshape(repmat(x13',n_int,1),1,[]),n_int,1),[],1);
y13 = reshape(repmat(reshape(repmat(y13',n_int,1),1,[]),n_int,1),[],1);
z13 = reshape(repmat(reshape(repmat(z13',n_int,1),1,[]),n_int,1),[],1);


ABC_12 = cross([x12-x11,y12-y11,z12-z11],[x13-x11,y13-y11,z13-z11]);
A12 = ABC_12(:,1);
B12 = ABC_12(:,2);
C12 = ABC_12(:,3);
D12 = A12.*x11 + B12.*y11 + C12.*z11;

x21 = reshape(nC_surf(elementNodesSurfInt{3}(:,1:n_int)',1),n_int,[])';
y21 = reshape(nC_surf(elementNodesSurfInt{3}(:,1:n_int)',2),n_int,[])';
z21 = reshape(nC_surf(elementNodesSurfInt{3}(:,1:n_int)',3),n_int,[])';

x21 = reshape(repmat(x21,1,n_int^2)',[],1);
y21 = reshape(repmat(y21,1,n_int^2)',[],1);
z21 = reshape(repmat(z21,1,n_int^2)',[],1);

x22 = reshape(nC_surf(elementNodesSurfInt{3}(:,n_int^2-n_int+1:n_int^2)',1),n_int,[])';
y22 = reshape(nC_surf(elementNodesSurfInt{3}(:,n_int^2-n_int+1:n_int^2)',2),n_int,[])';
z22 = reshape(nC_surf(elementNodesSurfInt{3}(:,n_int^2-n_int+1:n_int^2)',3),n_int,[])';

x22 = reshape(repmat(x22,1,n_int^2)',[],1);
y22 = reshape(repmat(y22,1,n_int^2)',[],1);
z22 = reshape(repmat(z22,1,n_int^2)',[],1);

x23 = reshape(nC_surf(elementNodesSurfInt{5}(:,1:n_int)',1),n_int,[])';
y23 = reshape(nC_surf(elementNodesSurfInt{5}(:,1:n_int)',2),n_int,[])';
z23 = reshape(nC_surf(elementNodesSurfInt{5}(:,1:n_int)',3),n_int,[])';

x23 = reshape(repmat(x23,1,n_int^2)',[],1);
y23 = reshape(repmat(y23,1,n_int^2)',[],1);
z23 = reshape(repmat(z23,1,n_int^2)',[],1);


ABC_35 = cross([x22-x21,y22-y21,z22-z21],[x23-x21,y23-y21,z23-z21]);
A35 = ABC_35(:,1);
B35 = ABC_35(:,2);
C35 = ABC_35(:,3);
D35 = A35.*x21 + B35.*y21 + C35.*z21;


x31 = reshape(nC_surf(elementNodesSurfInt{6}(:,1:n_int:n_int^2)',1),n_int,[])';
y31 = reshape(nC_surf(elementNodesSurfInt{6}(:,1:n_int:n_int^2)',2),n_int,[])';
z31 = reshape(nC_surf(elementNodesSurfInt{6}(:,1:n_int:n_int^2)',3),n_int,[])';

x31 = reshape(repmat(reshape(x31',1,[]),n_int^2,1),[],1);
y31 = reshape(repmat(reshape(y31',1,[]),n_int^2,1),[],1);
z31 = reshape(repmat(reshape(z31',1,[]),n_int^2,1),[],1);

x32 = reshape(nC_surf(elementNodesSurfInt{6}(:,n_int:n_int:n_int^2)',1),n_int,[])';
y32 = reshape(nC_surf(elementNodesSurfInt{6}(:,n_int:n_int:n_int^2)',2),n_int,[])';
z32 = reshape(nC_surf(elementNodesSurfInt{6}(:,n_int:n_int:n_int^2)',3),n_int,[])';

x32 = reshape(repmat(reshape(x32',1,[]),n_int^2,1),[],1);
y32 = reshape(repmat(reshape(y32',1,[]),n_int^2,1),[],1);
z32 = reshape(repmat(reshape(z32',1,[]),n_int^2,1),[],1);

x33 = reshape(nC_surf(elementNodesSurfInt{4}(:,1:n_int:n_int^2)',1),n_int,[])';
y33 = reshape(nC_surf(elementNodesSurfInt{4}(:,1:n_int:n_int^2)',2),n_int,[])';
z33 = reshape(nC_surf(elementNodesSurfInt{4}(:,1:n_int:n_int^2)',3),n_int,[])';

x33 = reshape(repmat(reshape(x33',1,[]),n_int^2,1),[],1);
y33 = reshape(repmat(reshape(y33',1,[]),n_int^2,1),[],1);
z33 = reshape(repmat(reshape(z33',1,[]),n_int^2,1),[],1);


ABC_64 = cross([x32-x31,y32-y31,z32-z31],[x33-x31,y33-y31,z33-z31]);
A64 = ABC_64(:,1);
B64 = ABC_64(:,2);
C64 = ABC_64(:,3);
D64 = A64.*x31 + B64.*y31 + C64.*z31;


W = A12.*B35.*C64 + A35.*B64.*C12 + A64.*B12.*C35 - A64.*B35.*C12 - A12.*B64.*C35 - A35.*B12.*C64;
Wx = D12.*B35.*C64 + D35.*B64.*C12 + D64.*B12.*C35 - D64.*B35.*C12 - D12.*B64.*C35 - D35.*B12.*C64;
Wy = A12.*D35.*C64 + A35.*D64.*C12 + A64.*D12.*C35 - A64.*D35.*C12 - A12.*D64.*C35 - A35.*D12.*C64;
Wz = A12.*B35.*D64 + A35.*B64.*D12 + A64.*B12.*D35 - A64.*B35.*D12 - A12.*B64.*D35 - A35.*B12.*D64;

X_vol = Wx./W; Y_vol = Wy./W; Z_vol = Wz./W;
nodeCoordinates3D = [nC_surf; X_vol, Y_vol, Z_vol];
nodesVolInt = reshape((1:numel(X_vol))' + size(nC_surf,1),n_int^3,[])';
elementNodesVolfInt = mat2cell(nodesVolInt,noElements,repmat(n_int^2,1,n_int));

%%
col_bot = cell(1,n_int);
% first layer
for i = 1 : n_int
    col_bot{i} = [elementNodesEdgeInt{4}(:,i) elementNodesSurfInt{1}(:,1+n_int*(i-1):n_int*i)...
        elementNodesEdgeInt{2}(:,i)];
end

col_bot = [elementNodes(:,1), elementNodesEdgeInt{1}, elementNodes(:,2), cell2mat(col_bot),...
    elementNodes(:,4), elementNodesEdgeInt{3}, elementNodes(:,3)];
% int layers

col_int = cell(n_int + 2, n_int);
for i = 1 : n_int
    for j = 1 : n_int + 2
        if j == 1
            col_int{j,i} = [elementNodesEdgeInt{9}(:,i), elementNodesSurfInt{3}(:,1+n_int*(i-1):n_int*i),...
                elementNodesEdgeInt{10}(:,i)];
        elseif j == n_int + 2
            col_int{j,i} = [elementNodesEdgeInt{12}(:,i), elementNodesSurfInt{5}(:,1+n_int*(i-1):n_int*i),...
                elementNodesEdgeInt{11}(:,i)];
        else
            col_int{j,i} = [elementNodesSurfInt{6}(:,(j-1)+n_int*(i-1)),...
                elementNodesVolfInt{i}(:,1+n_int*(j-2):n_int*(j-1)),...
                elementNodesSurfInt{4}(:,(j-1)+n_int*(i-1))];
        end
    end
end
col_int = cell2mat(reshape(col_int,1,[]));

col_top = cell(1,n_int);
% first layer
for i = 1 : n_int
    col_top{i} = [elementNodesEdgeInt{8}(:,i) elementNodesSurfInt{2}(:,1+n_int*(i-1):n_int*i)...
        elementNodesEdgeInt{6}(:,i)];
end

col_top = [elementNodes(:,5), elementNodesEdgeInt{5}, elementNodes(:,6), cell2mat(col_top),...
    elementNodes(:,8), elementNodesEdgeInt{7}, elementNodes(:,7)];

elementNodes3D = [col_bot col_int col_top];
%%
if plotOpt
    surfOnBoundary = cell(length(surfStruct),1);
    [surfOnBoundary{:}] = deal(surfStruct.boundary);
    surfOnBoundary = cellfun(@isempty, surfOnBoundary);
    surfOnBoundary = ~surfOnBoundary;
    boundaryElement = cell(size(find(surfOnBoundary)));
    boundaryNodes = cell(size(find(surfOnBoundary)));
    [boundaryElement{:}] = deal(surfStruct(surfOnBoundary).element);
    boundaryElement = cell2mat(boundaryElement);
    [boundaryNodes{:}] = deal(surfStruct(surfOnBoundary).boundary);
    boundaryNodes = cell2mat(boundaryNodes);
    corners = [1, n_x, n_x^2, n_x^2-n_x+1];
    col = {'r','r','r','r','r','b','b','b','b','b','g','g'};
    subSectBoundary =  subSect(boundaryElement);
    usubSect = unique(subSectBoundary);
    for iSection = 1 : length(usubSect)
        faces = elementNodes3D(sub2ind(size(elementNodes3D),...
            repmat(boundaryElement(subSectBoundary==(iSection)),1,4),...
            boundaryNodes(subSectBoundary==(iSection),corners)));
        patch('Faces',faces,'Vertices',nodeCoordinates3D,'FaceColor',col{iSection});
    end
    view(3)
    axis equal
    set(gca,'FontName','Times');
    set(gcf,'Color','w');
    set(gca,'Fontsize',10);
    axis off
    fig_width = 18; fig_height = 18; 
    fig = gcf;
    set(fig, 'Units','centimeters', 'Position',[10 5 fig_width fig_height]); 

    % remove unnecessary white space
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
    pause(3)
    fig = struct;
    fig.PaperPositionMode   = 'auto';
    print(figfilename,'-dpng', '-r300');
    pause(3)
    close all;
    pause(1)
end



% delete gmsh out m file
delete([mesh_output_path, mesh_filename,'.m']);
save([spec_mesh_output_path,mesh_filename,'.mat'],'elementNodes3D','nodeCoordinates3D');
%---------------------- END OF CODE---------------------- 

% ================ [hex2spec.m] ================  
