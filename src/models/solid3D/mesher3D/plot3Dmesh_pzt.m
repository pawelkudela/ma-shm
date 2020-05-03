function plot3Dmesh_pzt(nodes,coords,NofElNodesx,NofElNodesy,NofElNodesz,pztEl)
% PLOT3DMESH_PZT   One line description of what the function or script performs (H1 line) 
%    optional: more details about the function than in the H1 line 
%    optional: more details about the function than in the H1 line 
%    optional: more details about the function than in the H1 line 
% 
% Syntax: [output1,output2] = plot3Dmesh_pzt(input1,input2,input3) 
% 
% Inputs: 
%    input1 - Description, string, dimensions [m, n], Units: ms 
%    input2 - Description, logical, dimensions [m, n], Units: m 
%    input3 - Description, double, dimensions [m, n], Units: N 
% 
% Outputs: 
%    output1 - Description, integer, dimensions [m, n], Units: - 
%    output2 - Description, double, dimensions [m, n], Units: m/s^2 
% 
% Example: 
%    [output1,output2] = plot3Dmesh_pzt(input1,input2,input3) 
%    [output1,output2] = plot3Dmesh_pzt(input1,input2) 
%    [output1] = plot3Dmesh_pzt(input1,input2,input3) 
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

% corner nodes
mc=[1,NofElNodesx,NofElNodesx*NofElNodesy-NofElNodesx+1,NofElNodesx*NofElNodesy,...
    (NofElNodesz-1)*NofElNodesx*NofElNodesy+1,(NofElNodesz-1)*NofElNodesx*NofElNodesy+NofElNodesx,...
    (NofElNodesz-1)*NofElNodesx*NofElNodesy+NofElNodesx*NofElNodesy-NofElNodesx+1,...
    NofElNodesz*NofElNodesx*NofElNodesy];

host=setdiff([1:size(nodes,1)],pztEl);
% host 
BottomFaceQuadNodes = nodes(host,mc([1,2,4,3]));
TopFaceQuadNodes = nodes(host,mc([5,6,8,7]));
Side1FaceQuadNodes = nodes(host,mc([1,2,6,5]));
Side2FaceQuadNodes = nodes(host,mc([2,4,8,6]));
Side3FaceQuadNodes = nodes(host,mc([4,3,7,8]));
Side4FaceQuadNodes = nodes(host,mc([3,1,5,7]));

colour = 'rgybmcrgybmcrgybmcrgybmcrgybmcrgybmcrgybmcrgybmcrgybmcrgybmcrgybmc';

h=figure; hold all;

patch('faces',Side1FaceQuadNodes,'vertices',[coords(:,1) coords(:,2),coords(:,3)],'facecolor','y');
patch('faces',Side2FaceQuadNodes,'vertices',[coords(:,1) coords(:,2),coords(:,3)],'facecolor','y');
patch('faces',Side3FaceQuadNodes,'vertices',[coords(:,1) coords(:,2),coords(:,3)],'facecolor','y');
patch('faces',Side4FaceQuadNodes,'vertices',[coords(:,1) coords(:,2),coords(:,3)],'facecolor','y');

patch('faces',BottomFaceQuadNodes,'vertices',[coords(:,1) coords(:,2),coords(:,3)],'facecolor','y');
patch('faces',TopFaceQuadNodes,'vertices',[coords(:,1) coords(:,2),coords(:,3)],'facecolor','y');
%plot3(coords(:,1),coords(:,2),coords(:,3),'k.','MarkerSize',0.5);

% pzt
BottomFaceQuadNodes = nodes(pztEl,mc([1,2,4,3]));
TopFaceQuadNodes = nodes(pztEl,mc([5,6,8,7]));
Side1FaceQuadNodes = nodes(pztEl,mc([1,2,6,5]));
Side2FaceQuadNodes = nodes(pztEl,mc([2,4,8,6]));
Side3FaceQuadNodes = nodes(pztEl,mc([4,3,7,8]));
Side4FaceQuadNodes = nodes(pztEl,mc([3,1,5,7]));

patch('faces',Side1FaceQuadNodes,'vertices',[coords(:,1) coords(:,2),coords(:,3)],'facecolor','r');
patch('faces',Side2FaceQuadNodes,'vertices',[coords(:,1) coords(:,2),coords(:,3)],'facecolor','r');
patch('faces',Side3FaceQuadNodes,'vertices',[coords(:,1) coords(:,2),coords(:,3)],'facecolor','r');
patch('faces',Side4FaceQuadNodes,'vertices',[coords(:,1) coords(:,2),coords(:,3)],'facecolor','r');


patch('faces',BottomFaceQuadNodes,'vertices',[coords(:,1) coords(:,2),coords(:,3)],'facecolor','r');
patch('faces',TopFaceQuadNodes,'vertices',[coords(:,1) coords(:,2),coords(:,3)],'facecolor','r');
%plot3(coords(:,1),coords(:,2),coords(:,3),'k.','MarkerSize',0.5);

view(3);axis equal;
set(gcf,'Color','white');
%---------------------- END OF CODE---------------------- 

% ================ [plot3Dmesh_pzt.m] ================  
