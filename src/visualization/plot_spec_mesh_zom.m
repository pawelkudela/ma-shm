meshname = 'delam_Jochen_signals_D5_a_5mm_b_5mm_angle_0added_mass';
paper_fig_folder = 'E:\work\projects\nawa-bekker\ma-shm\reports\journal_papers\Elsevier\figs\';
meshfile = 'E:\work\projects\nawa-bekker\ma-shm\src\models\flat_shell\mesh\delam_Jochen_signals_D5_a_5mm_b_5mm_angle_0added_mass';
figfilename1 = [paper_fig_folder,meshname];
figfilename2 = [paper_fig_folder,meshname,'_spec_zoom'];
load(meshfile);
X=msh.POS(:,1);
Y=msh.POS(:,2);
nRegions = max(msh.QUADS(:,5));
nbElements = size(msh.QUADS,1);
%colour = 'rgybmcrgybmcrgybmcrgybmcrgybmcrgybmcrgybmcrgybmcrgybmcrgybmcrgybmc';
colour = 'rrrrrrrrrrrrgybmcrgybmcrgybmcrgybmcrgybmcrgybmcrgybmcrgybmcrgybmcrgybmcrgybmc';
h=figure; hold all;
c=1;
for k=1:nRegions
    c1=c;
    while(c <= nbElements && msh.QUADS(c,5)==k)
        c=c+1;
    end 
    nodes = msh.QUADS(c1:c-1,1:4);
    patch('faces',nodes,'vertices',[X Y],'facecolor',colour(k));
end
% add pzt numbering
pzt_x1=[0.45,0.37,0.29,0.21,0.13,0.05]+0.007;
pzt_x2=[0.45,0.37,0.29,0.21,0.13,0.05]+0.007;
pzt_y1 = 0.47;
pzt_y2 = 0.03;
for k=1:6
    text(pzt_x1(k),pzt_y1,num2str(k),'Color','red','FontSize',12);
end
for k=1:6
    text(pzt_x2(k),pzt_y2,num2str(6+k),'Color','red','FontSize',12);
end
% size 12cm by 8cm (1-column text)
fig_width = 12; fig_height = 12; 

axis equal;
axis([msh.MIN(1) msh.MAX(1) msh.MIN(2)  msh.MAX(2)]);
set(gca,'FontName','Times');
set(gcf,'Color','w');
set(gca,'Fontsize',10);
fig = gcf;
set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); 
% remove unnecessary white space
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
fig.PaperPositionMode   = 'auto';
print(figfilename1,'-dpng', '-r600'); 
plot(coords(:,1),coords(:,2),'k.');
for k=1:6
    text(pzt_x1(k),pzt_y1,num2str(k),'Color','red','FontSize',12);
end
for k=1:6
    text(pzt_x2(k),pzt_y2,num2str(6+k),'Color','red','FontSize',12);
end

axis([0.2 0.3 0.4 0.5]);
print(figfilename2,'-dpng', '-r600'); 
%saveas(h,figfilename,'png');
delete(h);
%---------------------- END OF CODE---------------------- 

% ================ [plot_mesh.m] ================  
