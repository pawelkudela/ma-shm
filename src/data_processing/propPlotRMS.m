function [h1] = propPlotRMS(RMS,ct,WL,PRT,n,Cmap,file)

h1 = figure('Name',file);
%h1 = figure('Position',[100 100 240+100 120+100],'Name',file) 

RMSc = RMS(1+ct:end-ct,1+ct:end-ct);

[rows cols] = size(RMS);

Xc = ct/(cols-1)*WL(1):WL(1)/(cols-1):WL(1)-ct/(cols-1)*WL(1);
Yc = ct/(rows-1)*WL(2):WL(2)/(rows-1):WL(2)-ct/(rows-1)*WL(2);

surf(Xc,Yc,RMSc);
max5 = max(max(RMSc));
min5 = min(min(RMSc));
shading interp; view(2); set(gcf,'Renderer','zbuffer'); axis equal; %axis off; 
caxis([min5 max5/n]);
%caxis([0.08e-16 max5/n]);
%caxis([0.2e-3 max5/n]);
ylim([ct/(cols-1) WL(2)-ct/(cols-1)]);
xlim([ct/(rows-1) WL(1)-ct/(rows-1)]);

xlabel('x (mm)','FontSize',16,'FontName','Arial');
ylabel('y (mm)','FontSize',16,'FontName','Arial');
set(gca,'FontSize',14,'FontName','Arial');
%set(gca,'YTick',[0 0.1 0.2 0.3 0.4])
%set(gca,'XTick',[0 0.1 0.2 0.3 0.4])

%xlabel('x [m]','FontSize',24,'FontName','Arial');
%ylabel('y [m]','FontSize',24,'FontName','Arial');
%set(gca,'FontSize',20,'FontName','Arial');
%set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5])
%set(gca,'XTick',[0 0.1 0.2 0.3 0.4 0.5])

%ZOOM
    %set(gca,'FontSize',50,'FontName','Arial');
    %ylim([0.229 0.27])
    %xlim([0.229 0.27])
    %set(gca,'YTick',[0.23 0.25 0.27])
    %set(gca,'XTick',[0.23 0.25 0.27]) 
    

colormap(Cmap)
if PRT == 1
    print('-dpng','-r300',['RMS\', file, '.png']);   
end

%{
X = 0:WL(1)/(cols-1):WL(1);
Y = 0:WL(2)/(rows-1):WL(2);
surf(X,Y,RMS);
max5 = max(max(RMS));
min5 = min(min(RMS));
shading interp; view(0,-90); set(gcf,'Renderer','zbuffer'); axis equal; %axis off; 
caxis([min5 max5/n]);
ylim([0 WL(1)]);
xlim([0 WL(2)]);
xlabel('x [m]');
ylabel('y [m]');
%}