% make png pictures from fig files

clc; clear all; close all;
disp('making figures ...');
surface='top'; % options: surface='top';surface='bottom'
% 827:835 damage 3 (crack) [48,64,96]
% 836:844 damage 4 (added_mass) [48,64,96]
% 845:851 triangle [48,64,96]
% 853:861 stiffener 0 deg [64,128,160]
% 862:870 stiffener 90 deg [64,128,160]

i=0;
for k_test=[845:851];
 i=i+1;
os='win';


%alpha=0.05;
%alpha=0.01;
%alpha=0.5;% no crack
alpha=0.5;% 
%alpha=0.3;


switch os 
    case 'win'
       
        filename=['outputs\output',num2str(k_test),'\plate_in_plane_',num2str(k_test),'_500x500',surface,'.mat']; 
    case 'linux'
        filename=['outputs/output',num2str(k_test),'/plate_in_plane_',num2str(k_test),'_500x250',surface,'.mat'];
    otherwise
        disp('Unknown system');
end
load(filename);
[s1,s2,s3]=size(Data);
c=0;
Smax=max(max(max(abs(Data))));
%Smax=3.2e-9;
Smin=-Smax;

for n=[48,64,96]%[64,128,160]%[48,64,96]

    disp('frame number ... total number of frames');
    [n,s3]
    c=c+1;
    switch os
    case 'win'
        
        
        figfilename=['outputs\output',num2str(k_test),'\A',num2str(i),'-',num2str(n),'_in_plane_triangle']; 
        figfilename2=['outputs\output',num2str(k_test),'\A',num2str(i),'-',num2str(n),'_in_plane_triangle.fig'];
       
    case 'linux'
       
         figfilename=['outputs\output',num2str(k_test),'\',num2str(n),'frame_in_plane_',surface];
         figfilename2=['outputs\output',num2str(k_test),'\',num2str(n),'frame_in_plane_',surface,'.fig'];
       
    end
     surf(squeeze(Data(:,:,n)));axis equal;view(2);set(gcf,'Renderer','zbuffer');shading interp;set(gcf,'color','white');axis off;
     %h=open(filename);set(h,'visible','on');
     %colorbar off;
     Smax=max(max(abs(squeeze(Data(:,:,n)))));
     Smin=-Smax;
     view(2);
     caxis([alpha*Smin,alpha*Smax]);
     saveas(gcf,figfilename2,'fig');
     close all;
     h=open(figfilename2);set(h,'visible','on');
     set(gca, 'Position',[0 0 1 1]); % figure without white border
     %set(h, 'Units','centimeters', 'Position',[10 10 10 5]); % size 10cm by 5cm
     %set(h, 'Units','centimeters', 'Position',[10 10 10 6]); % size 10cm by 6cm (triangle)
     set(h, 'Units','centimeters', 'Position',[10 10 5 5.04]); % size 5cm by 5cm (square)
     set(h,'PaperPositionMode','auto');
     colormap jet;
     %colorbar;
     %print -dtiff -r0 out.tiff
     print('-dpng', '-r600',figfilename); %
     delete(figfilename2);
     close all;
end
end
