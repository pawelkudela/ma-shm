% make png pictures from fig files

clc; clear all; close all;
disp('making figures ...');
surface='bottom'; % options: surface='top';surface='bottom'
extern_path='';
% 827:835 damage 3 (crack) frames 64 128 160
% 836:844 damage 4 (added_mass) frames 64 128 160
% 845:851 triangle frames 64 128 60
% 853:861 stiffener 0 deg frames 64 128 160
% 862:870 stiffener 90 deg 64 128 160


i=0;
for k_test=[862:870];
    i=i+1;
 
os='win';


%alpha=0.05;
%alpha=0.01;
%alpha=0.5;% no crack
%alpha=0.5;% 
%alpha=0.3;
alpha=0.2;

switch os 
    case 'win'
       
        filename=[extern_path,'outputs\output',num2str(k_test),'\plate_w_',num2str(k_test),'_500x500',surface,'.mat']; 
    case 'linux'
        filename=['outputs/output',num2str(k_test),'/plate_w_',num2str(k_test),'_500x250',surface,'.mat'];
    otherwise
        disp('Unknown system');
end
load(filename);
[s1,s2,s3]=size(Data);
c=0;
Smax=max(max(max(abs(Data))));
%Smax=3.2e-9;
Smin=-Smax;

for n=[64,128,160]%[32:32:512] 

    disp('frame number ... total number of frames');
    [n,s3]
    c=c+1;
    switch os
    case 'win'
        
%         figfilename=[extern_path,'outputs\output',num2str(k_test),'\',num2str(n),'frame_w_',surface];
%         figfilename2=[extern_path,'outputs\output',num2str(k_test),'\',num2str(n),'frame_w_',surface,'.fig'];
        figfilename=['outputs\output',num2str(k_test),'\A',num2str(i),'-',num2str(n),'_stiffener_90deg']; 
        figfilename2=['outputs\output',num2str(k_test),'\A',num2str(i),'-',num2str(n),'_stiffener_90deg.fig'];
       
    case 'linux'
       
        figfilename=[extern_path,'outputs\output',num2str(k_test),'\',num2str(n),'frame_w_',surface];
        figfilename2=[extern_path,'outputs\output',num2str(k_test),'\',num2str(n),'frame_w_',surface,'.fig'];
       
    end
     surf(squeeze(Data(:,:,n)));axis equal;view(2);set(gcf,'Renderer','zbuffer');shading interp;set(gcf,'color','white');axis off;
    
     %h=open(filename);set(h,'visible','on');
     %colorbar off;
     Smax=max(max(abs(squeeze(Data(:,:,n)))));
     Smin=-Smax;
     view(2);
     caxis([alpha*Smin,alpha*Smax]);
     %colormap jet;
     saveas(gcf,figfilename2,'fig');
     close all;
     h=open(figfilename2);set(h,'visible','on');
     set(gca, 'Position',[0 0 1 1]); % figure without white border
     %set(h, 'Units','centimeters', 'Position',[10 10 10 5]); % size 10cm by 5cm
     %set(h, 'Units','centimeters', 'Position',[10 10 10 6]); % size 10cm by 6cm (triangle)
     set(h, 'Units','centimeters', 'Position',[10 10 5 5.04]); % size 5cm by 5cm (square)
     set(h,'PaperPositionMode','auto');
     
     %colorbar;
     %print -dtiff -r0 out.tiff
     print('-dpng', '-r600',figfilename); %
     delete(figfilename2);
     close all;
end
end
