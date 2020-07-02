function A=frame2image_resized(frame, filename,desired_size_in_pixels)

% function converts frame of propagating waves (double values)
% to grayscale image: values [0, 255] and saves it to disc
% with the resizing to the desired size in pixels

[nx,ny]=size(frame);
[X,Y] = meshgrid(1:ny,1:nx);                                        % original value grid
Nx=desired_size_in_pixels;
Ny=desired_size_in_pixels;
[XI,YI] = meshgrid(1:(ny-1)/(Ny-1):ny,1:(nx-1)/(Nx-1):nx);          % new value grid
            
 frame_interp = interp2(X,Y,frame,XI,YI,'spline');

% image is flipped due to coordinates in the left upper corner
% whereas in the frame coordinates are in the left lower corner
% normalization 1
 A=uint8(((1+ frame_interp / max(max(abs(frame_interp))) )/2 )*255); 
 % normalization 2
%  f_mean=mean(mean(frame));
%  f_std=std(std(frame));
%  frame=(frame-f_mean)/f_std; % data centering at 0 is not necessary because displacement field is centered at zero
%  A=frame/max(max(abs(frame)));
%  A=uint8(((1+A)/2)*255);
 % correct flipping
 A=flipud(A);
 imwrite(A,[filename,'.png'],'png');