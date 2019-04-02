function plot_ellipse(xCenter,yCenter,xRadius,yRadius,rotAngle,N)
% function for plotting ellipse with center point (xCenter,yCenter)
% semi-major axis = xRadius
% semi-minor axis = yRadius
% rotation angle rotAngle [deg]
% N - numer of points in adjacent grid

alpha = [0 : 0.01 : 2*pi,0];
theta=rotAngle*pi/180;
% ellipse
x = xRadius * cos(alpha) ;
y = yRadius * sin(alpha) ;
% rotation
x1=x*cos(theta)-y*sin(theta);
y1=x*sin(theta)+y*cos(theta);
x1=x1+ xCenter;
y1=y1+ yCenter;
plot(x1, y1, 'LineWidth', 3);
axis square;
xlim([-5 N+5]);
ylim([-5 N+5]);
grid on;

% Draw lines along the axes
alpha=0;
x0 = xRadius * cos(alpha) ;
y0 = yRadius * sin(alpha) ;
% rotation
x0r=x0*cos(theta)-y0*sin(theta);
y0r=x0*sin(theta)+y0*cos(theta);
x0r=x0r+ xCenter;
y0r=y0r+ yCenter;
hold on;
line([xCenter, x0r], [yCenter, y0r],'LineWidth', 4, 'Color', [1,0,0]);
alpha=pi/2;
x0 = xRadius * cos(alpha) ;
y0 = yRadius * sin(alpha) ;
% rotation
x0r=x0*cos(theta)-y0*sin(theta);
y0r=x0*sin(theta)+y0*cos(theta);
x0r=x0r+ xCenter;
y0r=y0r+ yCenter;
hold on;
line([xCenter, x0r], [yCenter, y0r],'LineWidth', 4, 'Color', [1,0,0]);
% Enlarge figure to full screen.
set(gcf,'Color','White')
%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);