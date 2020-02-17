function [X,Y]=spec_nodes(ksi,x,y)
% generete coordinates of the nodes of spectral element
% taking into account 4 corner nodes specified by vectors x and y
% nodes distribution is specified by ksi (local coordinates of the 1D element)

c=0;
    n1=length(ksi);
    n2=length(ksi);
    
    x1=x(1);x2=x(2);x3=x(3);x4=x(4);
    y1=y(1);y2=y(2);y3=y(3);y4=y(4);
    bx=x4-x1;
    by=y4-y1;
    cx=x3-x2;
    cy=y3-y2;
    for m=1:n2
        x0=x1+bx/2*(1+ksi(m));
        y0=y1+by/2*(1+ksi(m));
        xc=x2+cx/2*(1+ksi(m));
        yc=y2+cy/2*(1+ksi(m));
        ax=xc-x0;
        ay=yc-y0;
        for n=1:n1
            c=c+1;
            X(c)=x0+ax/2*(1+ksi(n));
            Y(c)=y0+ay/2*(1+ksi(n));
        end
    end