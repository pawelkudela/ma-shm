function [invJ11,invJ12,invJ21,invJ22]=inv_jacp_2D(J11,J12,J21,J22)

a=1./(J11.*J22-J12.*J21);
invJ11=a.*J22;
invJ12=a.*(-J12);

invJ21=a.*(-J21);
invJ22=a.*J11;

