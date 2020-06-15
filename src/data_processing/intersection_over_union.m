% Intersection over union
 clear all;close all;clc
% Input
n=100;
A=zeros(n,n);
B=zeros(n,n);
% case 1
% delam_a=zeros(10,10)+1;
% delam_b=zeros(8,8)+1;
% A(20:20+10-1,20:20+10-1)=delam_a;
% B(18:18+8-1,18:18+8-1)=delam_b;
% case 2
% delam_a=zeros(10,10)+1;
% delam_b=delam_a;
% A(20:20+10-1,20:20+10-1)=delam_a;
% B(20:20+10-1,20:20+10-1)=delam_b;
% case 3
delam_a=zeros(10,10)+1;
delam_b=zeros(10,10)+1;
A(20:20+10-1,20:20+10-1)=delam_a;
B(40:40+10-1,20:20+10-1)=delam_b;
% case 4 multi delaminations
delam_a1=zeros(10,10)+1;
delam_b1=zeros(10,10)+1;
A(20:20+10-1,20:20+10-1)=delam_a1;
B(40:40+10-1,20:20+10-1)=delam_b1;

delam_a2=zeros(10,10)+1;
delam_b2=zeros(8,8)+1;
A(20:20+10-1,60:60+10-1)=delam_a2;
B(18:18+8-1,61:61+8-1)=delam_b2;
%%
spy(A); hold on; spy(B,'r');
% Intersection
I=A.*B;
spy(I,'g')
% Unity
U=A+B-A.*B;
figure;spy(U)
% Intersection over Union
if(nnz(U)~=0)
    IoU=nnz(I)/nnz(U);
else
    IoU=0;
end
IoU
IoU2=intersect_over_union_fun(A,B);
IoU2
