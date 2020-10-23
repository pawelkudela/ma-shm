function [Exx,Eyy,Ezz,Exy,Eyz,Exz]=strains_spec_p(Npx,Npy,Npz,invJ11,invJ12,invJ13,invJ21,invJ22,invJ23,invJ31,invJ32,invJ33,UX,UY,UZ)
    
% calculate strains at nodes of spectral elements
% parallel mode

% note - brackets are extremely important!
    
    Exx=(Npx*UX).*invJ11+(Npy*UX).*invJ21+(Npz*UX).*invJ31;
    
    Eyy=(Npx*UY).*invJ12+(Npy*UY).*invJ22+(Npz*UY).*invJ32;
    
    Ezz=(Npx*UZ).*invJ13+(Npy*UZ).*invJ23+(Npz*UZ).*invJ33;
    
    B41UX=(Npx*UX).*invJ12+(Npy*UX).*invJ22+(Npz*UX).*invJ32;
    B42UY=(Npx*UY).*invJ11+(Npy*UY).*invJ21+(Npz*UY).*invJ31;
    
    Exy=B41UX+B42UY;
    
    B52UY=(Npx*UY).*invJ13+(Npy*UY).*invJ23+(Npz*UY).*invJ33;
    B53UZ=(Npx*UZ).*invJ12+(Npy*UZ).*invJ22+(Npz*UZ).*invJ32;
    
    Eyz=B52UY+B53UZ;
    
    B61UX=(Npx*UX).*invJ13+(Npy*UX).*invJ23+(Npz*UX).*invJ33;
    B63UZ=(Npx*UZ).*invJ11+(Npy*UZ).*invJ21+(Npz*UZ).*invJ31;
    
    Exz=B61UX+B63UZ;