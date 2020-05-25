function [h,h1,h2,em,rhom,nim,vol,ef,rhof,nif,alpha]=lay_const(lay,lh,i_em,i_ef,i_rhom,i_rhof,i_nim,i_nif,lvol,lalpha,lmat,lfib)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate layer's constants
h2=zeros(1,lay);
h1=zeros(1,lay);
em=zeros(1,lay);
rhom=zeros(1,lay);
nim=zeros(1,lay);
vol=zeros(1,lay);
ef=zeros(1,lay);
rhof=zeros(1,lay);
nif=zeros(1,lay);
alpha=zeros(1,lay);
%%
h=sum(lh);
h2(1)=h/2; h1(1)=h/2-lh(1);
for i=2:lay;
  h2(i)=h2(i-1)-lh(i-1);
  h1(i)=h1(i-1)-lh(i);
end;
for i=1:lay;
% matrix properties %%%%%%fibres properties %%%%%%%%%%%%%%%%%
  em(i)=i_em(lmat(i));                ef(i)=i_ef(lfib(i));
  rhom(i)=i_rhom(lmat(i));            rhof(i)=i_rhof(lfib(i));
  nim(i)=i_nim(lmat(i));              nif(i)=i_nif(lfib(i));
  vol(i)=lvol(i);                     alpha(i)=lalpha(i);
end;
