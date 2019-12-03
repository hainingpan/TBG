function psi2=u2(vec,r,parameters)
Nvec=length(vec);
Nmax=(sqrt(Nvec/2)-1)/2;
Nrange=-Nmax:Nmax;
psi2=zeros(2*(2*Nmax+1)^2,1);
[h1index,h2index]=meshgrid(Nrange,Nrange);

[h1matX,h1matY]=meshgrid(h1index(:));
[h2matX,h2matY]=meshgrid(h2index(:));
h1mat=h1matX-h1matY;
h2mat=h2matX-h2matY;

for i=1:Nvec
   [ah1mat,ah2mat]=meshgrid(vec(1:(2*Nmax+1)^2,i));
   psit2=sum(exp(1i*dot((h1mat(:)*parameters.b1+h2mat(:)*parameters.b2),repmat(r,(2*Nmax+1)^4,1),2)).*ah1mat(:).*conj(ah2mat(:)));
   [ah1mat,ah2mat]=meshgrid(vec((2*Nmax+1)^2+1:2*(2*Nmax+1)^2,i));
   psib2=sum(exp(1i*dot((h1mat(:)*parameters.b1+h2mat(:)*parameters.b2),repmat(r,(2*Nmax+1)^4,1),2)).*ah1mat(:).*conj(ah2mat(:)));
   psi2(i)=psit2+psib2;
end
end