function psi2=u2(vec,r,parameters)
Nvec=length(vec);
Nmax=(sqrt(Nvec/2)-1)/2;
Nrange=-Nmax:Nmax;
[h1index,h2index]=meshgrid(Nrange,Nrange);

expo=exp(1i*dot((h1index(:)*parameters.bM1+h2index(:)*parameters.bM2),repmat(r,(2*Nmax+1)^2,1),2));
expo2=repmat(expo,[1,2*(2*Nmax+1)^2]);

psit=sum(expo2.*vec(1:(2*Nmax+1)^2,:),1);
psib=sum(expo2.*vec((2*Nmax+1)^2+1:2*(2*Nmax+1)^2,:),1);
psi2=abs(psit).^2+abs(psib).^2;
end