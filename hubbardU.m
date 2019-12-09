function re=hubbardU(wbgrid,wtgrid,rx,ry,parameters)
bM1=parameters.bM1;
bM2=parameters.bM2;
% [rx,ry]=meshgrid(linspace(-2*sqrt(3)*parameters.aM,2*sqrt(3)*parameters.aM,101));
[Nx,Ny]=size(rx);
% [wbgrid,wtgrid]=w(R,rx,ry,parameters);
alpha=0.00729735;
Nkx=200;
Nky=200;
kxlist=linspace(-3*bM2(1),3*bM2(1),Nkx);
kylist=linspace(-3*bM2(1),3*bM2(1),Nky);

[kxmap,kymap]=meshgrid(kxlist,kylist);
Fb=exp(1i*(repmat(reshape(kxmap,1,1,[]),Nx,Ny,1).*rx+repmat(reshape(kymap,1,1,[]),Nx,Ny,1).*ry)).*abs(wbgrid).^2;
intb=trapz(rx(1,:),trapz(ry(:,1),Fb,2));
Ft=exp(1i*(repmat(reshape(kxmap,1,1,[]),Nx,Ny,1).*rx+repmat(reshape(kymap,1,1,[]),Nx,Ny,1).*ry)).*abs(wtgrid).^2;
intt=trapz(rx(1,:),trapz(ry(:,1),Ft,2));

Fk=alpha/(2*pi)*1./(sqrt(kxmap.^2+kymap.^2)).*reshape(abs(intb).^2+abs(intt).^2,Nkx,Nky);
re=trapz(kxlist,trapz(kylist,Fk,2));
end