function [wbgrid,wtgrid]=w(R,rx,ry,parameters)
%R: center of wannier state
%r: scalar|array real space position

n=30;
state=1;

Nmax=parameters.Nmax;
bM1=parameters.bM1;
bM2=parameters.bM2;

% meshgrid on honeycomb
a1=-(2*bM1+bM2)/3/n;
a2=(bM1+2*bM2)/3/n;
% [rx,ry]=meshgrid(linspace(-sqrt(3).aM,sqrt(3).aM,Nrx));
[Nrx,Nry]=size(rx);
counter=1;
N=3*n^2+3*n+1;
kxmap=zeros(N,1);
kymap=zeros(N,1);
gauge=zeros(N,1);
expR=zeros(N,1);
psibline=zeros(N,Nrx*Nry);
psitline=zeros(N,Nrx*Nry);
for yindex=-n:n
    for xindex=max(-n,-n+yindex):min(n+yindex,n)
        k=xindex*a1+yindex*a2;
        kxmap(counter)=k(1);
        kymap(counter)=k(2);
        counter=counter+1;
    end
end

parfor index=1:N
    k=[kxmap(index),kymap(index)];
    [~,vec]=energyTMD(k(1),k(2),parameters);
    psib0=sum(vec(1:(2*Nmax+1)^2,state)); %bloch wf at r=(0,0) on the bottom layer
    gauge(index)=conj(abs(psib0)/psib0);
    [ubgrid,utgrid]=u(vec(:,state),rx,ry,parameters);
    expR(index)=exp(-1i*dot(k,R));
    psibgrid=ubgrid.*exp(1i*(k(1)*rx+k(2)*ry));
    psitgrid=utgrid.*exp(1i*(k(1)*rx+k(2)*ry));
    psibline(index,:)=psibgrid(:);
    psitline(index,:)=psitgrid(:);        
end

expRgauge=(gauge.*expR);

wbline=expRgauge.'*psibline/(N);
wtline=expRgauge.'*psitline/(N);

wbgrid=reshape(wbline,[Nrx,Nry]);
wtgrid=reshape(wtline,[Nrx,Nry]);

end


