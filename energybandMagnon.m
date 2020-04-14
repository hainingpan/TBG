param=mainMagnon('phi',pi);
Gamma=[0,0];
K=[4*pi/3,0];
M2=[0,2*pi/sqrt(3)];
Nx=50;
Ny=50;
kxlist=linspace(-K(1),K(1),Nx);
kylist=linspace(-M2(2),M2(2),Ny);

energymap=zeros(Nx,Ny,6);
for i=1:Nx
    for j=1:Ny
        kx=kxlist(i);
        ky=kylist(j);
        energymap(i,j,:)=energyMagnon(kx,ky,param);
    end
end

figure;
contourf(kxlist/K(1),kylist/M2(2),energymap(:,:,3),20);
xlabel('kx');
ylabel('ky');
