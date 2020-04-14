param=mainMagnon('phi',1.1*pi);
Gamma=[0,0];
K=[0,4*pi/3];
M2=[2*pi/sqrt(3),0];
Nx=100;
Ny=50;
kxlist=linspace(-M2(1),M2(1),Nx);
kylist=linspace(-K(2),K(2),Ny);

energymap=zeros(Nx,Ny);
for i=1:Nx
    for j=1:Ny
        kx=kxlist(i);
        ky=kylist(j);
        energymap(i,j)=energy2Magnon(kx,ky,param);
    end
end

figure;
contourf(kxlist,kylist,energymap',20);
xlabel('k_x');
ylabel('k_y');
daspect([1,1,1]);