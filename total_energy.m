function re=total_energy(energyall,wfall,energyall2,wfall2,Ux,parameters)
%alpha is up, beta is down
%energyall: n generation
%energyall2: n+1 generation

c=cc(energyall,wfall,parameters); 
c2=cc(energyall2,wfall2,parameters);


energyall2_sort=sort(energyall2(:));
mu=energyall2_sort(end/2);
occupied2=energyall2<=mu;

enH=mean(energyall2(occupied2));
K=enH-Ux*(c(1,1)*c2(2,2)+c(2,2)*c2(1,1)-c(1,2)*c2(2,1)-c(2,1)*c2(1,2));
re=K+Ux*(c2(1,1)*c2(2,2)-c2(1,2)*c2(2,1));
end