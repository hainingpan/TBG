gamma_kp_x=linspace(0,parameters.kn(1),20);
gamma_kp_y=linspace(0,parameters.kn(2),20);
kp_kn_x=linspace(parameters.kn(1),parameters.kp(1),20);
kp_kn_y=linspace(parameters.kn(2),parameters.kp(2),20);
kn_gamma_x=linspace(parameters.kp(1),0,20);
kn_gamma_y=linspace(parameters.kp(2),0,20);

kxlist=[gamma_kp_x,kp_kn_x,kn_gamma_x];
kylist=[gamma_kp_y,kp_kn_y,kn_gamma_y];
energylist=zeros(2*(2*parameters.Nmax+1)^2,length(kxlist)); %initialize
for i=1:length(kxlist)
    energylist(:,i)=energy(kxlist(i),kylist(i),parameters);
end

plot(1:length(kxlist),1000*energylist)
axis([1,40,0,25])
