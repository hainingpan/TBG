gamma_kp_x=linspace(0,parameters.kp(1),40);
gamma_kp_y=linspace(0,parameters.kp(2),40);
kp_kn_x=linspace(parameters.kp(1),parameters.kn(1),40);
kp_kn_y=linspace(parameters.kp(2),parameters.kn(2),40);
kn_gamma_x=linspace(parameters.kn(1),0,40);
kn_gamma_y=linspace(parameters.kn(2),0,40);

% kxlist=[gamma_kp_x,kp_kn_x,kn_gamma_x];
% kylist=[gamma_kp_y,kp_kn_y,kn_gamma_y];
kxlist=gamma_kp_x;
kylist=gamma_kp_y;
energylist=zeros(4*(2*parameters.Nmax+1)^2,length(kxlist)); %initialize
for i=1:length(kxlist)
    energylist(:,i)=energyTBG(kxlist(i),kylist(i),parameters);
end

plot(kylist,1000*energylist)
axis([kylist(1),kylist(end),-20,20])

ylim([-20,20])