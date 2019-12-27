thetalist=linspace(1,5,100);
clear tlist Ulist;
parfor i=1:length(thetalist)
    
    [tlist{i},Ulist{i}]=t_phase(thetalist(i));    
end
% figure;plot(Vzlist,phaselist/pi)
% xlabel('Vz (meV)');ylabel('\phi/\pi')
save('tU_vs_theta.mat','tlist','Ulist','thetalist');

function [t,U]=t_phase(theta)
parameters=mainTMD('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',theta,'Vz',40);

[t,~]=t_calc_func(1,parameters);
 U=U_calc_func(0,parameters);
% re=angle(t{2}(1));
end