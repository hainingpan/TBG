Vzlist=linspace(-100,100,100);
clear tlist Ulist;
for i=1:length(Vzlist)
    [tlist{i},Ulist{i}]=t_phase(Vzlist(i));    
end
% figure;plot(Vzlist,phaselist/pi)
% xlabel('Vz (meV)');ylabel('\phi/\pi')
save('tU_vs_Vz.mat','tlist','Ulist','Vzlist');

function [t,U]=t_phase(Vz)
parameters=mainTMD('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',4,'Vz',Vz);

[t,~]=t_calc_func(3,parameters);
 U=U_calc_func(3,parameters);
% re=angle(t{2}(1));
end