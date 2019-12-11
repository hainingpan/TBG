parameters.kpp=[-parameters.kp(1),parameters.kp(2)];

kpp_gamma_x=linspace(parameters.kpp(1),0,40);
kpp_gamma_y=linspace(parameters.kpp(2),0,40);
gamma_kn_x=linspace(0,parameters.kn(1),40);
gamma_kn_y=linspace(0,parameters.kn(2),40);
kn_kp_x=linspace(parameters.kn(1),parameters.kp(1),40);
kn_kp_y=linspace(parameters.kn(2),parameters.kp(2),40);
kp_kpp_x=linspace(parameters.kp(1),parameters.kpp(1),80);
kp_kpp_y=linspace(parameters.kp(2),parameters.kpp(2),80);



kxlist=[kpp_gamma_x,gamma_kn_x,kn_kp_x,kp_kpp_x];
kylist=[kpp_gamma_y,gamma_kn_y,kn_kp_y,kp_kpp_y];

segment=sqrt(diff(kxlist).^2+diff(kylist).^2);
klist=[0,cumsum(segment)];


energylist=real(tb([neighbor0list,neighbor1list,neighbor2list,neighbor3list,neighbor4list,neighbor5list],[t0,t1,t2,t3,t4,t5],kxlist,kylist,parameters));

figure;
plot(klist,1000*energylist)
hold on
line(klist([40,40]),[-200,30],'color','k','LineStyle','--','HandleVisibility','off');
line(klist([80,80]),[-200,30],'color','k','LineStyle','--','HandleVisibility','off');
line(klist([120,120]),[-200,30],'color','k','LineStyle','--','HandleVisibility','off');


xticks(klist([1,40,80,120,200]))
xticklabels({'\kappa_+^\prime','\gamma','\kappa_-','\kappa_+','\kappa_+^\prime'});
xlim([klist(1),klist(end)])
% ylim([min(1000*energylist(5,:)),1.2*max(1000*energylist(1,:))])