figure;
philist=2/3:1/6:1;
for i=1:length(philist)
    disp(philist(i));
end
legend

function disp(phi)
param=mainMagnon('phi',phi*pi);
Gamma=[0,0];
K=param.k;
M=param.m;
Gamma_M_x=linspace(Gamma(1),M(1),100);
Gamma_M_y=linspace(Gamma(2),M(2),100);
M_K_x=linspace(M(1),K(1),100);
M_K_y=linspace(M(2),K(2),100);
K_Gamma_x=linspace(K(1),Gamma(1),100);
K_Gamma_y=linspace(K(2),Gamma(2),100);

kxlist=[Gamma_M_x,M_K_x,K_Gamma_x];
kylist=[Gamma_M_y,M_K_y,K_Gamma_y];

segment=sqrt(diff(kxlist).^2+diff(kylist).^2);
klist=[0,cumsum(segment)];

energylist=zeros(1,length(kxlist)); %initialize
for i=1:length(kxlist)
    energylist(i)=energy2Magnon(kxlist(i),kylist(i),param);
end
% figure;
plot(klist,energylist,'DisplayName',strcat(num2str(phi),'\pi'));
hold on;

xticks(klist([1,100,200,300]))
xticklabels({'\Gamma','M','K','\Gamma'});
end
