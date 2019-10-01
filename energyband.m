kxlist=linspace(-parameters.kn(1),parameters.kn(1),100);
kylist=linspace(-parameters.kp(2),parameters.kp(2),100);
energymap=zeros(length(kxlist),length(kylist),4*(2*parameters.Nmax+1)^2);
for iky=1:length(kxlist)
    disp(iky)
    for ikx=1:length(kylist)
        energymap(ikx,iky,:)=energyTBG(kxlist(ikx),kylist(iky),parameters);
    end
end