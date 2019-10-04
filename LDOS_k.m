parameters=mainTBG();
parameters.a0=parameters.aM/50;
h=Htb(parameters);
st=500;
fprintf("diagonalizing...\n");
[vec,val]=eigs(h,st,'sm');
val=diag(val);
[val,order]=sort(val);
vec=vec(:,order);
vecbA=reshape(vec(1:parameters.NN^2,:),[parameters.NN,parameters.NN,st]);
vecbB=reshape(vec(1+parameters.NN^2:2*parameters.NN^2,:),[parameters.NN,parameters.NN,st]);
vectA=reshape(vec(1+2*parameters.NN^2:3*parameters.NN^2,:),[parameters.NN,parameters.NN,st]);
vectB=reshape(vec(1+3*parameters.NN^2:4*parameters.NN^2,:),[parameters.NN,parameters.NN,st]);
vecfbA=fft2(vecbA);
vecfbB=fft2(vecbB);
vecftA=fft2(vectA);
vecftB=fft2(vectB);

fprintf("FFT...\n");
for i=1:st
    vecfbA(:,:,i)=fftshift(vecfbA(:,:,i));
    vecfbB(:,:,i)=fftshift(vecfbB(:,:,i));
    vecftA(:,:,i)=fftshift(vecftA(:,:,i));
    vecftB(:,:,i)=fftshift(vecftB(:,:,i));
end

enlist=linspace(val(1),val(end),parameters.NN);
enmapk=zeros(parameters.NN,parameters.NN,length(enlist));
delta=1e-3;
fprintf("LDOS...\n");
for i=1:length(enlist)
%     fprintf("i_k=%d\n",i);
    deltaf=reshape(delta./((enlist(i)-val).^2+delta^2),[1,1,st]);
    psifbA=abs(vecfbA).^2;
    psifbB=abs(vecfbB).^2;
    psiftA=abs(vecftA).^2;
    psiftB=abs(vecftB).^2;
    psif=psifbA+psifbB+psiftA+psiftB;
%       psif=psifbA+psifbB;
    enmapk(:,:,i)=sum(deltaf.*psif,3);
end

save(sprintf("NN%da%d.mat",parameters.NN,parameters.aM/parameters.a0),'-v7.3');
% load("NN500a1.mat");
fig=figure;
xlength=parameters.NN*parameters.a0;
kxlist=2*pi/xlength*(-floor(parameters.NN/2):floor((parameters.NN-1)/2));
kylist=kxlist;
surf(kxlist,1000*enlist,squeeze((enmapk(parameters.NN/2,:,:)))','edgecolor','none');view(2);
title(sprintf("a0=aM/%d,NN=%d,Log(LDOS)",parameters.aM/parameters.a0,parameters.NN));
savefig(fig,sprintf("NN%da%d.fig",parameters.NN,parameters.a0/parameters.a))
fig2=figure;
surf(kxlist,1000*enlist,squeeze(log(enmapk(parameters.NN/2,:,:)))','edgecolor','none');view(2);
title(sprintf("a0=aM/%d,NN=%d,LDOS",parameters.aM/parameters.a0,parameters.NN));
savefig(fig2,sprintf("NN%da%d_log.fig",parameters.NN,parameters.a0/parameters.a))
fig3=figure;
surf(kxlist,kylist,log(squeeze(enmapk(:,:,parameters.NN/2))),'EdgeColor','none');view(2);
title(sprintf("a0=aM/%d,NN=%d,E=",parameters.aM/parameters.a0,parameters.NN,enlist(parameters.NN/2)));
savefig(fig3,sprintf("NN%da%d_kxky.fig",parameters.NN,parameters.a0/parameters.a))


