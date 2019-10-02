parameters=mainTBG();
h=Htb(parameters);
st=1000;
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

enlist=linspace(val(1),val(end),300);
enmapk=zeros(parameters.NN,parameters.NN,length(enlist));
delta=1e-5;
fprintf("LDOS...\n");
for i=1:length(enlist)
    fprintf("i_k=%d\n",i);
    deltaf=reshape(delta./((enlist(i)-val).^2+delta^2),[1,1,st]);
    psifbA=abs(vecfbA).^2;
    psifbB=abs(vecfbB).^2;
    psiftA=abs(vecftA).^2;
    psiftB=abs(vecftB).^2;
    psif=psifbA+psifbB+psiftA+psiftB;
%       psif=psifbA+psifbB;
    enmapk(:,:,i)=sum(deltaf.*psif,3);
end

save(sprintf("NN%da%d.mat",parameters.NN,parameters.a0/parameters.a),'-v7.3');
figure;
xlength=parameters.NN*parameters.a0;
kxlist=2*pi/xlength*(-floor(parameters.NN/2):floor((parameters.NN-1)/2));
kylist=kxlist;
surf(kxlist,1000*enlist,squeeze((enmapk(parameters.NN/2,:,:)))','edgecolor','none');view(2);
