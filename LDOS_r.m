enlist=linspace(val(1),val(end),300);
enmapr=zeros(parameters.NN,parameters.NN,length(enlist));
for i=1:length(enlist)
fprintf("i_r=%d\n",i);
deltaf=reshape(delta./((enlist(i)-val).^2+delta^2),[1,1,st]);
psibA=abs(vecbA).^2;
psibB=abs(vecbB).^2;
psi=psibA+psibB;
enmapr(:,:,i)=sum(deltaf.*psi,3);
end
