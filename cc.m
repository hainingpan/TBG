function c=cc(energyall,wfall,parameters)
energyall_sort=sort(energyall(:));
mu=energyall_sort(end/2);
occupied=energyall<=mu;
wf_occupied=[wfall{occupied}];
N=length(wf_occupied);
c=zeros(2);
Q=[0,0;parameters.kn-parameters.kp;parameters.kp-parameters.kn];
r=[0,0;sqrt(3)/2*parameters.aM,1/2*parameters.aM;sqrt(3)*parameters.aM,parameters.aM];

for alpha=1:2
    for beta=1:2
        for j=1:3 %j=A,B,C
            re=0;
            for k=1:N
                 wf=reshape(wf_occupied(:,k),3,2);
                 for q=1:3
                     bra=wf(:,alpha);
                     ketindex=mod([1,2,3]+q-1-1,3)+1;
                     ket=wf(ketindex,beta);
                     re=re+bra'*ket*exp(1i*(Q(q,:)*r(j,:)'));
                 end
            end
        end
        c(alpha,beta)=re/N/3;
    end
end

end