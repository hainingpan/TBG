function [S,mu]=S_calc(energyall,wfall)
%S(a,b,c), a,b=1(up) or 2(down) c=1,2,3 (0, Q, -Q)

energyall_sort=sort(energyall(:));
mu=energyall_sort(end/2);
occupied=(energyall<=mu);
wf_occupied=[wfall{occupied}];
S=zeros(2,2,3);
for a=1:2
    for b=1:2
        for c=1:2
            re=0;
            for i=1:length(wf_occupied)
                wf=reshape(wf_occupied(:,i),3,2);
                for q=1:3
                    re=re+conj(wf(q,(3-a)))*wf(mod(q+c-1-1,3)+1,(3-b))*((a==b)-(a~=b));   %3-a flips spin; 
                end
            end
            S(a,b,c)=re/length(wf_occupied);
        end
    end
end
S(:,:,3)=S(:,:,2)';
                