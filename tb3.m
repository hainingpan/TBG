function eigval=tb3(bond,tmat,kxlist,kylist,parameters)

Nk=length(kxlist);
AB=parameters.aM2-parameters.aM1;
AC=2*AB;
eigval=zeros(Nk,3);
parfor i=1:Nk
    kx=kxlist(i);
    ky=kylist(i);
    A=diag([1,exp(1i*([kx,ky]*AB')),exp(1i*([kx,ky]*AC'))]);
    n=length(bond);
    A1=parameters.A1;
    A2=parameters.A2;
    re=0;
    for j=1:n
        a=bond{j}(1)*A1+bond{j}(2)*A2;
        re=re+exp(-1i*([kx,ky]*a'))*tmat{j};
    end
    re=A'*re*A;
    eigval(i,:)=sort(eig(re));
end

