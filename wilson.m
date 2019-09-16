function re=wilson(n,k1,k2,k3,k4,parameters)
%berry curvature=arg(<u1|u2><u2|u3><u3|u4><u4|u1>)
[~,u1vec]=energy(k1(1),k1(2),parameters);
[~,u2vec]=energy(k2(1),k2(2),parameters);
[~,u3vec]=energy(k3(1),k3(2),parameters);
[~,u4vec]=energy(k4(1),k4(2),parameters);
u1=u1vec(:,n);
u2=u2vec(:,n);
u3=u3vec(:,n);
u4=u4vec(:,n);
re=angle(u1'*u2*u2'*u3*u3'*u4*u4'*u1);
end
