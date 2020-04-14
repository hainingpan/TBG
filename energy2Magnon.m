function re=energy2Magnon(kx,ky,param)
phi=param.phi;
S=param.S;
n=param.n;
dot=S{1}*S{2}';
cr=@(x,y) x(1)*y(2)-x(2)*y(1);
cross=cr(S{1},S{2});
sumcos=sum(arrayfun(@(x) cos([kx,ky]*n{x}),1:6));
% re=sumcos;
A=(1/2+cos(2*phi)*dot/2+sin(2*phi)*cross/2)*sumcos-6*(cos(2*phi)*dot+sin(2*phi)*cross);
B=(-1/2+cos(2*phi)*dot/2+sin(2*phi)*cross/2)*sumcos;
re=sqrt(A^2-B^2);
end