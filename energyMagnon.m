function val=energyMagnon(kx,ky,param)
N=param.N;
phi=param.phi;
S=param.S;
es=@(nlist) sum(arrayfun(@(x) exp(1i*[kx,ky]*N{x}),nlist));
A=[0,es([7,3,5]),es([3,5,4]);
    es([7,2,6]),0,es([7,3,5]);
    es([1,2,6]),es([7,2,6]),0];
SS=[0,S{1}*S{2}',S{1}*S{3}';
    S{2}*S{1}',0,S{2}*S{3}';
    S{3}*S{1}',S{3}*S{2}',0];
SA=SS.*A;
cr=@(x,y) x(1)*y(2)-x(2)*y(1);
S_DMI=[0,cr(S{1},S{2}),cr(S{3},S{1});
    cr(S{1},S{2}),0,cr(S{2},S{3});
    cr(S{3},S{1}),cr(S{2},S{3}),0];
S_DMIA=S_DMI.*A;

H=1/4*2*[A,-A;-A,A]+cos(2*phi)/4*2*[SA,SA;SA,SA]+cos(2*phi)/4*2*6*eye(6)+sin(2*phi)/4*2*[S_DMIA,S_DMIA;S_DMIA,S_DMIA];
H=diag([1,1,1,-1,-1,-1])*H;
val=sort(real(eig(H)),'descend');
