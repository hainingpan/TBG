function ds=ode(t,s,param)
n=param.n;
Jx=param.Jx;
Jy=param.Jy;
Jz=param.Jz;
Z=param.Z;
nn_lin=param.nn_lin;
N=length(n);
sx=s(1:N);
sy=s(N+1:2*N);
sz=s(2*N+1:3*N);

dsx=(sum(cell2mat(arrayfun(@(i) Jy(i)*sz.*sy(nn_lin{i})-Jz(i)*sy.*sz(nn_lin{i})-Z(i)*sz.*sx(nn_lin{i}),1:6,'UniformOutput',false)),2));
dsy=(sum(cell2mat(arrayfun(@(i) Jz(i)*sx.*sz(nn_lin{i})-Jx(i)*sz.*sx(nn_lin{i})-Z(i)*sz.*sy(nn_lin{i}),1:6,'UniformOutput',false)),2));
dsz=(sum(cell2mat(arrayfun(@(i) Jx(i)*sy.*sx(nn_lin{i})-Jy(i)*sx.*sy(nn_lin{i})+Z(i)*sx.*sx(nn_lin{i})+Z(i)*sy.*sy(nn_lin{i}),1:6,'UniformOutput',false)),2));

ds=[dsx;dsy;dsz];
end