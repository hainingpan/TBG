Vrange=linspace(0,20,50);
psirange=linspace(0,360,50);
int1=zeros(length(Vrange),length(psirange));
int2=zeros(length(Vrange),length(psirange));
NV=length(Vrange);
Npsi=length(psirange);
parfor Vindex=1:NV
    for psiindex=1:Npsi
        disp([Vindex,psiindex]);
        psi=psirange(psiindex);
        V=Vrange(Vindex);
        parameters=mainTMD('m',0.35,'psi',psi,'V',V,'w',10,'theta',3);
        [int1(Vindex,psiindex),int2(Vindex,psiindex)]=LDOS_TMD_r(parameters);
    end
end

save('int1.dat','int1','-ascii');
save('int2.dat','int2','-ascii');

% parameters=mainTMD('m',0.35,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',3);
% parameters=mainTMD('m',0.35,'psi',100,'V',20,'w',10);
% [ldosAA,ldosAB,intAA,intAB,enlist]=LDOS_TMD_r(parameters);
% [int1,int2]=LDOS_TMD_r(parameters);

% [ldos,enlist]=LDOS_TMD_rx(parameters);