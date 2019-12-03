Vrange=linspace(0,20,100);
psirange=linspace(0,360,100);
ldos=cell(length(Vrange),length(psirange));


for Vindex=1:length(Vrange)
    for psiindex=1:length(psirange)
        disp([Vindex,psiindex]);
        psi=psirange(psiindex);
        V=Vrange(Vindex);
        ldos{Vindex,psiindex}.parameters=mainTMD('m',0.35,'psi',psi,'V',V,'w',10);
        [ldos{Vindex,psiindex}.ldosAA,ldos{Vindex,psiindex}.ldosAB,...
            ldos{Vindex,psiindex}.enlist]=LDOS_TMD_r(ldos{Vindex,psiindex}.parameters);
    end
end