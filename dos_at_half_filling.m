function [energy_at_half,dos]=dos_at_half_filling(theta,Vz)
parameters=mainTMD('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',theta,'Vz',Vz);
neighborlist{1}={[0,0]};
neighborlist{2}={[-1,0],[0,-1],[1,-1],[1,0],[0,1],[-1,1]}; % direction on clock: 12, 10, 8, 6, 4, 2
neighborlist{3}={[-1,-1],[1,-2],[2,-1],[1,1],[-1,2],[-2,1]}; %directin on clock: 11, 9, 7, 5, 4, 1
neighborlist{4}={[-2,0],[0,-2],[2,-2],[2,0],[0,2],[-2,2]}; %direction on clock: 12, 10, 8, 6, 4, 2
neighborlist{5}={[-2,-1],[-1,-2],[1,-3],[2,-3],[3,-2],[3,-1],[2,1],[1,2],[-1,3],[-2,3],[-3,2],[-3,1]}; %counterclockwise from the first point in Quadrant II
neighborlist{6}={[-3,0],[0,-3],[3,-3],[3,0],[0,3],[-3,3]}; %direction on clock: 12, 10, 8, 6, 4, 2

[t,neighborlist]=t_calc_func(3,parameters);

n=50;
counter=1;
N=3*n^2+3*n+1;
kxlist=zeros(N,1);
kylist=zeros(N,1);
a1=-(2*parameters.bM1+parameters.bM2)/3/n;
a2=(parameters.bM1+2*parameters.bM2)/3/n;
for yindex=-n:n
    for xindex=max(-n,-n+yindex):min(n+yindex,n)
        k=xindex*a1+yindex*a2;
        kxlist(counter)=k(1);
        kylist(counter)=k(2);
        counter=counter+1;
    end
end
k=3;
energylist=real(tb([neighborlist{1:k+1}],[t{1:k+1}],kxlist,kylist,parameters));

eta=5e-4;

energylist_sort=sort(energylist);
energy_at_half=energylist_sort(floor(end/2));
deltaf=eta./((energy_at_half-energylist).^2+eta^2);
dos=sum(deltaf(:));
dos=dos/length(energylist);

end

