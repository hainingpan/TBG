parameters=mainTMD('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',4,'Vz',100);
neighborlist{1}={[0,0]};
neighborlist{2}={[-1,0],[0,-1],[1,-1],[1,0],[0,1],[-1,1]}; % direction on clock: 12, 10, 8, 6, 4, 2
neighborlist{3}={[-1,-1],[1,-2],[2,-1],[1,1],[-1,2],[-2,1]}; %directin on clock: 11, 9, 7, 5, 4, 1
neighborlist{4}={[-2,0],[0,-2],[2,-2],[2,0],[0,2],[-2,2]}; %direction on clock: 12, 10, 8, 6, 4, 2
neighborlist{5}={[-2,-1],[-1,-2],[1,-3],[2,-3],[3,-2],[3,-1],[2,1],[1,2],[-1,3],[-2,3],[-3,2],[-3,1]}; %counterclockwise from the first point in Quadrant II
neighborlist{6}={[-3,0],[0,-3],[3,-3],[3,0],[0,3],[-3,3]}; %direction on clock: 12, 10, 8, 6, 4, 2
% neighborlist={[-1,0]};

% t{1}= hoppingt(neighbor0list,parameters);
% t{2}= hoppingt(neighbor1list,parameters);
% t{3}= hoppingt(neighbor2list,parameters);
% t{4}= hoppingt(neighbor3list,parameters);
% t{5}= hoppingt(neighbor4list,parameters);
% t{6}= hoppingt(neighbor5list,parameters);
k=4;
t=cell(1,k);
for i=1:k
    t{i}=hoppingt(neighborlist{i},parameters);
end

% t_all= hoppingt([neighbor0list,neighbor1list,neighbor2list,neighbor3list,neighbor4list,neighbor5list],parameters);

% angle(t*1000)
% [rx,ry]=meshgrid(linspace(-2*sqrt(3)*parameters.aM,2*sqrt(3)*parameters.aM,101));
% [wbgrid,wtgrid]=w([0,0],rx,ry,parameters);
% uu=hubbardU(wbgrid,wtgrid,rx,ry,parameters);
% uu*1000

