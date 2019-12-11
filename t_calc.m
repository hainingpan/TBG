parameters=mainTMD('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',4,'Vz',0);
neighbor0list={[0,0]};
neighbor1list={[-1,0],[0,-1],[1,-1],[1,0],[0,1],[-1,1]};
neighbor2list={[-2,1],[-1,2],[1,1],[2,-1],[1,-2],[-1,-1]};
neighbor3list={[-2,0],[-2,2],[0,2],[2,0],[2,-2],[0,-2]};
neighbor4list={[-2,-1],[-1,-2],[1,-3],[2,-3],[3,-2],[3,-1],[2,1],[1,2],[-1,3],[-2,3],[-3,2],[-3,1]};
neighbor5list={[-3,0],[0,-3],[3,-3],[3,0],[0,3],[-3,3]};
% neighborlist={[-1,0]};

% t0= hoppingt(neighbor0list,parameters);
% t1= hoppingt(neighbor1list,parameters);
% t2= hoppingt(neighbor2list,parameters);
% t3= hoppingt(neighbor3list,parameters);
% t4= hoppingt(neighbor4list,parameters);
% t5= hoppingt(neighbor5list,parameters);

% t_all= hoppingt([neighbor0list,neighbor1list,neighbor2list,neighbor3list,neighbor4list,neighbor5list],parameters);

% angle(t*1000)
% [rx,ry]=meshgrid(linspace(-2*sqrt(3)*parameters.aM,2*sqrt(3)*parameters.aM,101));
% [wbgrid,wtgrid]=w([0,0],rx,ry,parameters);
% uu=hubbardU(wbgrid,wtgrid,rx,ry,parameters);
% uu*1000

