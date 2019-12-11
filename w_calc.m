parameters=mainTMD('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',4,'Vz',10);

neighbor0list={[0,0]};
neighbor1list={[-1,0],[0,-1],[1,-1],[1,0],[0,1],[-1,1]};
neighbor2list={[-2,1],[-1,2],[1,1],[2,-1],[1,-2],[-1,-1]};
neighbor3list={[-2,0],[-2,2],[0,2],[2,0],[2,-2],[0,-2]};
% neighbor4list={[-2,-1],[-1,-2],[1,-3],[2,-3],[3,-2],[3,-1],[2,1],[1,2],[-1,3],[-2,3],[-3,2],[-3,1]};
% neighbor5list={[-3,0],[0,-3],[3,-3],[3,0],[0,3],[-3,3]};

[wbgrid,wtgrid]=w_rec([neighbor0list,neighbor1list,neighbor2list,neighbor3list],rx,ry,parameters);

