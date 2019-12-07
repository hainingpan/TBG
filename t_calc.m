parameters=mainTMD('m',0.35,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',3,'Vz',10);
neighborlist={[-1,0],[0,-1],[1,-1],[1,0],[0,1],[-1,1]};
% neighborlist={[-1,0]};

t= hoppingt(neighborlist,parameters);

t*1000