% parameters=mainTMD('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',4,'Vz',0);

neighborlist{1}={[0,0]};
neighborlist{2}={[-1,0],[0,-1],[1,-1],[1,0],[0,1],[-1,1]}; % direction on clock: 12, 10, 8, 6, 4, 2
neighborlist{3}={[-1,-1],[1,-2],[2,-1],[1,1],[-1,2],[-2,1]}; %directin on clock: 11, 9, 7, 5, 4, 1
neighborlist{4}={[-2,0],[0,-2],[2,-2],[2,0],[0,2],[-2,2]}; %direction on clock: 12, 10, 8, 6, 4, 2
neighborlist{5}={[-2,-1],[-1,-2],[1,-3],[2,-3],[3,-2],[3,-1],[2,1],[1,2],[-1,3],[-2,3],[-3,2],[-3,1]}; %counterclockwise from the first point in Quadrant II
neighborlist{6}={[-3,0],[0,-3],[3,-3],[3,0],[0,3],[-3,3]}; %direction on clock: 12, 10, 8, 6, 4, 2

[rx,ry]=meshgrid(linspace(-3*sqrt(3)*parameters.aM,3*sqrt(3)*parameters.aM,101));
% neighbor0list={[0,0]};
% neighbor1list={[-1,0],[0,-1],[1,-1],[1,0],[0,1],[-1,1]};
% neighbor2list={[-2,1],[-1,2],[1,1],[2,-1],[1,-2],[-1,-1]};
% neighbor3list={[-2,0],[-2,2],[0,2],[2,0],[2,-2],[0,-2]};
% neighbor4list={[-2,-1],[-1,-2],[1,-3],[2,-3],[3,-2],[3,-1],[2,1],[1,2],[-1,3],[-2,3],[-3,2],[-3,1]};
% neighbor5list={[-3,0],[0,-3],[3,-3],[3,0],[0,3],[-3,3]};

[wbgrid,wtgrid]=w_rec([neighborlist{1:2}],rx,ry,parameters);
clear wb wt
counter=1;
for i=1:2
    for j=1:length(neighborlist{i})
        wb{i,j}=wbgrid(:,:,counter);
        wt{i,j}=wtgrid(:,:,counter);
        counter=counter+1;
    end
end

