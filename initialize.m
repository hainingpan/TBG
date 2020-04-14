function param=initialize(chirality,param)
%chirality=-1 -> phi<pi
%chirality=+1 -> phi>pi
n=param.n;
r=param.r;
N=length(n);
left=param.left;
right=param.right;
sx=zeros(N,1);
sy=zeros(N,1);
sz=zeros(N,1);
difn=mod(n(:,1)-n(:,2),3);
%(y-x)%3==0:(0) ..%3==1:(2pi*/3) ..%3==2:(4*pi/3)
sx(difn==0)=1;
sy(difn==0)=0;
sx(difn==1)=cos(-chirality*2*pi/3);
sy(difn==1)=sin(-chirality*2*pi/3);
sx(difn==2)=cos(-chirality*4*pi/3);
sy(difn==2)=sin(-chirality*4*pi/3);
%Sz=1 for middle region
down=40.5-10;
up=40.5+10;


filter=(((r(:,1)<=left | r(:,1) >= right)) & (r(:,2)>down | r(:,2)<up));
sx(filter)=0;
sy(filter)=0;
sz(filter)=0;

% filter=(r(:,2)>=down & r(:,2)<=up);
rx=8*sqrt(3);
ry=81/2;
radius=8*sqrt(3);
filter=((r(:,1)-rx).^2+(r(:,2)-ry).^2)<=radius^2;

sx(filter)=0;
sy(filter)=0;
sz(filter)=1;

param.s0=[sx;sy;sz];

end

