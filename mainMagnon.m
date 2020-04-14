function param=mainMagnon(varargin)
p=inputParser;
addParameter(p,'phi',3/4*pi);
addParameter(p,'aa',1);
parse(p,varargin{:});
phi=p.Results.phi;
aa=p.Results.aa;

a{1}=aa*[sqrt(3)/2;-1/2];
a{2}=aa*[sqrt(3)/2;1/2];
A{1}=aa*[sqrt(3)/2;3/2];
A{2}=aa*[-sqrt(3)/2;3/2];
b{1}=2*pi/aa*[1/sqrt(3);-1];
b{2}=2*pi/aa*[1/sqrt(3);1];
B{1}=2*pi/aa*[1/sqrt(3);1/3];
B{2}=2*pi/aa*[-1/sqrt(3);1/3];
m=(b{1}+b{2})/2;
M=(B{1}+B{2})/2;
k=(b{1}+2*b{2})/3;
kp=(2*b{1}+b{1})/3;
K=(B{1}+2*B{2})/3;
Kp=(2*B{1}+B{2})/3;

n=cell(1,7);
for i=1:6
    n{i}=[cos((i-1)*(pi/3)),-sin((i-1)*(pi/3));sin((i-1)*(pi/3)),cos((i-1)*(pi/3))]*[0;aa];
end
n{7}=[0;0];% reserved for [0,0]

N=cell(1,7);
for i=1:6
    N{i}=[cos((i-1)*(pi/3)),-sin((i-1)*(pi/3));sin((i-1)*(pi/3)),cos((i-1)*(pi/3))]*[0;3*aa];
end
N{7}=[0;0];% reserved for [0,0]

S=cell(1,3);
theta=2*pi/3*(2*(phi<=pi)-1);
for i=1:3
    S{i}=[cos((i-1)*theta),sin((i-1)*theta)];
end

param=struct('n',{n},'N',{N},'theta',theta,'S',{S},'phi',phi,'aa',aa,'a',a,'A',A,'b',b,'B',B,'m',m,'M',M,'k',k,'kp',kp,'K',K,'Kp',Kp);