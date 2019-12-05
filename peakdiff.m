% load('dos.mat')
[~,int1]=max(dos1,[],3);
[~,int2]=max(dos2,[],3);

for i=1:length(enlist1)
    for j=1:length(enlist1)
        peakd(i,j)=enlist1(i,j,int1(i,j))-enlist2(i,j,int2(i,j));
    end
end
% peakd=enlist1(int1)-enlist2(int2);