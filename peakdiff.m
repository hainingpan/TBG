% load('dos.mat')
[~,index1]=max(dos1,[],3);
[~,index2]=max(dos2,[],3);

for i=1:length(enlist1)
    for j=1:length(enlist1)
        peakd(i,j)=enlist1(i,j,index1(i,j))-enlist2(i,j,index2(i,j));
    end
end
