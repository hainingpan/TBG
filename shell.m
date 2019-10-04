function [centerlist,centercor]=shell(NN,rotmat)
%NN is defined as number of shell, total Moire unit cell within shell NN is
%3NN^2+3*NN+1
%Here unit is set to (-n..n)->(-1,1)
% aM1=rotmat*parameters.aM1';
% aM2=rotmat*parameters.aM2';
aM1=rotmat*[3/2;sqrt(3)/2];
aM2=rotmat*[-3/2;sqrt(3)/2];

centerlist=zeros(2,3*NN^2+3*NN+1);
counter=1;
for i=-NN:NN
    for j=max(-NN,i-NN):min(NN,i+NN)
        centerlist(:,counter)=[i;j];
        counter=counter+1;
    end
end
centercor=zeros(2,3*NN^2+3*NN+1);
counter=1;
for i=1:(3*NN^2+3*NN+1)
    centercor(:,counter)=[aM1,aM2]*centerlist(:,i);
    counter=counter+1;
end
end
