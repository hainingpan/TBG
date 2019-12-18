% [~,~,L]=size(wbgrid);
clear U
for i=1:2
    for j=1:length(neighborlist{i})
        U{i}(j)=hubbardU_fft(wb{1,1},wt{1,1},wb{i,j},wt{i,j},rx,ry,parameters);
    end
end
