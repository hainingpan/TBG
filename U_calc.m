[~,~,L]=size(wbgrid);
for i=1:L
    U(i)=hubbardU_fft(wbgrid(:,:,1),wtgrid(:,:,1),wbgrid(:,:,i),wtgrid(:,:,i),rx,ry,parameters);
end
