epoch=1;
% figure;scatter(param.r(:,1),param.r(:,2),[],s(epoch,1:length(s)/3),'.');colorbar;title('s_x');
% figure;scatter(param.r(:,1),param.r(:,2),[],s(epoch,length(s)/3+1:2*length(s)/3),'.');colorbar;title('s_y');
figure;scatter(param.r(:,1),param.r(:,2),[],s(epoch,2*length(s)/3+1:3*length(s)/3),'.');colorbar;title('s_z');
axis tight
daspect([1,1,1]);