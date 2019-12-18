% mm=m(:,:,end);
mm=m0;
% figure;
% hold on;
for i=1:3
p1 = [0,0];                  % First Point
p2 = mm(:,i);                       % Second Point
dp = p2-p1;                         % Difference
quiver(p1(1),p1(2),dp(1),dp(2),0);
end
daspect([1,1,1]);
pbaspect([1,1,1]);
legend('A','B','C')