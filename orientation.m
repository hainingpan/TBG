function [theta,phi]=orientation(m)
theta=zeros(1,3);
phi=zeros(1,3);
for i=1:3
    theta(i)=angle(m(1,i)+m(2,i)*1i)/pi*180;
    phi(i)=angle(m(3,i)+sqrt(m(1,i)^2+m(2,i)^2)*1i)/pi*180;
end
end