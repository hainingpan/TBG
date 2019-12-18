function theta=orientation(m)
theta=zeros(1,3);
for i=1:3
    theta(i)=angle(m(1,i)+m(2,i)*1i)/pi*180;
end
end