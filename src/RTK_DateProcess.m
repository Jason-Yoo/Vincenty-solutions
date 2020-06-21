clc 
clear
data = xlsread('date/RTK_Date.xlsx');
B=data(:,1);
L=data(:,2);
%H=data(:,3);
distance_True = data(:,4);
i=size(B);
date_out(i(1,1),3) = 0;
distance_Error(2,i(1,1)) = 0;
for n=1:i(1,1)
    
    point1 = [B(1,1) L(1,1)]; 
    point2 = [B(n,1) L(n,1)];
    [distance,Alpha1,Alpha2] = vincenty_inverse(point1,point2);
    date_out(n,1:3) = [distance,Alpha1,Alpha2];
    distance_Error(1,n) = (distance - distance_True(n));
    distance_Error(2,n) = abs(distance_Error(1,n))/distance_True(n)*100;
end
x = distance_True(:,1)';
y1 = distance_Error(1,:);
y2 = distance_Error(2,:);
y2(1,1) = 0;
plot(x,y1,'o-',x,y2);
legend('实际误差','误差百分比')

