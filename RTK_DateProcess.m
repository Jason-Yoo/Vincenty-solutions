clc 
clear
data = xlsread('RTK_Date.xlsx');
B=data(:,1);
L=data(:,2);
%H=data(:,3);
i=size(B);
date_out(i(1,1),3) = 0;
for n=1:i(1,1)
    
    point1 = [B(2,1) L(2,1)]; 
    point2 = [B(n,1) L(n,1)];
    [distance,Alpha1,Alpha2] = vincenty_inverse(point1,point2);
    date_out(n,1:3) = [distance,Alpha1,Alpha2];
    
end
