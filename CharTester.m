function [Q, V, Cvx, H] = CharTester(data, index, D, dz)
Q = double.empty;
V = double.empty;
Cvx = double.empty;
H = double.empty;
%%calculate transport rate, saltation velocity and airborne concentration
for i=1:size(data,1)
    Q(i,:) = getDimensionlessMassFlux(data(i),dz,D);
    V(i,:) = getSaltationVelocity(data(i), index, D);
    Cvx(i,:) = getAirborneConcentrationByVx(data(i),index, dz);
    H(i,:) = getH50(data(i), index, D);
end
end
 
function Q = getDimensionlessMassFlux(data, dz, D)
M = 0;
Q = double.empty;

[r,c] = size(data.z);
Density = data.Density;
VelocityX = data.VelocityX;
for t=1:c
for i=1:r
     M = M + Density(i,t)*VelocityX(i,t)*dz;   
end
Q(t) = M/sqrt((2650/1.225-1)*9.81*D^3)/2650;%dimensionless
M = 0;
end
end

function Cvx = getAirborneConcentrationByVx(data, index, dz)
M = 0;
Cvx = double.empty;
[r,c] = size(data.z);
Density = data.Density;
VelocityX = data.VelocityX;
ustar = sqrt(index*0.01*(2650-1.225)*9.81*dz/1.225);
ux_c = 0.1*ustar;
for t=1:c
    indexz0 = find(max(VelocityX(:,t)-ux_c,0),1);
    %disp(indexz0);
for i=indexz0:r
     M = M + Density(i,t)*dz;
end
Cvx(t) = M;%[kg/m^2]
M = 0;
end
end

function Vx = getSaltationVelocity(data, index, D)
[r,c] = size(data.z);
VelocityX=data.VelocityX;
VolumeFraction=data.VolumeFraction;
Vx=double.empty;

ustar = sqrt(index*0.01*(2650-1.225)*9.81*D/1.225);
ux_c = 0.1*ustar;

for t=1:c
    %indexz0 = find(min(VolumeFraction(:,t)-0.01,0),1);
    indexz0 = find(max(VelocityX(:,t)-ux_c,0),1);
    Cvel=0;Cvolume=0;
for i=indexz0:r
    if VelocityX(i,t)~=0
        Cvel = Cvel+VolumeFraction(i,t)*VelocityX(i,t);
        Cvolume = Cvolume+VolumeFraction(i,t);
    end
end
if Cvel == 0
    Vx(t)=0;
else
    Vx(t)=Cvel/Cvolume;%/sqrt(9.81*D*(2650/1.225-1));%non-dimensionalized velocity
end
end
end

function H50 = getH50(data, index, D)
[r,c] = size(data.z);
VolumeFraction=data.VolumeFraction;
VelocityX = data.VelocityX;
z=data.z;
H50=double.empty;

ustar = sqrt(index*0.01*(2650-1.225)*9.81*D/1.225);
ux_c = 0.1*ustar;

for t=1:c
    indexz0 = find(max(VelocityX(:,t)-ux_c,0),1);
    display(['indexz0=',num2str(indexz0)]);
    
    if mod(indexz0,1) == 0%判断Indexz0是不是整数(有无找到indexz0)
        cumulative_sum = cumsum(VolumeFraction(indexz0:end,t));
        % 计算总和
        total_sum = cumulative_sum(end);
        % 计算每个元素在总和中所占的比例
        percentages = cumulative_sum / total_sum;
        index50=find(max(percentages-0.5,0),1);
        %display(['index50=',num2str(index50)]);
        H50(t)=z(index50+indexz0-1)-z(indexz0-1);%/D;%non-dimensionalize
    else
        H50(t)=0;
    end
end
end
