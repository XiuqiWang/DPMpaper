function [Cz, Upz, Qz] = getZProfile(data, dz)
Cz = double.empty;
Upz = double.empty;
Qz = double.empty;

[r,c] = size(data.z);
VolumeFraction = data.VolumeFraction;
VelocityX = data.VelocityX;
Density = data.Density;
indexz0 = find(VolumeFraction<=0.4,1);
disp(indexz0);
for t=1:c
for i=1:r
     Cz(i,t) = Density(i,t)*dz;
     Upz(i,t) = VelocityX(i,t);
     Qz(i,t) = Cz(i,t)*Upz(i,t);
end
end
end