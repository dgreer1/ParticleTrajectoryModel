function p = updateParticles(p,dx,dy,settings)

% updateParticles
% -------------
%
% update particle position
% check for land and grid boundary
% returns p structure

dlon = m2lon(dx,p.lat);
dlat = m2lat(dy,p.lat);

lon_new = p.lon + dlon;
lat_new = p.lat + dlat;

% Greenwich meridian
lon_new = mod(lon_new,360);

% Check for coastlines 
id = getIndex(lon_new,settings.landmass.lon);
jd = getIndex(lat_new,settings.landmass.lat);

land  = zeros(1,p.np);

for k=1:p.np 
    land(k) = settings.landmass.data(id(k)+1,jd(k)+1);
end

lat_new(land==1) = p.lat(land==1);
lon_new(land==1) = p.lon(land==1);

p.lon = lon_new;
p.lat = lat_new;