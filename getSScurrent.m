function [u,v]=getSScurrent (p, settings)

% getSScurrent
% -------------
%
% reads in nc file of current from p.lon, p.lat and settings.date
% returns u,v vectors of length np of sea surface current

dateVec = datevec(settings.date);
y = dateVec(1);
m = dateVec(2);

ufile = [settings.SScurrentPath 'u_' num2str(y) '_' num2str(m) '.nc'];
vfile = [settings.SScurrentPath 'v_' num2str(y) '_' num2str(m) '.nc'];

ncidu = netcdf.open(ufile,'NOWRITE') ;
    time = netcdf.getVar(ncidu,0);
    time = datenum(2000,12,31,0,0,0) + time; % hycom time convention
    lon  = netcdf.getVar(ncidu,1);
    lat  = netcdf.getVar(ncidu,2);
netcdf.close(ncidu)

pLon = p.lon;
pLat = p.lat;


if min(lon)<0 % special case if longitude is referenced -180 to 180
    pLon(pLon>=180) = pLon(pLon>=180) - 360;
end

i = getIndex(pLon,lon);
j = getIndex(pLat,lat);
t = getIndex(settings.date,time) -1; % netcdf index starts at 0

is  =mod(i-1,length(lon))+1;
ie  =mod(i+1,length(lon))+1;
js  =mod(j-1,length(lat))+1;
je  =mod(j+1,length(lat))+1;
tdt =min(t+1,length(time) -1); % netcdf index starts at 0

uij=zeros(1,p.np);
vij=zeros(1,p.np);

uis=zeros(1,p.np);
vis=zeros(1,p.np);
uie=zeros(1,p.np);
vie=zeros(1,p.np);
ujs=zeros(1,p.np);
vjs=zeros(1,p.np);
uje=zeros(1,p.np);
vje=zeros(1,p.np);
udt=zeros(1,p.np);
vdt=zeros(1,p.np);
dx=zeros(1,p.np);
dy=zeros(1,p.np);

ncidu = netcdf.open(ufile,'NOWRITE') ;
ncidv = netcdf.open(vfile,'NOWRITE') ;

U   = netcdf.getVar( ncidu , 3 , [0,0,t] , [length(lon),length(lat),1] );
V   = netcdf.getVar( ncidv , 3 , [0,0,t] , [length(lon),length(lat),1] );
Udt = netcdf.getVar( ncidu , 3 , [0,0,tdt] , [length(lon),length(lat),1] );
Vdt = netcdf.getVar( ncidv , 3 , [0,0,tdt] , [length(lon),length(lat),1] );

netcdf.close(ncidu)
netcdf.close(ncidv)

% current scaling factor in Hycom files
U = double(U) /1000;
V = double(V) /1000;

% eliminate bad values
U ( U==-30 ) = 0;
V ( V==-30 ) = 0;

for k=1:p.np
    
%     uij(k) = netcdf.getVar( ncidu , 3 , [i(k),j(k),t] , [1,1,1] );
%     vij(k) = netcdf.getVar( ncidv , 3 , [i(k),j(k),t] , [1,1,1] );
%     
%     uis(k) = netcdf.getVar( ncidu , 3 , [is(k),j(k),t] , [1,1,1] );
%     vis(k) = netcdf.getVar( ncidv , 3 , [is(k),j(k),t] , [1,1,1] );
%     uie(k) = netcdf.getVar( ncidu , 3 , [ie(k),j(k),t] , [1,1,1] );
%     vie(k) = netcdf.getVar( ncidv , 3 , [ie(k),j(k),t] , [1,1,1] );
%     ujs(k) = netcdf.getVar( ncidu , 3 , [i(k),js(k),t] , [1,1,1] );
%     vjs(k) = netcdf.getVar( ncidv , 3 , [i(k),js(k),t] , [1,1,1] );
%     uje(k) = netcdf.getVar( ncidu , 3 , [i(k),je(k),t] , [1,1,1] );
%     vje(k) = netcdf.getVar( ncidv , 3 , [i(k),je(k),t] , [1,1,1] );
%     
%     udt(k) = netcdf.getVar( ncidu , 3 , [i(k),j(k),tdt] , [1,1,1] );
%     vdt(k) = netcdf.getVar( ncidv , 3 , [i(k),j(k),tdt] , [1,1,1] );
    
    
    uij(k) = U( i(k)  ,j(k) );
    vij(k) = V( i(k)  ,j(k) );
    
    uis(k) = U( is(k) ,j(k) );
    vis(k) = V( is(k) ,j(k) );
    uie(k) = U( ie(k) ,j(k) );
    vie(k) = V( ie(k) ,j(k) );
    ujs(k) = U( i(k)  ,js(k) );
    vjs(k) = V( i(k)  ,js(k) );
    uje(k) = U( i(k)  ,je(k) );
    vje(k) = V( i(k)  ,je(k) );

    udt(k) = Udt( i(k)  ,j(k) );
    vdt(k) = Vdt( i(k)  ,j(k) );
    
    
    % !! need to presave distances !
%   [dx(k) dx0]=lldistkm([ lat(j(k)) lon(is(k)) ]',[ lat(j(k)) lon(ie(k)) ]'); dx=dx * 1000; % in m
% 	[dy(k) dy0]=lldistkm([ lat(js(k)) lon(i(k)) ]',[ lat(je(k)) lon(i(k)) ]'); dy=dy * 1000; % in m
    dx=1;dy=1;

end


% calculate gradients

dt      = abs(time(tdt+1)-time(t+1)) * 24 *3600;


ux = (uie - uis) / dx;
uy = (uje - ujs) / dy;
vx = (vie - vis) / dx;
vy = (vje - vjs) / dy;
ut = (udt - uij) / dt;
vt = (vdt - vij) / dt;

% second order accurate advection scheme

u_ = uij + (dt*ut)/2;
v_ = vij + (dt*vt)/2;

u = ( u_ + ( uy.*v_ - vy.*u_ )*dt/2 ) ./ ( (1-ux*dt/2).*(1-vy*dt/2) - (uy.*vx*dt^2)/4 ) ;
v = ( v_ + ( vx.*u_ - ux.*v_ )*dt/2 ) ./ ( (1-uy*dt/2).*(1-vx*dt/2) - (ux.*vy*dt^2)/4 ) ;

u = uij;
v = vij;

if max(sqrt(u.^2+v.^2))>10 
    disp('velocities > 10m/s')
end


