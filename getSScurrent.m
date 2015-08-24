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
varid = netcdf.inqVarID(ncidu,'time');
time = netcdf.getVar(ncidu,varid);
time = datenum(2000,12,31,0,0,0) + time; % hycom time convention
varid = netcdf.inqVarID(ncidu,'lon');
lon  = netcdf.getVar(ncidu,varid);
varid = netcdf.inqVarID(ncidu,'lat');
lat  = netcdf.getVar(ncidu,varid);
netcdf.close(ncidu)

pLon = p.lon;
pLat = p.lat;


if min(lon)<0 % special case if longitude is referenced -180 to 180
    pLon(pLon>=180) = pLon(pLon>=180) - 360;
end

i = getIndex(pLon,lon);
j = getIndex(pLat,lat);
t = getIndex(settings.date,time) -1; % netcdf index starts at 0

is  =mod(i-1-1,length(lon))+1;
ie  =mod(i+1-1,length(lon))+1;
js  =mod(j-1-1,length(lat))+1;
je  =mod(j+1-1,length(lat))+1;
tdt =min(t+1,length(time) -1); % netcdf index starts at 0
if tdt==t; tdt = t-1; end % backward scheme for last time step

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

varidu = netcdf.inqVarID(ncidu,'water_u');
varidv = netcdf.inqVarID(ncidv,'water_v');
U   = netcdf.getVar( ncidu , varidu , [0,0,t] , [length(lon),length(lat),1] );
V   = netcdf.getVar( ncidv , varidv , [0,0,t] , [length(lon),length(lat),1] );
Udt = netcdf.getVar( ncidu , varidu , [0,0,tdt] , [length(lon),length(lat),1] );
Vdt = netcdf.getVar( ncidv , varidv , [0,0,tdt] , [length(lon),length(lat),1] );

netcdf.close(ncidu)
netcdf.close(ncidv)

% current scaling factor in Hycom files
U = double(U) /1000;
V = double(V) /1000;
Udt = double(Udt) /1000;
Vdt = double(Vdt) /1000;

% eliminate bad values
U ( U==-30 ) = 0;
V ( V==-30 ) = 0;
Udt ( Udt==-30 ) = 0;
Vdt ( Vdt==-30 ) = 0;

for k=1:p.np
    
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
    
    dx (k) = settings.grid.dx( i(k)  ,j(k) );
    dy (k) = settings.grid.dy( i(k)  ,j(k) );

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

if max(sqrt(u.^2+v.^2))>10 
    disp('velocities > 10m/s')
end


