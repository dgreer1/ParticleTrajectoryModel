function p = loadParticles(settings)


fileInfo=ncinfo(settings.SourceFilename);

if strcmp(fileInfo.Variables(1).Name,'time') % hotstart
    ncid = netcdf.open(settings.SourceFilename,'NOWRITE');
    time = netcdf.getVar(ncid,0);
    t = getIndex(settings.date,time) - 1; % NETCDF index starts at 0
    
    p.np = fileInfo.Dimensions(1).Length;
    p.id = netcdf.getVar(ncid,3);
    p.lon = netcdf.getVar(ncid,1,[t,0], [1,p.np]);
    p.lat = netcdf.getVar(ncid,2,[t,0], [1,p.np]);
    p.releaseDate = netcdf.getVar(ncid,4);
    p.UNSD = netcdf.getVar(ncid,5);
    
    netcdf.close(ncid)
    
else % source file
    ncid = netcdf.open(settings.SourceFilename,'NOWRITE');

    varid = netcdf.inqVarID(ncid,'id');
    p.id  = netcdf.getVar(ncid,varid)';
    varid = netcdf.inqVarID(ncid,'lon');
    p.lon = netcdf.getVar(ncid,varid)';
    varid = netcdf.inqVarID(ncid,'lat');
    p.lat = netcdf.getVar(ncid,varid)';
    varid = netcdf.inqVarID(ncid,'releaseDate');
    p.releaseDate = netcdf.getVar(ncid,varid)';
    varid = netcdf.inqVarID(ncid,'unsd');
    p.UNSD = netcdf.getVar(ncid,varid)';

    p.np = length(p.id);

    netcdf.close(ncid)
end


% 
% p.id(p.UNSD~=76)=[];
% p.lon(p.UNSD~=76)=[];
% p.lat(p.UNSD~=76)=[];
% p.releaseDate(p.UNSD~=76)=[];
% p.UNSD(p.UNSD~=76)=[];
% 
% p.id(p.lat<-10 | p.lat>-5)=[];
% p.lon(p.lat<-10 | p.lat>-5)=[];
% p.releaseDate(p.lat<-10 | p.lat>-5)=[];
% p.UNSD(p.lat<-10 | p.lat>-5)=[];
% p.lat(p.lat<-10 | p.lat>-5)=[];

% p.releaseDate=p.releaseDate*0;
% p.np=length(p.id);


p.LON=zeros( length(settings.outputDateList) ,p.np);
p.LAT=zeros( length(settings.outputDateList) ,p.np);



