% make waves data for NSW, Australia. SeanV
%clear all;
close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the transects data file (which is effectively the model grid/domain)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('transects','var')  % if the transect file does not exist in the workspace, then load it in ...
    % the transects file to load
    TRANSECTS_FILE='..\transects\transects.mat'; % path to the transects folder

    fprintf('loading transects file ...');
    load(TRANSECTS_FILE); % load the transects file ... and display progress
    fprintf(' done. \n');
end

Ntr=length(transects);

% waves start and stop time
tstart=datenum(1987,1,1);
%tstop=datenum(1988,1,1);
tstop=datenum(2023,6,1);

% daily wave time series
t=tstart:tstop;
len=length(t);

% initialize full wave data
Hs_max=NaN(Ntr,len);
Hs_mean=NaN(Ntr,len);
Tp=NaN(Ntr,len);
Tm=NaN(Ntr,len);
Dir=NaN(Ntr,len);
Dirp=NaN(Ntr,len);

% wave data directory
wave_data_dir='C:\Users\SeanV\Downloads\NSW_Regional\NSW_Regional\Waves';

drr=dir([wave_data_dir,filesep,'*.nc']); % get a list of the .nc files within the data directory

for j=1:length(drr) % for all wave files ...
    
    % file to load
    wave_data_file=[wave_data_dir,filesep,drr(j).name];

    fprintf('loading wave data from %s file (%d of %d) ... ',wave_data_file,j,length(drr)); % display progress

    %finfo=ncinfo([wave_data_dir,filesep,'wave_transects_byron.nc']);
    hh=ncread(wave_data_file,'time');      % load time        variable from .nc
    ID1=ncread(wave_data_file,'tran_id');  % load transect id variable from .nc
    Hs1=ncread(wave_data_file,'Hs');       % load Hs          variable from .nc
    Tp1=ncread(wave_data_file,'Tp');       % load Tp          variable from .nc
    Tm1=ncread(wave_data_file,'Tm02');     % load Tm          variable from .nc
    Dirp1=ncread(wave_data_file,'Dp');      % load Dp          variable from .nc
    Dir1=ncread(wave_data_file,'Dm');       % load Dm          variable from .nc

    fprintf('done.\n');

    % CONFIRM THAT THE START TIME IS INDEED Jan-1 1979!!!!!
    t1=datenum(1979,1,1)+double(hh)/24;

    Ntr_part=length(ID1); % get the length of the .nc transects

    for i=1:Ntr_part % for all transects in each file ...

        id=find(strcmp({transects.coastsat_tr_name}',ID1(i))); % compare transect name to current transect from .nc file

        fprintf('filling in waves for transect %d (%d of %d) ... ',id,i,Ntr_part); % display progress

        for n=1:len % for all days of the daily time series ...
            
            % get daily waves
            idw=find((t1-t(n))>=0 & (t1-t(n))<1);
            [Hs_max(id,n),idmax]=max(Hs1(idw,i));
            Hs_mean(id,n)=mean(Hs1(idw,i));
            Tp(id,n)=Tp1(idw(idmax),i);
            Tm(id,n)=Tm1(idw(idmax),i);
            Dir(id,n)=Dir1(idw(idmax),i);
            Dirp(id,n)=Dirp1(idw(idmax),i);
        end

        fprintf('done.\n');

    end

end

% save
save('waves.mat',"t","Hs_mean","Hs_max","Tp","Tm","Dir","Dirp");

[IDD,T]=meshgrid(1:Ntr,t); IDD=IDD'; T=T';

stridex=2;
stridey=2;

figure; pcolor(T(1:stridex:end,1:stridey:end),IDD(1:stridex:end,1:stridey:end),Hs_max(1:stridex:end,1:stridey:end)); caxis([0 5]); colormap(jet); shading flat; set(gca,'layer','top','Ydir','rev','Xtick',datenum(1979:2023,1,1)); colorbar; datetick('x','keeplimits','keepticks'); axis tight

figure; pcolor(T(1:stridex:end,1:stridey:end),IDD(1:stridex:end,1:stridey:end),Hs_mean(1:stridex:end,1:stridey:end)); caxis([0 5]); colormap(jet); shading flat; set(gca,'layer','top','Ydir','rev','Xtick',datenum(1979:2023,1,1)); colorbar; datetick('x','keeplimits','keepticks'); axis tight

figure; pcolor(T(1:stridex:end,1:stridey:end),IDD(1:stridex:end,1:stridey:end),Tp(1:stridex:end,1:stridey:end)); caxis([0 20]); colormap(jet); shading flat; set(gca,'layer','top','Ydir','rev','Xtick',datenum(1979:2023,1,1)); colorbar; datetick('x','keeplimits','keepticks'); axis tight

figure; pcolor(T(1:stridex:end,1:stridey:end),IDD(1:stridex:end,1:stridey:end),Tm(1:stridex:end,1:stridey:end)); caxis([0 20]); colormap(jet); shading flat; set(gca,'layer','top','Ydir','rev','Xtick',datenum(1979:2023,1,1)); colorbar; datetick('x','keeplimits','keepticks'); axis tight

figure; pcolor(T(1:stridex:end,1:stridey:end),IDD(1:stridex:end,1:stridey:end),Dir(1:stridex:end,1:stridey:end)); caxis([0 360]); colormap(jet); shading flat; set(gca,'layer','top','Ydir','rev','Xtick',datenum(1979:2023,1,1)); colorbar; datetick('x','keeplimits','keepticks'); axis tight

figure; pcolor(T(1:stridex:end,1:stridey:end),IDD(1:stridex:end,1:stridey:end),Dirp(1:stridex:end,1:stridey:end)); caxis([0 360]); colormap(jet); shading flat; set(gca,'layer','top','Ydir','rev','Xtick',datenum(1979:2023,1,1)); colorbar; datetick('x','keeplimits','keepticks'); axis tight


size(Dirp)
