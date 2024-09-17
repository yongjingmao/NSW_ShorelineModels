% make waves data for NSW, Australia. SeanV
clear all;
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
tstart=datenum(1980,1,1);
%tstop=datenum(1988,1,1);
tstop=datenum(2024,1,1);

% daily wave time series
t=tstart:tstop;
len=length(t);

% Beach state
d50 = 0.4; % d50 in mm
depth = 10; % depth of waves
g = 9.81; % Gravity acceleration

% Define phi range for ShoreFor calculation
phi_max = 2000; %days
phi_min = 5; %days


% initialize full wave data
Hs_max=NaN(Ntr,len);
Hs_mean=NaN(Ntr,len);
Tp=NaN(Ntr,len);
Tm=NaN(Ntr,len);
Dir=NaN(Ntr,len);
Dirp=NaN(Ntr,len);

% wave data directory
wave_data_dir='C:\Users\z3541792\OneDrive - UNSW\General - ARC FT220100009 - Splinter Shoreline Modelling\Shoreline Modelling\DATA\NSW_Regional\Waves_updated';

drr=dir([wave_data_dir,filesep,'*syd.nc']); % get a list of the .nc files within the data directory
%%
% Create the ranges and concatenate them
phis1 = phi_min:5:95;   % Equivalent to np.arange(phi_min,100,5)
phis2 = 100:10:490;    % Equivalent to np.arange(100,500,10)
phis3 = 500:20:phi_max; % Equivalent to np.arange(500,phi_max,20)

% Concatenate the ranges
phis = [phis1, phis2, phis3];

% Compute the logarithm base 10
phi_log = log10(phis);

% Compute fall velocity
w = fallvelocity(d50/1000, 15);

for j=1:length(drr) % for all wave files ...
    
    % file to load
    wave_data_file=[wave_data_dir,filesep,drr(j).name];

    fprintf('loading wave data from %s file (%d of %d) ... ',wave_data_file,j,length(drr)); % display progress

    %finfo=ncinfo([wave_data_dir,filesep,'wave_transects_byron.nc']);
    hh=ncread(wave_data_file,'time');      % load time        variable from .nc


    ID=ncread(wave_data_file,'tran_id');  % load transect id variable from .nc
    idx_in=find(ismember(ID, {transects.coastsat_tr_name})); % intersect transect id
    idx_range = min(idx_in):max(idx_in);
    ID1 = ID(idx_range);

    Hs1=ncread(wave_data_file,'Hs', [1, min(idx_in)], [inf, max(idx_in)-min(idx_in)+1]);  % load Hs          variable from .nc
    Tp1=ncread(wave_data_file,'Tp', [1, min(idx_in)], [inf, max(idx_in)-min(idx_in)+1]);     % load Tp          variable from .nc
    Tm1=ncread(wave_data_file,'Tm02', [1, min(idx_in)], [inf, max(idx_in)-min(idx_in)+1]);     % load Tm          variable from .nc
    Dirp1=ncread(wave_data_file,'Dp', [1, min(idx_in)], [inf, max(idx_in)-min(idx_in)+1]);      % load Dp          variable from .nc
    Dir1=ncread(wave_data_file,'Dm', [1, min(idx_in)], [inf, max(idx_in)-min(idx_in)]+1);       % load Dm          variable from .nc

    fprintf('done.\n');


    Ntr_part=length(ID1); % get the length of the .nc transects
    % CONFIRM THAT THE START TIME IS INDEED Jan-1 1979!!!!!
    t1=datenum(1979,1,1)+double(hh)/24;
    day1=unique(ceil(t1));
    dt = 1*24*3600; % dt=1day for wave data (in seconds)

    
    Hs = Hs1(:, ismember(ID1, {transects.coastsat_tr_name}));
    Tp = Tp1(:, ismember(ID1, {transects.coastsat_tr_name}));
    Tm = Tm1(:, ismember(ID1, {transects.coastsat_tr_name}));
    Dirp = Dirp1(:, ismember(ID1, {transects.coastsat_tr_name}));
    Dir = Dir1(:, ismember(ID1, {transects.coastsat_tr_name}));

    Hs = array2timetable(Hs,'RowTimes', datetime(t1,'ConvertFrom','datenum'));
    Hs_max = table2array(retime(Hs, 'daily', 'max'));
    Hs_mean = table2array(retime(Hs, 'daily', 'mean'));

    Tp = array2timetable(Tp,'RowTimes', datetime(t1,'ConvertFrom','datenum'));
    Tp = table2array(retime(Tp, 'daily', 'mean'));

    Tm = array2timetable(Tm,'RowTimes', datetime(t1,'ConvertFrom','datenum'));
    Tm = table2array(retime(Tm, 'daily', 'mean'));

    Dirp = array2timetable(Dirp,'RowTimes', datetime(t1,'ConvertFrom','datenum'));
    Dirp = table2array(retime(Dirp, 'daily', 'mean'));

    Dir = array2timetable(Dir,'RowTimes', datetime(t1,'ConvertFrom','datenum'));
    Dir = table2array(retime(Dir, 'daily', 'mean'));

    % for i=1:Ntr_part % for all transects in each file ...
    % 
    %     id=find(strcmp({transects.coastsat_tr_name}',ID1(i))); % compare transect name to current transect from .nc file
    % 
    %     if ~isempty(id)
    %         fprintf('filling in waves for transect %d (%d of %d) ... ',id,i,Ntr_part); % display progress
    % 
    %         for n=1:len % for all days of the daily time series ...
    % 
    %             % get daily waves
    %             idw=find((t1-t(n))>=0 & (t1-t(n))<1);
    %             [Hs_max(id,n),idmax]=max(Hs1(idw,i));
    %             Hs_mean(id,n)=mean(Hs1(idw,i));
    %             Tp(id,n)=Tp1(idw(idmax),i);
    %             Tm(id,n)=Tm1(idw(idmax),i);
    %             Dir(id,n)=Dir1(idw(idmax),i);
    %             Dirp(id,n)=Dirp1(idw(idmax),i);
    % 
    % 
    %         end
    % 
    %         fprintf('done.\n');
    %     end
    % 
    % 
    % end

    % Calculate offshore wave height
    y = 4.03*depth./(Tp.^2);
    kd2 = y.^2+y./(1+(0.666*y)+(0.355*y.^2)+(0.161*y.^3)+(0.0632*y.^4)+(0.0218*y.^5)+(0.00564*y.^6));
    kh = sqrt(kd2);
    Cg = (g.*Tp./(2 * pi)).*(tanh(kh)).*(0.5*(1+2*kh./sinh(2*kh)));
    Cgo = (1/4)*g*Tp/pi;
    Ks = sqrt(Cgo./Cg);
    Hs0 = Hs_max./Ks;

    % Calculate breaking wave height
    Hsb = 0.39*g^(1/5)*(Tp.*Hs0.^2).^(2/5);

    % Calculate diemnsionless fall velocity omega
    omega = Hsb./Tp./w;
    omega = fillmissing(omega, 'linear', 1);

    % Calculate wave power in shallow water depth
    P=Hsb.^2.5;

    % Calculate equilibrium omega
    % Fs = struct;
    % for i=1:length(phis)
    %     phi = phis(i);
    %     D = 2*phi;
    %     omega_eq = WS85FilterConv(omega, D, phi, dt);
    %     omega_delta = omega_eq-omega;
    %     F = omega_delta/std(omega_delta, "omitmissing").*P.^0.5;
    %     Fs.(['phi',num2str(phi)]) = F((day1>=tstart)&(day1<=tstop), :)';
    % end
    Fs = NaN(length(phis), Ntr, len); % Pre-calculate F for different sets of phi values
    for i=1:length(phis)
        phi = phis(i);
        D = 2*phi;
        omega_eq = WS85FilterConv(omega, D, phi, dt);
        omega_delta = omega_eq-omega;
        F = omega_delta/std(omega_delta, "omitmissing").*P.^0.5;
        Fs(i, :, :) = F((day1>=tstart)&(day1<=tstop), :)';
    end

    Hs_mean = Hs_mean((day1>=tstart)&(day1<=tstop), :)';
    Hs_max = Hs_max((day1>=tstart)&(day1<=tstop), :)';
    Hsb = Hsb((day1>=tstart)&(day1<=tstop), :)';
    omega = omega((day1>=tstart)&(day1<=tstop), :)';
    Tp = Tp((day1>=tstart)&(day1<=tstop), :)';
    Tm = Tm((day1>=tstart)&(day1<=tstop), :)';
    Dir = Dir((day1>=tstart)&(day1<=tstop), :)';
    Dirp = Dirp((day1>=tstart)&(day1<=tstop), :)';

end

% save
save('waves.mat',"t","Hs_mean","Hs_max","Tp","Tm","Dir","Dirp", "Fs", "Hsb", "omega");

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
