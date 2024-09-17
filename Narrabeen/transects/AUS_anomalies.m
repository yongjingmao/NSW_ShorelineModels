% Seasonality/anomalies of AUS beaches. SeanV + Kilian/Kristens's data.
% clear all;
close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add matlab helper routines to the path (the Google Earth .kml toolbox for model output is important here ... but others are needed too)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matlab_helpers_toolbox='K:\COSMOS\CoSMoS_COAST_California\matlab_helper_functions';
addpath(genpath(matlab_helpers_toolbox));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the transects data file (which is effectively the model grid/domain)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('transects','var')  % if the transect file does not exist in the workspace, then load it in ...
    % the transects file to load
    TRANSECTS_FILE='transects.mat'; % path to the transects folder

    fprintf('loading transects file ...');
    load(TRANSECTS_FILE); % load the transects file ... and display progress
    fprintf(' done. \n');
else
    clearvars -except transects t_obs Y_obs SAT Y_rms; % this clears a few variables except some that are time consuming to load in again
end

Ntr=length(transects);

% get the observations at each assimilation time
if ~exist('Y_obs')

    t_obs=[];
    for i=1:Ntr % for all 'full model', 'cross-shore only', and 'rate only' transects
        t_obs=union(t_obs,round(transects(i).t)); % get a union of ALL dates where observational data exist
    end

    t_obs=t_obs'; % transpose this (to make it into a row vector)

    t_obs=sort(t_obs,'ascend');

    % initial conditions
    t0=datenum(1985,3,1);
    tstop=datenum(2021,9,31);

    % remove assimilation times that are before or after the calibration period
    t_obs(t_obs<t0)=[];
    t_obs(t_obs>tstop)=[];

    % preallocate arrays for the observations at each assimilation time (for all transects)
    Y_obs=NaN(Ntr,length(t_obs));
    Y_rms=NaN(Ntr,length(t_obs));
    SAT=zeros(Ntr,length(t_obs));

    y_rms_gps=1;  % approximate RMS error [in meters] of lidar/gps shoreline observations
    y_rms_sat=14; % approximate RMS error [in meters] of the satellite-derived shoreline observations

    fprintf('preparing data for assimilation ... '); % display progress
    for n=1:length(t_obs)   % for all assimilation times
        for i=1:Ntr  % for all transects of interest
            [TF,id1]=ismember(t_obs(n),round(transects(i).t)); % determine if this transect contains data at time t_obs(n)
            if TF                                    % if it does, then save that data
                Y_obs(i,n)=transects(i).Y(id1);  % get the observation at that time
                if transects(i).SAT(id1)==1      % if that data source is from satellites
                    Y_rms(i,n)=y_rms_sat;        % get the satellite rms error at that time
                    SAT(i,n)=1;                  % and specify that that dat comes from a satellite
                else                                 % if not, then specify that the data comes from gps/lidar
                    Y_rms(i,n)=y_rms_gps;        % get the lidar/gps rms error at that time
                end
            end
        end
    end
    fprintf('done.\n');

    % if any assimilation times have no data as a result of the previous
    % removal step, then remove this assimilation time
    id_rm=find((sum(isnan(Y_obs))==Ntr));

    t_obs(id_rm)=[];   % remove elements identified with indices id_rm
    Y_obs(:,id_rm)=[];
    Y_rms(:,id_rm)=[];
    SAT(:,id_rm)=[];
end

Ya=detrend(Y_obs','linear','omitnan')';

Y_smooth=smoothn(Ya,100);

[ID,T]=meshgrid(1:Ntr,t_obs);ID=ID'; T=T';

figure; pcolor(T,ID,Y_smooth); caxis([-50 50]); colormap(jet); shading flat; set(gca,'layer','top','Ydir','rev','Xtick',datenum(1984:2021,1,1)); colorbar; datetick('x','keeplimits','keepticks'); axis tight

stop

Y_mean=nanmean(Ya,1);

Y_smooth=smoothn(Y_mean,100000);

figure; plot(t_obs,Y_mean,'r',t_obs,Y_smooth,'b'); datetick('x');

[YYYY,MM,DD]=datevec(t_obs);

t_cyc=datenum(1984,1:12,15);
Y_cyc=NaN(1,12)

for i=1:12;
    id=find(MM==i);
    Y_cyc(i)=nanmean(Y_smooth(id));
end

figure; plot(t_cyc,Y_cyc,'-bo'); datetick('x');

t_cyc=datenum(reshape(repmat(1984:2021,12,1),[],1)',repmat([1:12],1,length(1984:2021),1),15);
Y_cyc=repmat(Y_cyc',length(1984:2021),1)';

%% make plots
figure; hold on; box on;
plot(t_obs,Y_mean,'-b.');
han1=plot(t_obs,Y_smooth,'-r',t_cyc,Y_cyc,'m'); set(han1,'LineWidth',2); datetick('x','keeplimits','keepticks')

Yc=interp1(t_cyc,Y_cyc,t_obs); datetick('x');

figure; plot(t_cyc,Y_cyc,'b',t_obs,Yc,'r');  datetick('x');