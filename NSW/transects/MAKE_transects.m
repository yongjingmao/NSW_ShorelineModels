% Load in Killian's CoastSat data for Australia and make CoSMoS-COAST transects. SeanV.
close all; clc;
clearvars -except transects;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Add matlab helper routines toolbox to the path (the kml toolobx is the most important here ... but others might be too)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matlab_helpers_toolbox='K:\COSMOS\Matlab';
addpath(genpath(matlab_helpers_toolbox));

warning('off','all'); % turn off warnings, since I'm getting a wierd one with the matlab tables

% UTM zones for NSW Australia
UTM_ZONE=56;
UTM_HEM='S';

% initialize transects data structure (preallocate space)
transects=struct;

% use tide-corrected satellite shorelines (for raw shorelines set =0)
CORRECTED_SHORELINES=1;

% location of Sat-derived shoreline data folder
DATA_DIR='K:\COSMOS\Aus\NSW\shorelines\coastsat\csv_v6';

drr=dir(DATA_DIR); % get a list of the folders with data
drr(ismember({drr.name},{'.','..'}))=[]; % remove '.' and '..' from the list

all_sites={drr.name}';
all_sites1=strrep(all_sites,'s0','s'); % the spreadsheet name is missing a 0 from the site name, so I make a new variable to use when refering to the spreadsheet with the transect/littoral cell names

% the total number of sites
Nsites=length(drr);

% load in transect information table
T=readtable("TransectNames_NSW_KS_SV.xls");
sites=T.ID;
rank=T.RANK;
name=T.Name;
reorder=T.FLIPN_S; reorder=(~isnan(reorder));
remove=T.RemoveTransects; remove=(~isnan(remove));
id2rm=T.ID2Remove;

% transects counter variable
count=1;

for i=1:length(sites)  % for all data folders

    id1=find(rank==i); % get the row that is ranked in desired order 

    id=find(strcmp(sites{id1},all_sites1)); % find the name of the corresponding site

    fprintf('working on folder %s (%d of %d) ... ',all_sites{id},i,length(sites))
    
    % load transect information and beach slope
    data_file=[DATA_DIR,filesep,all_sites{id},filesep,'transect_coordinates_and_beach_slopes.csv']; % data file to load
    T=readtable(data_file); % load data as a matlab table

    tr_ID=T.TransectId; % get the transect id

    lon_on=T.Longitude_Origin; % get onshore transect coordinates
    lat_on=T.Latitude_Origin;

    lon_off=T.Longitude_SeawardsPoint; % get offshore transect coordinates
    lat_off=T.Latitude_SeawardsPoint;

    beach_slope=T.BeachFaceSlope; % get beach slope
    
    % convert transects lon,lat to UTM x,y
    [x_on,y_on]=wgs2utm(lat_on,lon_on,UTM_ZONE,UTM_HEM);
    [x_off,y_off]=wgs2utm(lat_off,lon_off,UTM_ZONE,UTM_HEM);
    
    phi=atan2d(y_off-y_on,x_off-x_on); % get the transect angle (pointing offshore)
    
    % fill in transects data in normal or reversed order
    if reorder(id1)
        fprintf(' ... REORDERING ... ');
        incr=length(x_on):-1:1;
    else
        incr=1:length(x_on);
    end

    for ii=incr
        
        transects(count,1).ID=count; % get a running counter

        transects(count,1).x_on=x_on(ii); % onshore coordinates
        transects(count,1).y_on=y_on(ii);

        transects(count,1).x_off=x_off(ii); % offshore coordinates
        transects(count,1).y_off=y_off(ii);
 
        transects(count,1).angle=phi(ii);   % transect angle

        transects(count,1).tanBeta=beach_slope(ii); % coastsat-derived beach slope

        transects(count,1).state='NSW';  % transect state/region/etc.

        transects(count,1).littoral_cell=name{id1}; % littoral cell name ... same as coastsat bounding box name for now

        transects(count,1).model_type='full model'; % set all model transects to be full model transects

        transects(count,1).Y0=NaN;          % the initial shoreline position (measured as distance from onshore end of transect)
        transects(count,1).Ymin=NaN;        % the minimum shoreline position (measured as distance from onshore end of transect) which represents the maximum eroded state or cliff toe or hard shoreline (i.e. the smallest distance the shoreline can be from the baseline)
       
        transects(count,1).t=NaN;           % times of the shoreline time series
        transects(count,1).Y=NaN;           % shoreline position data of the shoreline time series

        transects(count,1).SAT=NaN;         % boolean variable indicating shoreline position data point comes from satellites (=1) or not (=0)
        transects(count,1).UQ=NaN;          % boolean variable indicating the shoreline position uncertainty 

        transects(count,1).LTER=NaN;        % the long-term erosion rate (m/yr)

        transects(count,1).coastsat_tr_name=tr_ID{ii};                   % transect ID from coastsat
        transects(count,1).coastsat_bbox_name=tr_ID{ii}(1:7);            % coastsat bounding box name
        transects(count,1).coastsat_bbox_number=str2num(tr_ID{ii}(4:7)); % coastsat bounding box number

        count=count+1; % increment transect counter
    end

    if CORRECTED_SHORELINES % load tide-corrected satellite shorelines  
        % load corrected shoreline data folder
        data_file=[DATA_DIR,filesep,all_sites{id},filesep,'time_series_tidally_corrected.csv']; % point to current file in loop 
        T=readtable(data_file);                                                               % and load the data from a matlab table

        % get corrected shorelines
        transect_names={T.Properties.VariableNames{3:end}}'; % get all of the transect names in coastsat convention
        t=datenum(T.dates,'yyyy-mm-dd HH:MM:SS+00:00'); % get dates
        Y=table2array(T(:,3:end))';                   % get shoreline time series in 2-D array format         
    else % else load raw satellite shorelines
        % load raw shoreline data file
        data_file=[DATA_DIR,filesep,all_sites{id},filesep,'time_series_raw.csv']; % point to current file in loop 
        T=readtable(data_file);                                                 % and load the data from a matlab table

        % get corrected shorelines
        transect_names={T.Properties.VariableNames{3:end}}';
        t=datenum(T.dates,'yyyy-mm-dd HH:MM:SS+00:00'); % get dates
        Y=table2array(T(:,3:end))';                   % get shoreline time series in 2-D array format
    end

    transect_names=strrep(transect_names,'_','-'); % replace underscore with dash, since naming convention changed ever so sligtly

    % match transects name with Kilian's coastsat transect name
    [TF,id]=ismember(transect_names,{transects.coastsat_tr_name}');

    id=id(TF); % keep only the true entries

    if length(id)~=size(Y,1) % check to see if data arrays are the same size as expected
        error('Shouldn''t these be the same?. SeanV.')
    end

    % add shorelines time series to transects data structure
    for j=1:size(Y,1)
        id2=find(~isnan(Y(j,:)));
        transects(id(j)).t=t(id2);
        transects(id(j)).Y=Y(j,id2)';
        transects(id(j)).SAT=ones(size(Y(j,id2)'));
        transects(id(j)).UQ=NaN(size(Y(j,id2)'));
    end

    PLOT=0;
    if PLOT    % plots to show progress
        
        % shoreline positions 
        xs=x_on+Y.*cosd(phi);
        ys=y_on+Y.*sind(phi);
        
        figure(1);
        clf; hold on; box on;
        plot([x_on x_off]',[y_on y_off]','-r.');
        plot(x_on,y_on,'-g.');
        plot([[transects(id).x_on]; [transects(id).x_off]],[[transects(id).y_on]; [transects(id).y_off]],'--ko');
        plot(xs,ys,'-b.');
        axis equal;
        drawnow
        
    end

    fprintf('done.\n'); % display progress
    
end

% REMOVE TRANSECTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(remove)                                           % for all cells in excel sheet
    if remove(i)                                                 % if littoral cell is flagged for transects to be removed
        str=id2rm{i};                                            % get the transects to remove
        site2rm=strrep(sites{i},'s','s0');                       % add an extra zero to the string for compatibility with the transects format
        nums2rm=str2double(regexp(str,',','split'));             % get the numers to remove
        for j=1:length(nums2rm)                                  % for all of the transect numbers to remove
            site2rm1=site2rm+"-"+sprintf('%0.4d',nums2rm(j));         % get the coastsat transect name of the transect to remove
            [TF,id]=ismember(site2rm1,{transects.coastsat_tr_name}'); % find the transect to remove within the full list of transects
            transects(id).model_type='no prediction';                % set the transect to no prediction;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ntr=length(transects); % get total number of transects

% keep only unique data times (remove repeated observations)
for i=1:Ntr
    [~,ia,~]=unique(transects(i,1).t);
    transects(i,1).t=transects(i,1).t(ia);
    transects(i,1).Y=transects(i,1).Y(ia);
    transects(i,1).SAT=transects(i,1).SAT(ia);
    transects(i,1).UQ=transects(i,1).UQ(ia);
end

% sort the shoreline data by time (oldest to newest)
for i=1:Ntr
    [~,id]=sort(transects(i,1).t,'ascend');
    transects(i,1).t=transects(i,1).t(id);
    transects(i,1).Y=transects(i,1).Y(id);
    transects(i,1).SAT=transects(i,1).SAT(id);
    transects(i,1).UQ=transects(i,1).UQ(id);
    %transects(i,1).Y0=transects(i,1).Y(end);
end

% calculate long-term, linear erosion rate.
for i=1:Ntr
    
    DATA=transects(i,1).Y'; % get the shoreline position
    t_DATA=transects(i,1).t'; % and the time (matlab time)
    
    % remove the data this is NaN
    id=isnan(DATA);
    t_DATA(id)=[];
    DATA(id)=[];
    
    % remove the data that is prior to 1980
    id=t_DATA<datenum(1980,1,1);
    t_DATA(id)=[];
    DATA(id)=[];
    
    if length(DATA)>=3 % if we have enough data, do a robust fit
        % do a linear fit
        p1 = polyfit(t_DATA,DATA,1);
        transects(i,1).LTER=p1(1)*365.25; % convert to m/yr
    end
end

% REMOVE 'no prediction' TRANSECTS
if 1
    % remove no prediction transects
    idtrrm=find(strcmp({transects.model_type},'no prediction')); % find all 'no prediction' transects
    transects(idtrrm)=[]; % and remove them

    % renumber
    Ntr=length(transects); 
    for i=1:Ntr            % renumber each transect sequentially
        transects(i).ID=i;
    end
end

SAVE_TRANSECTS=1;
if SAVE_TRANSECTS
    save('transects','transects','-v7.3'); % save in current folder
end

% get littoral cell information from the transects
[littoral_cell_name_unique,littoral_cell_type,...
    littoral_cell_tr_start,littoral_cell_tr_end,...
    N_full_model_sections,N_cross_shore_only_sections,N_rate_only_sections,...
    id_full_model_start,id_full_model_end,...
    id_cross_shore_only_start,id_cross_shore_only_end,...
    id_rate_only_start,id_rate_only_end]=transects2cells(transects);

N_littoral_cells=length(littoral_cell_name_unique);

% write transects sctuct to kml file
transects2kml(transects,'transects.kmz')

% count the total number of shoreline data points for each transect
for i=1:Ntr
    NS(i)=length(transects(i).Y);
end

% plot the total number of shoreline data points across all transects (for good measure)
figure; plot(NS)
