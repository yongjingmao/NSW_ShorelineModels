% The CoSMoS-COAST New South Wales (NSW) Model. SeanV.
close all; clc;
clear all

% model name
Model_name='NSW';  fprintf('Starting CoSMoS-COAST %s model run ... \n',Model_name);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add matlab helper routines toolbox to the path (the kml toolobx is the most important here ... but others might be too)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matlab_helpers_toolbox='C:\Users\z3541792\OneDrive - UNSW\Documents\GitHub\MatlabToolbox';
addpath(genpath(matlab_helpers_toolbox));
% transects_toolbox='C:\Users\z3541792\OneDrive - UNSW\Research\CoSMoS-Coast\Australia\NSW\transects';
% addpath(genpath(transects_toolbox))

% surpress warnings
warning('off','all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seed random number generator so that the results are reproducible
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seed=1; rng(seed);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ensemble Model setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nens=200; % the # of members of the ensemble (typically 100-200, unless you can afford more)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model parameter setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ShoreFor_only = false; % Only consider ShoreFor
Longshore_var = true;
Nonstationary = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the transects data file (which is effectively the model grid/domain)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('transects','var')

    % the transects file to load
    TRANSECTS_FILE='..\transects\transects.mat';

    UTMZONE='56 H'; % the UTM zone that is used in the transects.mat file (this variable is used in the output/results scripts as well)

    fprintf('loading transects file ...');
    load(TRANSECTS_FILE); % load the transects file ... and display progress
    fprintf(' done. \n');

    % clear some clutter variables
    clear matlab_helpers_toolbox TRANSECTS_FILE;

else
    fprintf('using preivously loaded transects file.\n');
    clearvars -except transects Model_name UTMZONE Nens Ntr t Hs Tp Dir t_obs Y_obs Y_rms SAT wave_year1 wave_year2 ENSID WAVEID LOC R Hs_lt_trends;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup model time stepping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t0=datenum(1991,1,1);           % simulation start date

%tstop=datenum(2100,1,1);       % simulation stop date
tstop=datenum(2024,01,01);        % simulation stop date

tforecast=datenum(2017,01,01);   % the time when data asssimilation is turned off
tforecast2=datenum(2024,01,01);  % the time transition from a hindcast to a forecast (usually when no more data is available)

dt=1;                            % time step in [days] (1 day is almost always used for the time step)

Nsubcycles=1;                    % number of subcycles for longshore transport (increase this a bit if the model is unstable [due to the longshore transport term])

t1=(t0:dt:tstop);                % the simulation time
len=length(t1);                  % the number of time steps

nforecast2=find(t1>=tforecast2,1,'first'); % find the time step when of tforecast2 time

sim_time=NaN(1000,1); % a variable to help keep track of the computer time remaining

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model initial conditions part 1 (from transects file)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the model over specific transects
ids=[transects.ID];             % ALL transects to run the model over

% USE A LIMITED SUBSET OF TRANSECTS to run the model over
%ids=ids(ids>=500 & ids<=1000);         % some portion

% run the model for a certain littoral cell or coastsat bounding box only
% bbox=[transects.coastsat_bbox_number];
% ids=find(bbox==206);

[~,locb]=ismember(ids,[transects.ID]); % find the location of where the transects ID overlaps with the transects of interest specifed above

transects=transects(locb);   % overwrite transects file with transects of interst

Ntr=length(transects);       % transect length

ID=[transects.ID]';          % get the (overwritten) transect ID

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get transect type
if ShoreFor_only % Set all model_type to 'cross-shore only'
    for i = 1:length(transects)
        transects(i).model_type = 'cross-shore only';
    end
end

bool_full_model      =strcmp({transects.model_type},'full model')';        % transects that are 'full model'
bool_cross_shore_only=strcmp({transects.model_type},'cross-shore only')';  % transects that are 'cross-shore only'
bool_rate_only       =strcmp({transects.model_type},'rate only')';         % transects that are 'rate only'
bool_cliff_only      =strcmp({transects.model_type},'cliff only')';        % transects that are 'cliff only'
bool_no_prediction   =strcmp({transects.model_type},'no prediction')';     % transects that are 'no_prediction'

bool_full_model_mid=floor([0; 0.5*(bool_full_model(1:end-1)+bool_full_model(2:end)); 0]); % transects where the midpoints are sandwiched between two 'full model' transects

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE SOME TRANSECT ADJUSTMENTS
MAKE_transect_adjustments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% update the flag for the transects that might have become post-assigned 'no prediction'
bool_full_model_mid=bool_full_model_mid & ~ceil([0; 0.5*(bool_no_prediction(1:end-1)+bool_no_prediction(2:end)); 0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get some information contained within the transects data-struct (this is done to facilitate readability of the code, without constantly refering back to the large transects struct)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_on =[transects.x_on]';  % the x-coordinate of the onshore  baseline
y_on =[transects.y_on]';  % the y-corrdinate of the onshore  baseline
x_off=[transects.x_off]'; % the x-coordinate of the offshore baseline
y_off=[transects.y_off]'; % the y-corrdinate of the offshore baseline

% shoreline positions, angles, etc.
phi=[transects.angle]'; % the transect angle

Ymin=[transects.Ymin]'; % the "non-erodible" shoreline

%Ymin2=[transects.Ymin2]'; % the landward position of the transect crossing

% where Ymin does not exsit, set it to be the end of the transect.
Ymin(isnan(Ymin))=0;

% set the intial shoreline position
Y00=NaN(Ntr,1);  % the initial shoreline position (distance from offshore end of transect)
for i=1:Ntr
    [~,id]=min(abs(transects(i).t-t0));
    if ~isempty(id)
        Y00(i)=transects(i).Y(id);
    end
end

% set starting shoreline to be NaN for cliff only and no prediction shorelines
id=bool_cliff_only | bool_no_prediction;
Y00(id)=NaN;

% real world coordinates of the intitial shoreline position
x0=x_on+Y00.*cosd(phi);                 % the x-coordinate of the initial shoreline
y0=y_on+Y00.*sind(phi);                 % the y-coordinate of the initial shoreline

% set the mean shoreline position (only to calculate the closest wave observation points)
Ymean=Y00;  % the initial shoreline position (distance from offshore end of transect)
for i=1:Ntr
    if ~isnan(nanmean(transects(i).Y))  % if mean is not NaN, then get YMEAN as the mean shoreline data
        Ymean(i)=nanmean(transects(i).Y);
    end
end

% real world coordinates
xmean=x_on+Ymean.*cosd(phi); % the x-coordinate of the mean shoreline
ymean=y_on+Ymean.*sind(phi); % the y-coordinate of the mean shoreline

% if x/y-mean is still NaN then use the non-erodible position
id=isnan(xmean) | isnan(ymean);
xmean(id)=x_on(id)+Ymin(id).*cosd(phi(id)); % the x-coordinate of the mean shoreline
ymean(id)=y_on(id)+Ymin(id).*sind(phi(id)); % the y-coordinate of the mean shoreline

% if x/y-mean is still NaN then use the offshore location
id=isnan(xmean) | isnan(ymean);
xmean(id)=x_on(id); % the x-coordinate of the mean shoreline
ymean(id)=y_on(id); % the y-coordinate of the mean shoreline

% convert the x/y-mean positions to lat/lon, since the wave observation points are in lat/lon
[latmean, lonmean] = utm2deg(xmean,ymean,repmat(UTMZONE,Ntr,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model initial conditions part 2 (different shoreline change components)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set the ensemble intial shoreline position with NO perturbation ... if a
% perturbation is set here then that uncertainty will never go away, since
% Y0 is not changed via data assimilation
Y0=repmat(Y00,1,Nens); % the ensemble initial shoreline position

% set the ensemble short-term shoreline position with a perturbation of amp
amp=5; % perturbation to the short-term shoreline position in [m]
Yst=amp*repmat(randn(1,Nens),Ntr,1);   % the short-term shoreline position term (initialized as normally distributed ... not zero)

% set the ensemble shoreline anomaly term with a perturbation of amp
amp=0; % perturbation to the long-term shoreline anomaly term in [m] (I don't know if it's better to set this to uniformly zero or have some perturbation ... needs a bit more testing)
Ylst=amp*repmat(randn(1,Nens),Ntr,1);   % the (long-term) shoreline position term (initialized as normally distributed ... not zero)

Ylst_only=Ylst;       % the Ylst_ only is a variable to keep track of the perturbations to the long-term shoreline component WITHOUT positional adjustments due to data assimilation

Ybru=zeros(Ntr,Nens); % the Bruunian   recession component     (starts from zero)
Yvlt=zeros(Ntr,Nens); % the long-term rate       component     (starts from zero)

Ylt=Ylst+Ybru+Yvlt;       % the total long-term position, i.e., the sum of its components

Y  =Ylt+Yst+Y0;           % the total shoreline perturbation.  the sum of the short and long term parts, as well as the initial position

% real world coordinates of the current shoreline position (which includes preturbations)
x=repmat(x_on,1,Nens)+Y.*cosd(repmat(phi,1,Nens)); % the x-coordinate of the initial shoreline
y=repmat(y_on,1,Nens)+Y.*sind(repmat(phi,1,Nens)); % the y-coordinate of the initial shoreline

% repeat the minimum position
YMIN=repmat(Ymin,1,Nens);

% set all locations less than non-erodible position to the non-erodible position
Y(Y<YMIN)=YMIN(Y<YMIN);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid spacing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate the transect spacing, dx
DX00=NaN(Ntr,1);

DX00(2:end-1)=0.5*( sqrt( (x0(3:end  )-x0(2:end-1)).^2+(y0(3:end  )-y0(2:end-1)).^2 ) ...
    +sqrt( (x0(2:end-1)-x0(1:end-2)).^2+(y0(2:end-1)-y0(1:end-2)).^2 ) );

DX00(1  )=sqrt( (x0(2  )-x0(1    )).^2+(y0(2  )-y0(1    )).^2 );
DX00(end)=sqrt( (x0(end)-x0(end-1)).^2+(y0(end)-y0(end-1)).^2 );

% background/stock-standard transect spacing
DX0=100*ones(Ntr,1);

% ensemble transect spacing
DX=repmat(DX0,1,Nens);  % here I use stock standard grid spacing ... but the calculated values might also be used

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup littoral cells (mined from the transects data-struct)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get littoral cell information for all transect
littoral_cell_type_all={transects.model_type}';

% a function to take the transects and get the littoral cell information. SeanV.
[littoral_cell_names,littoral_cell_names_unique,littoral_cell_type,...
    littoral_cell_tr_start,littoral_cell_tr_end,...
    N_full_model_sections,N_cross_shore_only_sections,N_rate_only_sections,...
    id_full_model_start,id_full_model_end,...
    id_cross_shore_only_start,id_cross_shore_only_end,...
    id_rate_only_start,id_rate_only_end]=transects2cells(transects);

% get the length of each littoral cell
littoral_cell_tr_length=littoral_cell_tr_end-littoral_cell_tr_start+1;

% get the number of littoral cells
N_littoral_cells=length(littoral_cell_type);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load forcing conditions (waves and sea-level rise)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Wave conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (Re)load the wave forcing conditions
if ~exist('Hs','var')

    % the transects file to load
    WAVES_FILE='..\waves\waves.mat'; %load NEW bias corrected waves.

    fprintf('loading waves file ... ');
    w=load(WAVES_FILE);
    fprintf('done.\n');

    % associate wave conditions with transects one-to-one
    WAVEID=ids;

    % get wave conditions
    t=w.t;
    t_ID = find(ismember(t, t1));
    Hs=w.Hs_max(WAVEID,t_ID);
    Hsb=w.Hsb(WAVEID,t_ID);
    Tp=w.Tm(WAVEID,t_ID);
    Dir=w.Dir(WAVEID,t_ID);
    Fs=w.Fs(:, WAVEID, t_ID);
    omega=w.omega(WAVEID, t_ID);
    t=t(t_ID);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DEAL WITH NAN's IN THE EASIEST WAY POSSIBLE (for now)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Hs(isnan(Hs(:)))=nanmean(Hs(:));
    Tp(isnan(Tp(:)))=nanmean(Tp(:));
    Dir(isnan(Dir(:)))=nanmean(Dir(:));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CHANGE WAVE DIR to make it normal to the shoreline %%%%%%%%%%%%%%%%%%
    if 1
        % get (initial) shoreline angles
        shoreline_angle_mid=atan2d(diff(y0),diff(x0)); % the shoreline angle (Cartesian convention) at the transect midpoints
        alpha_shoreline_mid=shoreline_angle_mid-90;    % the shore-normal angle, pointing towards onshore location (-90 for West Coast setups when numbered south to north, +90 for East Coast setups)

        % get (initial) shoreline angle in natuical convention
        Dir_shoreline_mid=mod(270-alpha_shoreline_mid,360);

        % get (initial) shoreline angle at transects
        Dir_shoreline=cat(1,Dir_shoreline_mid(1),0.5*(Dir_shoreline_mid(1:end-1,:)+Dir_shoreline_mid(2:end,:)),Dir_shoreline_mid(end));

        % get mean wave direction
        Dir_mean=nanmean(Dir,2);

        % correct wave direction so that the mean wave direction is identical
        % to the (initial) shoreline direction
        Dir=Dir+repmat((Dir_shoreline-Dir_mean),1,len);

        % get (corrected) mean wave direction
        Dir_mean=nanmean(Dir,2);

        % compare shoreline direction and mean wave direction in Cartesian coordinates
        % plot(1:Ntr,Dir_shoreline,'-bo',1:Ntr,Dir_mean,'-r.');

        % get wave directions at midpoints between transects
        Dir_mid=0.5*(Dir(1:end-1,:)+Dir(2:end,:));

        alpha_wave_mid=270-Dir_mid;
        alpha_mid=alpha_wave_mid-alpha_shoreline_mid;   % the wave angle relative to the shoreline angle for each midpoint (alpha_wave_mid-alpha_shoreline_mid for West Coast setups, alpha_shoreline_mid-alpha_wave_mid for East Coast setups)

        % recenter alphas
        alpha_mid(alpha_mid>180)=alpha_mid(alpha_mid>180)-360;
        alpha_mid(alpha_mid<-180)=alpha_mid(alpha_mid<-180)+360;

    end


else
    fprintf('using preivously loaded waves files.\n');
end

if any(isnan(Hs(:))) || any(isnan(Tp(:))) || any(isnan(Dir(:)))
    error('NaN''s found in wave conditions (Hs, Tp, Dir)');
end

% set minimum wave conditions
Hs_min=0.25;
Tp_min=2;

% Sea-level rise% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SL=1.00;                    % Sea-level projection [meters] by the time tfuture (below).  This value is uniform for each transect.
SLRhist=0.003*1/365.25;     % the historical rate of sea level rise [meters/yr * yrs/day = meters/day]
tnow   =datenum(2000,1,1);  % start of the SLR curve
tfuture=datenum(2100,1,1);  % stop of SLR curve

t_SLR=t1;                   % the times to output the SLR scenario

% generate a quadratic SLR curve w.r.t. time
A=[tnow.^2 tnow 1; 2*tnow 1 0; tfuture.^2 tfuture 1]; % we assume the the sea-level curve is a 2nd-order polynomial and that the rate at t0 is SLRhist (=~2mm/yr) and that the sea level at tstop = SL and that the SL at t0 is 0
rhs=[0; SLRhist; SL];                                 % this gives us enough equations to determine the coefficients of the sea-level curve

xsl=A\rhs;                          % find the coefficients of the quadratic curve

asl=xsl(1); bsl=xsl(2); csl=xsl(3); % extract the SLR curve coefficients

S=asl*t_SLR.^2+bsl*t_SLR+csl;                   % the projected SLR curve S
S(t_SLR<tnow)=(t_SLR(t_SLR<tnow)-tnow)*SLRhist; % the historical part of the projected SLR curve is set to be the historical rate (given above)

if 0 % plots the SLR curve in case you want to check it
    plot(t1,S,'b',tnow,0,'ko',tfuture,SL,'ko'); han1=legend('SLR curve','t_{now}','t_{future}'); set(han1,'location','NorthWest'); datetick('x','keeplimits');
    stop  % my favorite matlab command is 'stop'. Matlab doesn't know what 'stop' means so it stops when it reaches the command 'stop'
end

% the mean sea-level scenario
Smean=S;    % if spatially varying sea-level rise simulations are used, then Smean is the mean of many SLR curves ...
% if only one SLR scenario exists then the mean scenario (Smean) is just replicated.

% the indices of the SLR scenarios used to set the SLR conditions on each model transect
SLRID=ones(Ntr,1); % this variable is analagous to WAVEID, but acts to assign the SLR projections (at their specific locations) to each transect
% If we are using any SLR scenarios at using any
% lat/longs from the SLR projections, then closest
% lat/lon of the transects can be used

% clear some clutter variables
clear han1 SLRhist tnow tfuture A rhs xsl asl bsl csl;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup model parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% state vector parameters used in data assimilation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initial background/mean parameter values (these initial values are uniform for all transects of the same type)

% equlibirum shoreline change parameters
DT00=50;   % equlibrium    time scale [days]
DY00=5;   % equlibrium erosion scale [meters]       (for full model transects)
DY01=5;   % equlibrium erosion scale [meters]       (for cross-shore only transects)

HSB0=nanmean(Hs(:,t>t0 & t<tforecast2,1),2); % the equlibrium background wave height parameter [meters] (variable for each transect)

% long-term parameters
c00  =1;    % Bruun coefficient               [-]
vlt00=0;    % long-term rate term             [m/yr]
K00  =100;   % longshore transport coefficient [-]    (for full model transects)

% noise parameter
sigma00=0.25; % adative noise [units vary but nominally meters]

% the spreads of initial background parameter values (these initial values are uniform for all transects of the same type)

% equlibirum shoreline change parameters
sig_DT=3.5;      % equlibrium    time scale spread [days]
sig_DY=3;        % equlibrium erosion scale spread [meters] (same for all transect types)

sig_HSB=0.075*nanmean(Hs(:,t>t0 & t<tforecast2,1),2);  % equlibrium wave height spread [meters] (variable for each transect)

% long-term parameters
sig_c  =0.05;  % Bruun coefficient spread        [-]
sig_vlt=0.025; % long-term rate term             [m/yr]
sig_K  =20;    % longshore transport coefficient [-]     (full model transects)

% noise parameter
sig_sigma=0.01;  % adative noise spread [units vary but nominally meters]

% upper and lower bounds of parameter values (these initial values are uniform for all transects of the same type)

% equlibirum shoreline change parameters
lb_DT =10;  ub_DT=120;  % equlibrium    time scale lower and upper bounds (respectively) [days]
lb_DY =4;  ub_DY=30;    % equlibrium erosion scale lower and upper bounds (respectively) [meters]  (same for all transect types)

lb_HSB=0.01; ub_HSB=4; % equlibrium wave height lower and upper bounds (respectively) [meters]  (variable for each transect)

% long-term parameters
lb_c  =0.01;   ub_c  =2;    % Bruun coefficient               lower and upper bounds (respectively)  [-]
lb_vlt=-5;     ub_vlt=5;    % long-term rate term             lower and upper bounds (respectively)  [m/yr]
lb_K  =0.01;   ub_K  =4000;  % longshore transport coefficient lower and upper bounds (respectively)  [-]     (full model transects)

% noise parameter
lb_sigma=0.0001;   ub_sigma=2;  % adative noise spread [units vary but nominally meters]

% make probability distribution objects (to generate the initial conditions ... generating truncated distributions is the most helpful thing about using these objects)
% these are not the acutal parameter arrays ... that comes later (and uses these objects)
fprintf('generating initial conditions ...');

% make probability distribution objects for the equlibirum shoreline change parameters
pdf_DT =truncate(makedist('Normal','mu',DT00,'sigma',sig_DT),lb_DT,ub_DT); % make a truncated normal distribution for DT,  the equlibrium    time scale [days]
pdf_DY =truncate(makedist('Normal','mu',DY00,'sigma',sig_DY),lb_DY,ub_DY); % make a truncated normal distribution for DY,  the equlibrium erosion scale [m]     (for full  model      transects)
pdf_DY1=truncate(makedist('Normal','mu',DY01,'sigma',sig_DY),lb_DY,ub_DY); % make a truncated normal distribution for DY,  the equlibrium erosion scale [m]     (for cross-shore only transects)


for i=Ntr:-1:1 % march in reverse to improve preallocation
    pdf_HSB(i,1)=truncate(makedist('Normal','mu',HSB0(i),'sigma',sig_HSB(i)),lb_HSB,ub_HSB); % make a truncated normal distribution for HSB, the equlibrium wave height scale [m]
end

% make probability distribution objects for the long-term parameters

%  generate the initial conditions for the Bruun coefficient
pdf_c  =truncate(makedist('Normal','mu',c00  ,'sigma',sig_c  ),lb_c  ,ub_c  ); % make a truncated normal distribution for c, the Bruun coefficient [-]

% long-term erosion rate (if available from the transects file,used to initialize the variable vlt)
LTER=[transects.LTER]';

% make some adjustments to the long-term erosion rate (i.e., set some to zero)
MAKE_rate_adjustments;

% use smoothed historical rates for initial conditions?
USE_SMOOTHED_RATES=1;

% calculate a smoothed long-term erosion rate (which is often used for rate only transects if data assimilation is turned off)
LTER_smooth=LTER;
LTER_smooth(bool_cliff_only | bool_no_prediction)=0;
LTER_smooth=smoothn(LTER_smooth,10); % smooth

% plot the smoothed rates
%plot(1:Ntr,LTER,'-b.',1:Ntr,LTER_smooth,'-r')

%  generate the initial conditions for the long-term erosion rate term (uniform alongshore)
% pdf_vlt=truncate(makedist('Normal','mu',vlt00,'sigma',sig_vlt),lb_vlt,ub_vlt); % make a truncated normal distribution for DY,  the equlibrium     erosion scale [m]

if USE_SMOOTHED_RATES
    for i=Ntr:-1:1 % march in reverse to improve preallocation (I can't quite remember why this was a good idea)
        if bool_rate_only(i)
            pdf_vlt(i,1)=truncate(makedist('Normal','mu',max(min(1.0*LTER_smooth(i),0.75*ub_vlt),0.75*lb_vlt),'sigma',sig_vlt),lb_vlt,max(1.5*LTER_smooth(i),1));
        else
            pdf_vlt(i,1)=truncate(makedist('Normal','mu',max(min(0.5*LTER_smooth(i),0.75*ub_vlt),0.75*lb_vlt),'sigma',sig_vlt),lb_vlt,max(1.5*LTER_smooth(i),1)); % make a truncated normal distribution for vlt
        end
    end
else
    for i=Ntr:-1:1 % march in reverse to improve preallocation (I can't quite remember why this was a good idea)
        if ~isnan(LTER(i))
            if bool_rate_only(i)
                pdf_vlt(i,1)=truncate(makedist('Normal','mu',max(min(1.0*LTER(i),0.75*ub_vlt),0.75*lb_vlt),'sigma',sig_vlt),lb_vlt,max(1.5*LTER(i),1)); % make a truncated normal distribution for vlt
            else
                pdf_vlt(i,1)=truncate(makedist('Normal','mu',max(min(0.5*LTER(i),0.75*ub_vlt),0.75*lb_vlt),'sigma',sig_vlt),lb_vlt,max(1.5*LTER(i),1)); % make a truncated normal distribution for vlt
            end
        else
            pdf_vlt(i,1)=truncate(makedist('Normal','mu',0,'sigma',sig_vlt),lb_vlt,ub_vlt); % make a truncated normal distribution for vlt
        end
    end
end

% generate the initial conditions for the long-term erosion rate term (uniform alongshore)
%pdf_K  =truncate(makedist('Normal','mu',K00  ,'sigma',sig_K  ),lb_K  ,ub_K  ); % make a truncated normal distribution for HSB, the equlibrium wave height scale [m]
pdf_K  =truncate(makedist('Uniform','Lower',lb_K,'Upper',ub_K),lb_K  ,ub_K  ); % make a truncated normal distribution for HSB, the equlibrium wave height scale [m]

% generate the initial conditions for thenoise parameter (which is assimilated ... I don't know if this is a good idea or not)
pdf_sigma=truncate(makedist('Normal','mu',sigma00,'sigma',sig_sigma),lb_sigma,ub_sigma); % make a truncated normal distribution for sigma, the addative noise parameter [units vary but nominally meters]

% preallocate initial parameter transect arrays

% preallocate arrays for the initial conditions of the equlibirum shoreline change parameters
DT0=NaN(Ntr,1);
DY0=NaN(Ntr,1);

% fill only full model or cross-shore only transects
id=(bool_full_model | bool_cross_shore_only ); DT0(id)=DT00;
id=(bool_full_model                         ); DY0(id)=DY00;
id=(                  bool_cross_shore_only ); DY0(id)=DY01;

% preallocate arrays for the long-term parameters
c0  =NaN(Ntr,1);
vlt0=NaN(Ntr,1);
K0  =NaN(Ntr,1);

% fill full model, cross-shore only, or rate only transects
id=(bool_full_model | bool_cross_shore_only | bool_rate_only);

c0(id)=c00;
vlt0(id)=LTER(id);

% fill only full model transects
id=bool_full_model;
K0(id)=K00;

% noise parameter
sigma0=NaN(Ntr,1);

id=(bool_full_model | bool_cross_shore_only | bool_rate_only);
sigma0(id)=sigma00;

% create initial parameter distribution ensembles, some of which are
% alongshore uniform within their littoral cell, and some have variable
% shapes to better assimilate the spatial patters of each model parameter

% create initial ensembles for the equlibirum shoreline change parameters
DTi =NaN(Ntr,Nens);
DYi =NaN(Ntr,Nens);
HSBi=NaN(Ntr,Nens);

%  create initial ensembles for the long-term parameters
ci  =NaN(Ntr,Nens);
vlti=NaN(Ntr,Nens);
Ki  =NaN(Ntr,Nens);

% create initial ensembles for the noise parameter
sigmai=NaN(Ntr,Nens);

% Make the parameter arrays with pdf objects

% initialize the background wave height parameter
for i=1:Ntr
    HSBi(i,:)=random(pdf_HSB(i),1,Nens);
end

% smoothing filter to smooth the longshore variability (only for HSBi)
[b_smooth,a_smooth] = butter(1,0.025,'low');     % design low-pass filter

% initialize spatially uniform parameter arrays for the Bruun and noise coefficients
id=(bool_full_model | bool_cross_shore_only | bool_rate_only);

% variable to mix random-by-transect and alongshore uniform random numbers
theta_rand=0.4;

ci(id,:)=theta_rand*random(pdf_c    ,sum(id),Nens)+(1-theta_rand)*repmat(random(pdf_c    ,1,Nens),sum(id),1);
sigmai(id,:)=theta_rand*random(pdf_sigma,sum(id),Nens)+(1-theta_rand)*repmat(random(pdf_sigma,1,Nens),sum(id),1);

% march thru each transect individually for the long term rate paraemter,
% vlt, since the pdf objects are different for each transect
for i=1:Ntr
    if id(i)
        vlti(i,:)=random(pdf_vlt(i),1,Nens);
    end
end

% hyperparameters for generating spatial variability
max_order=3;           % maximum order of the trig functions / Legendre polynomials
%coef=-0.5:0.01:0.5;  % the amplitude perturbations to these functions
%coef2=0.25:0.01:1.5;  % the amplitude perturbations to these functions
coef2=0.5:0.01:1.25;  % the amplitude perturbations to these functions
coef_DT=0.8:0.01:1.2;  % the amplitude perturbations to these functions

for i=1:N_littoral_cells  % for all littoral cells

    % get the id's of the 'full model' littoral cells
    id=(littoral_cell_tr_start(i):littoral_cell_tr_end(i))';

    % smooth the background wave height
    if length(id)>3
        HSBi(id,:)=filtfilt(b_smooth,a_smooth,HSBi(id,:));
    end

    if strcmp(littoral_cell_type(i),'full model') % for 'full model' littoral cells

        DTi(id,:)=repmat(random(pdf_DT,1,Nens),length(id),1); % make alongshore uniform distribution for DT
        DYi(id,:)=repmat(random(pdf_DY,1,Nens),length(id),1); % make alongshore uniform distribution for DY

        % adjust initial parameter ensemble to include some spatial variability based on sin/cos/Legendre Polynomials
        %LP=generate_alongshore_shape_functions_ensemble(length(id),Nens,max_order,coef);
        LP=generate_alongshore_shape_functions_ensemble2(length(id),Nens,coef_DT);
        DTi(id,:)=LP.*DTi(id,:);

        % adjust initial parameter ensemble to include some spatial variability based on sin/cos/Legendre Polynomials
        %LP=generate_alongshore_shape_functions_ensemble(length(id),Nens,max_order,coef);
        LP=generate_alongshore_shape_functions_ensemble2(length(id),Nens,coef2);
        DYi(id,:)=LP.*DYi(id,:);

        % initial parameter ensemble for longshore transport coefficient
        Ki(id,:)=repmat(random(pdf_K  ,1,Nens),length(id),1); % alongshore uniform distribution for K

        % adjust initial parameter ensemble to include some spatial variability based on sin/cos/Legendre Polynomials
        %LP=generate_alongshore_shape_functions_ensemble(length(id),Nens,max_order,coef);
        %LP=generate_alongshore_shape_functions_ensemble2(length(id),Nens,coef2);
        LP=generate_alongshore_shape_functions_ensemble3(length(id),Nens,coef2);
        Ki(id,:)=LP.*Ki(id,:);

    elseif strcmp(littoral_cell_type(i),'cross-shore only') % for 'cross-shore only' littoral cells

        DTi(id,:)=repmat(random(pdf_DT ,1,Nens),length(id),1); % make alongshore uniform distribution for DT
        DYi(id,:)=repmat(random(pdf_DY1,1,Nens),length(id),1); % make alongshore uniform distribution for DY (based on pdf_DY1)

        % adjust initial parameter ensemble to include some spatial variability based on sin/cos/Legendre Polynomials
        %LP=generate_alongshore_shape_functions_ensemble(length(id),Nens,max_order,coef);
        LP=generate_alongshore_shape_functions_ensemble2(length(id),Nens,coef_DT);
        DTi(id,:)=LP.*DTi(id,:);

        % adjust initial parameter ensemble to include some spatial variability based on sin/cos/Legendre Polynomials
        %LP=generate_alongshore_shape_functions_ensemble(length(id),Nens,max_order,coef);
        LP=generate_alongshore_shape_functions_ensemble2(length(id),Nens,coef2);
        DYi(id,:)=LP.*DYi(id,:);

    end

end

% bound the initial longshore transport parameter Ki
Ki(Ki>ub_K)=ub_K;
Ki(Ki<lb_K)=lb_K;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE SOME PARAMETER ADJUSTMENTS
MAKE_parameter_adjustments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% start the parameter arrays with the initial values
DT =DTi;
DY =DYi;
HSB=HSBi;
c  =ci;
vlt=vlti;
K  =Ki;
sigma=sigmai;

fprintf(' done. \n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% other parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% beach slopes
tanBeta=[transects.tanBeta]';                % get (subaerial) beach slopes from the model transects

%tanAlpha=[transects.tanAlpha]';              % get (average) beach slopes from the model transects
tanAlpha=1/50*ones(Ntr,1);  % constant average beach slope

tanBeta_Bruun=tanBeta;                       % set the Bruun rule beach slope to the passive flooding beach slope ... this is a bit ad hoc, but is somewhat consistent with Lidar profiles.  We apply the tanAlpha scaling later.

% depth of closure
Dc0=10; % in [m].  This is usually set to 10 m, for lack of a good alternative

Dc=Dc0*ones(Ntr,1); % uniform depth of closure for all transects

DC=repmat(Dc,1,Nens); % uniform depth of closure ensemble term for all ensembles

% smoother used in the main time-stepping loop
[b_smooth,a_smooth] = butter(2,0.05,'low');     % design low-pass filter to smooth longshore variability of parameters

% more wave conditions parameters

% wave angle taper function to deal with refraction/diffraction dur to long
% gryones at the end of littoral cells
f_alpha_taper=@(x,sig,p) tanh(x/sig).^p;
sig_taper=4;
pow_taper=1;
TAPER=ones(Ntr,1);

% smooth for each littoral cell
for i=1:N_littoral_cells

    % get id's of each cell
    id=(littoral_cell_tr_start(i):littoral_cell_tr_end(i))';

    % generate an alongshore coordinate helper function
    if mod(length(id),2)
        xtmp=[linspace(0,floor(length(id)/2),floor(length(id)/2)) floor(length(id)/2) linspace(floor(length(id)/2),0,floor(length(id)/2))]';
    else
        xtmp=[linspace(0,length(id)/2,length(id)/2) linspace(length(id)/2,0,length(id)/2) ]';
    end

    % apply taper function
    TAPER(id)=f_alpha_taper(xtmp,sig_taper,pow_taper);

end

% the taper function at the cell midpoints
TAPER_mid=0.5*(TAPER(1:end-1)+TAPER(2:end));

% specify a minimum of the taper function
TAPER_min=0.1;
TAPER_mid(TAPER_mid<TAPER_min)=TAPER_min;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup data assimilation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% turn on data assimilation?
DATA_ASSIMILATION=1;

% turn on data assimilation for rate only transects? (It is often perfererable to turn off data assimilation for rate only transects with very few data points)
DATA_ASSIMILATION_FOR_RATE_ONLY_TRANSECTS=0;

% additive noise
ADDITIVE_NOISE=1;
ADDITIVE_NOISE_LT=1;

y_rms_gps=1;  % rms error [in meters] of lidar/gps shoreline observations
y_rms_sat=10; % rms error [in meters] of the satellite-derived shoreline observations

sig=0.5; % widths of the modified state parameters

% the data assimilation variables that get assembled into the assimilated model state vector
SYMBOLIC_MODE=0; % turning SYMBOLIC_MODE is an in-development feature that allows one to change the variables that get assimilated in an automated way.
% it relies heavily on Matlab's 'eval' function, which can be a bit slow.
% turning this off is a bit faster ... for the time being.

da_vars           ={'Ylst'  'Yst'   'DT'   'DY'   'HSB'   'c'    'vlt'   'K'     'sigma'}; % the names of the variables used in data assimilation
da_var_is_positive=[   0      0      1      1      1       1      0       1       1     ]; % flag to indicate if the data assimilation variable should be forced to be positive
da_var_smooth     =[   0      0      1      1      1       1      1       1       1     ]; % flag to indicate if the data assimilation variable should be smoothed
da_var_alpha      =[   0      0     0.25    0.25   0.25    0.25   0.25   0.95     0.75  ]; % smoothing value
da_var_bound      =[   0      0      1      1      1       1      1       1       1     ]; % flag to indicate if the data assimilation variable should be bounded
da_var_add_noise  =[   1      1      1      1      1       1      1       1       1     ]; % flag to indicate if noise should be added to the data assimilation variable
%da_var_noise_fac  =[  0.5    2.0   0.025  0.025  0.0025  0.005    0.01   0.1     0.0005]; % noise factor coefficient
da_var_noise_fac  =[  1.2    0.5   0.025  0.025  0.0025  0.005    0.01   0.025     0.00025]; % noise factor coefficient

Nvars=length(da_vars); % the number of variables in the assimilated state vector

% make a state vector string (which is used later) that contains all of the (modified) variables in the state vector
state_vector_str='[';
for ii=1:Nvars
    state_vector_str=[state_vector_str,da_vars{ii},'1(:,k) '];
end
state_vector_str=[state_vector_str,']'];
% for example, [Ylst1(:,k) Yst1(:,k) DT1(:,k) DY1(:,k) HSB1(:,k) c1(:,k) vlt1(:,k) K1(:,k) sigma1(:,k)]

%%

if ~exist('Y_obs','var')

    % the times when data assimilation takes place on any (full model, cross-shore only, or rate only) transect
    t_obs=[];

    % get the id's for the transects of interest
    id=find(bool_full_model | bool_cross_shore_only | bool_rate_only);

    for i=1:length(id)
        t_obs=union(t_obs,round(transects(id(i)).t));
    end

    % remove assimilation times that are before or after the calibration period
    t_obs(t_obs<t0)=[];
    t_obs(t_obs>=tforecast2)=[];

    t_obs=t_obs';

    % preallocate arrays for the observations at each assimilation time (for all transects)
    Y_obs=NaN(Ntr,length(t_obs));
    Y_rms=NaN(Ntr,length(t_obs));
    SAT=zeros(Ntr,length(t_obs));

    % get the observations at each assimilation time
    fprintf('preparing data for assimilation ... ');
    for n=1:length(t_obs)   % for all assimilation times
        for i=1:length(id)  % for all transects of interest
            [TF,id1]=ismember(t_obs(n),round(transects(id(i)).t)); % determine if this transect contains data at time t_obs(n)
            if TF                                    % if it does, then save that data
                Y_obs(id(i),n)=transects(id(i)).Y(id1);  % get the observation at that time
                if transects(id(i)).SAT(id1)==1      % if that data source is from satellites
                    Y_rms(id(i),n)=y_rms_sat;        % get the satellite rms error at that time
                    SAT(id(i),n)=1;                  % and specify that that dat comes from a satellite
                else                                 % if not, then specify that the data comes from gps/lidar
                    Y_rms(id(i),n)=y_rms_gps;        % get the lidar/gps rms error at that time
                end
            end
        end
    end
    fprintf('done.\n');

    % if any assimilation times have no data as a result of the previous
    % removal step, then remove this assimilation time
    id_rm=find((sum(isnan(Y_obs))==Ntr));

    t_obs(id_rm)=[];
    Y_obs(:,id_rm)=[];
    Y_rms(:,id_rm)=[];
    SAT(:,id_rm)=[];

    % preform a wave setup correction of the observations on each transect
    setup_obs=NaN(Ntr,length(t_obs));
    setup_correction=NaN(Ntr,length(t_obs));
    for i=1:Ntr
        Hs_obs=interp1(t',Hs(i,:,1)',t_obs);                       % find      the wave height at the time of the observations
        Tp_obs=interp1(t',Tp(i,:,1)',t_obs);                       % find      the wave period at the time of the observations
        L0_obs=9.81*Tp_obs.^2/(2*pi);                            % calculate the wave period at the time of the observations
        setup_obs(i,:)=0.016*sqrt(Hs_obs.*L0_obs)+0.15;          % calculate the setup at the time of the observations
        setup_correction(i,:)=setup_obs(i,:)./nanmean(tanBeta(i,:)); % calculate the setup correction at the time of the observations
    end

    % for some reason we still have negative Hs, Tp generated from the
    % interpolation (extrapolation) ... and I have no idea why ... so we keep only the real part.
    setup_correction=real(setup_correction);

    if 0  % plot the setup correction
        figure; pcolor(setup_obs); shading flat; colormap(jet); colorbar;
        figure; pcolor(setup_correction); shading flat; colormap(jet); colorbar;
        stop
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MAKE SOME ADJUSTMENTS TO THE DATA-ASSIMILATED OBSERVATIONS
    MAKE_obs_adjustments
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Y_obs=Y_obs+SAT.*setup_correction; % add the setup correction to the actual observations

    % MHW vs MSL shoreline correction
    mhw=0.488;  % approximate mean high water level height above mean sea level for Duck, NC
    msl=0;
    mhw_correction=15*(mhw-msl); % median beach slope is ~1/15 based on previous CoastSat data

    Y_obs=Y_obs+(SAT==0).*mhw_correction;  % correct MHW (GPS/LIDAR) shorelines to MSL shorelines using beach slope

    % Use satellite data?
    USE_SAT_DATA=1;
    if ~USE_SAT_DATA
        % set the satellite derived obsevations to NaN
        Y_obs(SAT==1)=NaN;

        % and remove the satellite data
        id_rm=find((sum(isnan(Y_obs))==Ntr));

        t_obs(id_rm)=[];
        Y_obs(:,id_rm)=[];
        Y_rms(:,id_rm)=[];
        SAT(:,id_rm)=[];
    end

else
    fprintf('using preivously prepared assimilation data.\n');
end

% assess how much data comes from satellites
n_sat=0; n_all=0; n_obs=NaN(Ntr,1);
for i=1:Ntr
    n_all=n_all+length(transects(i).SAT);
    n_sat=n_sat+sum(transects(i).SAT);
    n_obs(i)=length(transects(i).Y);
end

% NEW INITIAL CONDITION (based on the first (corrected) observation on each transect where it exists)
for i=1:Ntr
    id=find(~isnan(Y_obs(i,:)),1,'first');
    if ~isempty(id)
        Y00(i)=Y_obs(i,id);
        Y0(i,:)=Y_obs(i,id);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE SOME ADJUSTMENTS TO THE INITIAL SHORELINE (after satellite corrections)
MAKE_initial_shoreline_adjustments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set starting shoreline to be NaN for cliff only and no prediction shorelines
id=bool_cliff_only | bool_no_prediction;
Y0(id,:)=NaN;
Y00(id)=NaN;

% real world coordinates of the intitial shoreline position (AGAIN)
x0=x_on+Y00.*cosd(phi);                 % the x-coordinate of the initial shoreline
y0=y_on+Y00.*sind(phi);                 % the y-coordinate of the initial shoreline

Y  =Ylt+Yst+Y0;                          % the total shoreline perturbation.  the sum of the short and long term parts, as well as the initial position

% real world coordinates of the current shoreline position (which includes preturbations)
x=repmat(x_on,1,Nens)+Y.*cosd(repmat(phi,1,Nens)); % the x-coordinate of the initial shoreline
y=repmat(y_on,1,Nens)+Y.*sind(repmat(phi,1,Nens)); % the y-coordinate of the initial shoreline

% set all locations less than non-erodible position to the non-erodible position
Y(Y<YMIN)=YMIN(Y<YMIN);
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Localization matrix (used for data assimilation across littoral cells)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('LOC','var') % if localization matrix variable exists

    fprintf('preparing localization matrix ... '); % display progress

    % the distance matrix that gives the distance from transect i to transect j
    TR_DIST=NaN(Ntr);
    for i=1:Ntr
        for j=1:Ntr
            TR_DIST(i,j)=abs(i-j);
        end
    end

    % calculate observation error matrix
    sigR1=1;   % the first exponential decay function coefficient           (units are # of transects: a value of 1=200 meters)
    sigR2=2;   % the tanh coefficient                                       (units are # of transects: a value of 2=400 meters)
    sigR3=25;  % the second (slower) exponential decay function coefficient (units are # of transects: a value of 25=5000 meters)

    % generate fake R matrix
    f=@(x) (exp(-x/sigR1)+0.5*tanh(x/sigR2)).*exp(-x/sigR3);

    % error covariance matrix
    R=f(TR_DIST);

    % number of transects which halves the influence of the observation point
    tr_decay_dist=2;

    % make the localization matrix (using a sparse matrix representation ... or not)
    % a sparse matrix representation may be helpful for LARGE model runs
    SPARSE=1;
    if SPARSE
        IN_LITTORAL_CELL=sparse(Ntr);  % a (sparse) boolean matrix that is true if transect i is in the same littoral cell as transect j
    else
        IN_LITTORAL_CELL=zeros(Ntr);   % a          boolean matrix that is true if transect i is in the same littoral cell as transect j
    end

    % make the boolean matrix to determine if transect i is in the same
    % littoral cell as transect j
    for j=1:Ntr
        id=littoral_cell_tr_start<=j & littoral_cell_tr_end>=j;
        IN_LITTORAL_CELL(littoral_cell_tr_start(id):littoral_cell_tr_end(id),j)=1; % this step is slow for sparse representation of large model domains ... but you only have to do it once
    end

    % make (possibly sparse) localization matrix
    LOC=2.^-(TR_DIST/tr_decay_dist).*IN_LITTORAL_CELL;

    fprintf('done.\n'); % display progress
else
    fprintf('using preivously localization matrices.\n');
end

% noise factor (reduce the addative noise for certain transects that seem to diverge)
NOISE_FAC=ones(Ntr,1);
NOISE_FAC_YLT=ones(Ntr,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % MAKE SOME ADJUSTMENTS TO THE ADDITIVE NOISE FACTOR (called "NOISE_FAC")
MAKE_additive_noise_adjustments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup model management scenarios & output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Management scenarios to run and directory to save results to
%HOLD_THE_LINE=1; CONTINUED_ACCRETION=0; OUTPUT_DIR='results\ManagmentCase1_hold_the_line_no_continued_accretion';
%HOLD_THE_LINE=1; CONTINUED_ACCRETION=1; OUTPUT_DIR='results\ManagmentCase2_hold_the_line_continued_accretion';
%HOLD_THE_LINE=0; CONTINUED_ACCRETION=0; OUTPUT_DIR='results\ManagmentCase3_no_hold_the_line_no_continued_accretion';
if ShoreFor_only
    HOLD_THE_LINE=0; CONTINUED_ACCRETION=1; OUTPUT_DIR='results\CS_only';
else
    HOLD_THE_LINE=0; CONTINUED_ACCRETION=1; OUTPUT_DIR='results\CS_LS';
end

SAVE_DATA=1;  % data is to be saved
if SAVE_DATA

    % make various output directories if they do not exist
    if ~exist(OUTPUT_DIR,'dir')
        mkdir(OUTPUT_DIR);
    end
    if ~exist([OUTPUT_DIR,filesep,'xlsx'],'dir')
        mkdir([OUTPUT_DIR,filesep,'xlsx']);
    end
    if ~exist([OUTPUT_DIR,filesep,'google_earth'],'dir')
        mkdir([OUTPUT_DIR,filesep,'google_earth']);
    end
    if ~exist([OUTPUT_DIR,filesep,'mat'],'dir')
        mkdir([OUTPUT_DIR,filesep,'mat']);
    end
    if ~exist([OUTPUT_DIR,filesep,'figures'],'dir')
        mkdir([OUTPUT_DIR,filesep,'figures']);
    end
    %     if ~exist([OUTPUT_DIR,filesep,'..',filesep,'LTER'],'dir')
    %         mkdir([OUTPUT_DIR,filesep,'..',filesep,'LTER']);
    %     end
    if ~exist([OUTPUT_DIR,filesep,'source'],'dir')
        mkdir([OUTPUT_DIR,filesep,'source']);
    end
    if ~exist([OUTPUT_DIR,filesep,'shp'],'dir')
        mkdir([OUTPUT_DIR,filesep,'shp']);
    end

    %copyfile([mfilename,'.m'],[OUTPUT_DIR,filesep,'source',filesep,'CoSMoS_COAST_code_for_current_run.m']);

    % copy correct legend to Google Earth output directory
    % legend_fig='Legend.png';
    % if HOLD_THE_LINE
    %     copyfile('Legend_Case1_2_FloSup.png',[OUTPUT_DIR,filesep,'google_earth',filesep,legend_fig]);
    % else
    %     copyfile('Legend_Case3_4_FloSup.png',[OUTPUT_DIR,filesep,'google_earth',filesep,legend_fig]);
    % end

    outputfilename=['CoSMoS_COAST_',Model_name,'_',num2str(SL*100,3),'cm_SLR_by_2100']; % name of output file

    % create a new kml object
    kmlout = kml(outputfilename);

    % times to save model results to kml
    t_save=[tstop ...
        t0+1 ...
        datenum(1990:5:2100,1,1) ...
        datenum(2044,1,1) ...
        datenum(2069,1,1) ...
        datenum(2089,1,1) ...
        t1(find(Smean>=0.25,1,'first')) ...
        t1(find(Smean>=0.50,1,'first')) ...
        t1(find(Smean>=0.75,1,'first')) ...
        t1(find(Smean>=1.0 ,1,'first')) ...
        t1(find(Smean>=1.25,1,'first')) ...
        t1(find(Smean>=1.5 ,1,'first')) ...
        t1(find(Smean>=1.75,1,'first')) ...
        t1(find(Smean>=2.0 ,1,'first')) ...
        t1(find(Smean>=2.5 ,1,'first')) ...
        t1(find(Smean>=3.0 ,1,'first')) ];  % results will always be saved at the end

    SAVE_XLSX_FILE=1; % save the excel file?

    % when to save xlsx file
    %t_save_xlsx=[datenum(2044,1,1) datenum(2069,1,1) datenum(2089,1,1) datenum(2100,1,1)];
    t_save_xlsx=datenum(2100,1,1);

    % standard output interval in days
    output_interval=14; % in days

    % model output times
    t_output_times=[t0+1 ...
        t0+1:output_interval:tstop ...
        tstop ...
        datenum(1990:5:2100,1,1) ...
        t1(find(Smean>=0.00,1,'first')) ...
        t1(find(Smean>=0.25,1,'first')) ...
        t1(find(Smean>=0.50,1,'first')) ...
        t1(find(Smean>=0.75,1,'first')) ...
        t1(find(Smean>=1.0 ,1,'first')) ...
        t1(find(Smean>=1.25,1,'first')) ...
        t1(find(Smean>=1.5 ,1,'first')) ...
        t1(find(Smean>=1.75,1,'first')) ...
        t1(find(Smean>=2.0 ,1,'first')) ...
        t1(find(Smean>=2.5 ,1,'first')) ...
        t1(find(Smean>=3.0 ,1,'first')) ];  % results will always be saved at the end

    t_output_times(t_output_times<t0+1)=[];
    t_output_times(t_output_times>tstop)=[];
    t_output_times=sort(unique(t_output_times));

    % preallocate space to save model output results

    len_interval=length(t_output_times)+1;

    SAVE_ASSIMILATION_ADJUSTMENTS=0;
    if SAVE_ASSIMILATION_ADJUSTMENTS  % This saves the pre- and post-assimilation model states. It may be better to turn this off for large "production" runs
        len_interval=len_interval+2*length(t_obs);
    end

    t_output=NaN(1,len_interval); t_output(1,1)=t1(1);

    YY=NaN(Ntr,len_interval); YY(:,1)=quantile(Y,0.5,2); % ensemble median

    YST=NaN(Ntr,len_interval); YST(:,1)=quantile(Yst,0.5,2);    % ensemble median of short-term shoreline position
    YLST=NaN(Ntr,len_interval); YLST(:,1)=quantile(Ylst,0.5,2); % ensemble median of longshore transport term
    YBRU=NaN(Ntr,len_interval); YBRU(:,1)=quantile(Ybru,0.5,2); % ensemble median of Bruunian recession term
    YVLT=NaN(Ntr,len_interval); YVLT(:,1)=quantile(Yvlt,0.5,2); % ensemble median of long-term rate term

    YLST_only=NaN(Ntr,len_interval); YLST_only(:,1)=quantile(Ylst,0.5,2); % ensemble median of longshore transport only term

    YCI1   =NaN(Ntr,len_interval); YCI1(:,1)   =quantile(Y   ,0.025,2)-0.0001;  % 95% confidence bands
    YCI2   =NaN(Ntr,len_interval); YCI2(:,1)   =quantile(Y   ,0.975,2)+0.0001;
    YSTCI1 =NaN(Ntr,len_interval); YSTCI1(:,1) =quantile(Yst ,0.025,2)-0.0001;  % 95% confidence bands
    YSTCI2 =NaN(Ntr,len_interval); YSTCI2(:,1) =quantile(Yst ,0.975,2)+0.0001;
    YLSTCI1=NaN(Ntr,len_interval); YLSTCI1(:,1)=quantile(Ylst,0.025,2)-0.0001;  % 95% confidence bands
    YLSTCI2=NaN(Ntr,len_interval); YLSTCI2(:,1)=quantile(Ylst,0.975,2)+0.0001;
    YBRUCI1=NaN(Ntr,len_interval); YBRUCI1(:,1)=quantile(Ybru,0.025,2)-0.0001;  % 95% confidence bands
    YBRUCI2=NaN(Ntr,len_interval); YBRUCI2(:,1)=quantile(Ybru,0.975,2)+0.0001;
    YVLTCI1=NaN(Ntr,len_interval); YVLTCI1(:,1)=quantile(Yvlt,0.025,2)-0.0001;  % 95% confidence bands
    YVLTCI2=NaN(Ntr,len_interval); YVLTCI2(:,1)=quantile(Yvlt,0.975,2)+0.0001;

    YLST_only_CI1=NaN(Ntr,len_interval); YLST_only_CI1(:,1)=quantile(Ylst_only,0.025,2)-0.0001;  % 95% confidence bands
    YLST_only_CI2=NaN(Ntr,len_interval); YLST_only_CI2(:,1)=quantile(Ylst_only,0.975,2)+0.0001;

    DTOUT=NaN(Ntr,len_interval); DTOUT(:,1)=quantile(DT,0.5,2);
    DTCI1=NaN(Ntr,len_interval); DTCI1(:,1)=quantile(DT,0.025,2)-1e-6; % 95% confidence bands
    DTCI2=NaN(Ntr,len_interval); DTCI2(:,1)=quantile(DT,0.975,2)+1e-6;

    DYOUT=NaN(Ntr,len_interval); DYOUT(:,1)=quantile(DY,0.5,2);
    DYCI1=NaN(Ntr,len_interval); DYCI1(:,1)=quantile(DY,0.025,2)-1e-6; % 95% confidence bands
    DYCI2=NaN(Ntr,len_interval); DYCI2(:,1)=quantile(DY,0.975,2)+1e-6;

    HSBOUT=NaN(Ntr,len_interval); HSBOUT(:,1)=quantile(HSB,0.5,2);
    HSBCI1=NaN(Ntr,len_interval); HSBCI1(:,1)=quantile(HSB,0.025,2)-1e-6; % 95% confidence bands
    HSBCI2=NaN(Ntr,len_interval); HSBCI2(:,1)=quantile(HSB,0.975,2)+1e-6;

    COUT=NaN(Ntr,len_interval); COUT(:,1)=quantile(c,0.5,2);
    CCI1=NaN(Ntr,len_interval); CCI1(:,1)=quantile(c,0.025,2)-1e-6; % 95% confidence bands
    CCI2=NaN(Ntr,len_interval); CCI2(:,1)=quantile(c,0.975,2)+1e-6;

    KOUT=NaN(Ntr,len_interval); KOUT(:,1)=quantile(K,0.5,2);
    KCI1=NaN(Ntr,len_interval); KCI1(:,1)=quantile(K,0.025,2)-1e-6; % 95% confidence bands
    KCI2=NaN(Ntr,len_interval); KCI2(:,1)=quantile(K,0.975,2)+1e-6;

    VLTOUT=NaN(Ntr,len_interval); VLTOUT(:,1)=quantile(vlt,0.5,2);
    VLTCI1=NaN(Ntr,len_interval); VLTCI1(:,1)=quantile(vlt,0.025,2)-1e-6; % 95% confidence bands
    VLTCI2=NaN(Ntr,len_interval); VLTCI2(:,1)=quantile(vlt,0.975,2)+1e-6;

    SIGMAOUT=NaN(Ntr,len_interval); SIGMAOUT(:,1)=quantile(sigma,0.5,2);
    SIGMACI1=NaN(Ntr,len_interval); SIGMACI1(:,1)=quantile(sigma,0.025,2)-1e-6; % 95% confidence bands
    SIGMACI2=NaN(Ntr,len_interval); SIGMACI2(:,1)=quantile(sigma,0.975,2)+1e-6;

    count_output=2;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up figure for plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PLOT=0;             % turn on plotting as simulation progresses?
PLOT_INTERVAL=28;   % the time interval (in days) for plotting

PLOT_SECTIONS=1; % plot the alongshore results in littoral cell sections? ... not for animation, but for the output

if PLOT
    % setup some figures for plotting
    figure(1); set(gcf,'PaperPositionMode','Auto','Position',[50 100 920 1000],'color','w');
    figure(2); set(gcf,'PaperPositionMode','Auto','Position',[985 100 920 1000],'color','w');
    figure(3); set(gcf,'PaperPositionMode','Auto','Position',[50 100 920 1000],'color','w');
    figure(4); set(gcf,'PaperPositionMode','Auto','Position',[100 100 920 700],'color','w');

    MAKE_ANIMATION=1;
    if MAKE_ANIMATION % make gif animation
        filename1=['CoSMoS_COAST_',Model_name,'_components.gif'];
        filename2=['CoSMoS_COAST_',Model_name,'_parameters.gif'];
        filename3=['CoSMoS_COAST_',Model_name,'_alongshore.gif'];
        count=1;
        count2=1;
        fps=35;
        fps2=24;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN simulation time-stepping loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

globalTimer=tic; % start global timer to record full simulation time

% main time stepping loop
for n=1:len-1 % for all time steps

    tic % start timer for each time step

    fprintf('COSMOS-COAST: timestep %s. %0.1f%% finished. ',datestr(t1(n)),n/len*100); % display progress

    % get year, month, and day of current timestep
    [YYYY,MM,DD]=datevec(t1(n));

    % get wave conditions at current time step
    Hs_now=repmat(Hs(:,n),1,Nens);
    Tp_now=repmat(Tp(:,n),1,Nens);
    Dir_now=repmat(Dir(:,n),1,Nens);

    % calculate wave height and direction at midpoints (by averaging adjacent transects)
    Hs_now_mid=0.5*(Hs_now(1:end-1,:)+Hs_now(2:end,:));
    Tp_now_mid=0.5*(Tp_now(1:end-1,:)+Tp_now(2:end,:));
    Dir_now_mid=0.5*(Dir_now(1:end-1,:)+Dir_now(2:end,:));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% The shorline model loop (the forward model)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % longshore transport term

    for nn=1:Nsubcycles  % for all subcycles

        % the long-term shoreline position (i.e., the sum of the longshore, Bruunian, long-term rate, and long-term anomaly components)
        Ylt=Ylst+Ybru+Yvlt;

        % the full shoreline position
        if ShoreFor_only
            Y=Yst+Y0; % Mao changes this to compare with ShoreFor
        else
            Y=Ylt+Yst+Y0;
        end

        x=repmat(x_on,1,Nens)+Y.*cosd(repmat(phi,1,Nens)); % the x-coordinate of the shoreline
        y=repmat(y_on,1,Nens)+Y.*sind(repmat(phi,1,Nens)); % the y-coordinate of the shoreline

        % get shoreline and wave angles
        shoreline_angle_mid=atan2d(diff(y),diff(x)); % the shoreline angle (Cartesian convention) at the transect midpoints
        alpha_shoreline_mid=shoreline_angle_mid-90;  % the shore-normal angle, pointing towards onshore location (-90 for West Coast setups when numbered south to north, +90 for East Coast setups)
        alpha_wave_mid=270-Dir_now_mid;
        alpha_mid=alpha_wave_mid-alpha_shoreline_mid;    % the wave angle relative to the shoreline angle for each midpoint (alpha_wave_mid-alpha_shoreline_mid for West Coast setups, alpha_shoreline_mid-alpha_wave_mid for East Coast setups)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ADDITIONAL_WAVE_REFRACTION=0; % apply additional refraction/diffraction calculations from the the loaded wave conditions to each transect?
        if ADDITIONAL_WAVE_REFRACTION

            % a simple refraction model (based on linear wave theory and preconstructed lookup-table-like functions)
            FAC1=real(0.543.*Hs_now_mid.^-0.398.*Tp_now_mid.^ 0.798.*(1+0.0002*alpha_mid.^2)); FAC1(FAC1<1.0)=1.0; FAC1(FAC1>12.42)=12.42;
            FAC2=real(1.318.*Hs_now_mid.^ 0.178.*Tp_now_mid.^-0.347.*(1+0.0002*alpha_mid.^2)); FAC2(FAC2<0.5)=0.5; FAC2(FAC2>1.557)=1.557;

            % refracted wave conditions
            alpha_mid=alpha_mid./FAC1;
            Hs_now_mid=Hs_now_mid./FAC2;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % alongshore smoothing of wave conditions
            SMOOTHING=1;
            if SMOOTHING
                SMOOTH_MAG=10;
                for i=1:N_littoral_cells % smooth wave height and wave angle along each littoral cell
                    % get id's of each cell
                    id=(littoral_cell_tr_start(i):littoral_cell_tr_end(i))'; id=min(id,Ntr-1);
                    id(isnan(alpha_mid(id,1)) | isnan(Hs_now_mid(id,1)))=[];
                    if length(id)>6
                        alpha_mid(id,:) =filtfilt(b_smooth,a_smooth,alpha_mid(id,:) );
                        Hs_now_mid(id,:) =filtfilt(b_smooth,a_smooth,Hs_now_mid(id,:) );
                        %alpha_mid(id,:) =smoothn(alpha_mid(id,:));
                        %Hs_now_mid(id,:)=smoothn(Hs_now_mid(id,:));
                        %alpha_mid(id,:) =smoothn(alpha_mid(id,:),SMOOTH_MAG);
                        %Hs_now_mid(id,:)=smoothn(Hs_now_mid(id,:),SMOOTH_MAG);
                    end
                end
            end

            % apply taper function to account for refraction/diffraction at the
            % end of littoral cells
            alpha_mid=alpha_mid.*repmat(TAPER_mid,1,Nens);

        end

        % remove NaN's
        alpha_mid(isnan(alpha_mid))=0;

        % recenter alphas
        alpha_mid(alpha_mid>180)=alpha_mid(alpha_mid>180)-360;

        % supress non-oblique wave approaches (relative wave angles greater than 90 deg)
        alpha_mid(abs(alpha_mid)>90)=0;

        % if angle exceeds critical angle, set the angle to the critical angle
        alpha_mid(alpha_mid>40)=40;
        alpha_mid(alpha_mid<-40)=-40;

        % Longshore transport limiter
        fac=0.25;
        p=10;
        limiter=@(y) tanh(fac*(y)).^p.*(y>0);
        %limiter=@(y) 1+0*y; % turn off limiter

        % calculate beach width to limit longshore transport
        beach_width=Y-YMIN; beach_width(isnan(beach_width))=Y(isnan(beach_width)); beach_width(repmat(~bool_full_model,1,Nens))=0;
        beach_width_mid=(alpha_mid>0).*beach_width(1:end-1,:)+(alpha_mid<0).*beach_width(2:end,:)+(alpha_mid==0).*(0.5*(beach_width(2:end,:)+beach_width(1:end-1,:)));

        K_mid=(alpha_mid>0).*K(1:end-1,:)+(alpha_mid<0).*K(2:end,:)+(alpha_mid==0).*(0.5*(K(2:end,:)+K(1:end-1,:)));

        % if any K's are NaN's they get set to zero
        K_mid(isnan(K_mid))=0;

        % calculate LST
        if HOLD_THE_LINE
            Q=K_mid.*Hs_now_mid.^2.*sind(2*alpha_mid).*limiter(beach_width_mid);  % calculate transport at transect midpoints, which is reduced when beach width becomes very narrow
        else
            Q=K_mid.*Hs_now_mid.^2.*sind(2*alpha_mid);                            % calculate transport at transect midpoints
        end

        % Boundary conditions on the longshore transport
        %Q=cat(1,zeros(1,Nens),Q,zeros(1,Nens)); % no flux bc's
        Q=cat(1,Q(1,:),Q,zeros(1,Nens)); % no gradient bc on south end

        % THIS PROBABLY WONT WORK AS IS ... NEED TO SET Q TO ZERO AT
        % MIDPOINTS BETWEEN LITTORAL CELLS

        % deal with surpressing Q for all non-'full model' transects
        % this also effectively sets boundary conditions between littoral cells
        Q(repmat(~bool_full_model_mid,1,Nens))=0;

        % update the long-term shoreline position
        Ylst=Ylst-dt/Nsubcycles*(1./DC.*diff(Q)./DX);
        Ylst_only=Ylst_only-dt/Nsubcycles*(1./DC.*diff(Q)./DX);
    end

    % STOP IF WE FIND NaN's
    if any(isnan(Q(:)))
        stop
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % short-term model
    Yeq=-DY.*((Hs_now.^2-HSB.^2)./HSB.^2); % the equlibirum (short-term) shoreline position
    tau=DT.*(HSB./Hs_now);                 % the equlibirum adjustment time scale
    Yst=Yst+dt*(Yeq-Yst)./tau;             % advance the short-term shoreline position

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Bruunian recession model
    S1=interp1(t_SLR',S',t1(n+1))';
    S0=interp1(t_SLR',S',t1(n))';
    dSdt=(S1-S0)/dt;                                                            % the rate of SLR at the current time step
    Ybru=Ybru-dt*c./repmat(tanBeta_Bruun,1,Nens).*repmat(dSdt(SLRID),1,Nens);   % the Bruunian recession position

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % long-term rate
    Yvlt=Yvlt+dt*vlt/365.25;           % convert vlt in m/yr to m/day

    % if any component is NaN, then it gets set to zero, because that
    % process is not included in the model for that particular transect
    Yst (isnan(Yst ))=0;
    Ybru(isnan(Ybru))=0;
    Yvlt(isnan(Yvlt))=0;

    % the long-term shoreline position (i.e., the sum of the longshore, Bruunian, and long-term rate components)
    Ylt=Ylst+Ybru+Yvlt;

    % the full shoreline position
    if ShoreFor_only
        Y=Yst+Y0; % Mao changes this to compare with ShoreFor
    else
        Y=Ylt+Yst+Y0;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Data assimilation loop (the inverse model)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if DATA_ASSIMILATION   % if data assimilation is performed

        [DATA_FOUND,nn]=ismember(t1(n+1),t_obs); % determine if the current time step is a data assimiation step

        if  DATA_FOUND % if the current time step is in t_obs (the set of observation times), then we do data assimilation

            if SAVE_DATA
                if SAVE_ASSIMILATION_ADJUSTMENTS
                    % save output pre assimilation
                    increment=-1/24;
                    save_CoSMoS_COAST_output
                end
            end

            % find the id's of transects with data for the current assimilation step
            id_tr_obs=find(~isnan(Y_obs(:,nn)));

            % the number of observations
            Ntr_obs=length(id_tr_obs);

            % get the names of all unique littoral cells that have observations at the current time step
            [littoral_cells_with_obs,id]=unique(littoral_cell_names(id_tr_obs));

            % get the littoral cell types of all unique littoral cells that have observations at the current time step
            littoral_cells_with_obs_type=littoral_cell_type_all(id_tr_obs(id));

            % The number of unique littoral cells that have observations at the current time step
            N_littoral_cells_with_obs=length(littoral_cells_with_obs);

            % initialize the post-assimilation model state variables (which
            % are initialized simply with the pre-assimilation values), but
            % the post-assimilation state variables have an '_new' post-fix
            for ii=1:Nvars                                     % for each assimilated variable ...
                eval([da_vars{ii},'_new=',da_vars{ii},';']);   % create a new, post-assimilation variable with an '_new' post-fix
            end
            % example: Ylst_new=Ylst;

            % for all littoral cells with data
            for i=1:N_littoral_cells_with_obs  % loop over all unique littoral cells

                if ~DATA_ASSIMILATION_FOR_RATE_ONLY_TRANSECTS && strcmp(littoral_cells_with_obs_type{i},'rate only') % don't do data assimilation for 'rate only' model types if it's turned off
                    continue
                end

                % get the id's of the transects with observations only for each unique littoral cell
                id_tr_obs_cell=id_tr_obs(strcmp(littoral_cells_with_obs{i},littoral_cell_names(id_tr_obs)));

                % get the number of transects with observations on a given littoral cell
                Ntr_obs_cell=length(id_tr_obs_cell);

                % get the id's of the transects to assimilate.  i.e., all transects within the same littoral cell
                % (Note: we only assimilate data one littoral cell at a time to improve scalability)
                id_tr_assim=find(strcmp(littoral_cells_with_obs{i},littoral_cell_names));

                % get the number of transects that are assimilated on a given littoral cell (i.e., all of the transects on that littoral cell)
                Ntr_assim=length(id_tr_assim);

                % do data assimilation on each littoral cell seperately
                y_obs=Y_obs(id_tr_obs_cell,nn)-Y0(id_tr_obs_cell); % the observed shoreline positions on the transects with data (minus the initial shoreline position)
                y_rms=nanmean(Y_rms(id_tr_obs_cell,nn));           % this kindof assumes that on a given day there are not BOTH satellite data and lidar/gps data for a given littoral cell

                Y_mod=Y(id_tr_obs_cell,:)-repmat(Y0(id_tr_obs_cell),1,Nens);  % the (ensemble of) modeled  shoreline positions on the transects with data (minus the initial shoreline position)
                y_mod_bar=nanmean(Y_mod,2);                                   % the ensemble mean of the modeled shoreline positions on the transects with data

                % (inverse) error covariance matrix of observations
                Rloc=y_rms^2*R(id_tr_obs_cell,id_tr_obs_cell);

                % mistmatch / error between observerations and model
                % Y_error=repmat(y_obs,1,Nens)-Y_mod;   % error across the ensemble
                Y_error=mvnrnd(y_obs,Rloc,Nens)'-Y_mod; % error across the ensemble w/ perturbed observations
                y_error=nanmean(Y_error,2);             % mean error across the ensemble

                % initialize the modified model state variables (which
                % are initialized as a portion of the full state vector for a given littoral cell), but
                % the post-assimilation state variables have an '1' post-fix
                for ii=1:Nvars                   % for each assimilated variable ...
                    if ~da_var_is_positive(ii)   % if variable is not positive ...
                        eval([da_vars{ii},'1=',da_vars{ii},'(id_tr_assim,:);']); % extract portion of that (unmodified) assimilated state/paramter for a given littoral cell
                    else
                        eval([da_vars{ii},'1=1/sig*log(max(',da_vars{ii},'(id_tr_assim,:),lb_',da_vars{ii},')./repmat(',da_vars{ii},'0(id_tr_assim,:),1,Nens));']); % extract portion of that (modified) assimilated state/paramter for a given littoral cell
                    end
                end
                % example: Ylst1=Ylst(id_tr_assim,:); (for unmodified state/parameter)
                %     or   DT1  =1/sig*log(max(DT(id_tr_assim,:)   ,lb_DT     )./repmat(   DT0(id_tr_assim,:),1,Nens)); (for modified state/parameter)

                % assemble (modified) paremeters into state vector for data assimilation (we use the full model state for all transects)
                E=NaN(Nvars*Ntr_assim,Nens);  % initialize full ensemble state vector E
                if SYMBOLIC_MODE
                    for k=1:Nens
                        E(:,k)=reshape(eval(state_vector_str)',[],1);  % group all of the (semi-symbolic) transect variables (i.e., the states and parameters) togeather because of the localization part of the Kalman filter operation
                    end
                else
                    for k=1:Nens
                        E(:,k)=reshape([Ylst1(:,k) Yst1(:,k) DT1(:,k) DY1(:,k) HSB1(:,k) c1(:,k) vlt1(:,k) K1(:,k) sigma1(:,k)]',[],1);  % group all of the transect variables (i.e., the states and parameters) togeather because of the localization part of the Kalman filter operation
                    end
                end

                if ~isreal(E)
                    error('Imaginary numbers found in the model state vector!?! SeanV.');
                end

                % perform data assimilation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % calculate ensemble mean (x_bar)
                x_bar=mean(E,2);

                % calculate the matrix of normalized anomalies (X).  Note that the model covariance matrix P = X*X^T.
                X=E-repmat(x_bar,1,Nens);

                % calculate the deviations from the observations in ensemble space (Y_obs)
                Y_anom=Y_mod-repmat(y_mod_bar,1,Nens);

                % Ensemble Kalman Filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % calculate the error covariance matrix
                rho=1.1; % covariance inflation factor
                P=rho./(Nens-1)*Y_anom*Y_anom'+Rloc;

                % calculate the Kalman gain
                %KAL=1./(Nens-1)*X*Y_anom'*inv(P);  % the old way (slow)
                KAL=1./(Nens-1)*X*(Y_anom'/P);      % the new way (fast)

                KAL(isnan(KAL))=0;                  % remove the NaN values from KAL, NaNs in the analysis ensemble will still be retained

                % calculate the observation perturbations
                y_perturb=mvnrnd(zeros(length(Rloc),1),Rloc,Nens)';
                y_perturb=y_perturb-repmat(nanmean(y_perturb,2),1,Nens); % we just remove the mean for good measure

                % tile the localization matrix (LOC) for the transects with observations for Nvars variables
                LOC_KRON=kron(LOC(id_tr_assim,id_tr_obs_cell),ones(Nvars,1));

                % calculate the ensemble analysis step
                Enew=E+(LOC_KRON.*KAL)*(repmat(y_obs,1,Nens)+y_perturb-Y_mod);

                % convert the state vector ensemble back into individual
                % state / parameter vectors (and apply positive variable modification/conversion if necessary)
                for ii=1:Nvars                    % for each assimilated variable ...
                    if ~da_var_is_positive(ii)    % if variable is not positive ...
                        eval([da_vars{ii},'_new(id_tr_assim,:)=Enew(',num2str(ii,'%d'),':Nvars:end,:);']); % convert Enew back into state/parameters
                    else
                        eval([da_vars{ii},'_new(id_tr_assim,:)=repmat(',da_vars{ii},'0(id_tr_assim,:),1,Nens).*exp(sig*Enew(',num2str(ii,'%d'),':Nvars:end,:));']); % else convert Enew back into modified state/parameters
                    end
                end
                % example: Ylst_new(id_tr_assim,:)=Enew(1:Nvars:end,:);                                               (standard)
                %    or    DT_new(id_tr_assim,:)  =repmat(DT0(id_tr_assim,:),1,Nens).*exp(sig*Enew(3:Nvars:end,:));   (if parameter is positive)

            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % alongshore smoothing
            SMOOTHING=1;
            if SMOOTHING

                % smooth for each littoral cell that was assimilated
                for i=1:N_littoral_cells_with_obs

                    % get the id's of the transects that were assimilated.  i.e., all transects within the same littoral cell
                    % (Note: we only assimilate data one littoral cell at a time to improve scalability)
                    id_tr_assim=find(strcmp(littoral_cells_with_obs{i},littoral_cell_names));

                    % only smooth if the littoral cell has more than 6 transects
                    if length(id_tr_assim)>6

                        if strcmp(littoral_cells_with_obs{i},'cliff only') || strcmp(littoral_cells_with_obs{i},'no prediction')
                            continue % no smoothing for 'cliff only' and 'no prediction' cells
                        end

                        if SYMBOLIC_MODE % smoothing using the semi-symbolic variables that are assimilated

                            id_smooth_vars=find(da_var_smooth);                                        % for 'full model' cells keep all variables needing smoothing
                            if strcmp(littoral_cells_with_obs{i},'cross-shore only')                   % for 'cross-shore only' cells
                                id_smooth_vars(strcmp('K',da_vars(id_smooth_vars)))=[];                % remove 'K' from needed any smoothing
                            elseif strcmp(littoral_cells_with_obs{i},'rate only')                      % for 'rate only' cells
                                if DATA_ASSIMILATION_FOR_RATE_ONLY_TRANSECTS
                                    [~,locb]=ismember({'DT','DY','HSB','K'},da_vars(id_smooth_vars));                    % remove 'DT','DY','HSB', and 'K' from needed any smoothing
                                else
                                    [~,locb]=ismember({'DT','DY','HSB','K','vlt','c','sigma'},da_vars(id_smooth_vars));  % remove 'DT','DY','HSB', 'K', 'vlt', 'c', and 'sigma' from needed any smoothing
                                end
                                id_smooth_vars(locb)=[];
                            end

                            for ii=id_smooth_vars  % for each (remaining) assimilated variable needing smoothing...
                                eval([da_vars{ii},'_new(id_tr_assim,:)=da_var_alpha(',num2str(ii,'%d'),')*movmean(',da_vars{ii},'_new(id_tr_assim,:),5,1)+(1-da_var_alpha(',num2str(ii,'%d'),'))*',da_vars{ii},'_new(id_tr_assim,:);']);
                            end
                            % example: DT_new(id,:) =da_var_alpha(ii)*movmean(DT_new(id,:),5,1)+(1-da_var_alpha)*DT_new(id,:);

                        else % traditional/manual smoothing NOT using semi-symbolic variables (i.e., not using 'eval's)

                            % Smoothing for short-term and long-term parameters
                            if (strcmp(littoral_cells_with_obs{i},'cross-shore only') || strcmp(littoral_cells_with_obs{i},'full model'))

                                %DT_new(id_tr_assim,:)   =filtfilt(b_smooth,a_smooth,DT_new(id_tr_assim,:)   );
                                %DY_new(id_tr_assim,:)   =filtfilt(b_smooth,a_smooth,DY_new(id_tr_assim,:)   );
                                %HSB_new(id_tr_assim,:)  =filtfilt(b_smooth,a_smooth,HSB_new(id_tr_assim,:)  );

                                DT_new(id_tr_assim,:) =da_var_alpha(strcmp(da_vars,'DT'))*movmean(DT_new(id_tr_assim,:) ,5,1)+(1-da_var_alpha(strcmp(da_vars,'DT')))*DT_new(id_tr_assim,:);
                                DY_new(id_tr_assim,:) =da_var_alpha(strcmp(da_vars,'DY'))*movmean(DY_new(id_tr_assim,:) ,5,1)+(1-da_var_alpha(strcmp(da_vars,'DY')))*DY_new(id_tr_assim,:);
                                HSB_new(id_tr_assim,:)=da_var_alpha(strcmp(da_vars,'HSB'))*movmean(HSB_new(id_tr_assim,:),5,1)+(1-da_var_alpha(strcmp(da_vars,'HSB')))*HSB_new(id_tr_assim,:);

                            end

                            % Smoothing for short-term and long-term parameters
                            if strcmp(littoral_cells_with_obs{i},'full model')

                                %K_new(id_tr_assim,:)    =filtfilt(b_smooth,a_smooth,K_new(id_tr_assim,:)    );

                                K_new(id_tr_assim,:)=da_var_alpha(strcmp(da_vars,'K'))*movmean(K_new(id_tr_assim,:),5,1)+(1-da_var_alpha(strcmp(da_vars,'K')))*K_new(id_tr_assim,:);

                            end

                            if ~DATA_ASSIMILATION_FOR_RATE_ONLY_TRANSECTS && strcmp(littoral_cells_with_obs_type{i},'rate only') % don't do smoothing for 'rate only' model types if it's turned off
                                continue
                            end

                            % assimilate the rate, Bruun, and noise variables
                            vlt_new(id_tr_assim,:)  =da_var_alpha(strcmp(da_vars,'vlt'  ))*movmean(  vlt_new(id_tr_assim,:),5,1)+(1-da_var_alpha(strcmp(da_vars,'vlt'  )))*  vlt_new(id_tr_assim,:);
                            c_new(id_tr_assim,:)    =da_var_alpha(strcmp(da_vars,'c'    ))*movmean(    c_new(id_tr_assim,:),5,1)+(1-da_var_alpha(strcmp(da_vars,'c'    )))*    c_new(id_tr_assim,:);
                            sigma_new(id_tr_assim,:)=da_var_alpha(strcmp(da_vars,'sigma'))*movmean(sigma_new(id_tr_assim,:),5,1)+(1-da_var_alpha(strcmp(da_vars,'sigma')))*sigma_new(id_tr_assim,:);

                        end

                    end
                end
            end

            % check if upper or lower bounds are exceeded, and clip/resample if they are
            BOUND=1;
            if BOUND

                % only do bounding of DT, DY, and HSB for 'full model' or 'cross-shore only' transects
                id=bool_full_model | bool_cross_shore_only;

                for ii=1:Nvars            % for each assimilated variable ...
                    if da_var_bound(ii)   % ... if bounding is turned on
                        if strcmp(da_vars{ii},'DT') || strcmp(da_vars{ii},'DY') || strcmp(da_vars{ii},'HSB')
                            eval([da_vars{ii},'_new(id,:)=clip2bounds(',da_vars{ii},'_new(id,:),lb_',da_vars{ii},',ub_',da_vars{ii},');']);
                        elseif strcmp(da_vars{ii},'K')
                            K_new(K_new>ub_K)=ub_K; K_new(K_new<lb_K)=lb_K;
                        else
                            eval([da_vars{ii},'_new=clip2bounds(',da_vars{ii},'_new,lb_',da_vars{ii},',ub_',da_vars{ii},');']);
                        end
                    end
                end
                % for example: DT_new(id,:)=clip2bounds(DT_new(id,:) ,lb_DT ,ub_DT );
                %     or       vlt_new     =clip2bounds(vlt_new      ,lb_vlt,ub_vlt);
            end

            % re-calculate the long-term shoreline position (i.e., the sum of the longshore, Bruunian, and long-term rate components)
            if any(strcmp(da_vars,'Ylst'))    % if 'Ylst' variable was assimilated,
                Ylt_new=Ylst_new+Ybru+Yvlt;   % account for the assimilated variable (i.e., 'Ylst_new') in 'Ylt_new'
            else                              % if 'Ylst' variable was not assimilated,
                Ylt_new=Ylst+Ybru+Yvlt;       % then we just use the old value for 'Ylst'
            end

            % recalculate the full shoreline position
            if any(strcmp(da_vars,'Yst'))   % if 'Yst' variable was assimilated,
                Y_new=Ylt_new+Yst_new+Y0;   % account for the assimilated variable (i.e., 'Yst_new') in 'Ylt_new'
            else                            % if 'Yst' variable was not assimilated,
                Y_new=Ylt_new+Yst+Y0;       % then we just use the old value for 'Yst'
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % re-insert the assimilated variables
            Y   =Y_new;
            Ylt =Ylt_new;

            for ii=1:Nvars                                     % for each assimilated variable ...
                eval([da_vars{ii},'=',da_vars{ii},'_new;']);   % reinsert the new post-assimilationvariable back into the prior-state variable
            end
            % example: Ylst=Ylst_new;
            %    or    DT  =DT_new  ;

            if SAVE_DATA
                if SAVE_ASSIMILATION_ADJUSTMENTS
                    % save output post assimilation
                    increment=+1/24;
                    save_CoSMoS_COAST_output
                end
            end

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % alongshore smoothing
    LST_SMOOTHING=1;
    if LST_SMOOTHING

        % smooth for each littoral cell
        for i=1:N_littoral_cells

            % for 'full model' littoral cells
            if strcmp(littoral_cell_type{i},'full model')

                % get id's of each cell
                id=(littoral_cell_tr_start(i):littoral_cell_tr_end(i))';

                % smoothing paramter (larger alpha means more smoothing)
                alpha=0.01;
                if length(id)>=3; % if
                    % do a convex combination-type smoothing
                    Ylst(id(1),:)=(1-alpha)*Ylst(id(1),:)+alpha*Ylst(id(2),:);
                    Ylst(id(end),:)=(1-alpha)*Ylst(id(end),:)+alpha*Ylst(id(end-1),:);
                    Ylst(id(2:end-1),:)=alpha*Ylst(id(1:end-2),:)+(1-2*alpha)*Ylst(id(2:end-1),:)+alpha*Ylst(id(3:end),:);
                end
            end

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Additive noise loop
    if ADDITIVE_NOISE

        NOISE=sigma;

        % In general, we add noise to assimilated variables, which
        % counteracts the data assimilation step where noise/variance of
        % the variable is surpressed.  Getting a balance is tricky, but
        % that is the goal to get a "good" (i.e., converged and accurate)
        % parameter estimate. Sadly, this is more of an 'art' than a
        % 'science' at this point ...

        % add noise Yst in full model or cross-shore only transects ...
        % this is always done, even if data assimilation is turned off.
        % The balance between noise and damping sets the uncertainty
        % vis-a-vis me et al. 2021 - JGR-ES.

        if da_var_add_noise(strcmp(da_vars,'Yst'))
            id=bool_full_model | bool_cross_shore_only;
            Yst(id,:)  =Yst(id,:)+da_var_noise_fac(strcmp(da_vars,'Yst'))*NOISE(id,:).*randn(sum(id),Nens); % additive noise on the short-term component
        end

        if ADDITIVE_NOISE_LT             % additive noise on the long-term components ... This is turned off when data assimilation is finished

            % we apply different rules for adding noise to each assimilated variable...

            % for the equlibirum components
            id=bool_full_model | bool_cross_shore_only;
            if SYMBOLIC_MODE
                [~,locb]=ismember({'DT','DY','HSB'},da_vars);
                for ii=locb              % for each assimilated variable ...
                    if da_var_add_noise(ii)  % ... if additive nose is turned on
                        eval([da_vars{ii},'(id,:)=',da_vars{ii},'(id,:)+da_var_noise_fac(ii)*repmat(NOISE_FAC(id,:),1,Nens).*NOISE(id,:).*randn(sum(id),Nens);']);
                    end
                end
                %for example: DT(id,:) =DT(id,:)+da_var_noise_fac(ii)*repmat(NOISE_FAC(id,:),1,Nens).*NOISE(id,:).*randn(sum(id),Nens);
            else
                if da_var_add_noise(strcmp(da_vars,'DT'));  DT(id,:) =DT(id,:) +da_var_noise_fac(strcmp(da_vars,'DT'))*repmat(NOISE_FAC(id,:),1,Nens).*NOISE(id,:).*randn(sum(id),Nens);  end
                if da_var_add_noise(strcmp(da_vars,'DY'));  DY(id,:) =DY(id,:) +da_var_noise_fac(strcmp(da_vars,'DY'))*repmat(NOISE_FAC(id,:),1,Nens).*NOISE(id,:).*randn(sum(id),Nens);  end
                if da_var_add_noise(strcmp(da_vars,'HSB')); HSB(id,:)=HSB(id,:)+da_var_noise_fac(strcmp(da_vars,'HSB'))*repmat(NOISE_FAC(id,:),1,Nens).*NOISE(id,:).*randn(sum(id),Nens); end
            end

            % add noise to Ylst for the full model transects ... we don't need any 'evals' here since there is only one, known variable (i.e., Ylst) that we must deal with here
            if da_var_add_noise(strcmp(da_vars,'Ylst'))
                for i=1:N_full_model_sections
                    id=id_full_model_start(i):id_full_model_end(i);
                    Ylst(id,:)=Ylst(id,:)+repmat(NOISE_FAC_YLT(id,:),1,Nens).*(da_var_noise_fac(strcmp(da_vars,'Ylst'))*NOISE(id,:).*repmat(randn(1,Nens),length(id),1)+0.2*da_var_noise_fac(strcmp(da_vars,'Ylst'))*NOISE(id,:).*randn(length(id),Nens));
                end

                % SHOULD I ADD NOISE to Ylst in cross-shore only transects... or should
                % they stay at zero? LET's TRY IT ...
                for i=1:N_cross_shore_only_sections
                    id=id_cross_shore_only_start(i):id_cross_shore_only_end(i);
                    Ylst(id,:)=Ylst(id,:)+repmat(NOISE_FAC_YLT(id,:),1,Nens).*(da_var_noise_fac(strcmp(da_vars,'Ylst'))*NOISE(id,:).*repmat(randn(1,Nens),length(id),1)+0.2*da_var_noise_fac(strcmp(da_vars,'Ylst'))*NOISE(id,:).*randn(length(id),Nens));
                end
            end

            % add noise to K for the full model transects ... we don't need any 'evals' here since there is only one, known variable (i.e., K) that we must deal with here
            if da_var_add_noise(strcmp(da_vars,'K'))
                for i=1:N_full_model_sections
                    id=id_full_model_start(i):id_full_model_end(i);
                    K(id,:)   =K(id,:)   +da_var_noise_fac(strcmp(da_vars,'K'))*repmat(NOISE_FAC(id,:),1,Nens).*NOISE(id,:).*repmat(randn(1,Nens),length(id),1)+da_var_noise_fac(strcmp(da_vars,'K'))*repmat(NOISE_FAC(id,:),1,Nens).*NOISE(id,:).*randn(length(id),Nens);
                end
            end

            % noise component
            sigma=sigma+da_var_noise_fac(strcmp(da_vars,'sigma'))*repmat(NOISE_FAC,1,Nens).*NOISE.*randn(Ntr,Nens);  % Yo dog, I heard you like noise, so I'm adding noise to your noise

            % for all other components (e.g., vlt, c)
            if SYMBOLIC_MODE
                [~,locb]=ismember({'vlt','c'},da_vars);
                for ii=locb              % for each assimilated variable ...
                    if da_var_add_noise(ii)  % ... if additive nose is turned on
                        eval([da_vars{ii},'=',da_vars{ii},'+da_var_noise_fac(ii)*repmat(NOISE_FAC,1,Nens).*NOISE.*randn(Ntr,Nens);']);
                    end
                end
                %for example: vlt=vlt+da_var_noise_fac(ii)*NOISE.*randn(Ntr,Nens);
                %    or       var=var+da_var_noise_fac(ii)*NOISE.*randn(Ntr,Nens);
            else
                vlt=vlt+da_var_noise_fac(strcmp(da_vars,'vlt'))*repmat(NOISE_FAC,1,Nens).*NOISE.*randn(Ntr,Nens);
                c  =c  +da_var_noise_fac(strcmp(da_vars,'c'  ))*repmat(NOISE_FAC,1,Nens).*NOISE.*randn(Ntr,Nens);
            end

        end
    end

    % TURN OFF DATA ASSIMILATION (and the addative noise) during the forecast period
    if DATA_ASSIMILATION && t1(n+1)>=tforecast
        DATA_ASSIMILATION=0;
        ADDITIVE_NOISE_LT=0;
        PLOT=0;
        PLOT_INTERVAL=30;
        vlt_assim=vlt;
    end

    if ~CONTINUED_ACCRETION && n+1==nforecast2
        vlt(vlt>0)=0;  % suppress accretion in the 'no continued accretion' scenario
    end

    % since we've modified the components from data-assimilation, additive noise, and smoothing, we update the total shoreline position to reflect those changes

    % update the long-term shoreline position (i.e., the sum of the longshore, Bruunian, and long-term rate components)
    Ylt=Ylst+Ybru+Yvlt;

    % update the full shoreline position
    if ShoreFor_only
        Y=Yst+Y0; % Mao changes this to compare with ShoreFor
    else
        Y=Ylt+Yst+Y0;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Model output loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if SAVE_DATA
        % periodically save the model output (i.e., save the model state into 'OUT' variables)
        if ismember(t1(n+1),t_output_times)
            increment=0;
            save_CoSMoS_COAST_output; % save the model data into the matlab output arrays
        end
        if ismember(t1(n+1),t_save)
            %save_CoSMoS_COAST_results; % save the results to .kml
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot stuff
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if PLOT
        if mod(n,PLOT_INTERVAL)==1
            ID2PLOT=500;
            [~,id_tr]=ismember(ID2PLOT,[transects.ID]);
            plot_CoSMoS_COAST_results; % plot the results
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % keep track of model timing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sim_time(mod(n,1000)+1)=toc;
    min_remaining=nanmean(sim_time)*(len-1-n)/60; % calculate minutes remaining in the simulation

    if min_remaining>1
        if min_remaining>60
            fprintf(' (~%0.1f hours remaining)\n',min_remaining/60); % display progress
        else
            fprintf(' (~%0.1f min remaining)\n',min_remaining); % display progress
        end
    else
        fprintf(' (~%0.0f sec remaining)\n',ceil(min_remaining*60)); % display progress
    end

end

total_sim_time=toc(globalTimer);

if total_sim_time>60
    if total_sim_time>(60*60)
        fprintf('total model run time = ~%0.1f hours.\n',total_sim_time/(60*60)); % display total run time
    else
        fprintf('total model run time = ~%0.1f min.\n',total_sim_time/60); % display total run time
    end
else
    fprintf('total model run time = ~%0.1f sec.\n',total_sim_time); % display total run time
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model post-processing (make figs, calculate skill metrics, etc.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save_CoSMoS_COAST_skill_metrics

% plot skill results every so often and for transects with high RMS error
% ids2plot=ID(1:1:Ntr);
% ids2plot=cat(1,ids2plot,ID(find(RMSE_sat>50)));
% ids2plot=sort(unique(ids2plot));

save_CoSMoS_COAST_figures % (for the production runs we will often not save the figures)

% Save model results to .mat files
if SAVE_DATA
    fprintf('saving output file ... ');
    matname_tr=[OUTPUT_DIR,filesep,'mat',filesep,'CoSMoS_COAST_results_',datestr(t1(n+1)+1),'_SL=',num2str(100*Smean(n+1),'%2.0f'),'cm_transects.mat'];
    matname1=[OUTPUT_DIR,filesep,'mat',filesep,'CoSMoS_COAST_results_',datestr(t1(n+1)+1),'_SL=',num2str(100*Smean(n+1),'%2.0f'),'cm_state.mat'];
    matname2=[OUTPUT_DIR,filesep,'mat',filesep,'CoSMoS_COAST_results_',datestr(t1(n+1)+1),'_SL=',num2str(100*Smean(n+1),'%2.0f'),'cm_params.mat'];

    save(matname_tr,'transects','ids','x_on','y_on','x_off','y_off','phi','Ymin','bool_full_model','bool_cross_shore_only','bool_rate_only','bool_cliff_only','bool_no_prediction','bool_full_model_mid','littoral_cell_names','littoral_cell_names_unique','littoral_cell_type','littoral_cell_tr_start','littoral_cell_tr_end','N_full_model_sections','N_cross_shore_only_sections','N_rate_only_sections','id_full_model_start','id_full_model_end','id_cross_shore_only_start','id_cross_shore_only_end','id_rate_only_start','id_rate_only_end','littoral_cell_tr_length','N_littoral_cells','-v7.3');
    save(matname1,'x0','y0','Y00','Ymin','Ntr','ID','HOLD_THE_LINE','CONTINUED_ACCRETION','Model_name','UTMZONE','SL','t1','t0','tstop','tforecast','tforecast2','dt','Nsubcycles','WAVEID','SLRID','y_rms_sat','y_rms_gps','t_obs','Y_obs','Y_rms','SAT','S','t_SLR','Smean','t_output','Y','Ylst','Ylt','Yst','Ybru','Yvlt','Y0','YY','YST','YLST','YBRU','YVLT','YCI1','YCI2','YSTCI1','YSTCI2','YLSTCI1','YLSTCI2','YBRUCI1','YBRUCI2','YVLTCI1','YVLTCI2','YLST_only','YLST_only_CI1','YLST_only_CI2','NOBS','NOBS_cal','NOBS_val','R','MSE','RMSE','RMSM','MAE','SKILL','IDX_AGREEMENT','LAMBDA','NOBS_cal_sat','NOBS_val_sat','NOBS_sat','R_sat','MSE_sat','RMSE_sat','RMSM_sat','MAE_sat','SKILL_sat','IDX_AGREEMENT_sat','LAMBDA_sat','PCT_sat','-v7.3');
    save(matname2,'t_output','DT','DY','HSB','c','K','sigma','vlt','vlt_assim','DTOUT','DTCI1','DTCI2','DYOUT','DYCI1','DYCI2','HSB0','HSBOUT','HSBCI1','HSBCI2','KOUT','KCI1','KCI2','COUT','CCI1','CCI2','VLTOUT','VLTCI1','VLTCI2','LTER','SIGMAOUT','SIGMACI1','SIGMACI2','sigma0','sigma00','da_vars','da_var_is_positive','da_var_smooth','da_var_alpha','da_var_bound','da_var_add_noise','da_var_noise_fac','-v7.3');

    fprintf('done.\n');
end

%% some error analysis
[max_er,id]=sort(RMSE_sat,'descend');

id(isnan(max_er))=[];
max_er(isnan(max_er))=[];


