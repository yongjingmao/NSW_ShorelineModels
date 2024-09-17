% analyze waves data for NSW, Australia. SeanV
%clear all;
close all; clc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Add matlab helper routines toolbox to the path (the kml toolobx is the most important here ... but others might be too)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matlab_helpers_toolbox='C:\Users\z3541792\OneDrive - UNSW\Documents\GitHub\MatlabToolbox';
addpath(genpath(matlab_helpers_toolbox));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the transects data file (which is effectively the model grid/domain)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('transects','var')  % if the transect file does not exist in the workspace, then load it in ...
    % the transects file to load
    TRANSECTS_FILE='..\transects\transects.mat'; % path to the transects folder

    fprintf('loading transects file ...');
    load(TRANSECTS_FILE); % load the transects file ... and display progress
    fprintf(' done. \n');
else
    fprintf('using previously loaded transects file.\n');
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Setup littoral cells (mined from the transects data-struct)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
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

faceColors = makesymbolspec('Polygon',{'INDEX',[1 N_littoral_cells],'FaceColor',polcmap(N_littoral_cells)});
cmap2=faceColors.FaceColor{1,3};


% Number of transects
Ntr=length(transects);

% get transect coordinates
ID=[transects.ID]';
x_on=[transects.x_on]';
y_on=[transects.y_on]';
x_off=[transects.x_off]';
y_off=[transects.y_off]';
phi=[transects.angle]';

% get mean shoreline position
Y0=NaN(Ntr,1);
for i=1:Ntr
    Y0(i,1)=nanmean(transects(i).Y);
end

% plot starting shoreline position
plot(Y0);

% get mean shoreline coordinates
x0=x_on+Y0.*cosd(phi);
y0=y_on+Y0.*sind(phi);

dx=x0(2:end)-x0(1:end-1);
dy=y0(2:end)-y0(1:end-1);

% convert to Lat/Lon from meters UTM zone
UTMZONE='56 H';
[lat0, lon0] = utm2deg(x0,y0,repmat(UTMZONE,Ntr,1));
[lat_on, lon_on] = utm2deg(x_on,y_on,repmat(UTMZONE,Ntr,1));
[lat_off, lon_off] = utm2deg(x_off,y_off,repmat(UTMZONE,Ntr,1));

geoplot(lat0,lon0,'-b.');

shoreline_angle=atan2d(dy,dx)-90;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the waves data file (which is effectively the model grid/domain)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('Hs_mean','var')  % if the transect file does not exist in the workspace, then load it in ...
    % the transects file to load
    WAVES_FILE='waves.mat'; % path to the transects folder

    fprintf('loading waves file ...');
    load(WAVES_FILE); % load the transects file ... and display progress
    fprintf(' done. \n');
else
    fprintf('using previously loaded waves file.\n');
end

figure; set(gcf,'PaperPositionMode','auto','Position',[100 100 1600 800],'color','w');

YMAX_SCALE=1.7;
FACE_ALPHA=0.2;
FONT_SIZE=12;
FONT_SIZE_TEXT=12;

XMIN=-6;
XMAX=6;
XTEXT=2;

for i=1:N_littoral_cells
    fill([XMIN,XMIN,XMAX,XMAX],[littoral_cell_tr_start(i),littoral_cell_tr_end(i),littoral_cell_tr_end(i),littoral_cell_tr_start(i)],cmap2(i,:),'FaceAlpha',FACE_ALPHA,'EdgeColor',cmap2(i,:),'EdgeAlpha',0.25,'HandleVisibility','off');
end

for i=1:N_littoral_cells
    if littoral_cell_tr_length(i)>40
        cellname_1=text(XTEXT,0.5*(littoral_cell_tr_start(i)+littoral_cell_tr_end(i)),littoral_cell_names_unique{i});
        set(cellname_1,'Rotation',0,'FontSize',FONT_SIZE_TEXT);
    end
end

Dir_mean=nanmean(Dir,2);
Dir_mid_cart=270-0.5*(Dir_mean(1:end-1)+Dir_mean(2:end))-360;
ID_mid=0.5*(ID(1:end-1)+ID(2:end));
plot(shoreline_angle,ID_mid,'-r',Dir_mid_cart,ID_mid,'-b'); set(gca,'YDir','reverse')
ylabel('Transect #')
xlabel('wave/shoreline angle (Cartesican convention)');
legend('shoreline angle','wave angle');
