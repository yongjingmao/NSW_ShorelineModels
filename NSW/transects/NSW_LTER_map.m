% Make a nice domain figure with long-term erosion rates. SeanV.
%clear all; 
close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add matlab helper routines toolbox to the path (the kml toolobx is the most important here)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matlab_helpers_toolbox='K:\COSMOS\Matlab';
addpath(genpath(matlab_helpers_toolbox));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the transects data file (which is effectively the model grid/domain)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('transects','var')
    
    % the transects file to load
    TRANSECTS_FILE='transects.mat';

    % UTM zones for NSW Australia
    UTM_ZONE=56;
    UTM_HEM='S';
    UTMZONE='56 H';

    fprintf('loading transects file ...');
    load(TRANSECTS_FILE); % load the transects file ... and display progress
    fprintf(' done. \n');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Setup littoral cells (mined from the transects data-struct)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % get littoral cell information for all transect
    littoral_cell_type_all={transects.model_type}';

    % a function to take the transects and get the littoral cell information. SeanV.
    [littoral_cell_names_unique,littoral_cell_type,...
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

else
    fprintf('using preivously loaded transects file.\n');
end

% number of transects
Ntr=length(transects);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get some information contained within the transects data-struct (this is done to facilitate readability of the code, without constantly refering back to the large transects struct)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_on =[transects.x_on]';  % the x-coordinate of the onshore  baseline
y_on =[transects.y_on]';  % the y-corrdinate of the onshore  baseline
x_off=[transects.x_off]'; % the x-coordinate of the offshore baseline
y_off=[transects.y_off]'; % the y-corrdinate of the offshore baseline

phi=[transects.angle]'; % the transect angle

% convert the x/y-mean positions to lat/lon, since the wave observation points are in lat/lon
[lat_on, lon_on] = utm2deg(x_on,y_on,repmat(UTMZONE,Ntr,1));
[lat_off, lon_off] = utm2deg(x_off,y_off,repmat(UTMZONE,Ntr,1));

ID=[transects.ID];
LTER=[transects.LTER];
STATE={transects.state};

Ntr=length(ID);

id=strcmp(STATE,'NSW');

ID_NSW=ID(id);
lat_on_NSW=lat_on(id);
lon_on_NSW=lon_on(id);
LTER_NSW=LTER(id);

id=find(id);

N_shorelines=NaN(length(id),1);
t_min=NaN(length(id),1);
t_max=NaN(length(id),1);

for i=1:length(id)
    N_shorelines(i)=length(transects(id(i)).Y);
    t_min(i)=min(transects(id(i)).t);
    t_max(i)=max(transects(id(i)).t);
end

Ntr_NSW=length(LTER_NSW);

cmap=flipud(redblue(64));

figure; set(gcf,'PaperPositionMode','auto','Position',[100 100 1600 800],'color','w');

xp0=0.05;
yp0=0.09;
width=0.4;
height=0.85;
xsep=0.02;

% PLOT 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nplot=1; subplot('Position',[xp0+(nplot-1)*(width+xsep) yp0 width height]); box on;

geoscatter(lat_on_NSW,lon_on_NSW,20,LTER_NSW,'filled'); colormap(cmap); colorbar; caxis([-2 2]);

geobasemap('satellite'); set(gca,'FontSize',14);
gx=gca;

geolimits([min(lat_on_NSW)-0.1 max(lat_on_NSW)+0.1],[min(lon_on_NSW)-0.1 max(lon_on_NSW)+0.1]);

geotickformat -dd;

hold on;

for i=round(min(find(strcmp(STATE,'NSW')))):1000:max(find(strcmp(STATE,'NSW')));
    han1=geoplot(lat_off(i),lon_off(i),'wo'); set(han1,'LineWidth',1,'MarkerFaceColor','w','MarkerEdgeColor','k');
    han1=text(lat_off(i),lon_off(i)+0.05,num2str(i)); set(han1,'color','y');
end

% PLOT 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
width2=0.47;

nplot=2; subplot('Position',[xp0+(nplot-1)*(width+xsep) yp0 width2 height]); hold on; box on;

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

han1=plot([0 0],[1 Ntr],'k--','HandleVisibility','off');

for j=1:N_littoral_cells
        [~,id]=ismember(littoral_cell_tr_start(j):littoral_cell_tr_end(j),ID);
        han1=plot(LTER(id),ID(id),'r'); set(han1,'LineWidth',2);
end

for i=1:N_littoral_cells
    if littoral_cell_tr_length(i)>40
        cellname_1=text(XTEXT,0.5*(littoral_cell_tr_start(i)+littoral_cell_tr_end(i)),littoral_cell_names_unique{i});
        set(cellname_1,'Rotation',0,'FontSize',FONT_SIZE_TEXT);
    end
end

axis([XMIN XMAX min(find(strcmp(STATE,'NSW'))) max(find(strcmp(STATE,'NSW')))]);

set(gca,'FontSize',14,'layer','top','Xtick',XMIN:5:XMAX,'Ytick',round(min(find(strcmp(STATE,'NSW'))),-3):1000:max(find(strcmp(STATE,'NSW'))),'YAxisLocation','right','YDir','reverse');

xlabel({'long-term erosion rate [m/yr]'},'FontSize',18);
ylabel('transect #','FontSize',18);

ytickformat('%0g')

ax = gca;
ax.YAxis.Exponent=0;
ax.YAxis.TickLabelFormat='%.0f';

annotation(gcf,'textbox', [0.0539395227860321 0.587921383647801 0.241060477213968 0.12105031446541],'Color',[1 1 0],...
    'String',['New South Wales,',sprintf('\n'),'Australia'],'LineStyle','none','FontWeight','bold','FontSize',30,'FitBoxToText','off');

han1=text(-8.52127659574468,14422.2161764706,'erosion');set(han1,'Color',[1 0 0],'LineStyle','none','FontWeight','bold','FontSize',20);
han1=text(5.52127659574468,14422.2161764706,'accretion');set(han1,'Color',[0 0 1],'LineStyle','none','FontWeight','bold','FontSize',20);
han1=text(-25.9313725490196,30500.1630892256,'Long-term erosion rate [m/yr]');set(han1,'LineStyle','none','FontWeight','bold','FontSize',20,'Rotation',90);

annotation(gcf,'textarrow',[0.313680555555556 0.299236111111111],[0.796416666666667 0.841416666666667],'String',{'transect #'},'Color',[1 1 0],'FontWeight','bold','FontSize',20);

annotation(gcf,'textbox',[0.086601851851851 0.964666666666667 0.879648148148149 0.0326666666666688],'String',{'Long-term shoreline erosion rates [m/yr] derived from satellite-based shoreline observations (1984-2022)'},'FitBoxToText','off','LineStyle','none','FontWeight','bold','FontSize',20);

print(gcf,'Fig_CoastSat_Rates_NSW.jpg','-djpeg','-r300');
