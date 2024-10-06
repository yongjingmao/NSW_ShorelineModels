%create summary error figures of CoSMos_COAST outputs
clear all
close all
addpath(genpath('/Users/z3392050/Library/CloudStorage/OneDrive-UNSW/Research/ARC_FutureFellowship/Shoreline Modelling/CosMoS-Coast/MatlabToolbox'))
%load in files
pn = '/Users/z3392050/Library/CloudStorage/OneDrive-UNSW/Research/ARC_FutureFellowship/Shoreline Modelling/CosMoS-Coast/model/results/ManagmentCase4_no_hold_the_line_continued_accretion_SydMeshOnly_CorrectedWaves/mat/'
cd(pn)
load('CoSMoS_COAST_results_02-Jun-2023_SL=11cm_transects.mat')
load('CoSMoS_COAST_results_02-Jun-2023_SL=11cm_params.mat')
load('CoSMoS_COAST_results_02-Jun-2023_SL=11cm_state.mat')
dxmap = 100;
lit_size = 15;

%default values that keep everything
tran_min = 1;
tran_max = length(transects);
stLit = 1;
endLit = N_littoral_cells;

%% start turn off if looking at certain areas 
%turn these on if you want to do a subplot of a certain area
for ii=1:length(transects)
    %tran_keep = 
tmpMin(ii) = strcmp(transects(ii).coastsat_tr_name,'aus0169-0000');
tmpMax(ii) = strcmp(transects(ii).coastsat_tr_name,'aus0237-0026');
end
tran_min = find(tmpMin==1,1,'first')
tran_max = find(tmpMax==1,1,'last')
tmp3 = strcmp(littoral_cell_names_unique,littoral_cell_names(tran_min))
tmp4 = strcmp(littoral_cell_names_unique,littoral_cell_names(tran_max))
stLit = find(tmp3==1,1,'first')
endLit = find(tmp4==1,1,'first')
%%end turn of shorten of data set

%figure of RMSE and cells

UTMZONE='56 H';
[lat0, lon0] = utm2deg(x0,y0,repmat(UTMZONE,Ntr,1));
[lat_on, lon_on] = utm2deg(x_on,y_on,repmat(UTMZONE,Ntr,1));
[lat_off, lon_off] = utm2deg(x_off,y_off,repmat(UTMZONE,Ntr,1));
ID=[transects.ID];
STATE={transects.state};



id=strcmp(STATE,'NSW');

ID_NSW=ID(id);
lat_on_NSW=lat_on(id);
lon_on_NSW=lon_on(id);
figure; set(gcf,'PaperPositionMode','auto','Position',[100 100 1600 800],'color','w');

xp0=0.05;
yp0=0.09;
width=0.5;
height=0.85;
xsep=0.02;
nplot=1; subplot('Position',[xp0+(nplot-1)*(width+xsep) yp0 width height]); box on;

geoplot(lat0,lon0,'-b.');
geobasemap('satellite'); set(gca,'FontSize',14);
gx=gca;

geolimits([min(lat_on_NSW(tran_min:tran_max))-0.1 max(lat_on_NSW(tran_min:tran_max))+0.1],[min(lon_on_NSW(tran_min:tran_max))-0.1 max(lon_on_NSW(tran_min:tran_max))+0.1]);

geotickformat -dd;
hold on;
geoplot(lat0,lon0,'-b.');
geobasemap('satellite'); set(gca,'FontSize',14);
gx=gca;

geolimits([min(lat_on_NSW(tran_min:tran_max))-0.1 max(lat_on_NSW(tran_min:tran_max))+0.1],[min(lon_on_NSW(tran_min:tran_max))-0.1 max(lon_on_NSW(tran_min:tran_max))+0.1]);
for i=round(min(find(strcmp(STATE,'NSW')))):dxmap:max(find(strcmp(STATE,'NSW')));
    han1=geoplot(lat_off(i),lon_off(i),'wo'); set(han1,'LineWidth',1,'MarkerFaceColor','w','MarkerEdgeColor','k');
    han1=text(lat_off(i),lon_off(i)+0.05,num2str(i)); set(han1,'color','y');
end

N_littoral_cells=length(littoral_cell_type);

faceColors = makesymbolspec('Polygon',{'INDEX',[1 N_littoral_cells],'FaceColor',polcmap(N_littoral_cells)});
cmap2=faceColors.FaceColor{1,3};
cmap=flipud(redblue(64));
geoscatter(lat_on_NSW,lon_on_NSW,20,RMSE_sat,'filled'); colormap(cmap); colorbar; caxis([0 50]);
geolimits([min(lat_on_NSW(tran_min:tran_max))-0.1 max(lat_on_NSW(tran_min:tran_max))+0.1],[min(lon_on_NSW(tran_min:tran_max))-0.1 max(lon_on_NSW(tran_min:tran_max))+0.1]);

% PLOT 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
width2=0.4;

nplot=2; subplot('Position',[xp0+(nplot-1)*(width+xsep) yp0 width2 height]); hold on; box on;

YMAX_SCALE=1.7;
FACE_ALPHA=0.2;
FONT_SIZE=12;
FONT_SIZE_TEXT=12;

XMIN=0;
XMAX=250;
XTEXT=125;

for i=stLit:endLit
    fill([XMIN,XMIN,XMAX,XMAX],[littoral_cell_tr_start(i),littoral_cell_tr_end(i),littoral_cell_tr_end(i),littoral_cell_tr_start(i)],cmap2(i,:),'FaceAlpha',FACE_ALPHA,'EdgeColor',cmap2(i,:),'EdgeAlpha',0.25,'HandleVisibility','off');
end

han1=plot([0 0],[stLit endLit],'k--','HandleVisibility','off');

for j=stLit:endLit
        [~,id]=ismember(littoral_cell_tr_start(j):littoral_cell_tr_end(j),ID);
        han1=plot(RMSE_sat(id),ID(id),'r'); set(han1,'LineWidth',2);
end

for i=stLit:endLit
    if littoral_cell_tr_length(i)>lit_size
        cellname_1=text(XTEXT,0.5*(littoral_cell_tr_start(i)+littoral_cell_tr_end(i)),littoral_cell_names_unique{i});
        set(cellname_1,'Rotation',0,'FontSize',FONT_SIZE_TEXT);
    end
end

%axis([XMIN XMAX min(find(strcmp(STATE,'NSW'))) max(find(strcmp(STATE,'NSW')))]);
axis([XMIN XMAX ID(tran_min) ID(tran_max)]);

set(gca,'FontSize',14,'layer','top','Xtick',XMIN:20:XMAX,'Ytick',ID(tran_min):100:ID(tran_max),'YAxisLocation','right','YDir','reverse');
xlabel('RMSE SAT (m)')

