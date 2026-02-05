% This is a script (indended to be run within another script ... [who knew])
% that saves the results of a CoSMoC COAST simulation as a kml file.

fprintf(' Saving CoSMoS-COAST results ... ');

Yout=Y; % Yout is the main shoreline results variable

if HOLD_THE_LINE
    % clip to minimum shoreline
    Yout=max(Yout,Ymin);
end

% get +/- uncertainty bands
Y_uq_min=quantile(Yout,0.025,2)-0.0001;
Y_uq_max=quantile(Yout,0.975,2)+0.0001;

% take ensemble median
Yout=quantile(Yout,0.5,2);

% get the eroded shoreline representing 95% CI bands for the seasonal erosion
SIGMA=1.959964;
Y_eroded=SIGMA*nanstd(YST,[],2); % 95% CI
Y_eroded(bool_cliff_only | bool_rate_only | bool_no_prediction)=0;

% Y_eroded_max=Y_uq_min;
% Y_eroded_min=min(Y_eroded_max,Yout-Y_eroded);

% sum the uncertainties in quadrature
% Y_eroded_max=Y_uq_min;
% Y_eroded_quadrature=sqrt((Yout-Y_uq_min).^2+(Y_eroded).^2);
% Y_eroded_min=min(Y_eroded_max,Yout-Y_eroded_quadrature);

Yst_1yr=NaN(Ntr,1);
Yst_20yr=NaN(Ntr,1);
Yst_100yr=NaN(Ntr,1);

id=find(t_output>tforecast2);

if length(t_output(id))>75

    for i=1:Ntr

        %fprintf('working on %d of %d ... \n',i,Ntr);

        % 1-yr storm-driven erosion
        [YTR,mu,~,~]=GEV_TR(t_output(id),-YST(i,id),[1.1 20 100]);

        if isnan(YTR(1))
            Yst_1yr(i)=mu; % if the GEV fit fails then set the 1-yr to Delta Y
        else
            Yst_1yr(i)=YTR(1);
        end

        if isnan(YTR(2))
            Yst_20yr(i)=1.5*mu; % if the GEV fit fails then set the 20-yr to 1.5 x Delta Y
        else
            Yst_20yr(i)=YTR(2);
        end

        if isnan(YTR(3))
            Yst_100yr(i)=2*mu; % if the GEV fit fails then set the 100-yr to 2 x Delta Y
        else
            Yst_100yr(i)=YTR(3);
        end

    end

else
    for i=1:Ntr
        Yst_1yr(i)=1*nanmean(DY(i,:),2); % calculated from DY scale
        Yst_20yr(i)=1.5*nanmean(DY(i,:),2); % calculated from DY scale
        Yst_100yr(i)=2*nanmean(DY(i,:),2); % calculated from DY scale
    end
end

Yst_100yr=min(max(0,Yst_100yr),2.0*nanmean(DY,2));
Yst_20yr=min(max(0,Yst_20yr),Yst_100yr);
Yst_1yr=min(max(0,Yst_1yr),Yst_20yr);

[b_smooth,a_smooth] = butter(2,0.01,'low');     % design low-pass filter to smooth longshore variability
Yst_1yr=filtfilt(b_smooth,a_smooth,Yst_1yr);
Yst_20yr=filtfilt(b_smooth,a_smooth,Yst_20yr);
Yst_100yr=filtfilt(b_smooth,a_smooth,Yst_100yr);

%figure; plot(1:Ntr,Yst_1yr,'g',1:Ntr,Yst_20yr,'y',1:Ntr,Yst_100yr,'r');

Y_eroded_max=Y_uq_min;
Y_eroded_min=Y_eroded_max-Yst_1yr;
Y_eroded_min20=Y_eroded_max-Yst_20yr;
Y_eroded_min100=Y_eroded_max-Yst_100yr;

if HOLD_THE_LINE
    
    Y_eroded_min=max(Y_eroded_min,Ymin);
    
    % Here I am NOT clipping the erosion at the 20-year and 100-year to the
    % non-erodible shoreline.  This is intended (in an ad hoc way) to capture undermining due
    % to large wave events.  If undermining is not sought, then the results
    % can be clipped to the non-erodible shoreline for hold-the-line-type
    % scenarios
    
    %Y_eroded_min=max(Y_eroded_min20,Ymin);
    %Y_eroded_min=max(Y_eroded_min100,Ymin);
    
end

% remove output for any 'no prediction' transects
x0p=x0;
y0p=y0;

x0p(bool_no_prediction | bool_cliff_only)=NaN;
y0p(bool_no_prediction | bool_cliff_only)=NaN;

Yout(bool_no_prediction | bool_cliff_only)=NaN;

Y_uq_min(bool_no_prediction | bool_cliff_only)=NaN;
Y_uq_max(bool_no_prediction | bool_cliff_only)=NaN;
Yout(bool_no_prediction | bool_cliff_only)=NaN;

Y_eroded_min(bool_no_prediction | bool_cliff_only | bool_rate_only)=NaN;
Y_eroded_min20(bool_no_prediction | bool_cliff_only | bool_rate_only)=NaN;
Y_eroded_min100(bool_no_prediction | bool_cliff_only | bool_rate_only)=NaN;
Y_eroded_max(bool_no_prediction | bool_cliff_only | bool_rate_only)=NaN;

x=x_on+Yout.*cosd(phi);            % final shoreline position
y=y_on+Yout.*sind(phi);

x_uq_min=x_on+Y_uq_min.*cosd(phi); % uncertainty bands min shoreline position
y_uq_min=y_on+Y_uq_min.*sind(phi); % uncertainty bands min shoreline position

x_uq_max=x_on+Y_uq_max.*cosd(phi); % uncertainty bands max shoreline position
y_uq_max=y_on+Y_uq_max.*sind(phi); % uncertainty bands max shoreline position

x_ero_min=x_on+Y_eroded_min.*cosd(phi); % uncertainty bands eroded min shoreline position
y_ero_min=y_on+Y_eroded_min.*sind(phi); % uncertainty bands eroded min shoreline position

x_ero_min20=x_on+Y_eroded_min20.*cosd(phi); % uncertainty bands eroded min shoreline position
y_ero_min20=y_on+Y_eroded_min20.*sind(phi); % uncertainty bands eroded min shoreline position

x_ero_min100=x_on+Y_eroded_min100.*cosd(phi); % uncertainty bands eroded min shoreline position
y_ero_min100=y_on+Y_eroded_min100.*sind(phi); % uncertainty bands eroded min shoreline position

x_ero_max=x_on+Y_eroded_max.*cosd(phi); % uncertainty bands eroded min shoreline position
y_ero_max=y_on+Y_eroded_max.*sind(phi); % uncertainty bands eroded min shoreline position

x_MIN=x_on+Ymin.*cosd(phi);             % non-eroidble shoreline
y_MIN=y_on+Ymin.*sind(phi);

% convert to lat-long
[lat0, lon0] = utm2deg(x0p,y0p,repmat(UTMZONE,size(x0)));
[lat, lon] = utm2deg(x,y,repmat(UTMZONE,size(x)));

[lat_uq_min, lon_uq_min] = utm2deg(x_uq_min,y_uq_min,repmat(UTMZONE,size(x_uq_min)));
[lat_uq_max, lon_uq_max] = utm2deg(x_uq_max,y_uq_max,repmat(UTMZONE,size(x_uq_max)));

[lat_ero_max, lon_ero_max] = utm2deg(x_ero_max,y_ero_max,repmat(UTMZONE,size(x_ero_max)));
[lat_ero_min, lon_ero_min] = utm2deg(x_ero_min,y_ero_min,repmat(UTMZONE,size(x_ero_min)));
[lat_ero_min20, lon_ero_min20]   = utm2deg(x_ero_min20 ,y_ero_min20 ,repmat(UTMZONE,size(x_ero_min20 )));
[lat_ero_min100, lon_ero_min100] = utm2deg(x_ero_min100,y_ero_min100,repmat(UTMZONE,size(x_ero_min100)));

[lat_MIN, lon_MIN] = utm2deg(x_MIN,y_MIN,repmat(UTMZONE,size(x_MIN)));

[lat_tr_on, lon_tr_on] = utm2deg(x_on,y_on,repmat(UTMZONE,size(x_on)));
[lat_tr_off, lon_tr_off] = utm2deg(x_off,y_off,repmat(UTMZONE,size(x_off)));

PLOT_CONTINUOUS_RESULTS=0;
if PLOT_CONTINUOUS_RESULTS
    % split plotted results by NaN's (results will NOT have a gap between adjacent transects (with predictions) belonging to different littoral cells)
    [lat_cells0,lon_cells0]=polysplit(lat0,lon0);
    [lat_cells,lon_cells]=polysplit(lat,lon);
    [lat_cells_min,lon_cells_min]=polysplit(lat_uq_min,lon_uq_min);
    [lat_cells_max,lon_cells_max]=polysplit(lat_uq_max,lon_uq_max);
    [lat_cells_ero_min,lon_cells_ero_min]=polysplit(lat_ero_min,lon_ero_min);
    [lat_cells_ero_min20,lon_cells_ero_min20]=polysplit(lat_ero_min20,lon_ero_min20);
    [lat_cells_ero_min100,lon_cells_ero_min100]=polysplit(lat_ero_min100,lon_ero_min100);
    [lat_cells_ero_max,lon_cells_ero_max]=polysplit(lat_ero_max,lon_ero_max);
    [lat_cells_MIN,lon_cells_MIN]=polysplit(lat_MIN,lon_MIN);
else
    % split plotted results by cells (results will have a gap between adjacent transects belonging to different littoral cells)
    lat_cells0=mat2cell(lat0,littoral_cell_tr_length);
    lon_cells0=mat2cell(lon0,littoral_cell_tr_length);
    lat_cells=mat2cell(lat,littoral_cell_tr_length);
    lon_cells=mat2cell(lon,littoral_cell_tr_length);
    lat_cells_min=mat2cell(lat_uq_min,littoral_cell_tr_length);
    lon_cells_min=mat2cell(lon_uq_min,littoral_cell_tr_length);
    lat_cells_max=mat2cell(lat_uq_max,littoral_cell_tr_length);
    lon_cells_max=mat2cell(lon_uq_max,littoral_cell_tr_length);
    lat_cells_ero_min=mat2cell(lat_ero_min,littoral_cell_tr_length);
    lon_cells_ero_min=mat2cell(lon_ero_min,littoral_cell_tr_length);
    lat_cells_ero_min20=mat2cell(lat_ero_min20,littoral_cell_tr_length);
    lon_cells_ero_min20=mat2cell(lon_ero_min20,littoral_cell_tr_length);
    lat_cells_ero_min100=mat2cell(lat_ero_min100,littoral_cell_tr_length);
    lon_cells_ero_min100=mat2cell(lon_ero_min100,littoral_cell_tr_length);
    lat_cells_ero_max=mat2cell(lat_ero_max,littoral_cell_tr_length);
    lon_cells_ero_max=mat2cell(lon_ero_max,littoral_cell_tr_length);
    lat_cells_MIN=mat2cell(lat_MIN,littoral_cell_tr_length);
    lon_cells_MIN=mat2cell(lon_MIN,littoral_cell_tr_length);
end

if SAVE_XLSX_FILE
    
    if ismember(t1(n+1),t_save_xlsx)
        
        % OUTPUT XLSX FILE
        header = {'transect number','onshore transect position','','offshore transect position','','initial shoreline','','final shoreline','',...
            'final shoreline uncertainty band (min)','','final shoreline uncertainty band (med)','','final shoreline uncertainty band (max)',''...
            'final shoreline + potential winter erosion uncertainty band (1-yr)','','final shoreline + potential winter erosion uncertainty band (20-yr)','','final shoreline + potential winter erosion uncertainty band (100-yr)','','final shoreline + potential winter erosion uncertainty band (max)','','non-erodible shoreline',''};
        header_data1 = {'transect number','x','y','x','y','x','y','x','y','x','y','x','y','x','y','x','y','x','y','x','y'};
        header_data2 = {'transect number','lat','lon','lat','lon','lat','lon','lat','lon','lat','lon','lat','lon','lat','lon','lat','lon','lat','lon','lat','lon'};
        
        DATA1=[ID x_on y_on x_off y_off x0 y0 x y x_uq_min y_uq_min x_uq_max y_uq_max x_ero_min y_ero_min x_ero_min20 y_ero_min20 x_ero_min100 y_ero_min100 x_ero_max y_ero_max x_MIN y_MIN];
        DATA2=[ID lat_tr_on lon_tr_on lat_tr_off lon_tr_off lat0 lon0 lat lon lat_uq_min lon_uq_min lat lon lat_uq_max lon_uq_max lat_ero_min lon_ero_min lat_ero_min20 lon_ero_min20 lat_ero_min100 lon_ero_min100 lat_ero_max lon_ero_max lat_MIN lon_MIN];
        
        xlsfilename=[OUTPUT_DIR,filesep,'xlsx',filesep,outputfilename,'.xlsx'];
        
        xlswrite(xlsfilename,header,['CoSMoS_',datestr(t1(n+1)),'_SL_',num2str(100*Smean(n+1),'%2.0f'),'cm UTM'],'A1');           % Write column header 1
        xlswrite(xlsfilename,header_data1,['CoSMoS_',datestr(t1(n+1)),'_SL_',num2str(100*Smean(n+1),'%2.0f'),'cm UTM'],'A2');     % Write column header 2
        xlswrite(xlsfilename,DATA1,['CoSMoS_',datestr(t1(n+1)),'_SL_',num2str(100*Smean(n+1),'%2.0f'),'cm UTM'],'A3');            % Write DATA 1
        
        xlswrite(xlsfilename,header,['CoSMoS_',datestr(t1(n+1)),'_SL_',num2str(100*Smean(n+1),'%2.0f'),'cm'],'A1');           % Write column header 1
        xlswrite(xlsfilename,header_data2,['CoSMoS_',datestr(t1(n+1)),'_SL_',num2str(100*Smean(n+1),'%2.0f'),'cm'],'A2');     % Write column header 3
        xlswrite(xlsfilename,DATA2,['CoSMoS_',datestr(t1(n+1)),'_SL_',num2str(100*Smean(n+1),'%2.0f'),'cm'],'A3');            % Write DATA 2
        
    end
    
end

% write out initial conditions file if this is the very first output time
if t1(n+1)==t0+1
    
     % create a folder in the kml output file named transects.
    f1 = kmlout.createFolder('initial conditions');
    
    % INITIAL shoreline
    
    [r1,g1,b1,a1]= deal(0,255,0,255);
    [rhex, ghex, bhex, ahex ]= deal(dec2hex(r1),dec2hex(g1),dec2hex(b1),dec2hex(a1));
    if length(rhex)==1,rhex=['0' rhex];end
    if length(ghex)==1,ghex=['0' ghex];end
    if length(bhex)==1,bhex=['0' bhex];end
    if length(ahex)==1,ahex=['0' ahex];end
    
    colorHex = [ahex bhex ghex rhex]; colorHexg=colorHex;
    
    f2 = f1.createFolder('initial shoreline');
    
    for i=1:length(lat_cells0)
        [lattmp,lontmp] = polysplit(lat_cells0{i},lon_cells0{i});
        for j=1:length(lattmp)
            f2.plot(lontmp{j},lattmp{j},'name',['initial shoreline #',num2str(i),' (segment #',num2str(j),')'],'tessellate',false,'lineWidth',1,'lineColor',colorHex,'description','initial shoreline');
        end
    end
    
    % non-erodible shoreline
    
    % create a folder in the kml output for the non-erodible shoreline
    f2 = f1.createFolder('non-erodible shoreline');

    [r1,g1,b1,a1] = deal(0,0,0,255);
    [rhex, ghex, bhex, ahex ]= deal(dec2hex(r1),dec2hex(g1),dec2hex(b1),dec2hex(a1));
    if length(rhex)==1,rhex=['0' rhex];end
    if length(ghex)==1,ghex=['0' ghex];end
    if length(bhex)==1,bhex=['0' bhex];end
    if length(ahex)==1,ahex=['0' ahex];end
    
    colorHex = [ahex bhex ghex rhex]; colorHexk=colorHex;
    
    for i=1:length(lat_cells_MIN)
        [lattmp,lontmp] = polysplit(lat_cells_MIN{i},lon_cells_MIN{i});
        for j=1:length(lattmp)
            f2.plot(lontmp{j},lattmp{j},'name',['non-erodible shoreline #',num2str(i),' (segment #',num2str(j),')'],'tessellate',false,'lineWidth',2,'lineColor',colorHex,'description','non-erodible shoreline');
        end
    end
    
end

% Write results to .kml file

f1 = kmlout.createFolder(['CoSMoS_COAST_results_',datestr(t1(n+1)),'_SL=',num2str(100*Smean(n+1),'%2.0f'),'cm']);

% FINAL shoreline
[r1,g1,b1,a1]=deal(255,0,0,255);
[rhex, ghex, bhex, ahex ]= deal(dec2hex(r1),dec2hex(g1),dec2hex(b1),dec2hex(a1));
if length(rhex)==1,rhex=['0' rhex];end
if length(ghex)==1,ghex=['0' ghex];end
if length(bhex)==1,bhex=['0' bhex];end
if length(ahex)==1,ahex=['0' ahex];end

colorHex = [ahex bhex ghex rhex]; colorHexr=colorHex;

f2 = f1.createFolder('final shoreline');

for i=1:length(lat_cells)
    [lattmp,lontmp] = polysplit(lat_cells{i},lon_cells{i});
    for j=1:length(lattmp)
        f2.plot(lontmp{j},lattmp{j},'name',['final shoreline #',num2str(i),' (segment #',num2str(j),')'],'tessellate',false,'lineWidth',1,'lineColor',colorHex,'description','final shoreline');
    end  
end

% uncertainty bands model

[r1,g1,b1,a1] = deal(255,255,0,120);
[rhex, ghex, bhex, ahex ]= deal(dec2hex(r1),dec2hex(g1),dec2hex(b1),dec2hex(a1));
if length(rhex)==1,rhex=['0' rhex];end
if length(ghex)==1,ghex=['0' ghex];end
if length(bhex)==1,bhex=['0' bhex];end
if length(ahex)==1,ahex=['0' ahex];end

colorHex = [ahex bhex ghex rhex]; colorHexy=colorHex;

f2 = f1.createFolder('final shoreline uncertainty');

for i=1:length(lat_cells_min)
    [lattmp_min,lontmp_min] = polysplit(lat_cells_min{i},lon_cells_min{i});
    [lattmp_med,lontmp_med] = polysplit(lat_cells{i},lon_cells{i});
    %[lattmp_max,lontmp_max] = polysplit(lat_cells_max{i},lon_cells_max{i});
    for j=1:length(lattmp_min)
        f2.poly([lontmp_min{j}; flipud(lontmp_med{j})],[lattmp_min{j}; flipud(lattmp_med{j})],'altitudeMode','clampToGround','tessellate',false,'name',['(high-end) shoreline position uncertainty #',num2str(i),' (segment #',num2str(j),')'],'lineWidth',0,'polyColor',colorHex,'description','(high-end) uncertainty in final shoreline position');
    end
end

for i=1:length(lat_cells_max)
    %[lattmp_min,lontmp_min] = polysplit(lat_cells_min{i},lon_cells_min{i});
    [lattmp_med,lontmp_med] = polysplit(lat_cells{i},lon_cells{i});
    [lattmp_max,lontmp_max] = polysplit(lat_cells_max{i},lon_cells_max{i});
    for j=1:length(lattmp_max)
        f2.poly([lontmp_med{j}; flipud(lontmp_max{j})],[lattmp_med{j}; flipud(lattmp_max{j})],'altitudeMode','clampToGround','tessellate',false,'name',['(low-end) shoreline position uncertainty #',num2str(i),' (segment #',num2str(j),')'],'lineWidth',0,'polyColor',colorHex,'description','(low-end) uncertainty in final shoreline position');
    end
end

% uncertainty bands erosion

[r1,g1,b1,a1] = deal(255,140,0,120);
[rhex, ghex, bhex, ahex ]= deal(dec2hex(r1),dec2hex(g1),dec2hex(b1),dec2hex(a1));
if length(rhex)==1,rhex=['0' rhex];end
if length(ghex)==1,ghex=['0' ghex];end
if length(bhex)==1,bhex=['0' bhex];end
if length(ahex)==1,ahex=['0' ahex];end

colorHex = [ahex bhex ghex rhex];

f2 = f1.createFolder('final shoreline + potential erosion uncertainty');

for i=1:length(lat_cells_ero_min)
    [lattmp_min,lontmp_min] = polysplit(lat_cells_ero_min{i},lon_cells_ero_min{i});
    [lattmp_max,lontmp_max] = polysplit(lat_cells_ero_max{i},lon_cells_ero_max{i});
    for j=1:length(lattmp_min)
        f2.poly([lontmp_min{j}; flipud(lontmp_max{j})],[lattmp_min{j}; flipud(lattmp_max{j})],'altitudeMode','clampToGround','tessellate',false,'name',['shoreline position + winter erosion #',num2str(i),' (segment #',num2str(j),') (1-yr storm)'],'lineWidth',0,'polyColor',colorHex,'description','uncertainty in final shoreline position + potential winter erosion (1-yr storm)');
    end
end

[r1,g1,b1,a1] = deal(255,0,0,120);
[rhex, ghex, bhex, ahex ]= deal(dec2hex(r1),dec2hex(g1),dec2hex(b1),dec2hex(a1));
if length(rhex)==1,rhex=['0' rhex];end
if length(ghex)==1,ghex=['0' ghex];end
if length(bhex)==1,bhex=['0' bhex];end
if length(ahex)==1,ahex=['0' ahex];end

colorHex = [ahex bhex ghex rhex];

for i=1:length(lat_cells_ero_min)
    [lattmp_min,lontmp_min] = polysplit(lat_cells_ero_min{i},lon_cells_ero_min{i});
    [lattmp_min20,lontmp_min20] = polysplit(lat_cells_ero_min20{i},lon_cells_ero_min20{i});
    for j=1:length(lattmp_min)
        f2.poly([lontmp_min20{j}; flipud(lontmp_min{j})],[lattmp_min20{j}; flipud(lattmp_min{j})],'altitudeMode','clampToGround','tessellate',false,'name',['shoreline position + winter erosion #',num2str(i),' (segment #',num2str(j),') (20-yr storm)'],'lineWidth',0,'polyColor',colorHex,'description','uncertainty in final shoreline position + potential winter erosion (20-yr storm)');
    end
end

[r1,g1,b1,a1] = deal(102,0,0,120);
[rhex, ghex, bhex, ahex ]= deal(dec2hex(r1),dec2hex(g1),dec2hex(b1),dec2hex(a1));
if length(rhex)==1,rhex=['0' rhex];end
if length(ghex)==1,ghex=['0' ghex];end
if length(bhex)==1,bhex=['0' bhex];end
if length(ahex)==1,ahex=['0' ahex];end

colorHex = [ahex bhex ghex rhex];

for i=1:length(lat_cells_ero_min20)
    [lattmp_min20,lontmp_min20] = polysplit(lat_cells_ero_min20{i},lon_cells_ero_min20{i});
    [lattmp_min100,lontmp_min100] = polysplit(lat_cells_ero_min100{i},lon_cells_ero_min100{i});
    for j=1:length(lattmp_min20)
        f2.poly([lontmp_min100{j}; flipud(lontmp_min20{j})],[lattmp_min100{j}; flipud(lattmp_min20{j})],'altitudeMode','clampToGround','tessellate',false,'name',['shoreline position + winter erosion #',num2str(i),' (segment #',num2str(j),') (100-yr storm)'],'lineWidth',0,'polyColor',colorHex,'description','uncertainty in final shoreline position + potential winter erosion (100-yr storm)');
    end
end

% write transects for the very last time step
if n+1>=len

    % create a folder in the kml output file named transects.
    f1 = kmlout.createFolder('transects');
    
    [lat_tr_on , lon_tr_on ] = utm2deg(x_on ,y_on ,repmat(UTMZONE,size(x_on)));
    [lat_tr_off, lon_tr_off] = utm2deg(x_off,y_off,repmat(UTMZONE,size(x_off)));
    
    colorHexg='FF00FF00'; % colors
    colorHexy='7800FFFF';
    colorHexr='FF0000FF';
    colorHexm='FFCC33BA';
    colorHexk='FF000000';
    
    for i=1:Ntr % for all transects
        
        if strcmp(transects(i).model_type,'full model') || strcmp(transects(i).model_type,'cross-shore only') || strcmp(transects(i).model_type,'rate only')
            
            description=['<html><table border="1">', ...
                '<tr><td>Transect ID                              [-]    </td><td>',num2str(ID(i))                 ,'</td></tr>', ...
                '<tr><td>Transect number                          [-]    </td><td>',num2str(i)                     ,'</td></tr>', ...
                ...%'<tr><td>Transect State                           [-]    </td><td>','CA'                        ,'</td></tr>', ...
                ...%'<tr><td>Littoral cell name                       [-]    </td><td>',littoral_cell_names(i)             ,'</td></tr>', ...
                '<tr><td>Shoreline type                           [-]    </td><td>',char(transects(i).model_type)             ,'</td></tr>', ...
                '<tr><td>Number of shorelines                     [-]    </td><td>',num2str(length(transects(i).Y))                ,'</td></tr>', ...
                ...%'<tr><td>CoastSat transect name                   [-]    </td><td>',transects(i).coastsat_tr_name             ,'</td></tr>', ...
                ...%'<tr><td>CoastSat bounding box name               [-]    </td><td>',transects(i).coastsat_bbox_name           ,'</td></tr>', ...
                ...%'<tr><td>CoastSat bounding box number             [-]    </td><td>',num2str(transects(i).coastsat_bbox_number),'</td></tr>', ...
                '<tr><td>Shoreline change rate (+=accre./-=ero.)  [m/yr] </td><td>',num2str(transects(i).LTER)         ,'</td></tr>', ...
                '<tr><td>v_lt (assimilated)                       [m/yr] </td><td>',num2str(nanmean(vlt_assim(i,:)),4)  ,'</td></tr>', ...
                '<tr><td>v_lt (projected)                         [m/yr] </td><td>',num2str(nanmean(vlt(i,:)),4)        ,'</td></tr>', ...
                '<tr><td>DT  (eq. shoreline change time scale)    [days] </td><td>',num2str(nanmean(DT(i,:)),4)         ,'</td></tr>', ...
                '<tr><td>DY  (eq. shoreline change erosion scale) [m]    </td><td>',num2str(nanmean(DY(i,:)),4)         ,'</td></tr>', ...
                '<tr><td>Hsb (eq. background wave height)         [m]    </td><td>',num2str(nanmean(HSB(i,:)),4)        ,'</td></tr>', ...
                '<tr><td>c (Bruun coefficient)                    [-]    </td><td>',num2str(nanmean(c(i,:)),4)          ,'</td></tr>', ...
                '<tr><td>K (longshore transport coefficient)      [-]    </td><td>',num2str(nanmean(K(i,:)),4)          ,'</td></tr>', ...
                '<tr><td>sigma (additive noise paramater)         [m]    </td><td>',num2str(nanmean(sigma(i,:)),4)      ,'</td></tr>', ...
                '</table> '];
            
            INCLUDE_FIGS=1;
            if INCLUDE_FIGS
                description_figs=['<br></br> <body> <img src="..\figures\CoSMoS-COAST_calibration_validation_tr_ID_',num2str(ID(i),'%0.5d'),'.png" alt="CoSMoS-COAST results: model validation" width="920" height="690"> </body>' ...
                    '<br></br> <body> <img src="..\figures\CoSMoS-COAST_results_components_tr_ID_',num2str(ID(i),'%0.5d'),'.png" alt="CoSMoS-COAST results: model components" width="920" height="1000"> </body>' ...
                    '<br></br> <body> <img src="..\figures\CoSMoS-COAST_results_parameters_tr_ID_',num2str(ID(i),'%0.5d'),'.png" alt="CoSMoS-COAST results: model parameters" width="920" height="1000"> </body>' ...
                    '<br></br> <body> <img src="..\figures\CoSMoS-COAST_results_alongshore_tr_ID_',num2str(ID(i),'%0.5d'),'.png" alt=""CoSMoS-COAST results: model parameters (alongshore)" width="920" height="1000"> </body>' ...
                    '<br></br> <body> <img src="..\..\LTER\LTER_tr_',num2str(ID(i),'%0.5d'),'.png" alt="Long-term erosion rate" width="900" height="400"> </body>'];
                
                    description=strcat(description,description_figs);
            end
    
            description=strcat(description,'</html>');
            
        else
            
            description=['<html><table border="1">', ...
                '<tr><td>Transect ID                              [-]          </td><td>',num2str(ID(i))                    ,'</td></tr>', ...
                '<tr><td>Transect number                          [-]          </td><td>',num2str(i)                        ,'</td></tr>', ...
                '<tr><td>Shoreline type                           [-]          </td><td>',char(transects(i).model_type)     ,'</td></tr>', ...
                '<tr><td>Shoreline change rate (+=accre./-=ero.)  [m/yr]       </td><td>',num2str(transects(i).LTER)        ,'</td></tr>', ...
                '</table> ' ...
                '</html>'];
            
        end
        
        if strcmp(transects(i).model_type,'full model')
            color=colorHexg;
        elseif strcmp(transects(i).model_type,'cross-shore only')
            color=colorHexy;
        elseif strcmp(transects(i).model_type,'rate only')
            color=colorHexr;
        elseif  strcmp(transects(i).model_type,'cliff only')
            color=colorHexm;
        elseif strcmp(transects(i).model_type,'no prediction')
            color=colorHexk;
        else
            error('model type not found')
        end
        
        % seperate "no prediction" transects since the no prediction can be
        % turned on later to surpress output for certain transects
        if strcmp(transects(i).model_type,'no prediction') || bool_no_prediction(i)
            color=colorHexk;
        end
        
        f1.plot([lon_tr_on(i) lon_tr_off(i)],[lat_tr_on(i) lat_tr_off(i)],'name',['transect ID: ',num2str(ID(i))],'tessellate',false,'lineWidth',2,'lineColor',color,'description',description);
        
    end
    
    if SAVE_XLSX_FILE && exist('xlsfilename','var')
        
        excelFileName = [OUTPUT_DIR,filesep,'xlsx',filesep,xlsfilename];
        
        if exist(excelFileName,'file')  % remove sheet1, sheet2 and sheet3 from data file
            
            sheetName = 'Sheet'; % EN: Sheet, DE: Tabelle, etc. (Lang. dependent)
            % Open Excel file.
            objExcel = actxserver('Excel.Application');
            objExcel.Workbooks.Open(excelFileName); % Full path is necessary!
            
            % Delete sheets.
            try
                % Throws an error if the sheets do not exist.
                objExcel.ActiveWorkbook.Worksheets.Item([sheetName '1']).Delete;
                objExcel.ActiveWorkbook.Worksheets.Item([sheetName '2']).Delete;
                objExcel.ActiveWorkbook.Worksheets.Item([sheetName '3']).Delete;
            catch
                % Do nothing.
            end
            
            % Save, close and clean up.
            objExcel.ActiveWorkbook.Save;
            objExcel.ActiveWorkbook.Close;
            objExcel.Quit;
            objExcel.delete;
            
            disp('Sheets 1,2,3 deleted.');
        end
        
    end
    
end

pwd1=pwd;
cd([OUTPUT_DIR,filesep,'google_earth'])
kmlout.save;
cd(pwd1);

fprintf('done.\n');