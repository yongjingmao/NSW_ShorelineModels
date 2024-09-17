function [] = transects2kml(transects,kml_output,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write CoSMoS-COAST transects to a .kml output. SeanV
% input: "transects" is the CoSMoS-COAST transects input
% input: "kml_output" is the CoSMoS-COAST transects input
% optional input: the optional 3rd input argument is the bounding box struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(varargin) % get bounding box as the optional 3rd argument
    BB=varargin{:};
end

if isa(class(kml_output),'kml')                % if kmlout is already kml object
    fprintf('adding transects folder to current .kml object ... ');
    f1 = kml_output.createFolder('transects'); % create a folder called transects
elseif ischar(kml_output)                      % else if kmlout is just a sring of a filename to write to
    fprintf('saving transects to .kml file ... ');
    kmlout = kml(kml_output);              % create a new kml object
    f1 = kmlout.createFolder('transects'); % and make a folder within it
else
    error('input variable ''kmlout'' should be a string/char or a kml object. SeanV.');
end

% number of transects
Ntr=length(transects);

ID=[transects.ID]';
x_on=[transects.x_on]';
y_on=[transects.y_on]';
x_off=[transects.x_off]';
y_off=[transects.y_off]';

% convert to Lat/Lon from meters UTM zone
UTMZONE='56 H';
[lat_tr_on, lon_tr_on] = utm2deg(x_on,y_on,repmat(UTMZONE,Ntr,1));
[lat_tr_off, lon_tr_off] = utm2deg(x_off,y_off,repmat(UTMZONE,Ntr,1));

% colorHexg='FF00FF00'; % colors
% colorHexy='7800FFFF';
% colorHexr='FF0000FF';
% colorHexm='FFCC33BA';
% colorHexk='FF000000';

colorHexr=rgba2hex(255,0  ,0  ,255); % colors for different transects
colorHexg=rgba2hex(0  ,255,0  ,255);
colorHexy=rgba2hex(255,255,0  ,255);
colorHexm=rgba2hex(255,0  ,255,255);
colorHexk=rgba2hex(0  ,0  ,0  ,255);
colorHexc=rgba2hex(0,255,255,128);
colorHexb=rgba2hex(0  ,0  ,255,255);

% colorHexr=rgb2hex(uint8([255,0  ,0  ])); % colors for different transects
% colorHexg=rgb2hex(uint8([0  ,255,0  ]));
% colorHexy=rgb2hex(uint8([255,255,0  ]));
% colorHexm=rgb2hex(uint8([255,0  ,255]));
% colorHexk=rgb2hex(uint8([0  ,0  ,0  ]));
% colorHexc=rgb2hex(uint8([0,255,255]));
% colorHexb=rgb2hex(uint8([0  ,0  ,255]));


for i=1:Ntr % for all transects
    
    description=['<html><table border="1">', ...
        '<tr><td>Transect ID                  [-]    </td><td>',num2str(ID(i))                      ,'</td></tr>', ...
        '<tr><td>Transect number              [-]    </td><td>',num2str(i)                          ,'</td></tr>', ...
...%        '<tr><td>State                        [-]    </td><td>',transects(i).state                  ,'</td></tr>', ...
        '<tr><td>Littoral cell name           [-]    </td><td>',transects(i).littoral_cell          ,'</td></tr>', ...
        '<tr><td>Model type                   [-]    </td><td>',transects(i).model_type             ,'</td></tr>', ...
        '<tr><td>Number of shorelines         [-]    </td><td>',num2str(sum(~isnan(transects(i).Y))),'</td></tr>', ...
        '<tr><td>CoastSat transect name       [-]    </td><td>',transects(i).coastsat_tr_name             ,'</td></tr>', ...
        '<tr><td>CoastSat bounding box name   [-]    </td><td>',transects(i).coastsat_bbox_name           ,'</td></tr>', ...
        '<tr><td>CoastSat bounding box number [-]    </td><td>',num2str(transects(i).coastsat_bbox_number),'</td></tr>', ...
        ...'<tr><td>Shoreline change rate (+=accre./-=ero.)  [m/yr] </td><td>',num2str(transects(i).LTER)         ,'</td></tr>', ...
        ...'<tr><td>v_lt (assimilated)                       [m/yr] </td><td>',num2str(nanmean(vlt_assim(i,:)),4)  ,'</td></tr>', ...
        ...'<tr><td>v_lt (projected)                         [m/yr] </td><td>',num2str(nanmean(vlt(i,:)),4)        ,'</td></tr>', ...
        ...'<tr><td>DT  (eq. shoreline change time scale)    [days] </td><td>',num2str(nanmean(DT(i,:)),4)         ,'</td></tr>', ...
        ...'<tr><td>DY  (eq. shoreline change erosion scale) [m]    </td><td>',num2str(nanmean(DY(i,:)),4)         ,'</td></tr>', ...
        ...'<tr><td>Hsb (eq. background wave height)         [m]    </td><td>',num2str(nanmean(HSB(i,:)),4)        ,'</td></tr>', ...
        ...'<tr><td>c (Bruun coefficient)                    [-]    </td><td>',num2str(nanmean(c(i,:)),4)          ,'</td></tr>', ...
        ...'<tr><td>K (longshore transport coefficient)      [-]    </td><td>',num2str(nanmean(K(i,:)),4)          ,'</td></tr>', ...
        ...'<tr><td>sigma (additive noise paramater)         [m]    </td><td>',num2str(nanmean(sigma(i,:)),4)      ,'</td></tr>', ...
        ...'<tr><td>notes:                              </td><td>',transects(i).notes,'</td></tr>', ...
        '</table> '];
    
    %description=['is broken'];
    
    INCLUDE_FIGS=0;
    if INCLUDE_FIGS
        description_figs=['<br></br> <body> <img src="..\figures\CoSMoS-COAST_calibration_validation_tr_ID_',num2str(ID(i),'%0.5d'),'.png" alt="CoSMoS-COAST results: model validation" width="920" height="690"> </body>' ...
            '<br></br> <body> <img src="..\figures\CoSMoS-COAST_results_components_tr_ID_',num2str(ID(i),'%0.5d'),'.png" alt="CoSMoS-COAST results: model components" width="920" height="1000"> </body>' ...
            '<br></br> <body> <img src="..\figures\CoSMoS-COAST_results_parameters_tr_ID_',num2str(ID(i),'%0.5d'),'.png" alt="CoSMoS-COAST results: model parameters" width="920" height="1000"> </body>' ...
            '<br></br> <body> <img src="..\figures\CoSMoS-COAST_results_alongshore_tr_ID_',num2str(ID(i),'%0.5d'),'.png" alt=""CoSMoS-COAST results: model parameters (alongshore)" width="920" height="1000"> </body>' ...
            '<br></br> <body> <img src="..\..\LTER\LTER_tr_',num2str(ID(i),'%0.5d'),'.png" alt="Long-term erosion rate" width="900" height="400"> </body>'];
        
        description=strcat(description,description_figs);
    end
    
    description=strcat(description,'</html>');

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
    
    % write individual transect to kml
    color
    f1.plot([lon_tr_on(i) lon_tr_off(i)],[lat_tr_on(i) lat_tr_off(i)],'name',['transect ID: ',num2str(ID(i))],'tessellate',false,'lineWidth',2,'lineColor',color,'description',description);
end

% get unique littoral cell names
littoral_cell_name_unique=unique({transects.littoral_cell});

f2 = kmlout.createFolder('littoral_cells'); % and make a folder within it

% get littoral cell information from the transects
[littoral_cell_name_unique,littoral_cell_type,...
    littoral_cell_tr_start,littoral_cell_tr_end,...
    N_full_model_sections,N_cross_shore_only_sections,N_rate_only_sections,...
    id_full_model_start,id_full_model_end,...
    id_cross_shore_only_start,id_cross_shore_only_end,...
    id_rate_only_start,id_rate_only_end]=transects2cells(transects);

N_littoral_cells=length(littoral_cell_name_unique);

cmap=round(255*colormap(prism(N_littoral_cells)));

for i=1:N_littoral_cells % for each cell

    % color
    colorHex=rgba2hex(cmap(i,1),cmap(i,2),cmap(i,3),255); % colors for different transects
    % colorHex=rgb2hex(uint8([cmap(i,1),cmap(i,2),cmap(i,3)])); % colors for different transects

    % description the kml data
    description=['littoral cell indicator for the ''',littoral_cell_name_unique{i},''' cell (#',num2str(i,'%d'),' of ',num2str(N_littoral_cells,'%d'),') with model type=''',littoral_cell_type{i},''''];

    % which cells to plot
    id=littoral_cell_tr_start(i):littoral_cell_tr_end(i);

    for j=1:length(id)
       f2.scatter(lon_tr_off(id(j)),lat_tr_off(id(j)),'name',['littoral cell: ',littoral_cell_name_unique{i}],'iconColor',colorHex,'iconScale',0.5,'description',description);
    end
  
end

if exist('BB','var') % if the bounding box was given as an input

    f3 = kmlout.createFolder('bounding_boxes'); % and make a folder within it

    for i=1:length(BB) % for all bounding boxes

        % convert to Lat/Lon from UTM
        [latbb, lonbb] = utm2deg(BB(i).bb(1,:),BB(i).bb(2,:),repmat(UTMZONE,4,1));
        [lats, lons] = utm2deg(BB(i).shoreline.x,BB(i).shoreline.y,repmat(UTMZONE,length(BB(i).shoreline.x),1));

        f4 = f3.createFolder(['BB_',num2str(i,'%03d')]); % and make a folder within it

        description=['<html><table border="1">', ...
            '<tr><td>Bounding Box #             [-]     </td><td>',num2str(BB(i).ID)         ,'</td></tr>', ...
            '<tr><td>Bounding Box name          [-]     </td><td>',BB(i).name                ,'</td></tr>', ...
            '<tr><td>Bounding Box State         [-]     </td><td>',BB(i).state               ,'</td></tr>', ...
            '<tr><td>CoastSat Bounding Box name [-]     </td><td>',BB(i).coastsat_bbox_name  ,'</td></tr>', ...
            '<tr><td>Transect start #           [-]     </td><td>',num2str(BB(i).tr_start)   ,'</td></tr>', ...
            '<tr><td>Transect stop #            [-]     </td><td>',num2str(BB(i).tr_end)     ,'</td></tr>', ...
            '<tr><td>Area                       [km^2]  </td><td>',num2str(BB(i).area)       ,'</td></tr>', ...
            '<tr><td>L1                         [km]  </td><td>',num2str(BB(i).L1)           ,'</td></tr>', ...
            '<tr><td>L2                         [km]  </td><td>',num2str(BB(i).L2)           ,'</td></tr>', ...
            '</table><br></br></html>'];

        f4.poly(lonbb,latbb,'name',['Bounding Box # ',num2str(i),' ',BB(i).name],'tessellate',false,'polyColor',colorHexc,'lineColor',colorHexk,'lineWidth',3,'description',description,'altitudeMode','clampToGround');

        description=['shoreline for Bounding Box # ',num2str(i),' ',BB(i).name];

        f4.plot(lons,lats,'name',['shoreline for Bounding Box # ',num2str(i)],'tessellate',false,'lineWidth',4,'lineColor',colorHexb,'description',description);

    end

end
    
f1 = kmlout.createFolder('shorelines');

% initial shoreline
x0=[transects.x_on]+[transects.Y0].*cosd([transects.angle]);
y0=[transects.y_on]+[transects.Y0].*sind([transects.angle]);

[lat0      , lon0     ]=utm2deg(x0,y0,repmat(UTMZONE,Ntr,1));
[lat_cells0,lon_cells0]=polysplit(lat0,lon0);

f2 = f1.createFolder('initial shoreline');
for i=1:length(lat_cells0)
    f2.plot(lon_cells0{i},lat_cells0{i},'name',['non-erodible shoreline #',num2str(i)],'tessellate',false,'lineWidth',1,'lineColor',colorHexg,'description','non-erodible shoreline');
end

% non-eroidble shoreline
x_MIN=[transects.x_on]+[transects.Ymin].*cosd([transects.angle]);
y_MIN=[transects.y_on]+[transects.Ymin].*sind([transects.angle]);

[lat_MIN      , lon_MIN     ]=utm2deg(x_MIN,y_MIN,repmat(UTMZONE,Ntr,1));
[lat_cells_MIN,lon_cells_MIN]=polysplit(lat_MIN,lon_MIN);

f2 = f1.createFolder('non-erodible shoreline');
for i=1:length(lat_cells_MIN)
    f2.plot(lon_cells_MIN{i},lat_cells_MIN{i},'name',['non-erodible shoreline #',num2str(i)],'tessellate',false,'lineWidth',2,'lineColor',colorHexk,'description','non-erodible shoreline');
end

f1 = kmlout.createFolder('misc.');

f2 = f1.createFolder('onshore points');
f2.scatter(lon_tr_on,lat_tr_on,'name','onshore_points','iconColor',colorHexg,'iconScale',0.5,'description','onshore points location','visibility',false);     

f2 = f1.createFolder('offshore points');
f2.scatter(lon_tr_off,lat_tr_off,'name','offshore_points','iconColor',colorHexr,'iconScale',0.5,'description','offshore points location','visibility',false);     

f2 = f1.createFolder('offshore points (ordering)');
f2.plot(lon_tr_off,lat_tr_off,'name','offshore_points_ordering','tessellate',false,'lineWidth',4,'lineColor',colorHexk,'description','offshore points location','visibility',false);

if ischar(kml_output) % if kml_output is just a sring of a filename to write to
    kmlout.save;      % write the .kml file
end  

fprintf('done.\n');