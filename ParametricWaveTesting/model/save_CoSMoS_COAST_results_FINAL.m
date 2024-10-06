% This is a script (indended to be run within another script ... [who knew])
% that saves the FINAL results of a CoSMoS-COAST simulation as a "focus-group-approved" kml file.

%outputfilename=['CoSMoS_COAST_',Model_name,'_FINAL']; % name of output file

% create a new kml object
kmlout = kml(outputfilename);

% zipped or not? save kmz (zip=1) or kml (zip=0)
kmlout.zip=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load in "focus-group-approved" language (from word document) describing each model ouput

product_info_doc='product_info_NEW2.docx';     % the word doc to load in
str = extractFileText(product_info_doc);      % load in the specified word document

% get "Product Information" text %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id1=strfind(str,'Product Information');                           id1=id1(1);  % find the (first) occurance of the phrase "Product Information"
id2=strfind(str,'Amplifying information for specific variables'); id2=id2(1);  % find the (first) occurance of the phrase "Amplifying information for specific variables"

product_information=str{1}(id1+21:id2-4); % the "Product Information" text is sandwiched in between these two

product_information=strrep(product_information,'%','%%'); % replace single percent signs with double percent signs, which are needed to display properly in the .kml

str{1}(1:id2+82)=[]; % get rid of the previously extracted text, which is not needed anymore

% get insertable text %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id1=strfind(product_information,'TEXT 1 CALIFORNIA');    id1=id1(1);  % find the (first) occurance of the phrase "TEXT 1 CALIFORNIA"
id2=strfind(product_information,'(Sweet et al. 2022).'); id2=id2(1);  % find the (first) occurance of the phrase "(Sweet et al. 2022)."

scenario_txt=product_information(id1:id2+19);

product_information(id1-1:id2+19)=[]; % get rid of the previously extracted text, which is not needed anymore

id1=strfind(scenario_txt,'TEXT 1 CALIFORNIA');id1=id1(1);  % find the (first) occurance of the phrase "TEXT 1 CALIFORNIA"
id2=strfind(scenario_txt,'TEXT 1 FLOSUP');    id2=id2(1);  % find the (first) occurance of the phrase "TEXT 1 FLOSUP"

text1_ca=scenario_txt(id1+19:id2-3);

id1=strfind(scenario_txt,'TEXT 1 FLOSUP');     id1=id1(1);  % find the (first) occurance of the phrase "TEXT 1 FLOSUP"
id2=strfind(scenario_txt,'TEXT 2 CALIFORNIA'); id2=id2(1);  % find the (first) occurance of the phrase "TEXT 2 CALIFORNIA"

text1_flosup=scenario_txt(id1+15:id2-3);

id1=strfind(scenario_txt,'TEXT 2 CALIFORNIA'); id1=id1(1);  % find the (first) occurance of the phrase "TEXT 2 CALIFORNIA"
id2=strfind(scenario_txt,'TEXT 2 FLOSUP');     id2=id2(1);  % find the (first) occurance of the phrase "TEXT 2 FLOSUP"

text2_ca=scenario_txt(id1+19:id2-3);

id1=strfind(scenario_txt,'TEXT 2 FLOSUP'); id1=id1(1);  % find the (first) occurance of the phrase "TEXT 2 FLOSUP"
                                           id2=length(scenario_txt);  % find the end of the string

text2_flosup=scenario_txt(id1+15:id2-1);

if strcmp(Model_name,'California')
    product_information=strrep(product_information,'[INSERT TEXT 1]',text1_ca);
    product_information=strrep(product_information,'[INSERT TEXT 2]',text2_ca);
elseif strcmp(Model_name,'FloSup')
    product_information=strrep(product_information,'[INSERT TEXT 1]',text1_flosup);
    product_information=strrep(product_information,'[INSERT TEXT 2]',text2_flosup);
end

% get slope scenario text %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[id1,id2]=regexp(product_information,'Slope A.*?Additional'); % find the (first) occurance of the phrase "TEXT 1 CALIFORNIA"

scenario_txt=product_information(id1:id2-12);

product_information(id1:id2-12)=[];

% get text for case 1
[id1,id2]=regexp(scenario_txt,'Slope A.*?Slope B'); case1_text=scenario_txt(id1+16:id2-9);
[id1,id2]=regexp(scenario_txt,'Slope B.*?Slope C'); case2_text=scenario_txt(id1+23:id2-9);
[id1,~]  =regexp(scenario_txt,'Slope C'); id2=length(scenario_txt); case3_text=scenario_txt(id1+17:id2);

if     BRUUN_CASE==1
    product_information=strrep(product_information,'[INSERT LANGUAGE FOR CORRESPONDING SLOPE CASE BELOW]',case1_text);
elseif BRUUN_CASE==2
    product_information=strrep(product_information,'[INSERT LANGUAGE FOR CORRESPONDING SLOPE CASE BELOW]',case2_text);
elseif BRUUN_CASE==3
    product_information=strrep(product_information,'[INSERT LANGUAGE FOR CORRESPONDING SLOPE CASE BELOW]',case3_text);
else
    error('Some strange case is appearing. SeanV.\n');
end

% get hold the line / continued accretion scenario text %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[id1,id2]=regexp(product_information,'CASE 1.*?Uncertainty'); % find the (first) occurance of the phrase "TEXT 1 CALIFORNIA"

scenario_txt=product_information(id1:id2-19);

product_information(id1-1:id2-18)=[];

% get text for case 1
[id1,id2]=regexp(scenario_txt,'CASE 1.*?CASE 2');                  case1_text=scenario_txt(id1+39:id2-8);
[id1,id2]=regexp(scenario_txt,'CASE 2.*?CASE 3');                  case2_text=scenario_txt(id1+46:id2-8);
[id1,id2]=regexp(scenario_txt,'CASE 3.*?CASE 4');                  case3_text=scenario_txt(id1+42:id2-8);
[id1,~]  =regexp(scenario_txt,'CASE 4'); id2=length(scenario_txt); case4_text=scenario_txt(id1+49:id2);

if     HOLD_THE_LINE==1 && CONTINUED_ACCRETION==0
    product_information=strrep(product_information,'[INSERT LANGUAGE FOR CORRESPONDING CASE FROM BELOW]',case1_text);
elseif HOLD_THE_LINE==1 && CONTINUED_ACCRETION==1
    product_information=strrep(product_information,'[INSERT LANGUAGE FOR CORRESPONDING CASE FROM BELOW]',case2_text);
elseif HOLD_THE_LINE==0 && CONTINUED_ACCRETION==0
    product_information=strrep(product_information,'[INSERT LANGUAGE FOR CORRESPONDING CASE FROM BELOW]',case3_text);
elseif HOLD_THE_LINE==0 && CONTINUED_ACCRETION==1
    product_information=strrep(product_information,'[INSERT LANGUAGE FOR CORRESPONDING CASE FROM BELOW]',case4_text);
else
    error('Some strange case is appearing. SeanV.\n');
end

% add some html formatting

product_information=strrep(product_information,'[NEWLINE]','<br> <br>');
product_information=strrep(product_information,'Model Scenarios','<b>Model Scenarios</b><br><br>');
product_information=strrep(product_information,'Model Uncertainty','<b>Model Uncertainty</b><br><br>');

product_information=['<![CDATA[',product_information,']]>']; % add necessary html brackets

% get "Initial shoreline" text %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id1=strfind(str,'Initial shoreline'); id1=id1(1);  % find the (first) occurance of the phrase "Initial shoreline"
id2=strfind(str,'Modeled shoreline'); id2=id2(1);  % find the (first) occurance of the phrase "Modeled shoreline"

init_shoreline_folder_description=str{1}(id1+19:id2-3); % the "Initial shoreline" text is sandwiched in between these two

init_shoreline_folder_description=strrep(init_shoreline_folder_description,'%','%%'); % replace single percent signs with double percent signs, which are needed to display properly in the .kml

str{1}(1:id2-1)=[]; % get rid of the previously extracted text, which is not needed anymore

% get "Modeled shoreline" text %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id1=strfind(str,'Modeled shoreline');             id1=id1(1);  % find the (first) occurance of the phrase "Modeled shoreline"
id2=strfind(str,'Modeled shoreline uncertainty'); id2=id2(1);  % find the (first) occurance of the phrase "Modeled shoreline uncertainty"

modeled_shoreline_folder_description=str{1}(id1+19:id2-3); % the "Modeled shoreline" text is sandwiched in between these two

modeled_shoreline_folder_description=strrep(modeled_shoreline_folder_description,'%','%%'); % replace single percent signs with double percent signs, which are needed to display properly in the .kml

str{1}(1:id2-1)=[]; % get rid of the previously extracted text, which is not needed anymore

% get "Modeled shoreline uncertainty" text %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id1=strfind(str,'Modeled shoreline uncertainty');       id1=id1(1);  % find the (first) occurance of the phrase "Modeled shoreline uncertainty"
id2=strfind(str,'Potential storm erosion uncertainty'); id2=id2(1);  % find the (first) occurance of the phrase "Potential storm erosion uncertainty"

modeled_shoreline_uncertainty_folder_description=str{1}(id1+31:id2-4); % the "Modeled shoreline uncertainty" text is sandwiched in between these two

modeled_shoreline_uncertainty_folder_description=strrep(modeled_shoreline_uncertainty_folder_description,'%','%%'); % replace single percent signs with double percent signs, which are needed to display properly in the .kml

str{1}(1:id2-1)=[]; % get rid of the previously extracted text, which is not needed anymore

% get "Potential storm erosion uncertainty" text %%%%%%%%%%%%%%%%%%%%%%%%%%
id1=strfind(str,'Potential storm erosion uncertainty'); id1=id1(1);  % find the (first) occurance of the phrase "Potential storm erosion uncertainty"
id2=strfind(str,'Unresolved process uncertainty');      id2=id2(1);  % find the (first) occurance of the phrase "Unresolved process uncertainty"

potential_storm_erosion_uncertainty_folder_description=str{1}(id1+37:id2-3); % the "Potential storm erosion uncertainty" text is sandwiched in between these two

potential_storm_erosion_uncertainty_folder_description=strrep(potential_storm_erosion_uncertainty_folder_description,'%','%%'); % replace single percent signs with double percent signs, which are needed to display properly in the .kml

str{1}(1:id2-1)=[]; % get rid of the previously extracted text, which is not needed anymore

% get "Unresolved process uncertainty" text %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id1=strfind(str,'Unresolved process uncertainty'); id1=id1(1);  % find the (first) occurance of the phrase "Unresolved process uncertainty"
id2=strfind(str,'Transects');                      id2=id2(1);  % find the (first) occurance of the phrase "Transects"

unresolved_process_uncertainty_folder_description=str{1}(id1+32:id2-3); % the "Modeled shoreline uncertainty" text is sandwiched in between these two

unresolved_process_uncertainty_folder_description=strrep(unresolved_process_uncertainty_folder_description,'%','%%'); % replace single percent signs with double percent signs, which are needed to display properly in the .kml

str{1}(1:id2-1)=[]; % get rid of the previously extracted text, which is not needed anymore

% get "Transects" text %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id1=strfind(str,'Transects');               id1=id1(1);  % find the (first) occurance of the phrase "Transects"
id2=strfind(str,'Landward Model Boundary'); id2=id2(1);  % find the (first) occurance of the phrase "Landward Model Boundary"

transects_folder_description=str{1}(id1+11:id2-4); % the "Transects" text is sandwiched in between these two

transects_folder_description=strrep(transects_folder_description,'%','%%'); % replace single percent signs with double percent signs, which are needed to display properly in the .kml

str{1}(1:id2-1)=[]; % get rid of the previously extracted text, which is not needed anymore

% get "Landward Model Boundary" text %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id1=strfind(str,'Landward Model Boundary');      id1=id1(1);  % find the (first) occurance of the phrase "Landward Model Boundary"
id2=strfind(str,'Shoreline change hazard zone'); id2=id2(1);  % find the (first) occurance of the phrase "Shoreline change hazard zone"

landward_bnd_folder_description=str{1}(id1+25:id2-5); % the "Landward Model Boundary" text is sandwiched in between these two

landward_bnd_folder_description=strrep(landward_bnd_folder_description,'%','%%'); % replace single percent signs with double percent signs, which are needed to display properly in the .kml

str{1}(1:id2-1)=[]; % get rid of the previously extracted text, which is not needed anymore

% get "Shoreline change hazard zone" text %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id1=strfind(str,'Shoreline change hazard zone'); id1=id1(1);  % find the (first) occurance of the phrase "Shoreline change hazard zone"
id2=strfind(str,'Extreme storm hazard zone');    id2=id2(1);  % find the (first) occurance of the phrase "Extreme storm shoreline change hazard zone"

shoreline_change_hazard_zone_folder_description=str{1}(id1+30:id2-3); % the "Shoreline change hazard zone" text is sandwiched in between these two

shoreline_change_hazard_zone_folder_description=strrep(shoreline_change_hazard_zone_folder_description,'%','%%'); % replace single percent signs with double percent signs, which are needed to display properly in the .kml

str{1}(1:id2-1)=[]; % get rid of the previously extracted text, which is not needed anymore

% get "Extreme storm shoreline change hazard zone" text %%%%%%%%%%%%%%%%%%%
id1=strfind(str,'Extreme storm hazard zone');             id1=id1(1);  % find the (first) occurance of the phrase "Extreme storm shoreline change hazard zone"
id2=strfind(str,'Definitions for lines/polygons/points'); id2=id2(1);  % find the end of the string

extreme_storm_hazard_zone_folder_description=str{1}(id1+27:id2-3); % the "Extreme storm shoreline change hazard zone" text is sandwiched in between these two

extreme_storm_hazard_zone_folder_description=strrep(extreme_storm_hazard_zone_folder_description,'%','%%'); % replace single percent signs with double percent signs, which are needed to display properly in the .kml

str{1}(1:id2-1)=[]; % get rid of the previously extracted text, which is not needed anymore


% remove next header line

id1=strfind(str,'Definitions for lines/polygons/points');

str{1}(1:id1+38)=[];

% get "initial shoreline" definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id1=strfind(str,'initial shoreline'); id1=id1(1);  % find the (first) occurance of the phrase "initial shoreline"
id2=strfind(str,'modeled shoreline'); id2=id2(1);  % find the (first) occurance of the phrase "modeled shoreline"

initial_shoreline_definition0=str{1}(id1+19:id2-3); % the "initial shoreline" text is sandwiched in between these two

initial_shoreline_definition0=strrep(initial_shoreline_definition0,'%','%%'); % replace single percent signs with double percent signs, which are needed to display properly in the .kml

str{1}(1:id2-1)=[]; % get rid of the previously extracted text, which is not needed anymore

% get "modeled shoreline" definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id1=strfind(str,'modeled shoreline');                         id1=id1(1);  % find the (first) occurance of the phrase "modeled shoreline"
id2=strfind(str,'upper bound modeled shoreline uncertainty'); id2=id2(1);  % find the (first) occurance of the phrase "upper bound modeled shoreline uncertainty"

modeled_shoreline_definition0=str{1}(id1+19:id2-3); % the "Modeled shoreline" text is sandwiched in between these two

modeled_shoreline_definition0=strrep(modeled_shoreline_definition0,'%','%%'); % replace single percent signs with double percent signs, which are needed to display properly in the .kml

str{1}(1:id2-1)=[]; % get rid of the previously extracted text, which is not needed anymore

% get "upper bound modeled shoreline uncertainty" definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id1=strfind(str,'upper bound modeled shoreline uncertainty'); id1=id1(1);  % find the (first) occurance of the phrase "upper bound modeled shoreline uncertainty"
id2=strfind(str,'lower bound modeled shoreline uncertainty'); id2=id2(1);  % find the (first) occurance of the phrase "lower bound modeled shoreline uncertainty"

ub_modeled_shoreline_uncertainty_definition0=str{1}(id1+43:id2-3); % the "Modeled shoreline" text is sandwiched in between these two

ub_modeled_shoreline_uncertainty_definition0=strrep(ub_modeled_shoreline_uncertainty_definition0,'%','%%'); % replace single percent signs with double percent signs, which are needed to display properly in the .kml

str{1}(1:id2-1)=[]; % get rid of the previously extracted text, which is not needed anymore

% get "lower bound modeled shoreline uncertainty" definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id1=strfind(str,'lower bound modeled shoreline uncertainty'); id1=id1(1);  % find the (first) occurance of the phrase "upper bound modeled shoreline uncertainty"
id2=strfind(str,'potential storm erosion uncertainty');       id2=id2(1);  % find the (first) occurance of the phrase "potential storm erosion uncertainty"

lb_modeled_shoreline_uncertainty_definition0=str{1}(id1+43:id2-3); % the "Modeled shoreline" text is sandwiched in between these two

lb_modeled_shoreline_uncertainty_definition0=strrep(lb_modeled_shoreline_uncertainty_definition0,'%','%%'); % replace single percent signs with double percent signs, which are needed to display properly in the .kml

str{1}(1:id2-1)=[]; % get rid of the previously extracted text, which is not needed anymore

% get "potential storm erosion uncertainty" definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id1=strfind(str,'potential storm erosion uncertainty'); id1=id1(1);  % find the (first) occurance of the phrase "potential storm erosion uncertainty"
id2=strfind(str,'Landward Model Boundary');       id2=id2(1);  % find the (first) occurance of the phrase "Landward Model Boundary"

potential_storm_erosion_definition0=str{1}(id1+53:id2-3); % the "Modeled shoreline" text is sandwiched in between these two

potential_storm_erosion_definition0=strrep(potential_storm_erosion_definition0,'%','%%'); % replace single percent signs with double percent signs, which are needed to display properly in the .kml

str{1}(1:id2-1)=[]; % get rid of the previously extracted text, which is not needed anymore

% get "Landward Model Boundary" definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id1=strfind(str,'Landward Model Boundary');      id1=id1(1);  % find the (first) occurance of the phrase "Landward Model Boundary"
id2=strfind(str,'shoreline change hazard zone'); id2=id2(1);  % find the (first) occurance of the phrase "shoreline change hazard zone"

landward_model_boundary_definition0=str{1}(id1+69:id2-3); % the "Modeled shoreline" text is sandwiched in between these two

landward_model_boundary_definition0=strrep(landward_model_boundary_definition0,'%','%%'); % replace single percent signs with double percent signs, which are needed to display properly in the .kml

str{1}(1:id2-1)=[]; % get rid of the previously extracted text, which is not needed anymore

% get "shoreline change hazard zone" definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id1=strfind(str,'shoreline change hazard zone');      id1=id1(1);  % find the (first) occurance of the phrase "shoreline change hazard zone"
id2=strfind(str,'extreme storm hazard zone'); id2=id2(1);  % find the (first) occurance of the phrase "extreme storm hazard zone"

shoreline_change_hazard_zone_definition0=str{1}(id1+30:id2-3); % the "Modeled shoreline" text is sandwiched in between these two

shoreline_change_hazard_zone_definition0=strrep(shoreline_change_hazard_zone_definition0,'%','%%'); % replace single percent signs with double percent signs, which are needed to display properly in the .kml

str{1}(1:id2-1)=[]; % get rid of the previously extracted text, which is not needed anymore

% get "extreme storm hazard zone" definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id1=strfind(str,'extreme storm hazard zone');      id1=id1(1);  % find the (first) occurance of the phrase "extreme storm hazard zone"
id2=strfind(str,'unresolved process uncertainty'); id2=id2(1);  % find the (first) occurance of the phrase "unresolved process uncertainty"

extreme_storm_hazard_zone_definition0=str{1}(id1+27:id2-3); % the "Modeled shoreline" text is sandwiched in between these two

extreme_storm_hazard_zone_definition0=strrep(extreme_storm_hazard_zone_definition0,'%','%%'); % replace single percent signs with double percent signs, which are needed to display properly in the .kml

str{1}(1:id2-1)=[]; % get rid of the previously extracted text, which is not needed anymore

% get "unresolved process uncertainty" definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id1=strfind(str,'unresolved process uncertainty');      id1=id1(1);  % find the (first) occurance of the phrase "extreme storm hazard zone"
                                                        id2=length(str{1});

unresolved_process_uncertainty_definition0=str{1}(id1+32:id2); % the "Modeled shoreline" text is sandwiched in between these two

unresolved_process_uncertainty_definition0=strrep(unresolved_process_uncertainty_definition0,'%','%%'); % replace single percent signs with double percent signs, which are needed to display properly in the .kml

% old product information
init_conditions_folder_description='Initial conditions (such as the starting shoreline position or model parameter estimates) are needed to run the the model.';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Legend code
legend_code=['\n\t<Folder id="legend0">',...
		'\n\t\t<name>Legend</name>',...
		'\n\t\t<Style id="LegendRadioFolder0">',...
	'\n\t\t		<IconStyle>',...
	'\n\t\t			<scale>0</scale>',...
	'\n\t\t		</IconStyle>',...
	'\n\t\t		<ListStyle>',...
	'\n\t\t			<listItemType>radioFolder</listItemType>',...
	'\n\t\t			<bgColor>00ffffff</bgColor>',...
	'\n\t\t			<maxSnippetLines>2</maxSnippetLines>',...
	'\n\t\t		</ListStyle>',...
	'\n\t\t	</Style>',...
	'\n\t\t	<ScreenOverlay id="A0">',...
	'\n\t\t		<name>Upper Left</name>',...
	'\n\t\t		<color>b0ffffff</color>',...
	'\n\t\t		<Icon>',...
	'\n\t\t			<href>',legend_fig,'</href>',...
	'\n\t\t		</Icon>',...
	'\n\t\t		<overlayXY x="0.01" y="0.99" xunits="fraction" yunits="fraction"/>',...
	'\n\t\t		<screenXY x="0.01" y="0.99" xunits="fraction" yunits="fraction"/>',...
	'\n\t\t		<rotationXY x="0.5" y="0.5" xunits="fraction" yunits="fraction"/>',...
	'\n\t\t		<size x="-1" y="-1" xunits="pixels" yunits="pixels"/>',...
	'\n\t\t	</ScreenOverlay>',...
	'\n\t\t	<ScreenOverlay id="B0">',...
	'\n\t\t		<name>Upper Center</name>',...
	'\n\t\t		<visibility>0</visibility>',...
	'\n\t\t		<color>b0ffffff</color>',...
	'\n\t\t		<Icon>',...
	'\n\t\t			<href>',legend_fig,'</href>',...
	'\n\t\t		</Icon>',...
	'\n\t\t		<overlayXY x="0.5" y="0.99" xunits="fraction" yunits="fraction"/>',...
	'\n\t\t		<screenXY x="0.5" y="0.99" xunits="fraction" yunits="fraction"/>',...
	'\n\t\t		<rotationXY x="0.5" y="0.5" xunits="fraction" yunits="fraction"/>',...
	'\n\t\t		<size x="-1" y="-1" xunits="pixels" yunits="pixels"/>',...
	'\n\t\t	</ScreenOverlay>',...
	'\n\t\t	<ScreenOverlay id="C0">',...
	'\n\t\t		<name>Upper Right</name>',...
	'\n\t\t	<visibility>0</visibility>',...
	'\n\t\t		<color>b0ffffff</color>',...
	'\n\t\t		<Icon>',...
	'\n\t\t			<href>',legend_fig,'</href>',...
	'\n\t\t		</Icon>',...
	'\n\t\t		<overlayXY x="0.99" y="0.99" xunits="fraction" yunits="fraction"/>',...
	'\n\t\t		<screenXY x="0.99" y="0.99" xunits="fraction" yunits="fraction"/>',...
	'\n\t\t		<rotationXY x="0.5" y="0.5" xunits="fraction" yunits="fraction"/>',...
	'\n\t\t		<size x="-1" y="-1" xunits="pixels" yunits="pixels"/>',...
	'\n\t\t	</ScreenOverlay>',...
	'\n\t\t	<ScreenOverlay id="D0">',...
	'\n\t\t		<name>Center Left</name>',...
	'\n\t\t		<visibility>0</visibility>',...
	'\n\t\t		<color>b0ffffff</color>',...
	'\n\t\t		<Icon>',...
	'\n\t\t			<href>',legend_fig,'</href>',...
	'\n\t\t		</Icon>',...
	'\n\t\t		<overlayXY x="0.01" y="0.5" xunits="fraction" yunits="fraction"/>',...
	'\n\t\t		<screenXY x="0.01" y="0.5" xunits="fraction" yunits="fraction"/>',...
	'\n\t\t		<rotationXY x="0.5" y="0.5" xunits="fraction" yunits="fraction"/>',...
	'\n\t\t		<size x="-1" y="-1" xunits="pixels" yunits="pixels"/>',...
	'\n\t\t	</ScreenOverlay>',...
	'\n\t\t	<ScreenOverlay id="E0">',...
	'\n\t\t		<name>Center</name>',...
	'\n\t\t		<visibility>0</visibility>',...
	'\n\t\t		<color>b0ffffff</color>',...
	'\n\t\t		<Icon>',...
	'\n\t\t		<href>',legend_fig,'</href>',...
	'\n\t\t		</Icon>',...
	'\n\t\t		<overlayXY x="0.5" y="0.5" xunits="fraction" yunits="fraction"/>',...
	'\n\t\t		<screenXY x="0.5" y="0.5" xunits="fraction" yunits="fraction"/>',...
	'\n\t\t		<rotationXY x="0.5" y="0.5" xunits="fraction" yunits="fraction"/>',...
	'\n\t\t		<size x="-1" y="-1" xunits="pixels" yunits="pixels"/>',...
	'\n\t\t	</ScreenOverlay>',...
	'\n\t\t	<ScreenOverlay id="F0">',...
	'\n\t\t		<name>Center Right</name>',...
	'\n\t\t		<visibility>0</visibility>',...
	'\n\t\t		<color>b0ffffff</color>',...
	'\n\t\t		<Icon>',...
	'\n\t\t			<href>',legend_fig,'</href>',...
	'\n\t\t		</Icon>',...
	'\n\t\t		<overlayXY x="0.99" y="0.5" xunits="fraction" yunits="fraction"/>',...
	'\n\t\t		<screenXY x="0.99" y="0.5" xunits="fraction" yunits="fraction"/>',...
	'\n\t\t		<rotationXY x="0.5" y="0.5" xunits="fraction" yunits="fraction"/>',...
	'\n\t\t		<size x="-1" y="-1" xunits="pixels" yunits="pixels"/>',...
	'\n\t\t	</ScreenOverlay>',...
	'\n\t\t	<ScreenOverlay id="G0">',...
	'\n\t\t		<name>Lower Left</name>',...
	'\n\t\t		<visibility>0</visibility>',...
	'\n\t\t		<color>b0ffffff</color>',...
	'\n\t\t		<Icon>',...
	'\n\t\t			<href>',legend_fig,'</href>',...
	'\n\t\t		</Icon>',...
	'\n\t\t		<overlayXY x="0.01" y="0.01" xunits="fraction" yunits="fraction"/>',...
	'\n\t\t		<screenXY x="0.01" y="0.01" xunits="fraction" yunits="fraction"/>',...
	'\n\t\t		<rotationXY x="0.5" y="0.5" xunits="fraction" yunits="fraction"/>',...
	'\n\t\t		<size x="-1" y="-1" xunits="pixels" yunits="pixels"/>',...
	'\n\t\t	</ScreenOverlay>',...
	'\n\t\t	<ScreenOverlay id="H0">',...
	'\n\t\t	<name>Lower Center</name>',...
	'\n\t\t		<visibility>0</visibility>',...
	'\n\t\t		<color>b0ffffff</color>',...
	'\n\t\t		<Icon>',...
	'\n\t\t			<href>',legend_fig,'</href>',...
	'\n\t\t		</Icon>',...
	'\n\t\t		<overlayXY x="0.5" y="0.01" xunits="fraction" yunits="fraction"/>',...
	'\n\t\t		<screenXY x="0.5" y="0.01" xunits="fraction" yunits="fraction"/>',...
	'\n\t\t		<rotationXY x="0.5" y="0.5" xunits="fraction" yunits="fraction"/>',...
	'\n\t\t		<size x="-1" y="-1" xunits="pixels" yunits="pixels"/>',...
	'\n\t\t	</ScreenOverlay>',...
	'\n\t\t	<ScreenOverlay id="I0">',...
	'\n\t\t		<name>Lower Right</name>',...
	'\n\t\t		<visibility>0</visibility>',...
	'\n\t\t		<color>b0ffffff</color>',...
	'\n\t\t		<Icon>',...
	'\n\t\t		    <href>',legend_fig,'</href>',...
	'\n\t\t		</Icon>',...
	'\n\t\t	<overlayXY x="0.99" y="0.01" xunits="fraction" yunits="fraction"/>',...
	'\n\t\t		<screenXY x="0.99" y="0.01" xunits="fraction" yunits="fraction"/>',...
	'\n\t\t		<rotationXY x="0.5" y="0.5" xunits="fraction" yunits="fraction"/>',...
	'\n\t\t		<size x="-1" y="-1" xunits="pixels" yunits="pixels"/>',...
	'\n\t\t	</ScreenOverlay>',...
	'\n\t\t	<Placemark id="ID_0">',...
	'\n\t\t		<name>None</name>',...
	'\n\t\t		<visibility>0</visibility>',...
	'\n\t\t		<styleUrl>#LegendRadioFolder0</styleUrl>',...
	'\n\t\t	</Placemark>',...
	'\n\t\t</Folder>'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with all the various SLR scenarios

SL_SCENARIOS=[0 0.25 0.5 0.75 1.0 1.5 2.0 3.0]; % sea-level rise scenarios in meters [m]

% id's of the transects to include
%ids2output=1:size(YY,1); Ntr_output=length(ids2output);
%ids2output=find(ID>=19321 & ID<=28899); Ntr_output=length(ids2output); %  (NC MODEL) NC only, cut off VA
ids2output=find(ID>=28900 & ID<=32271); Ntr_output=length(ids2output); %  (VA MODEL) VA only, cut off NC & MD portion
%ids2output=find(ID>=32272 & ID<=33273); Ntr_output=length(ids2output); % (MD MODEL) MD only, cut off VA & MD portion
%ids2output=find(ID>=33274 & ID<=34067); Ntr_output=length(ids2output); % (DE MODEL) DE only, cut off MD portion

% the colormap for the different SLR predictions
% cmap=hsv(length(SL_SCENARIOS));

% focus-group-approved color palette
cmap=[255 0 0 ;...
      255 191 0 ; ...
      128 255 0 ; ...
      0 145 3; ...
      0 197 255; ...
      0 92 230; ...
      169 0 230; ...
      255 0 191];

for ii=1:length(SL_SCENARIOS)

    SL_scenario=SL_SCENARIOS(ii);

    if strcmp(Model_name,'FloSup')
        t_SL_scenario=interp1(Smean,t1,SL_scenario,'nearest','extrap');
    else
        t_SL_scenario=interp1(S,t_SLR,SL_scenario,'nearest','extrap');
    end

    if SL_scenario<=0  % the name of the folder is slightly different for the zero SLR scenario
        fprintf('Saving       CoSMoS-COAST results kml for SLR=%2.2f m scenario ... ',SL_scenario);
    else
        fprintf('Saving FINAL CoSMoS-COAST results kml for SLR=%2.2f m scenario ... ',SL_scenario);
    end

    if t_SL_scenario<datenum(2100,1,1)
        
        % find the output time that best matches the time of the desired SLR scenario
        [~,id_output]=min(abs(t_output-t_SL_scenario));

        % get the median shoreline and upper and lower confidence intervals associated with this output time
        Yout=YY(ids2output,id_output);       % Yout is the main shoreline results variable

        dY_uq_min=Yout-YCI1(ids2output,id_output); % get size of +/- uncertainty bands
        dY_uq_max=YCI2(ids2output,id_output)-Yout;

        % scale the results by the applied Bruun factor
        Yout=YY(ids2output,id_output)-YBRU(ids2output,id_output) + BRUUN_FACTOR(ids2output).*YBRU(ids2output,id_output);
        
        Y_uq_min=Yout-dY_uq_min; % get +/- uncertainty bands
        Y_uq_max=Yout+dY_uq_max;

    else
            
        % get the factor of increase for the Bruunian recession component (in my estimation, this is an acceptable alternative to re-running the entire model for the new SLR scenario, since term is purely deterministic based on Ybru=(c/tanAlpha)*SLR)
        SLR_FAC=SL_scenario/SL;

        %if exist('Ylst','var') % if we have the state variable Ylst (and  Y, Ybru, Yvlt), use one method ...

        % update the long-term shoreline position (i.e., the sum of the longshore, Bruunian, and long-term rate components)
        Ylt=Ylst+(SLR_FAC)*repmat(BRUUN_FACTOR,1,Nens).*Ybru+Yvlt;

        % update the full shoreline position
        Y=Ylt+Yst+Y0;

        Yout=Y(ids2output,:); % Yout is the main shoreline results variable

        % get +/- uncertainty bands
        Y_uq_min=quantile(Yout,0.025,2)-0.0001;
        Y_uq_max=quantile(Yout,0.975,2)+0.0001;

        % take ensemble median
        Yout=quantile(Yout,0.5,2);

    end

    if HOLD_THE_LINE
        % clip to minimum shoreline
        Yout=max(Yout,Ymin(ids2output));
        Y_uq_min=max(Y_uq_min,Ymin(ids2output));
        Y_uq_max=max(Y_uq_max,Ymin(ids2output));
    end

    % get the hazard zone bands
    Y_haz_min=min(Y00(ids2output),Y_uq_min);
    Y_haz_max=max(Y00(ids2output),Y_uq_max);

    % get +/- uncertainty bands
    SIGMA=1.959964;
    Y_unresolved_min=Yout-SIGMA*RMSE_sat(ids2output);
    Y_unresolved_med=Yout-0*RMSE_sat(ids2output);
    Y_unresolved_max=Yout+SIGMA*RMSE_sat(ids2output);

    % get the eroded shoreline representing 95% CI bands for the seasonal erosion
    SIGMA=1.959964;
    Y_eroded=SIGMA*nanstd(YST(ids2output,:),[],2); % 95% CI
    Y_eroded(bool_cliff_only(ids2output) | bool_rate_only(ids2output) | bool_no_prediction(ids2output))=0;

    % Y_eroded_max=Y_uq_min;
    % Y_eroded_min=min(Y_eroded_max,Yout-Y_eroded);

    % sum the uncertainties in quadrature
    % Y_eroded_max=Y_uq_min;
    % Y_eroded_quadrature=sqrt((Yout-Y_uq_min).^2+(Y_eroded).^2);
    % Y_eroded_min=min(Y_eroded_max,Yout-Y_eroded_quadrature);

    Yst_1yr=NaN(Ntr_output,1);
    Yst_20yr=NaN(Ntr_output,1);
    Yst_100yr=NaN(Ntr_output,1);

    id=find(t_output>tforecast2);

    if length(t_output(id))>75

        for i=1:Ntr_output

            % 1-yr storm-driven erosion
            [YTR,mu,~,~]=GEV_TR(t_output(id),-YST(ids2output(i),id),[1.1 20 100]);

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
        for i=1:Ntr_output
            Yst_1yr(i)=1*nanmean(DY(ids2output(i),:),2); % calculated from DY scale
            Yst_20yr(i)=1.5*nanmean(DY(ids2output(i),:),2); % calculated from DY scale
            Yst_100yr(i)=2*nanmean(DY(ids2output(i),:),2); % calculated from DY scale
        end
    end

    Yst_100yr=min(max(0,Yst_100yr),2.0*nanmean(DY(ids2output,:),2));
    Yst_20yr=min(max(0,Yst_20yr),Yst_100yr);
    Yst_1yr=min(max(0,Yst_1yr),Yst_20yr);

    [b_smooth,a_smooth] = butter(2,0.01,'low');     % design low-pass filter to smooth longshore variability
    Yst_1yr=filtfilt(b_smooth,a_smooth,Yst_1yr);
    Yst_20yr=filtfilt(b_smooth,a_smooth,Yst_20yr);
    Yst_100yr=filtfilt(b_smooth,a_smooth,Yst_100yr);

    Yst_100yr=min(max(0,Yst_100yr),2.0*nanmean(DY(ids2output,:),2));
    Yst_20yr=min(max(0,Yst_20yr),Yst_100yr);
    Yst_1yr=min(max(0,Yst_1yr),Yst_20yr);

    Y_eroded_max=Y_uq_min;
    Y_eroded_min=Y_eroded_max-Yst_1yr;
    Y_eroded_min20=Y_eroded_max-Yst_20yr;
    Y_eroded_min100=Y_eroded_max-Yst_100yr;

    if HOLD_THE_LINE

        % Here I am NOT clipping the erosion at the 1-year, 20-year, and 100-year to the
        % non-erodible shoreline.  
        % 
        % sometimes clipping only the high frequency extremes is intended (in an ad hoc way) to capture undermining due
        % to low frequency extreme wave events.  If undermining is not sought, then the results
        % can be clipped to the non-erodible shoreline for hold-the-line-type
        % scenarios

        % Y_eroded_min=max(Y_eroded_min,Ymin(ids2output));
        % Y_eroded_min=max(Y_eroded_min20,Ymin(ids2output));
        % Y_eroded_min=max(Y_eroded_min100,Ymin(ids2output));

    end

    if REMOVE_TRANSECT_CROSSING
        
        % identify transects where model results cross on the landward side (during extraordinary erosion extents)

        %id_x=Y_unresolved_min<Ymin2(ids2output) | Y_eroded_min100<Ymin2(ids2output) | Y_uq_min<Ymin2(ids2output);     % to include unresolved processes here or not?
        id_x=Y_eroded_min100<Ymin2(ids2output) | Y_uq_min<Ymin2(ids2output); % they are not included here
    else
        id_x=zeros(size(Y_uq_min)); % set all crossing id's to zeros
    end

    Yout(bool_no_prediction(ids2output) | bool_cliff_only(ids2output))=NaN;

    Y_uq_min(bool_no_prediction(ids2output) | bool_cliff_only(ids2output))=NaN;
    Y_uq_max(bool_no_prediction(ids2output) | bool_cliff_only(ids2output))=NaN;
    Yout(bool_no_prediction(ids2output)| bool_cliff_only(ids2output))=NaN;

    Y_haz_min(bool_no_prediction(ids2output) | bool_cliff_only(ids2output))=NaN;
    Y_haz_max(bool_no_prediction(ids2output) | bool_cliff_only(ids2output))=NaN;

    %Y_storm_haz_min(bool_no_prediction(ids2output) | bool_cliff_only(ids2output))=NaN;
    %Y_storm_haz_max(bool_no_prediction(ids2output) | bool_cliff_only(ids2output))=NaN;

    % I'm not sure about keeping the unresolved process uncertainty for 'rate only' shorelines

%     Y_unresolved_min(bool_no_prediction(ids2output) | bool_cliff_only(ids2output) | bool_rate_only(ids2output))=NaN;
%     Y_unresolved_med(bool_no_prediction(ids2output) | bool_cliff_only(ids2output) | bool_rate_only(ids2output))=NaN;
%     Y_unresolved_max(bool_no_prediction(ids2output) | bool_cliff_only(ids2output) | bool_rate_only(ids2output))=NaN;

    Y_unresolved_min(bool_no_prediction(ids2output) | bool_cliff_only(ids2output) )=NaN;
    Y_unresolved_med(bool_no_prediction(ids2output) | bool_cliff_only(ids2output) )=NaN;
    Y_unresolved_max(bool_no_prediction(ids2output) | bool_cliff_only(ids2output) )=NaN;

    Y_eroded_min(bool_no_prediction(ids2output) | bool_cliff_only(ids2output) | bool_rate_only(ids2output))=NaN;
    Y_eroded_min20(bool_no_prediction(ids2output) | bool_cliff_only(ids2output) | bool_rate_only(ids2output))=NaN;
    Y_eroded_min100(bool_no_prediction(ids2output) | bool_cliff_only(ids2output) | bool_rate_only(ids2output))=NaN;
    Y_eroded_max(bool_no_prediction(ids2output) | bool_cliff_only(ids2output) | bool_rate_only(ids2output))=NaN;

    % get the storm hazard zone bands
    %Y_storm_haz_min=min(Y00(ids2output),Y_eroded_min100);
    %Y_storm_haz_max=max(Y00(ids2output),Y_uq_max);

    % get the storm hazard zone bands (doing this a slightly different way)
    Y_storm_haz_min=min(Y_haz_min,Y_eroded_min100);
    Y_storm_haz_max=max(Y00(ids2output),Y_uq_max);

    % remove output for any 'no prediction' transects
    x0p=x0(ids2output);
    y0p=y0(ids2output);

    x0p(bool_no_prediction(ids2output) | bool_cliff_only(ids2output))=NaN;
    y0p(bool_no_prediction(ids2output) | bool_cliff_only(ids2output))=NaN;

    x_on_p=x_on(ids2output);
    y_on_p=y_on(ids2output);

    x_on_p(bool_no_prediction(ids2output) | bool_cliff_only(ids2output))=NaN;
    y_on_p(bool_no_prediction(ids2output) | bool_cliff_only(ids2output))=NaN;

    x_off_p=x_off(ids2output);
    y_off_p=y_off(ids2output);

    x_off_p(bool_no_prediction(ids2output) | bool_cliff_only(ids2output))=NaN;
    y_off_p(bool_no_prediction(ids2output) | bool_cliff_only(ids2output))=NaN;

    x=x_on(ids2output)+Yout.*cosd(phi(ids2output));            % final shoreline position
    y=y_on(ids2output)+Yout.*sind(phi(ids2output));

    x_uq_min=x_on(ids2output)+Y_uq_min.*cosd(phi(ids2output)); % uncertainty bands min shoreline position
    y_uq_min=y_on(ids2output)+Y_uq_min.*sind(phi(ids2output)); % uncertainty bands min shoreline position

    x_uq_max=x_on(ids2output)+Y_uq_max.*cosd(phi(ids2output)); % uncertainty bands max shoreline position
    y_uq_max=y_on(ids2output)+Y_uq_max.*sind(phi(ids2output)); % uncertainty bands max shoreline position
    
    x_haz_min=x_on(ids2output)+Y_haz_min.*cosd(phi(ids2output)); % uncertainty bands min shoreline position
    y_haz_min=y_on(ids2output)+Y_haz_min.*sind(phi(ids2output)); % uncertainty bands min shoreline position

    x_haz_max=x_on(ids2output)+Y_haz_max.*cosd(phi(ids2output)); % uncertainty bands max shoreline position
    y_haz_max=y_on(ids2output)+Y_haz_max.*sind(phi(ids2output)); % uncertainty bands max shoreline position

    x_storm_haz_min=x_on(ids2output)+Y_storm_haz_min.*cosd(phi(ids2output)); % uncertainty bands min shoreline position
    y_storm_haz_min=y_on(ids2output)+Y_storm_haz_min.*sind(phi(ids2output)); % uncertainty bands min shoreline position

    x_storm_haz_max=x_on(ids2output)+Y_storm_haz_max.*cosd(phi(ids2output)); % uncertainty bands max shoreline position
    y_storm_haz_max=y_on(ids2output)+Y_storm_haz_max.*sind(phi(ids2output)); % uncertainty bands max shoreline position

    x_ero_min=x_on(ids2output)+Y_eroded_min.*cosd(phi(ids2output)); % uncertainty bands eroded min shoreline position
    y_ero_min=y_on(ids2output)+Y_eroded_min.*sind(phi(ids2output)); % uncertainty bands eroded min shoreline position

    x_ero_min20=x_on(ids2output)+Y_eroded_min20.*cosd(phi(ids2output)); % uncertainty bands eroded min shoreline position
    y_ero_min20=y_on(ids2output)+Y_eroded_min20.*sind(phi(ids2output)); % uncertainty bands eroded min shoreline position

    x_ero_min100=x_on(ids2output)+Y_eroded_min100.*cosd(phi(ids2output)); % uncertainty bands eroded min shoreline position
    y_ero_min100=y_on(ids2output)+Y_eroded_min100.*sind(phi(ids2output)); % uncertainty bands eroded min shoreline position

    x_ero_max=x_on(ids2output)+Y_eroded_max.*cosd(phi(ids2output)); % uncertainty bands eroded min shoreline position
    y_ero_max=y_on(ids2output)+Y_eroded_max.*sind(phi(ids2output)); % uncertainty bands eroded min shoreline position

    x_unresolved_min=x_on(ids2output)+Y_unresolved_min.*cosd(phi(ids2output)); % uncertainty bands min shoreline position
    y_unresolved_min=y_on(ids2output)+Y_unresolved_min.*sind(phi(ids2output)); % uncertainty bands min shoreline position

    x_unresolved_med=x_on(ids2output)+Y_unresolved_med.*cosd(phi(ids2output)); % uncertainty bands min shoreline position
    y_unresolved_med=y_on(ids2output)+Y_unresolved_med.*sind(phi(ids2output)); % uncertainty bands min shoreline position

    x_unresolved_max=x_on(ids2output)+Y_unresolved_max.*cosd(phi(ids2output)); % uncertainty bands max shoreline position
    y_unresolved_max=y_on(ids2output)+Y_unresolved_max.*sind(phi(ids2output)); % uncertainty bands max shoreline position

    x_MIN=x_on(ids2output)+Ymin(ids2output).*cosd(phi(ids2output));            % non-eroidble shoreline
    y_MIN=y_on(ids2output)+Ymin(ids2output).*sind(phi(ids2output));

    INFILL_ISOLATED_NANS=0;
    if INFILL_ISOLATED_NANS
        % deal with isolated non-NaN's in the results
        % this just helps to deal with isolated transects neighboored by "no prediction" (i.e., NaNs)
        Y_nan_buf=25;

        [x0p,y0p]=fill_isolated_nonNaNs(x0p,y0p,phi,Y_nan_buf);

        [x,y]=fill_isolated_nonNaNs(x,y,phi,Y_nan_buf);

        [x_uq_min,y_uq_min]=fill_isolated_nonNaNs(x_uq_min,y_uq_min,phi,Y_nan_buf);
        [x_uq_max,y_uq_max]=fill_isolated_nonNaNs(x_uq_max,y_uq_max,phi,Y_nan_buf);

        [x_haz_min,y_haz_min]=fill_isolated_nonNaNs(x_haz_min,y_haz_min,phi,Y_nan_buf);
        [x_haz_max,y_haz_max]=fill_isolated_nonNaNs(x_haz_max,y_haz_max,phi,Y_nan_buf);

        [x_storm_haz_min,y_storm_haz_min]=fill_isolated_nonNaNs(x_storm_haz_min,y_storm_haz_min,phi,Y_nan_buf);
        [x_storm_haz_max,y_storm_haz_max]=fill_isolated_nonNaNs(x_storm_haz_max,y_storm_haz_max,phi,Y_nan_buf);

        [x_unresolved_min,y_unresolved_min]=fill_isolated_nonNaNs(x_unresolved_min,y_unresolved_min,phi,Y_nan_buf);
        [x_unresolved_med,y_unresolved_med]=fill_isolated_nonNaNs(x_unresolved_med,y_unresolved_med,phi,Y_nan_buf);
        [x_unresolved_max,y_unresolved_max]=fill_isolated_nonNaNs(x_unresolved_max,y_unresolved_max,phi,Y_nan_buf);

        [x_ero_max   ,y_ero_max   ]=fill_isolated_nonNaNs(x_ero_max   ,y_ero_max   ,phi,Y_nan_buf);
        [x_ero_min   ,y_ero_min   ]=fill_isolated_nonNaNs(x_ero_min   ,y_ero_min   ,phi,Y_nan_buf);
        [x_ero_min20 ,y_ero_min20 ]=fill_isolated_nonNaNs(x_ero_min20 ,y_ero_min20 ,phi,Y_nan_buf);
        [x_ero_min100,y_ero_min100]=fill_isolated_nonNaNs(x_ero_min100,y_ero_min100,phi,Y_nan_buf);

        % this plus dealing with the transect crossing is not implemented
    end

    % convert to lat-long
    [lat0, lon0] = utm2deg(x0p,y0p,repmat([num2str(UTMZONE,'%d'),' S'],size(x0p)));
    [lat, lon] = utm2deg(x,y,repmat([num2str(UTMZONE,'%d'),' S'],size(x)));

    [lat_on, lon_on] = utm2deg(x_on_p,y_on_p,repmat([num2str(UTMZONE,'%d'),' S'],size(x_on_p)));
    [lat_off, lon_off] = utm2deg(x_off_p,y_off_p,repmat([num2str(UTMZONE,'%d'),' S'],size(x_off_p)));

    [lat_uq_min, lon_uq_min] = utm2deg(x_uq_min,y_uq_min,repmat([num2str(UTMZONE,'%d'),' S'],size(x_uq_min)));
    [lat_uq_max, lon_uq_max] = utm2deg(x_uq_max,y_uq_max,repmat([num2str(UTMZONE,'%d'),' S'],size(x_uq_max)));

    [lat_haz_min, lon_haz_min] = utm2deg(x_haz_min,y_haz_min,repmat([num2str(UTMZONE,'%d'),' S'],size(x_haz_min)));
    [lat_haz_max, lon_haz_max] = utm2deg(x_haz_max,y_haz_max,repmat([num2str(UTMZONE,'%d'),' S'],size(x_haz_max)));

    [lat_storm_haz_min, lon_storm_haz_min] = utm2deg(x_storm_haz_min,y_storm_haz_min,repmat([num2str(UTMZONE,'%d'),' S'],size(x_storm_haz_min)));
    [lat_storm_haz_max, lon_storm_haz_max] = utm2deg(x_storm_haz_max,y_storm_haz_max,repmat([num2str(UTMZONE,'%d'),' S'],size(x_storm_haz_max)));

    [lat_unresolved_min, lon_unresolved_min] = utm2deg(x_unresolved_min,y_unresolved_min,repmat([num2str(UTMZONE,'%d'),' S'],size(x_unresolved_min)));
    [lat_unresolved_med, lon_unresolved_med] = utm2deg(x_unresolved_med,y_unresolved_med,repmat([num2str(UTMZONE,'%d'),' S'],size(x_unresolved_med)));
    [lat_unresolved_max, lon_unresolved_max] = utm2deg(x_unresolved_max,y_unresolved_max,repmat([num2str(UTMZONE,'%d'),' S'],size(x_unresolved_max)));

    [lat_ero_max, lon_ero_max] = utm2deg(x_ero_max,y_ero_max,repmat([num2str(UTMZONE,'%d'),' S'],size(x_ero_max)));
    [lat_ero_min, lon_ero_min] = utm2deg(x_ero_min,y_ero_min,repmat([num2str(UTMZONE,'%d'),' S'],size(x_ero_min)));
    [lat_ero_min20, lon_ero_min20]   = utm2deg(x_ero_min20 ,y_ero_min20 ,repmat([num2str(UTMZONE,'%d'),' S'],size(x_ero_min20 )));
    [lat_ero_min100, lon_ero_min100] = utm2deg(x_ero_min100,y_ero_min100,repmat([num2str(UTMZONE,'%d'),' S'],size(x_ero_min100)));

    [lat_MIN, lon_MIN] = utm2deg(x_MIN,y_MIN,repmat([num2str(UTMZONE,'%d'),' S'],size(x_MIN)));

    [lat_tr_on, lon_tr_on] = utm2deg(x_on_p,y_on_p,repmat([num2str(UTMZONE,'%d'),' S'],size(x_on_p)));
    [lat_tr_off, lon_tr_off] = utm2deg(x_off_p,y_off_p,repmat([num2str(UTMZONE,'%d'),' S'],size(x_off_p)));

    PLOT_CONTINUOUS_RESULTS=0;
    if PLOT_CONTINUOUS_RESULTS
        % split plotted results by NaN's (results will NOT have a gap between adjacent transects (with predictions) belonging to different littoral cells)
        [lat_cells0,lon_cells0]=polysplit(lat0,lon0);
        [lat_cells,lon_cells]=polysplit(lat,lon);
        [lat_cells_min,lon_cells_min]=polysplit(lat_uq_min,lon_uq_min);
        [lat_cells_max,lon_cells_max]=polysplit(lat_uq_max,lon_uq_max);
        [lat_cells_haz_min,lon_cells_haz_min]=polysplit(lat_haz_min,lon_haz_min);
        [lat_cells_haz_max,lon_cells_haz_max]=polysplit(lat_haz_max,lon_haz_max);
        [lat_cells_storm_haz_min,lon_cells_storm_haz_min]=polysplit(lat_storm_haz_min,lon_storm_haz_min);
        [lat_cells_storm_haz_max,lon_cells_storm_haz_max]=polysplit(lat_storm_haz_max,lon_storm_haz_max);
        [lat_cells_unresolved_min,lon_cells_unresolved_min]=polysplit(lat_unresolved_min,lon_unresolved_min);
        [lat_cells_unresolved_med,lon_cells_unresolved_med]=polysplit(lat_unresolved_med,lon_unresolved_med);
        [lat_cells_unresolved_max,lon_cells_unresolved_max]=polysplit(lat_unresolved_max,lon_unresolved_max);
        [lat_cells_ero_min,lon_cells_ero_min]=polysplit(lat_ero_min,lon_ero_min);
        [lat_cells_ero_min20,lon_cells_ero_min20]=polysplit(lat_ero_min20,lon_ero_min20);
        [lat_cells_ero_min100,lon_cells_ero_min100]=polysplit(lat_ero_min100,lon_ero_min100);
        [lat_cells_ero_max,lon_cells_ero_max]=polysplit(lat_ero_max,lon_ero_max);
        [lat_cells_MIN,lon_cells_MIN]=polysplit(lat_MIN,lon_MIN);
    else

        littoral_cell_tr_output_start=littoral_cell_tr_start;
        littoral_cell_tr_output_end=littoral_cell_tr_end;

        id=find(littoral_cell_tr_output_end<min(ids2output));

        littoral_cell_tr_output_start(id)=[];
        littoral_cell_tr_output_end(id)=[];

        id=find(littoral_cell_tr_output_start>max(ids2output));

        littoral_cell_tr_output_start(id)=[];
        littoral_cell_tr_output_end(id)=[];

        littoral_cell_tr_output_start(1)=min(ids2output);
        littoral_cell_tr_output_end(end)=max(ids2output);

        littoral_cell_tr_output_length=littoral_cell_tr_output_end-littoral_cell_tr_output_start+1;
        
        % split plotted results by cells (results will have a gap between adjacent transects belonging to different littoral cells)
        lat_cells0=mat2cell(lat0,littoral_cell_tr_output_length);
        lon_cells0=mat2cell(lon0,littoral_cell_tr_output_length);
        lat_cells=mat2cell(lat,littoral_cell_tr_output_length);
        lon_cells=mat2cell(lon,littoral_cell_tr_output_length);
        lat_cells_min=mat2cell(lat_uq_min,littoral_cell_tr_output_length);
        lon_cells_min=mat2cell(lon_uq_min,littoral_cell_tr_output_length);
        lat_cells_max=mat2cell(lat_uq_max,littoral_cell_tr_output_length);
        lon_cells_max=mat2cell(lon_uq_max,littoral_cell_tr_output_length);

        lat_cells_haz_min=mat2cell(lat_haz_min,littoral_cell_tr_output_length);
        lon_cells_haz_min=mat2cell(lon_haz_min,littoral_cell_tr_output_length);
        lat_cells_haz_max=mat2cell(lat_haz_max,littoral_cell_tr_output_length);
        lon_cells_haz_max=mat2cell(lon_haz_max,littoral_cell_tr_output_length);

        lat_cells_storm_haz_min=mat2cell(lat_storm_haz_min,littoral_cell_tr_output_length);
        lon_cells_storm_haz_min=mat2cell(lon_storm_haz_min,littoral_cell_tr_output_length);
        lat_cells_storm_haz_max=mat2cell(lat_storm_haz_max,littoral_cell_tr_output_length);
        lon_cells_storm_haz_max=mat2cell(lon_storm_haz_max,littoral_cell_tr_output_length);

        lat_cells_unresolved_min=mat2cell(lat_unresolved_min,littoral_cell_tr_output_length);
        lon_cells_unresolved_min=mat2cell(lon_unresolved_min,littoral_cell_tr_output_length);
        lat_cells_unresolved_med=mat2cell(lat_unresolved_med,littoral_cell_tr_output_length);
        lon_cells_unresolved_med=mat2cell(lon_unresolved_med,littoral_cell_tr_output_length);
        lat_cells_unresolved_max=mat2cell(lat_unresolved_max,littoral_cell_tr_output_length);
        lon_cells_unresolved_max=mat2cell(lon_unresolved_max,littoral_cell_tr_output_length);

        lat_cells_ero_min=mat2cell(lat_ero_min,littoral_cell_tr_output_length);
        lon_cells_ero_min=mat2cell(lon_ero_min,littoral_cell_tr_output_length);
        lat_cells_ero_min20=mat2cell(lat_ero_min20,littoral_cell_tr_output_length);
        lon_cells_ero_min20=mat2cell(lon_ero_min20,littoral_cell_tr_output_length);
        lat_cells_ero_min100=mat2cell(lat_ero_min100,littoral_cell_tr_output_length);
        lon_cells_ero_min100=mat2cell(lon_ero_min100,littoral_cell_tr_output_length);
        lat_cells_ero_max=mat2cell(lat_ero_max,littoral_cell_tr_output_length);
        lon_cells_ero_max=mat2cell(lon_ero_max,littoral_cell_tr_output_length);
        lat_cells_MIN=mat2cell(lat_MIN,littoral_cell_tr_output_length);
        lon_cells_MIN=mat2cell(lon_MIN,littoral_cell_tr_output_length);

        % deal with the transect crossing
        id_x_cells=mat2cell(id_x,littoral_cell_tr_output_length);
    end

    if ii==1 % the model initial conditions and the initial shoreline are only written out once

        % create a folder in the kml output file named transects.
        f1 = kmlout.createFolder('model_initial_conditions');

        % INITIAL shoreline

        f2 = f1.createFolder('initial_shoreline');

        colorHex=rgba2hex(0,0,0,255);

        init_shoreline_description=sprintf('<html><table border="1"> <tr><td>SLR_cm  </td><td>NA  </td></tr> <tr><td>Definition  </td><td>%s </td></tr> </table>  </html>',initial_shoreline_definition0);
    
        count_init_shoreline=1;
        for i=1:length(lat_cells0)
            [lattmp,lontmp] = polysplit(lat_cells0{i},lon_cells0{i});
            for j=1:length(lattmp)
                f2.plot(lontmp{j},lattmp{j},'name',['initial_shoreline_#',num2str(count_init_shoreline)],'tessellate',false,'lineWidth',3,'lineColor',colorHex,'description',init_shoreline_description,'altitudeMode','clampToGround');
                count_init_shoreline=count_init_shoreline+1;
            end
        end

    end

    % make a folder for the modeled shoreline position and uncertainty bands

    if SL_scenario<=0  % the name of the folder is slightly different for the zero SLR scenario
        f1 = kmlout.createFolder(['model_position_SLR_',num2str(100*SL_scenario,'%2.0f'),'cm']);
    else
        f1 = kmlout.createFolder(['model_FINAL_position_SLR_',num2str(100*SL_scenario,'%2.0f'),'cm']);
    end

    % FINAL shoreline

    f2 = f1.createFolder('modeled_shoreline');

    % original red
    % colorHex=rgba2hex(255,0,0,255);

    % colormap color
    colorHex=rgba2hex(cmap(ii,1),cmap(ii,2),cmap(ii,3),255);

    modeled_shoreline_definition=strrep(modeled_shoreline_definition0,'[INSERT DATE]',datestr(t_SL_scenario,'dd-mmm-yyyy'));

    if t_SL_scenario<tforecast2
        modeled_shoreline_definition=strrep(modeled_shoreline_definition,'projection/hindcast','hindcast');
    else
        modeled_shoreline_definition=strrep(modeled_shoreline_definition,'projection/hindcast','projection');
    end

    modeled_shoreline_description=sprintf('<html><table border="1"> <tr><td>SLR_cm  </td><td>%s  </td></tr> <tr><td>Definition  </td><td>%s </td></tr> </table>  </html>',num2str(SL_scenario*100,3),modeled_shoreline_definition);

    count_mod_shoreline=1;
    for i=1:length(lat_cells)
        [lattmp,lontmp] = polysplit(lat_cells{i},lon_cells{i});

        % code to handle model output transect crossings
        id_x_cells_nan=double(id_x_cells{i});                     % get the transect crossing id vector for each cell and convert it to an int, so that NaN's can be inserted
        id_x_cells_nan(isnan(lat_cells{i}))=NaN;                 % impose the same NaN structure as the model output
        [id_x_tmp,~] = polysplit(id_x_cells_nan,id_x_cells_nan); % and (thus) split the id's in the same form as above

        for j=1:length(lattmp)
            id_not_x=~logical(id_x_tmp{j});
            if ~any(id_not_x)
                continue;
            end

            f2.plot(lontmp{j}(id_not_x),lattmp{j}(id_not_x),'name',['modeled_shoreline_#',num2str(count_mod_shoreline)],'tessellate',false,'lineWidth',3,'lineColor',colorHex,'description',modeled_shoreline_description);
            count_mod_shoreline=count_mod_shoreline+1;
        end
    end

    % uncertainty bands model

    f2 = f1.createFolder('modeled_shoreline_uncertainty');

    colorHex=rgba2hex(255,255,0,120);

    ub_modeled_shoreline_uncertainty_definition=strrep(ub_modeled_shoreline_uncertainty_definition0,'[INSERT DATE]',datestr(t_SL_scenario,'dd-mmm-yyyy'));

    if t_SL_scenario<tforecast2
        ub_modeled_shoreline_uncertainty_definition=strrep(ub_modeled_shoreline_uncertainty_definition,'projection/hindcast','hindcast');    
    else
        ub_modeled_shoreline_uncertainty_definition=strrep(ub_modeled_shoreline_uncertainty_definition,'projection/hindcast','projection');
    end

    ub_modeled_shoreline_unceratinty_description=sprintf('<html><table border="1"> <tr><td>SLR_cm  </td><td>%s  </td></tr> <tr><td>Definition  </td><td>%s </td></tr> </table>  </html>',num2str(SL_scenario*100,3),ub_modeled_shoreline_uncertainty_definition);

    count_mod_shoreline_uncert1=1;
    for i=1:length(lat_cells_min)
        [lattmp_min,lontmp_min] = polysplit(lat_cells_min{i},lon_cells_min{i});
        [lattmp_med,lontmp_med] = polysplit(lat_cells{i},lon_cells{i});

        % code to handle model output transect crossings
        id_x_cells_nan=double(id_x_cells{i});                     % get the transect crossing id vector for each cell and convert it to an int, so that NaN's can be inserted
        id_x_cells_nan(isnan(lat_cells{i}))=NaN;                 % impose the same NaN structure as the model output
        [id_x_tmp,~] = polysplit(id_x_cells_nan,id_x_cells_nan); % and (thus) split the id's in the same form as above

        for j=1:length(lattmp_min)
            id_not_x=~logical(id_x_tmp{j});
            if ~any(id_not_x)
                continue;
            end

            f2.poly([lontmp_min{j}(id_not_x); flipud(lontmp_med{j}(id_not_x))],[lattmp_min{j}(id_not_x); flipud(lattmp_med{j}(id_not_x))],'altitudeMode','clampToGround','tessellate',false,'name',['(upper_bound)_modeled_shoreline_uncertainty_#',num2str(count_mod_shoreline_uncert1)],'lineWidth',0,'polyColor',colorHex,'description',ub_modeled_shoreline_unceratinty_description,'visibility',false,'drawOrder',3);
            count_mod_shoreline_uncert1=count_mod_shoreline_uncert1+1;
        end
    end

    lb_modeled_shoreline_uncertainty_definition=strrep(lb_modeled_shoreline_uncertainty_definition0,'[INSERT DATE]',datestr(t_SL_scenario,'dd-mmm-yyyy'));

    if t_SL_scenario<tforecast2
        lb_modeled_shoreline_uncertainty_definition=strrep(lb_modeled_shoreline_uncertainty_definition,'projection/hindcast','hindcast');    
    else
        lb_modeled_shoreline_uncertainty_definition=strrep(lb_modeled_shoreline_uncertainty_definition,'projection/hindcast','projection');
    end

    lb_modeled_shoreline_unceratinty_description=sprintf('<html><table border="1"> <tr><td>SLR_cm  </td><td>%s  </td></tr> <tr><td>Definition  </td><td>%s </td></tr> </table>  </html>',num2str(SL_scenario*100,3),lb_modeled_shoreline_uncertainty_definition);

    count_mod_shoreline_uncert2=1;
    for i=1:length(lat_cells_max)
        [lattmp_med,lontmp_med] = polysplit(lat_cells{i},lon_cells{i});
        [lattmp_max,lontmp_max] = polysplit(lat_cells_max{i},lon_cells_max{i});

        % code to handle model output transect crossings
        id_x_cells_nan=double(id_x_cells{i});                     % get the transect crossing id vector for each cell and convert it to an int, so that NaN's can be inserted
        id_x_cells_nan(isnan(lat_cells{i}))=NaN;                 % impose the same NaN structure as the model output
        [id_x_tmp,~] = polysplit(id_x_cells_nan,id_x_cells_nan); % and (thus) split the id's in the same form as above

        for j=1:length(lattmp_max)
            id_not_x=~logical(id_x_tmp{j});
            if ~any(id_not_x)
                continue;
            end

            f2.poly([lontmp_med{j}(id_not_x); flipud(lontmp_max{j}(id_not_x))],[lattmp_med{j}(id_not_x); flipud(lattmp_max{j}(id_not_x))],'altitudeMode','clampToGround','tessellate',false,'name',['(lower_bound)_modeled_shoreline_uncertainty_#',num2str(count_mod_shoreline_uncert2)],'lineWidth',0,'polyColor',colorHex,'description',lb_modeled_shoreline_unceratinty_description,'visibility',false,'drawOrder',3);
            count_mod_shoreline_uncert2=count_mod_shoreline_uncert2+1;
        end
    end

    % uncertainty bands erosion

    f2 = f1.createFolder('potential_storm_erosion_uncertainty');

    colorHex=rgba2hex(255,140,0,120);

    potential_storm_erosion_definition=strrep(potential_storm_erosion_definition0,'[INSERT DATE]',datestr(t_SL_scenario,'dd-mmm-yyyy'));

    if t_SL_scenario<tforecast2
        potential_storm_erosion_definition=strrep(potential_storm_erosion_definition,'projection/hindcast','hindcast');    
    else
        potential_storm_erosion_definition=strrep(potential_storm_erosion_definition,'projection/hindcast','projection');
    end

    potential_storm_erosion_unceratinty_description=sprintf('<html><table border="1"> <tr><td>SLR_cm  </td><td>%s  </td></tr> <tr><td>Definition  </td><td>%s </td></tr> </table>  </html>',num2str(SL_scenario*100,3),potential_storm_erosion_definition);

    count_mod_shoreline_ero_uncert1=1;
    for i=1:length(lat_cells_ero_min)
        [lattmp_min,lontmp_min] = polysplit(lat_cells_ero_min{i},lon_cells_ero_min{i});
        [lattmp_max,lontmp_max] = polysplit(lat_cells_ero_max{i},lon_cells_ero_max{i});

        % code to handle model output transect crossings
        id_x_cells_nan=double(id_x_cells{i});                     % get the transect crossing id vector for each cell and convert it to an int, so that NaN's can be inserted
        id_x_cells_nan(isnan(lat_cells_ero_min{i}))=NaN;         % impose the same NaN structure as the model output
        [id_x_tmp,~] = polysplit(id_x_cells_nan,id_x_cells_nan); % and (thus) split the id's in the same form as above

        for j=1:length(lattmp_min)
            id_not_x=~logical(id_x_tmp{j});
            if ~any(id_not_x)
                continue;
            end

            f2.poly([lontmp_min{j}(id_not_x); flipud(lontmp_max{j}(id_not_x))],[lattmp_min{j}(id_not_x); flipud(lattmp_max{j}(id_not_x))],'altitudeMode','clampToGround','tessellate',false,'name',['potential_storm_erosion_uncertainty_(1-yr_storm)_#',num2str(count_mod_shoreline_ero_uncert1)],'lineWidth',0,'polyColor',colorHex,'description',potential_storm_erosion_unceratinty_description,'visibility',false,'drawOrder',3);
            count_mod_shoreline_ero_uncert1=count_mod_shoreline_ero_uncert1+1;
        end
    end

    colorHex=rgba2hex(255,0,0,120);

    count_mod_shoreline_ero_uncert2=1;
    for i=1:length(lat_cells_ero_min)
        [lattmp_min,lontmp_min] = polysplit(lat_cells_ero_min{i},lon_cells_ero_min{i});
        [lattmp_min20,lontmp_min20] = polysplit(lat_cells_ero_min20{i},lon_cells_ero_min20{i});

        % code to handle model output transect crossings
        id_x_cells_nan=double(id_x_cells{i});                     % get the transect crossing id vector for each cell and convert it to an int, so that NaN's can be inserted
        id_x_cells_nan(isnan(lat_cells_ero_min{i}))=NaN;         % impose the same NaN structure as the model output
        [id_x_tmp,~] = polysplit(id_x_cells_nan,id_x_cells_nan); % and (thus) split the id's in the same form as above

        for j=1:length(lattmp_min)
            id_not_x=~logical(id_x_tmp{j});
            if ~any(id_not_x)
                continue;
            end

            f2.poly([lontmp_min20{j}(id_not_x); flipud(lontmp_min{j}(id_not_x))],[lattmp_min20{j}(id_not_x); flipud(lattmp_min{j}(id_not_x))],'altitudeMode','clampToGround','tessellate',false,'name',['potential_storm_erosion_uncertainty_(20-yr_storm)_#',num2str(count_mod_shoreline_ero_uncert2)],'lineWidth',0,'polyColor',colorHex,'description',potential_storm_erosion_unceratinty_description,'visibility',false,'drawOrder',3);
            count_mod_shoreline_ero_uncert2=count_mod_shoreline_ero_uncert2+1;
        end
    end

    colorHex=rgba2hex(102,0,0,120);

    count_mod_shoreline_ero_uncert3=1;
    for i=1:length(lat_cells_ero_min20)
        [lattmp_min20,lontmp_min20] = polysplit(lat_cells_ero_min20{i},lon_cells_ero_min20{i});
        [lattmp_min100,lontmp_min100] = polysplit(lat_cells_ero_min100{i},lon_cells_ero_min100{i});

        % code to handle model output transect crossings
        id_x_cells_nan=double(id_x_cells{i});                     % get the transect crossing id vector for each cell and convert it to an int, so that NaN's can be inserted
        id_x_cells_nan(isnan(lat_cells_ero_min20{i}))=NaN;       % impose the same NaN structure as the model output
        [id_x_tmp,~] = polysplit(id_x_cells_nan,id_x_cells_nan); % and (thus) split the id's in the same form as above

        for j=1:length(lattmp_min20)
            id_not_x=~logical(id_x_tmp{j});
            if ~any(id_not_x)
                continue;
            end

            f2.poly([lontmp_min100{j}(id_not_x); flipud(lontmp_min20{j}(id_not_x))],[lattmp_min100{j}(id_not_x); flipud(lattmp_min20{j}(id_not_x))],'altitudeMode','clampToGround','tessellate',false,'name',['potential_storm_erosion_uncertainty_(100-yr_storm)_#',num2str(count_mod_shoreline_ero_uncert3)],'lineWidth',0,'polyColor',colorHex,'description',potential_storm_erosion_unceratinty_description,'visibility',false,'drawOrder',3);
            count_mod_shoreline_ero_uncert3=count_mod_shoreline_ero_uncert3+1;
        end
    end

    % we don't output these hazard zones for FloSup ... for some reason
    if ~strcmp(Model_name,'FloSup')

        f2 = f1.createFolder('shoreline_change_hazard_zones');

        colorHex=rgba2hex(170,255,255,120);

        shoreline_change_hazard_zone_definition=strrep(shoreline_change_hazard_zone_definition0,'[INSERT DATE]',datestr(t_SL_scenario,'dd-mmm-yyyy'));

        if t_SL_scenario<tforecast2
            shoreline_change_hazard_zone_definition=strrep(shoreline_change_hazard_zone_definition,'projection/hindcast','hindcast');
        else
            shoreline_change_hazard_zone_definition=strrep(shoreline_change_hazard_zone_definition,'projection/hindcast','projection');
        end

        shoreline_change_hazard_zone_description=sprintf('<html><table border="1"> <tr><td>SLR_cm  </td><td>%s  </td></tr> <tr><td>Definition  </td><td>%s </td></tr> </table>  </html>',num2str(SL_scenario*100,3),shoreline_change_hazard_zone_definition);

        count_mod_shoreline_haz1=1;
        for i=1:length(lat_cells_haz_min)
            [lattmp_min,lontmp_min] = polysplit(lat_cells_haz_min{i},lon_cells_haz_min{i});
            [lattmp_max,lontmp_max] = polysplit(lat_cells_haz_max{i},lon_cells_haz_max{i});

            % code to handle model output transect crossings
            id_x_cells_nan=double(id_x_cells{i});                     % get the transect crossing id vector for each cell and convert it to an int, so that NaN's can be inserted
            id_x_cells_nan(isnan(lat_cells_haz_min{i}))=NaN;         % impose the same NaN structure as the model output
            [id_x_tmp,~] = polysplit(id_x_cells_nan,id_x_cells_nan); % and (thus) split the id's in the same form as above

            for j=1:length(lattmp_min)
                id_not_x=~logical(id_x_tmp{j});
                if ~any(id_not_x)
                    continue;
                end

                f2.poly([lontmp_min{j}(id_not_x); flipud(lontmp_max{j}(id_not_x))],[lattmp_min{j}(id_not_x); flipud(lattmp_max{j}(id_not_x))],'altitudeMode','clampToGround','tessellate',false,'name',['shoreline_change_hazard_zone_#',num2str(count_mod_shoreline_haz1)],'lineWidth',0,'polyColor',colorHex,'description',shoreline_change_hazard_zone_description,'visibility',false,'drawOrder',2);
                count_mod_shoreline_haz1=count_mod_shoreline_haz1+1;
            end
        end

        f2 = f1.createFolder('extreme_storm_hazard_zones');

        colorHex=rgba2hex(0,170,255,120);

        extreme_storm_hazard_zone_definition=strrep(extreme_storm_hazard_zone_definition0,'[INSERT DATE]',datestr(t_SL_scenario,'dd-mmm-yyyy'));

        if t_SL_scenario<tforecast2
            extreme_storm_hazard_zone_definition=strrep(extreme_storm_hazard_zone_definition,'projection/hindcast','hindcast');
        else
            extreme_storm_hazard_zone_definition=strrep(extreme_storm_hazard_zone_definition,'projection/hindcast','projection');
        end

        extreme_storm_hazard_zone_description=sprintf('<html><table border="1"> <tr><td>SLR_cm  </td><td>%s  </td></tr> <tr><td>Definition  </td><td>%s </td></tr> </table>  </html>',num2str(SL_scenario*100,3),extreme_storm_hazard_zone_definition);

        count_mod_shoreline_haz2=1;
        for i=1:length(lat_cells_storm_haz_min)
            [lattmp_min,lontmp_min] = polysplit(lat_cells_storm_haz_min{i},lon_cells_storm_haz_min{i});
            [lattmp_max,lontmp_max] = polysplit(lat_cells_storm_haz_max{i},lon_cells_storm_haz_max{i});

            % code to handle model output transect crossings
            id_x_cells_nan=double(id_x_cells{i});                     % get the transect crossing id vector for each cell and convert it to an int, so that NaN's can be inserted
            id_x_cells_nan(isnan(lat_cells_storm_haz_min{i}))=NaN;   % impose the same NaN structure as the model output
            [id_x_tmp,~] = polysplit(id_x_cells_nan,id_x_cells_nan); % and (thus) split the id's in the same form as above

            for j=1:length(lattmp_min)
                id_not_x=~logical(id_x_tmp{j});
                if ~any(id_not_x)
                    continue;
                end

                f2.poly([lontmp_min{j}(id_not_x); flipud(lontmp_max{j}(id_not_x))],[lattmp_min{j}(id_not_x); flipud(lattmp_max{j}(id_not_x))],'altitudeMode','clampToGround','tessellate',false,'name',['extreme_storm_hazard_zone_#',num2str(count_mod_shoreline_haz2)],'lineWidth',0,'polyColor',colorHex,'description',extreme_storm_hazard_zone_description,'visibility',false,'drawOrder',1);
                count_mod_shoreline_haz2=count_mod_shoreline_haz2+1;
            end
        end

    end

    if SL_scenario>=0.25  % unresolved process uncertainty is only shown for SLR >= 0.25 m scenario

        f2 = f1.createFolder('unresolved_process_uncertainty');

        % uncertainty bands color
        colorHex=rgba2hex(146,125,147,120);

        unresolved_process_uncertainty_definition=strrep(unresolved_process_uncertainty_definition0,'[INSERT DATE]',datestr(t_SL_scenario,'dd-mmm-yyyy'));

        if t_SL_scenario<tforecast2
            unresolved_process_uncertainty_definition=strrep(unresolved_process_uncertainty_definition,'projection/hindcast','hindcast');
        else
            unresolved_process_uncertainty_definition=strrep(unresolved_process_uncertainty_definition,'projection/hindcast','projection');
        end

        unresoved_process_unceratinty_description=sprintf('<html><table border="1"> <tr><td>SLR_cm  </td><td>%s  </td></tr> <tr><td>Definition  </td><td>%s </td></tr> </table>  </html>',num2str(SL_scenario*100,3),unresolved_process_uncertainty_definition);

        count_unresolved_process_uncert1=1;
        for i=1:length(lat_cells_unresolved_min)
            [lattmp_min,lontmp_min] = polysplit(lat_cells_unresolved_min{i},lon_cells_unresolved_min{i});
            [lattmp_med,lontmp_med] = polysplit(lat_cells_unresolved_med{i},lon_cells_unresolved_med{i});

            % code to handle model output transect crossings
            id_x_cells_nan=double(id_x_cells{i});                     % get the transect crossing id vector for each cell and convert it to an int, so that NaN's can be inserted
            id_x_cells_nan(isnan(lat_cells_unresolved_min{i}))=NaN;  % impose the same NaN structure as the model output
            [id_x_tmp,~] = polysplit(id_x_cells_nan,id_x_cells_nan); % and (thus) split the id's in the same form as above

            for j=1:length(lattmp_min)
                id_not_x=~logical(id_x_tmp{j});
                if ~any(id_not_x)
                    continue;
                end

                f2.poly([lontmp_min{j}(id_not_x); flipud(lontmp_med{j}(id_not_x))],[lattmp_min{j}(id_not_x); flipud(lattmp_med{j}(id_not_x))],'altitudeMode','clampToGround','tessellate',false,'name',['(upper_bound)_unresolved_process_uncertainty_#',num2str(count_unresolved_process_uncert1)],'lineWidth',0,'polyColor',colorHex,'description',unresoved_process_unceratinty_description,'visibility',false,'drawOrder',0);
                count_unresolved_process_uncert1=count_unresolved_process_uncert1+1;
            end
        end

        count_unresolved_process_uncert2=1;
        for i=1:length(lat_cells_unresolved_max)
            [lattmp_med,lontmp_med] = polysplit(lat_cells_unresolved_med{i},lon_cells_unresolved_med{i});
            [lattmp_max,lontmp_max] = polysplit(lat_cells_unresolved_max{i},lon_cells_unresolved_max{i});

            % code to handle model output transect crossings
            id_x_cells_nan=double(id_x_cells{i});                     % get the transect crossing id vector for each cell and convert it to an int, so that NaN's can be inserted
            id_x_cells_nan(isnan(lat_cells_unresolved_med{i}))=NaN;  % impose the same NaN structure as the model output
            [id_x_tmp,~] = polysplit(id_x_cells_nan,id_x_cells_nan); % and (thus) split the id's in the same form as above

            for j=1:length(lattmp_max)
                id_not_x=~logical(id_x_tmp{j});
                if ~any(id_not_x)
                    continue;
                end

                f2.poly([lontmp_med{j}(id_not_x); flipud(lontmp_max{j}(id_not_x))],[lattmp_med{j}(id_not_x); flipud(lattmp_max{j}(id_not_x))],'altitudeMode','clampToGround','tessellate',false,'name',['(lower_bound)_unresolved_process_uncertainty_#',num2str(count_unresolved_process_uncert2)],'lineWidth',0,'polyColor',colorHex,'description',unresoved_process_unceratinty_description,'visibility',false,'drawOrder',0);
                count_unresolved_process_uncert2=count_unresolved_process_uncert2+1;
            end
        end

    end

    if ii==length(SL_SCENARIOS)

        % create a folder in the kml output file named transects.
        f1 = kmlout.createFolder('transects');

        colorHexg='FF00FF00'; % colors
        colorHexy='7800FFFF';
        colorHexr='FF0000FF';
        colorHexm='FFCC33BA';
        colorHexk='FF000000';

        for i=1:Ntr_output % for all transects

            if strcmp(transects(ids2output(i)).model_type,'full model') || strcmp(transects(i).model_type,'cross-shore only') || strcmp(transects(i).model_type,'rate only') || ~bool_no_prediction(i)

                % seperate "no prediction" transects since the no prediction can be
                % turned on later to surpress output for certain transects
                if bool_no_prediction(ids2output(i))
                    continue
                end

                description=['<html><table border="1">', ...
                    '<tr><td>Trans_ID    [-]   </td><td>',num2str(ID(ids2output(i)))                    ,'</td></tr>', ...
                    '<tr><td>ShrType     [-]   </td><td>',char(transects(ids2output(i)).model_type)     ,'</td></tr>', ...
                    '<tr><td>ChgRate     [m/yr]</td><td>',num2str(transects(ids2output(i)).LTER)        ,'</td></tr>', ...
                    '<tr><td>TrgSlope^-1 [m/m] </td><td>',num2str(BRUUN_FACTOR(ids2output(i))./tanBeta_Bruun(ids2output(i))),'</td></tr>', ...
                    '<tr><td>v_lt_assim  [m/yr]</td><td>',num2str(nanmean(vlt_assim(ids2output(i),:)),4),'</td></tr>', ...
                    '<tr><td>v_lt_proj   [m/yr]</td><td>',num2str(nanmean(vlt(ids2output(i),:)),4)      ,'</td></tr>', ...
                    '<tr><td>DT_days     [days]</td><td>',num2str(nanmean(DT(ids2output(i),:)),4)       ,'</td></tr>', ...
                    '<tr><td>DY_m        [m]   </td><td>',num2str(nanmean(DY(ids2output(i),:)),4)       ,'</td></tr>', ...
                    '<tr><td>Hsb_m       [m]   </td><td>',num2str(nanmean(HSB(ids2output(i),:)),4)      ,'</td></tr>', ...
                    '<tr><td>c_BrunnCo   [-]   </td><td>',num2str(nanmean(c(ids2output(i),:)),4)        ,'</td></tr>', ...
                    '<tr><td>K_LongShrT  [-]   </td><td>',num2str(nanmean(K(ids2output(i),:)),4)        ,'</td></tr>', ...
                    '<tr><td>sigma_m     [m]   </td><td>',num2str(nanmean(sigma(ids2output(i),:)),4)    ,'</td></tr>', ...
                    '</table>  </html>'];

                if strcmp(transects(ids2output(i)).model_type,'full model')
                    color=colorHexg;
                elseif strcmp(transects(ids2output(i)).model_type,'cross-shore only')
                    color=colorHexy;
                elseif strcmp(transects(ids2output(i)).model_type,'rate only')
                    color=colorHexr;
                elseif  strcmp(transects(ids2output(i)).model_type,'cliff only')
                    color=colorHexm;
                else
                    error('model type not found')
                end

                f1.plot([lon_tr_on(i) lon_tr_off(i)],[lat_tr_on(i) lat_tr_off(i)],'name',['transect_ID:',num2str(ID(ids2output(i)))],'tessellate',false,'lineWidth',2,'lineColor',color,'description',description,'altitudeMode','clampToGround');

            end

        end

        % non-erodible shoreline

        if HOLD_THE_LINE

            % create a folder in the kml output for the non-erodible shoreline
            f1 = kmlout.createFolder('landward_model_boundary');

            landward_bnd_description=sprintf('<html><table border="1"> <tr><td>SLR_cm  </td><td> NA </td></tr> <tr><td>Definition  </td><td>%s</td></tr> </table>  </html>',landward_model_boundary_definition0);

            colorHex=rgba2hex(0,0,0,255);

            for i=1:Ntr_output
                if ~isnan(lat_MIN(i)) && ~bool_no_prediction(i)   
                    f1.point(lon_MIN(i),lat_MIN(i),0,'name',['landward_model_boundary_(transect_ID:',num2str(ID(i)),')'],'iconColor',colorHex,'iconScale',0.25,'description',landward_bnd_description,'visibility',true,'altitudeMode','clampToGround','labelScale',0);
                end
            end

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % writing the shapefile output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if WRITE_SHAPEFILES

        % make the transects, initial shoreline, non-erodible shoreline .shp files?
        % (which are uniform for all SLR and Management Scenarios)

        MAKE_TRANSECTS_INIT_SL_NEL_SHP=1;
        if MAKE_TRANSECTS_INIT_SL_NEL_SHP && ii==1

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Make transects .shp
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Make a shape file struct
            SHP=struct;
            for i=1:Ntr_output

                % Geometry
                SHP(i,1).Geometry='Line';

                % coordiantes
                SHP(i,1).X=[lon_on(i) lon_off(i) NaN];
                SHP(i,1).Y=[lat_on(i) lat_off(i) NaN];

                % bounding box
                xmin=min([lon_on(i) lon_off(i)]);
                xmax=max([lon_on(i) lon_off(i)]);

                ymin=min([lat_on(i) lat_off(i)]);
                ymax=max([lat_on(i) lat_off(i)]);

                % get bounding box
                SHP(i,1).BoundingBox=[xmin ymin; xmax ymax];

                % ID
                SHP(i,1).ID=i;
            end

            % Make the dbf specification object
            dbfspec=struct;
            dbfspec.ID.FieldName='ID';
            dbfspec.ID.FieldType='N';
            dbfspec.ID.FieldLength=5;
            dbfspec.ID.FieldDecimalCount=0;

            % write a transects shapefile & display progress
            shapewrite(SHP,[OUTPUT_DIR,filesep,'shp',filesep,'transects.shp'],'DbfSpec',dbfspec);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Make initial shoreline .shp
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Make a shape file struct
            SHP=struct;
            for i=1:length(lat_cells0)

                % Geometry
                SHP(i,1).Geometry='Line';

                % coordiantes
                SHP(i,1).X=[lon_cells0{i}; NaN];
                SHP(i,1).Y=[lat_cells0{i}; NaN];

                % bounding box
                xmin=min(lon_cells0{i});
                xmax=max(lon_cells0{i});

                ymin=min(lat_cells0{i});
                ymax=max(lat_cells0{i});

                % get bounding box
                SHP(i,1).BoundingBox=[xmin ymin; xmax ymax];

                % ID
                SHP(i,1).ID=i;
            end

            % Make the dbf specification object
            dbfspec=struct;
            dbfspec.ID.FieldName='ID';
            dbfspec.ID.FieldType='N';
            dbfspec.ID.FieldLength=5;
            dbfspec.ID.FieldDecimalCount=0;

            % write a transects shapefile & display progress
            shapewrite(SHP,[OUTPUT_DIR,filesep,'shp',filesep,'initial_shoreline.shp'],'DbfSpec',dbfspec);

            if HOLD_THE_LINE
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Make non-erodible shoreline .shp
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % Make a shape file struct
                SHP=struct;
                for i=1:length(lon_cells_MIN)

                    % Geometry
                    SHP(i,1).Geometry='Line';

                    % coordiantes
                    SHP(i,1).X=[lon_cells_MIN{i}; NaN];
                    SHP(i,1).Y=[lat_cells_MIN{i}; NaN];

                    % bounding box
                    xmin=min(lon_cells_MIN{i});
                    xmax=max(lon_cells_MIN{i});

                    ymin=min(lat_cells_MIN{i});
                    ymax=max(lat_cells_MIN{i});

                    % get bounding box
                    SHP(i,1).BoundingBox=[xmin ymin; xmax ymax];

                    % ID
                    SHP(i,1).ID=i;
                end

                % Make the dbf specification object
                dbfspec=struct;
                dbfspec.ID.FieldName='ID';
                dbfspec.ID.FieldType='N';
                dbfspec.ID.FieldLength=5;
                dbfspec.ID.FieldDecimalCount=0;

                % write a transects shapefile & display progress
                shapewrite(SHP,[OUTPUT_DIR,filesep,'shp',filesep,'nonerodible_shoreline.shp'],'DbfSpec',dbfspec);
            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Make final shoreline .shp
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Make a shape file struct
        SHP=struct;
        for i=1:length(lat_cells)

            % Geometry
            SHP(i,1).Geometry='Line';

            % coordiantes
            SHP(i,1).X=[lon_cells{i}(~id_x_cells{i}); NaN];
            SHP(i,1).Y=[lat_cells{i}(~id_x_cells{i}); NaN];

            % bounding box
            xmin=min(lon_cells{i});
            xmax=max(lon_cells{i});

            ymin=min(lat_cells{i});
            ymax=max(lat_cells{i});

            % get bounding box
            SHP(i,1).BoundingBox=[xmin ymin; xmax ymax];

            % ID
            SHP(i,1).ID=i;
        end

        % Make the dbf specification object
        dbfspec=struct;
        dbfspec.ID.FieldName='ID';
        dbfspec.ID.FieldType='N';
        dbfspec.ID.FieldLength=5;
        dbfspec.ID.FieldDecimalCount=0;

        % write a transects shapefile & display progress
        shapewrite(SHP,[OUTPUT_DIR,filesep,'shp',filesep,'final_shoreline_SLR_',num2str(SL_scenario*100,3),'cm_TrgxSlope',num2str(BRUUN_CASE),'.shp'],'DbfSpec',dbfspec);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Make final shoreline + uncertainty .shp files
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Make a shape file struct
        SHP=struct;

        count_shp=1;
        for i=1:length(lat_cells_min)

            [lattmp_min,lontmp_min] = polysplit(lat_cells_min{i},lon_cells_min{i});
            [lattmp_max,lontmp_max] = polysplit(lat_cells_max{i},lon_cells_max{i});

            % code to handle model output transect crossings
            id_x_cells_nan=double(id_x_cells{i});                     % get the transect crossing id vector for each cell and convert it to an int, so that NaN's can be inserted
            id_x_cells_nan(isnan(lat_cells_min{i}))=NaN;             % impose the same NaN structure as the model output
            [id_x_tmp,~] = polysplit(id_x_cells_nan,id_x_cells_nan); % and (thus) split the id's in the same form as above

            for j=1:length(lattmp_min)

                id_not_x=~logical(id_x_tmp{j});
                if ~any(id_not_x)
                    continue;
                end

                % Geometry
                SHP(count_shp,1).Geometry='Polygon';

                %lon_shp=[lontmp_min{j}; flipud(lontmp_max{j}); lontmp_min{1}; NaN];
                %lat_shp=[lattmp_min{j}; flipud(lattmp_max{j}); lattmp_min{1}; NaN];

                lon_shp=[lontmp_min{j}(id_not_x); flipud(lontmp_max{j}(id_not_x)); NaN];
                lat_shp=[lattmp_min{j}(id_not_x); flipud(lattmp_max{j}(id_not_x)); NaN];

                % coordiantes
                SHP(count_shp,1).X=lon_shp;
                SHP(count_shp,1).Y=lat_shp;

                % bounding box
                xmin=min(lon_shp);
                xmax=max(lon_shp);

                ymin=min(lat_shp);
                ymax=max(lat_shp);

                % get bounding box
                SHP(count_shp,1).BoundingBox=[xmin ymin; xmax ymax];

                % ID
                SHP(count_shp,1).ID=count_shp;

                count_shp=count_shp+1;

            end
        end

        % Make the dbf specification object
        dbfspec=struct;
        dbfspec.ID.FieldName='ID';
        dbfspec.ID.FieldType='N';
        dbfspec.ID.FieldLength=5;
        dbfspec.ID.FieldDecimalCount=0;

        % write a transects shapefile & display progress
        shapewrite(SHP,[OUTPUT_DIR,filesep,'shp',filesep,'final_shoreline_uncertainty_SLR_',num2str(SL_scenario*100,3),'cm_TrgxSlope',num2str(BRUUN_CASE),'.shp'],'DbfSpec',dbfspec);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Make potential storm erosion uncertainty bands
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Make a shape file struct
        SHP=struct;

        count_shp=1;
        for i=1:length(lat_cells_ero_max)

            [lattmp_max,lontmp_max] = polysplit(lat_cells_ero_max{i},lon_cells_ero_max{i});
            [lattmp_min,lontmp_min] = polysplit(lat_cells_ero_min100{i},lon_cells_ero_min100{i});

            % code to handle model output transect crossings
            id_x_cells_nan=double(id_x_cells{i});                     % get the transect crossing id vector for each cell and convert it to an int, so that NaN's can be inserted
            id_x_cells_nan(isnan(lat_cells_ero_max{i}))=NaN;         % impose the same NaN structure as the model output
            [id_x_tmp,~] = polysplit(id_x_cells_nan,id_x_cells_nan); % and (thus) split the id's in the same form as above

            for j=1:length(lattmp_min)

                id_not_x=~logical(id_x_tmp{j});
                if ~any(id_not_x)
                    continue;
                end

                % Geometry
                SHP(count_shp,1).Geometry='Polygon';

                %lon_shp=[lontmp_min{j}; flipud(lontmp_max{j}); lontmp_min{1}; NaN];
                %lat_shp=[lattmp_min{j}; flipud(lattmp_max{j}); lattmp_min{1}; NaN];

                lon_shp=[lontmp_min{j}(id_not_x); flipud(lontmp_max{j}(id_not_x)); NaN];
                lat_shp=[lattmp_min{j}(id_not_x); flipud(lattmp_max{j}(id_not_x)); NaN];

                % coordiantes
                SHP(count_shp,1).X=lon_shp;
                SHP(count_shp,1).Y=lat_shp;

                % bounding box
                xmin=min(lon_shp);
                xmax=max(lon_shp);

                ymin=min(lat_shp);
                ymax=max(lat_shp);

                % get bounding box
                SHP(count_shp,1).BoundingBox=[xmin ymin; xmax ymax];

                % ID
                SHP(count_shp,1).ID=count_shp;

                count_shp=count_shp+1;

            end
        end

        % Make the dbf specification object
        dbfspec=struct;
        dbfspec.ID.FieldName='ID';
        dbfspec.ID.FieldType='N';
        dbfspec.ID.FieldLength=5;
        dbfspec.ID.FieldDecimalCount=0;

        % write a transects shapefile & display progress
        shapewrite(SHP,[OUTPUT_DIR,filesep,'shp',filesep,'potential_storm_erosion_uncertainty_SLR_',num2str(SL_scenario*100,3),'cm_TrgxSlope',num2str(BRUUN_CASE),'.shp'],'DbfSpec',dbfspec);

        % Make shapefile of shoreline change hazard zone  %%%%%%%%%%%%%%%%%%%%%

        % Make a shape file struct
        SHP=struct;

        count_shp=1;
        for i=1:length(lat_cells_haz_min)

            [lattmp_min,lontmp_min] = polysplit(lat_cells_haz_min{i},lon_cells_haz_min{i});
            [lattmp_max,lontmp_max] = polysplit(lat_cells_haz_max{i},lon_cells_haz_max{i});

            % code to handle model output transect crossings
            id_x_cells_nan=double(id_x_cells{i});                     % get the transect crossing id vector for each cell and convert it to an int, so that NaN's can be inserted
            id_x_cells_nan(isnan(lat_cells_haz_min{i}))=NaN;         % impose the same NaN structure as the model output
            [id_x_tmp,~] = polysplit(id_x_cells_nan,id_x_cells_nan); % and (thus) split the id's in the same form as above

            for j=1:length(lattmp_min)

                id_not_x=~logical(id_x_tmp{j});
                if ~any(id_not_x)
                    continue;
                end

                % Geometry
                SHP(count_shp,1).Geometry='Polygon';

                %lon_shp=[lontmp_min{j}; flipud(lontmp_max{j}); lontmp_min{1}; NaN];
                %lat_shp=[lattmp_min{j}; flipud(lattmp_max{j}); lattmp_min{1}; NaN];

                lon_shp=[lontmp_min{j}(id_not_x); flipud(lontmp_max{j}(id_not_x)); NaN];
                lat_shp=[lattmp_min{j}(id_not_x); flipud(lattmp_max{j}(id_not_x)); NaN];

                % coordiantes
                SHP(count_shp,1).X=lon_shp;
                SHP(count_shp,1).Y=lat_shp;

                % bounding box
                xmin=min(lon_shp);
                xmax=max(lon_shp);

                ymin=min(lat_shp);
                ymax=max(lat_shp);

                % get bounding box
                SHP(count_shp,1).BoundingBox=[xmin ymin; xmax ymax];

                % ID
                SHP(count_shp,1).ID=count_shp;

                count_shp=count_shp+1;

            end
        end

        % Make the dbf specification object
        dbfspec=struct;
        dbfspec.ID.FieldName='ID';
        dbfspec.ID.FieldType='N';
        dbfspec.ID.FieldLength=5;
        dbfspec.ID.FieldDecimalCount=0;

        % write a transects shapefile & display progress
        shapewrite(SHP,[OUTPUT_DIR,filesep,'shp',filesep,'shoreline_change_hazard_zone_SLR_',num2str(SL_scenario*100,3),'cm_TrgxSlope',num2str(BRUUN_CASE),'.shp'],'DbfSpec',dbfspec);

        % Make shapefile of the extreme storm hazard zone  %%%%%%%%%%%%%%%%%%%%

        % Make a shape file struct
        SHP=struct;

        count_shp=1;
        for i=1:length(lat_cells_storm_haz_min)

            [lattmp_min,lontmp_min] = polysplit(lat_cells_storm_haz_min{i},lon_cells_storm_haz_min{i});
            [lattmp_max,lontmp_max] = polysplit(lat_cells_storm_haz_max{i},lon_cells_storm_haz_max{i});

            % code to handle model output transect crossings
            id_x_cells_nan=double(id_x_cells{i});                     % get the transect crossing id vector for each cell and convert it to an int, so that NaN's can be inserted
            id_x_cells_nan(isnan(lat_cells_storm_haz_min{i}))=NaN;   % impose the same NaN structure as the model output
            [id_x_tmp,~] = polysplit(id_x_cells_nan,id_x_cells_nan); % and (thus) split the id's in the same form as above

            for j=1:length(lattmp_min)

                id_not_x=~logical(id_x_tmp{j});
                if ~any(id_not_x)
                    continue;
                end

                % Geometry
                SHP(count_shp,1).Geometry='Polygon';

                %lon_shp=[lontmp_min{j}; flipud(lontmp_max{j}); lontmp_min{1}; NaN];
                %lat_shp=[lattmp_min{j}; flipud(lattmp_max{j}); lattmp_min{1}; NaN];

                lon_shp=[lontmp_min{j}(id_not_x); flipud(lontmp_max{j}(id_not_x)); NaN];
                lat_shp=[lattmp_min{j}(id_not_x); flipud(lattmp_max{j}(id_not_x)); NaN];

                % coordiantes
                SHP(count_shp,1).X=lon_shp;
                SHP(count_shp,1).Y=lat_shp;

                % bounding box
                xmin=min(lon_shp);
                xmax=max(lon_shp);

                ymin=min(lat_shp);
                ymax=max(lat_shp);

                % get bounding box
                SHP(count_shp,1).BoundingBox=[xmin ymin; xmax ymax];

                % ID
                SHP(count_shp,1).ID=count_shp;

                count_shp=count_shp+1;

            end
        end

        % Make the dbf specification object
        dbfspec=struct;
        dbfspec.ID.FieldName='ID';
        dbfspec.ID.FieldType='N';
        dbfspec.ID.FieldLength=5;
        dbfspec.ID.FieldDecimalCount=0;

        % write a transects shapefile & display progress
        shapewrite(SHP,[OUTPUT_DIR,filesep,'shp',filesep,'extreme_storm_hazard_zone_SLR_',num2str(SL_scenario*100,3),'cm_TrgxSlope',num2str(BRUUN_CASE),'.shp'],'DbfSpec',dbfspec);

        % Make shapefile of the unresolved process uncertainty  %%%%%%%%%%%%%%%

        % Make a shape file struct
        SHP=struct;

        count_shp=1;
        for i=1:length(lat_cells_unresolved_min)

            [lattmp_min,lontmp_min] = polysplit(lat_cells_unresolved_min{i},lon_cells_unresolved_min{i});
            [lattmp_max,lontmp_max] = polysplit(lat_cells_unresolved_max{i},lon_cells_unresolved_max{i});

            % code to handle model output transect crossings
            id_x_cells_nan=double(id_x_cells{i});                     % get the transect crossing id vector for each cell and convert it to an int, so that NaN's can be inserted
            id_x_cells_nan(isnan(lat_cells_unresolved_min{i}))=NaN;   % impose the same NaN structure as the model output
            [id_x_tmp,~] = polysplit(id_x_cells_nan,id_x_cells_nan); % and (thus) split the id's in the same form as above

            for j=1:length(lattmp_min)

                id_not_x=~logical(id_x_tmp{j});
                if ~any(id_not_x)
                    continue;
                end                

                % Geometry
                SHP(count_shp,1).Geometry='Polygon';

                %lon_shp=[lontmp_min{j}; flipud(lontmp_max{j}); lontmp_min{1}; NaN];
                %lat_shp=[lattmp_min{j}; flipud(lattmp_max{j}); lattmp_min{1}; NaN];

                lon_shp=[lontmp_min{j}(id_not_x); flipud(lontmp_max{j}(id_not_x)); NaN];
                lat_shp=[lattmp_min{j}(id_not_x); flipud(lattmp_max{j}(id_not_x)); NaN];

                % coordiantes
                SHP(count_shp,1).X=lon_shp;
                SHP(count_shp,1).Y=lat_shp;

                % bounding box
                xmin=min(lon_shp);
                xmax=max(lon_shp);

                ymin=min(lat_shp);
                ymax=max(lat_shp);

                % get bounding box
                SHP(count_shp,1).BoundingBox=[xmin ymin; xmax ymax];

                % ID
                SHP(count_shp,1).ID=count_shp;

                count_shp=count_shp+1;

            end
        end

        % Make the dbf specification object
        dbfspec=struct;
        dbfspec.ID.FieldName='ID';
        dbfspec.ID.FieldType='N';
        dbfspec.ID.FieldLength=5;
        dbfspec.ID.FieldDecimalCount=0;

        % write a transects shapefile & display progress
        shapewrite(SHP,[OUTPUT_DIR,filesep,'shp',filesep,'unresolved_process_uncertainty_SLR_',num2str(SL_scenario*100,3),'cm_TrgxSlope',num2str(BRUUN_CASE),'.shp'],'DbfSpec',dbfspec);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    fprintf('done.\n');

end

tic
fprintf('writing .kml to file ... ');

pwd1=pwd;
cd([OUTPUT_DIR,filesep,'google_earth'])
kmlout.save;
cd(pwd1);

fprintf('done. ');
toc

% get a handle to the .kml text file
tic
fprintf('reading in saved .kml (text) file ... ');
[FID, msg] = fopen([OUTPUT_DIR,filesep,'google_earth',filesep,outputfilename,'.kml'],'rt');
if FID<0
    error(msg)
end
txt = fread(FID,'uint8=>char')';  % read the txt
fclose(FID);                      % close the file
fprintf('done. ');
toc

new_txt=txt;

tic
fprintf('modifying output kml ...');

% add legend code to kml
txtold=['<name>',outputfilename,'</name>'];
txtnew=[txtold legend_code];
new_txt=strrep(new_txt,txtold,txtnew);

% add product information
txtold=['<name>',outputfilename,'</name>'];
txtnew=[txtold ['\n\t<Placemark>\n\t<name>Product Information</name>\n\t\t<Snippet maxLines="0"></Snippet>\n\t\t<description>',product_information,'\n\t\t</description><styleUrl>#m_ylw-pushpin</styleUrl>\n\t</Placemark>'] ];
new_txt=strrep(new_txt,txtold,txtnew);

% add line to default the kml to open
txtold=['<name>',outputfilename,'</name>'];
txtnew=[txtold '\n\t<open>1</open>'];
new_txt=strrep(new_txt,txtold,txtnew);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add folder information ...

% ... to the 'model_initial_conditions' folder
txtold='<name>model_initial_conditions</name>';
txtnew=[txtold ['<Snippet maxLines="0"></Snippet><description>',init_conditions_folder_description,'</description>'] ];
new_txt=strrep(new_txt,txtold,txtnew);

% ... to the 'initial_shoreline' folder
txtold='<name>initial_shoreline</name>';
txtnew=[txtold ['<Snippet maxLines="0"></Snippet><description>',init_shoreline_folder_description,'</description>'] ];
new_txt=strrep(new_txt,txtold,txtnew);

% ... to the 'landward_model_boundary' folder
if HOLD_THE_LINE
    txtold='<name>landward_model_boundary</name>';
    txtnew=[txtold ['<Snippet maxLines="0"></Snippet><description>',landward_bnd_folder_description,'</description>'] ];
    new_txt=strrep(new_txt,txtold,txtnew);
end

% ... to the 'transects' folder
txtold='<name>transects</name>';
txtnew=[txtold ['<Snippet maxLines="0"></Snippet><description>',transects_folder_description,'</description>'] ];
new_txt=strrep(new_txt,txtold,txtnew);

% ... to the 'modeled_shoreline' folder
txtold='<name>modeled_shoreline</name>';
txtnew=[txtold ['<Snippet maxLines="0"></Snippet><description>',modeled_shoreline_folder_description,'</description>'] ];
new_txt=strrep(new_txt,txtold,txtnew);

% ... to the 'modeled_shoreline_uncertainty' folder
txtold='<name>modeled_shoreline_uncertainty</name>';
txtnew=[txtold ['<Snippet maxLines="0"></Snippet><description>',modeled_shoreline_uncertainty_folder_description,'</description>'] ];
new_txt=strrep(new_txt,txtold,txtnew);

% ... to the 'potential_storm_erosion_uncertainty' folder
txtold='<name>potential_storm_erosion_uncertainty</name>';
txtnew=[txtold ['<Snippet maxLines="0"></Snippet><description>',potential_storm_erosion_uncertainty_folder_description,'</description>'] ];
new_txt=strrep(new_txt,txtold,txtnew);

% we don't output these hazard zones for FloSup ... for some reason
if ~strcmp(Model_name,'FloSup')

    % ... to the 'shoreline_change_hazard_zone' folder
    txtold='<name>shoreline_change_hazard_zones</name>';
    txtnew=[txtold ['<Snippet maxLines="0"></Snippet><description>',shoreline_change_hazard_zone_folder_description,'</description>'] ];
    new_txt=strrep(new_txt,txtold,txtnew);

    % ... to the 'shoreline_change_hazard_zone' folder
    txtold='<name>extreme_storm_hazard_zones</name>';
    txtnew=[txtold ['<Snippet maxLines="0"></Snippet><description>',extreme_storm_hazard_zone_folder_description,'</description>'] ];
    new_txt=strrep(new_txt,txtold,txtnew);

end

% ... to the 'unresolved_process_uncertainty' folder
txtold='<name>unresolved_process_uncertainty</name>';
txtnew=[txtold ['<Snippet maxLines="0"></Snippet><description>',unresolved_process_uncertainty_folder_description,'</description>'] ];
new_txt=strrep(new_txt,txtold,txtnew);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add the 'Snippet maxLines=0' piece to tighten up display text for different kml placemarks ...

% ... for initial shoreline
for i=1:count_init_shoreline
    txtold=['<name>initial_shoreline_#',num2str(i),'</name>'];
    txtnew=[txtold '\n\t\t\t\t<Snippet maxLines="0"></Snippet>'];
    new_txt=strrep(new_txt,txtold,txtnew);
end

% ... for modeled shoreline
for i=1:count_mod_shoreline
    txtold=['<name>modeled_shoreline_#',num2str(i),'</name>'];
    txtnew=[txtold '\n\t\t\t\t<Snippet maxLines="0"></Snippet>'];
    new_txt=strrep(new_txt,txtold,txtnew);
end

% ... for modeled shoreline uncertainty (upper)
for i=1:count_mod_shoreline_uncert1
    txtold=['<name>(upper_bound)_modeled_shoreline_uncertainty_#',num2str(i),'</name>'];
    txtnew=[txtold '\n\t\t\t\t<Snippet maxLines="0"></Snippet>'];
    new_txt=strrep(new_txt,txtold,txtnew);
end

% ... for modeled shoreline uncertainty (lower)
for i=1:count_mod_shoreline_uncert2
    txtold=['<name>(lower_bound)_modeled_shoreline_uncertainty_#',num2str(i),'</name>'];
    txtnew=[txtold '\n\t\t\t\t<Snippet maxLines="0"></Snippet>'];
    new_txt=strrep(new_txt,txtold,txtnew);
end

% ... for potential storm erosion uncerainty (1-yr storm)
for i=1:count_mod_shoreline_ero_uncert1
    txtold=['<name>potential_storm_erosion_uncertainty_(1-yr_storm)_#',num2str(i),'</name>'];
    txtnew=[txtold '\n\t\t\t\t<Snippet maxLines="0"></Snippet>'];
    new_txt=strrep(new_txt,txtold,txtnew);
end

% ... for potential storm erosion uncerainty (20-yr storm)
for i=1:count_mod_shoreline_ero_uncert2
    txtold=['<name>potential_storm_erosion_uncertainty_(20-yr_storm)_#',num2str(i),'</name>'];
    txtnew=[txtold '\n\t\t\t\t<Snippet maxLines="0"></Snippet>'];
    new_txt=strrep(new_txt,txtold,txtnew);
end

% ... for potential storm erosion uncerainty (100-yr storm)
for i=1:count_mod_shoreline_ero_uncert3
    txtold=['<name>potential_storm_erosion_uncertainty_(100-yr_storm)_#',num2str(i),'</name>'];
    txtnew=[txtold '\n\t\t\t\t<Snippet maxLines="0"></Snippet>'];
    new_txt=strrep(new_txt,txtold,txtnew);
end

% we don't output these hazard zones for FloSup ... for some reason
if ~strcmp(Model_name,'FloSup')

    % ... for shoreline change hazard zone uncertainty
    for i=1:count_mod_shoreline_haz1
        txtold=['<name>shoreline_change_hazard_zone_#',num2str(i),'</name>'];
        txtnew=[txtold '\n\t\t\t\t<Snippet maxLines="0"></Snippet>'];
        new_txt=strrep(new_txt,txtold,txtnew);
    end

    % ... for shoreline change hazard zone uncertainty
    for i=1:count_mod_shoreline_haz2
        txtold=['<name>extreme_storm_hazard_zone_#',num2str(i),'</name>'];
        txtnew=[txtold '\n\t\t\t\t<Snippet maxLines="0"></Snippet>'];
        new_txt=strrep(new_txt,txtold,txtnew);
    end

end

% ... for unresolved process uncertainty (upper)
for i=1:count_unresolved_process_uncert1
    txtold=['<name>(upper_bound)_unresolved_process_uncertainty_#',num2str(i),'</name>'];
    txtnew=[txtold '\n\t\t\t\t<Snippet maxLines="0"></Snippet>'];
    new_txt=strrep(new_txt,txtold,txtnew);
end

% ... for unresolved process uncertainty (lower)
for i=1:count_unresolved_process_uncert2
    txtold=['<name>(lower_bound)_unresolved_process_uncertainty_#',num2str(i),'</name>'];
    txtnew=[txtold '\n\t\t\t\t<Snippet maxLines="0"></Snippet>'];
    new_txt=strrep(new_txt,txtold,txtnew);
end

% ... for the landward model boundary
if HOLD_THE_LINE

    % find the landward model boundary text
    expr='<Folder>\s+<name>landward_model_boundary</name>.*?</Folder>';
    [tr_id1,tr_id2]=regexp(new_txt,expr);
    tr_txt=new_txt(tr_id1:tr_id2);

    for i=1:Ntr_output
        if ~isnan(lat_MIN(i))
            txtold=['<name>landward_model_boundary_(transect_ID:',num2str(ID(i)),')</name>'];
            txtnew=[txtold '\n\t\t\t<Snippet maxLines="0"></Snippet>'];
            tr_txt=strrep(tr_txt,txtold,txtnew);
        end
    end

    new_txt=strcat(new_txt(1:tr_id1-1),tr_txt,new_txt(tr_id2+1:end));
end

% ... for the transects

% find the transect text
expr='<Folder>\s+<name>transects</name>.*?</Folder>';
[tr_id1,tr_id2]=regexp(new_txt,expr);
tr_txt=new_txt(tr_id1:tr_id2);

for i=1:Ntr_output % for all transects
    if strcmp(transects(ids2output(i)).model_type,'full model') || strcmp(transects(ids2output(i)).model_type,'cross-shore only') || strcmp(transects(ids2output(i)).model_type,'rate only') || ~bool_no_prediction(ids2output(i))

        % seperate "no prediction" transects since the no prediction can be
        % turned on later to surpress output for certain transects
        if bool_no_prediction(ids2output(i))
            continue
        end

        txtold=['<name>transect_ID:',num2str(ID(ids2output(i))),'</name>'];      
        txtnew=[txtold '\n\t\t\t<Snippet maxLines="0"></Snippet>'];
        tr_txt=strrep(tr_txt,txtold,txtnew);
    end
end

new_txt=strcat(new_txt(1:tr_id1-1),tr_txt,new_txt(tr_id2+1:end));

fprintf('done. ');
toc

tic 
fprintf('(over)writing new .kml (text) file ... ');
fid = fopen([OUTPUT_DIR,filesep,'google_earth',filesep,outputfilename,'.kml'],'wt');
fprintf(fid,new_txt);
fclose(fid);
fprintf('done. ');
toc