function [littoral_cell_name,littoral_cell_name_unique,littoral_cell_type,littoral_cell_tr_start,littoral_cell_tr_end,N_full_model_sections,N_cross_shore_only_sections,N_rate_only_sections,id_full_model_start,id_full_model_end,id_cross_shore_only_start,id_cross_shore_only_end,id_rate_only_start,id_rate_only_end]=transects2cells(transects)
% a function to take the transects and get the littoral cell information. SeanV. 

% number of transects
Ntr=length(transects);

% get start and end indices of the full model transects
littoral_cell_name={transects.littoral_cell}';                 % get littoral cell names
[littoral_cell_name_unique,ia,~]=unique(littoral_cell_name); % get unique littoral cell names

[~,id]=sort(ia);  % sort the cell names by transect id #
ia_sort=ia(id);
littoral_cell_name_unique=littoral_cell_name_unique(id); % ordered (south to north) unique littoral cell names

id=find(cellfun(@isempty,littoral_cell_name_unique));
littoral_cell_name_unique(id)=[];
ia_sort(id)=[];
littoral_cell_type={transects(ia_sort).model_type}';

N_littoral_cells=length(littoral_cell_name_unique);

% count up the littoral cell sections in the following loop (starting from zero)
N_full_model_sections=0; N_cross_shore_only_sections=0; N_rate_only_sections=0;

% inizialize the indices as empty and then they are filled in the following loop
id_full_model_start=[];        id_full_model_end=[];
id_cross_shore_only_start=[];  id_cross_shore_only_end=[];
id_rate_only_start=[];         id_rate_only_end=[];

littoral_cell_tr_start=NaN(N_littoral_cells,1);
littoral_cell_tr_end  =NaN(N_littoral_cells,1);

for i=1:N_littoral_cells
    
    id=find(strcmp(littoral_cell_name,littoral_cell_name_unique(i)));
    
    littoral_cell_tr_start(i)=min(id);
    littoral_cell_tr_end(i)  =max(id);
    
    if strcmp(littoral_cell_type{i},'full model') % if the current unique cell is "full model"
        
        id_full_model_start=cat(1,id_full_model_start,min(id));
        id_full_model_end  =cat(1,id_full_model_end,max(id));
        N_full_model_sections=N_full_model_sections+1;
        
    elseif strcmp(littoral_cell_type{i},'cross-shore only') % if the current unique cell is "cross-shore only"
        
        id_cross_shore_only_start=cat(1,id_cross_shore_only_start,min(id));
        id_cross_shore_only_end  =cat(1,id_cross_shore_only_end,max(id));
        N_cross_shore_only_sections=N_cross_shore_only_sections+1;
        
    elseif strcmp(littoral_cell_type{i},'rate only')  % if the current unique cell is "rate only"
        
        id_rate_only_start=cat(1,id_rate_only_start,min(id));
        id_rate_only_end  =cat(1,id_rate_only_end,max(id));
        N_rate_only_sections=N_rate_only_sections+1;
        
    else
        % do nothing?
    end
    
end