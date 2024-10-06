% MAKE SOME ADJUSTMENTS TO THE ADDITIVE NOISE

% reduce additive noise for transects with very little data
id=n_obs<25;
NOISE_FAC(id)=0.5;

% e.g., reduce additive noise for transects where model uncertainty seems to blow up

% % some ID's with initially large RMSE
% ID_RMSE=[]';
% 
% % reduce the addative noise to transects with initially large RMSE
% [id,~]=ismember(ID,ID_RMSE); NOISE_FAC(id)=0.01;