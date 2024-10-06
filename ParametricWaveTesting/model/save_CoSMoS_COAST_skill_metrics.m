%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save CoSMoS-COAST skill metrics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% additional performance metric variables to save
NOBS=NaN(Ntr,1);
NOBS_cal=NaN(Ntr,1);
NOBS_val=NaN(Ntr,1);
R=NaN(Ntr,1);
MSE=NaN(Ntr,1);
RMSE=NaN(Ntr,1);
RMSM=NaN(Ntr,1);
MAE=NaN(Ntr,1);
SKILL=NaN(Ntr,1);
IDX_AGREEMENT=NaN(Ntr,1);
LAMBDA=NaN(Ntr,1);

NOBS_sat=NaN(Ntr,1);
NOBS_cal_sat=NaN(Ntr,1);
NOBS_val_sat=NaN(Ntr,1);
R_sat=NaN(Ntr,1);
MSE_sat=NaN(Ntr,1);
RMSE_sat=NaN(Ntr,1);
RMSM_sat=NaN(Ntr,1);
MAE_sat=NaN(Ntr,1);
SKILL_sat=NaN(Ntr,1);
IDX_AGREEMENT_sat=NaN(Ntr,1);
LAMBDA_sat=NaN(Ntr,1);

PCT_sat=NaN(Ntr,1);

for i=1:Ntr
    
    % get the indices of the non-satellite observations
    nn=~isnan(Y_obs(i,:)) & ~SAT(i,:);
    
    % RMS error during forecast period
    tobs=t_obs(nn);
    obs=Y_obs(i,nn)-Y00(i);
    
    NOBS(i)=length(obs);

    % # of obs during calibration period
    id=(tobs>t0 | tobs<tforecast);
    NOBS_cal(i)=length(obs(id));

    % remove data that is outside of validation period
    
    tobs(id)=[];
    obs(id)=[];
    
    NOBS_val(i)=length(obs);
    
    % if observations exist
    if NOBS_val(i)>=2
        
        % model output
        id=~isnan(t_output);
        tmod=t_output(id);
        model=YY(i,id)-Y00(i);
        
        % interp onto observations
        if length(model)>=3
            mod_obs=interp1(tmod,model,tobs);
        else
            mod_obs=NaN*obs;
        end
        
        % the error associated with one data point
        y_error=obs-mod_obs;
        
        % error/skill metrics
        RM=corrcoef(mod_obs,obs);
        R(i)=RM(2,1);
        MSE(i)=nanmean((y_error).^2);
        RMSE(i)=sqrt(MSE(i));
        MAE(i)=nanmean(abs(y_error));
        %RMSM(i)=sqrt(nanmean((obs-obs(1)).^2));
        RMSM(i)=sqrt(nanmean((obs).^2));
        
        SKILL(i)=1-RMSE(i)/RMSM(i);
        
        IDX_AGREEMENT(i)=1-sum((obs-mod_obs).^2)/sum( (abs(mod_obs-nanmean(obs))+abs(obs-nanmean(obs))).^2 );
        
        LAMBDA(i)=1-MSE(i)/(var(obs)+var(mod_obs)+(nanmean(obs)-nanmean(mod_obs)).^2);
  
    end
    
    % get the indices of the satellite observations
    nn=~isnan(Y_obs(i,:)) & SAT(i,:);

    % RMS error during forecast period
    tobs=t_obs(nn);
    obs=Y_obs(i,nn)-Y00(i);
    
    NOBS_sat(i)=length(obs);

    % # of obs during calibration period
    id=(tobs>t0 & tobs<tforecast);
    NOBS_cal_sat(i)=length(obs(id));

    % remove data that is outside of validation period
    id=(tobs<tforecast | tobs>tforecast2);
    tobs(id)=[];
    obs(id)=[];
    
    NOBS_val_sat(i)=length(obs);
    
    % if observations exist
    if NOBS_val_sat(i)>=2
        
        % model output
        id=~isnan(t_output);
        tmod=t_output(id);
        model=YY(i,id)-Y00(i);
        obs_CI1=obs-2*y_rms_sat;
        obs_CI2=obs+2*y_rms_sat;
        
        % interp onto observations
        if length(model)>=3
            mod_obs=interp1(tmod,model,tobs);
        else
            mod_obs=NaN*obs;
        end
        
        % the error associated with one data point
        y_error=obs-mod_obs;
        
        % error/skill metrics
        RM_sat=corrcoef(mod_obs,obs);
        R(i)=RM_sat(2,1);
        MSE_sat(i)=nanmean((y_error).^2);
        RMSE_sat(i)=sqrt(MSE_sat(i));
        MAE_sat(i)=nanmean(abs(y_error));
        %RMSM_sat(i)=sqrt(nanmean((obs-obs(1)).^2));
        RMSM_sat(i)=sqrt(nanmean((obs).^2));
        
        SKILL_sat(i)=1-RMSE_sat(i)/RMSM_sat(i);
        
        IDX_AGREEMENT_sat(i)=1-sum((obs-mod_obs).^2)/sum( (abs(mod_obs-nanmean(obs))+abs(obs-nanmean(obs))).^2 );
        
        LAMBDA_sat(i)=1-MSE_sat(i)/(var(obs)+var(mod_obs)+(nanmean(obs)-nanmean(mod_obs)).^2);
        
        within_CI=(mod_obs>=obs_CI1 & mod_obs<=obs_CI2);
        PCT_sat(i)=(sum(within_CI)./length(mod_obs))*100;
        
    end
    
end