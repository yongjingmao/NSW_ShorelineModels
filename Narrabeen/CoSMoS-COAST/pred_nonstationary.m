rho = 1025; % Density of water
g = 9.81; %
gamma = 0.78; % Shallow water Hsb/h

w.Hsb2 = w.Hsb.^2; % Hsb square
w.power = (1/16)*rho*g*w.Hsb2*sqrt(g)*(1/sqrt(gamma)).*(w.Hsb).^0.5; % Wave power

dt_wave = t1(2) - t1(1); % dt of wave data in days
% dummyC = sqrt(2229)*24*3600;
% PHI = log10(PHIOUT);
% CA = CAOUT*dummyC;
% CE = CEOUT*dummyC;
DT = DTOUT;
DY = DYOUT;
HSB = HSBOUT;



%=========================================================================
% Create covariates (the moving std and average of wave variables)
%=========================================================================
N_months = 6: 3: 120; % Number of months used for rolling
roll_windows = int64(N_months/12*365./dt_wave); % Lenght of rolling window in days
num_window = length(roll_windows); % Number of rolling windows
[num_tran, num_step] = size(Hsb); % Number of transects and number of simulation time steps
num_step0 = size(w.Hsb, 2); % Number of raw wave data time steps; row wave data starts ealier than the simulation

var_names = {'Hsb', 'Hsb2', 'omega'};
for var_idx=1:length(var_names)                                     % for each assimilated variable ...
    eval([var_names{var_idx},'_avg = NaN(num_window, num_tran, num_step0);']); % Averaged Hsb with different rolling windows [num_windows, num_transects, timesteps]
    % example Hsb_avg = NaN(num_window, num_tran, num_step);
    eval([var_names{var_idx},'_std = NaN(num_window, num_tran, num_step0);']);  % Std Hsb with different rolling windows [num_windows, num_transects, timesteps]
    % example Hsb_std = NaN(num_window, num_tran, num_step);   
end


for window=1:length(roll_windows)
    roll_window = roll_windows(window);
    for var_idx=1:length(var_names) 
        varname = var_names{var_idx};
        eval([varname '_avg(window, :, :)=[NaN(size(w.' varname ', 1), roll_window-1), movmean(w.' varname ', roll_window, 2, "Endpoints", "discard")];']); % Calculate rolling average
        % example Hsb_avg(i, :, :) = [NaN(size(w.Hsb, 1), roll_window-1), movmean(w.Hsb, roll_window, 2, "Endpoints", "discard")];
        eval([varname '_std(window, :, :)=[NaN(size(w.' varname ', 1), roll_window-1), movstd(w.' varname ', roll_window, 0, 2, "Endpoints", "discard")];']); % Calculate rolling std
        % example Hsb_std(i, :, :) = [NaN(size(w.Hsb, 1), roll_window-1), movstd(w.Hsb, roll_window, 0, 2, "Endpoints", "discard")];
    end
end

%%
%=========================================================================
% Calculate correlation coefficients
%=========================================================================
% Starting time of calibration. Nonstationary model also requires sometime to stabilize
n_stable = 1; % This value is arbitrary so far but can work to determine a trend changing point
t_cali = t_obs((t_obs>=t_obs(n_stable))&(t_obs<tforecast)); % Time step for calibration

targ_vars = {'DY', 'DT', 'HSB'}; % Target variables
pred_vars = {}; % Pred variables
for var_idx=1:length(var_names) 
    varname = var_names{var_idx};
    pred_vars{2*(var_idx-1)+1} = [varname, '_avg'];
    pred_vars{2*(var_idx-1)+2} = [varname, '_std'];
end

num_preds = length(pred_vars); %num of predictors

for t_var=1:length(targ_vars)
    targ_var = targ_vars{t_var};
    eval([targ_var '_corrs= NaN(num_preds, num_window, num_tran);']) % Save correlations [num_preds, num_windows, num_transects]
    % example CA_corrs = NaN(num_preds, num_window, num_tran);
    for p_var=1:length(pred_vars)
        pred_var = pred_vars{p_var};

        % Loop over each window and transect to calculate corrcoef
        for window = 1:num_window
            for tran = 1:num_tran
                predictor = eval(['squeeze(' pred_var '(window, tran, ismember(w.t, t_cali)))']);
                %example: pred = squeeze(Hsb_avg(window, tran, ismember(w.t, t_cali)))';  % A_window is now of size [N_transects, time]
                targ = eval([targ_var '(tran, ismember(t_output, t_cali))']);
                %example: targ = CA(tran, ismember(t_output, t_cali));
                predictor = predictor';

                % Calculate the correlation coefficients between A_window and B
                corr = corrcoef(predictor(~isnan(predictor)), targ(~isnan(predictor)));
                eval([targ_var '_corrs(p_var, window, tran) = corr(1, 2);'])
                %example CA_corrs(p_var, window, tran) = corr(1, 2);
            end
        end
    end
end


%% Find combinations of CA, CE and PHI
DT_cali = NaN(num_tran, Nens, num_step); % CA predictions
DY_cali = NaN(num_tran, Nens, num_step); % CE predictions
HSB_cali = NaN(num_tran, Nens, num_step); % PHI predictions

for tran = 1:num_tran
    count_scen = 0;
    corr_thd = 0.7;
    while count_scen<Nens
        DT_filter = abs(DT_corrs(:, :, tran))>=corr_thd;
        DT_linearIndices = find(DT_filter);
        [DT_pred_idx, DT_window_idx] = ind2sub(size(DT_corrs), DT_linearIndices);

        DY_filter = abs(DY_corrs(:, :, tran))>=corr_thd;
        DY_linearIndices = find(DY_filter);
        [DY_pred_idx, DY_window_idx] = ind2sub(size(DY_corrs), DY_linearIndices);

        HSB_filter = abs(HSB_corrs(:, :, tran))>=corr_thd;
        HSB_linearIndices = find(HSB_filter);
        [HSB_pred_idx, HSB_window_idx] = ind2sub(size(HSB_corrs), HSB_linearIndices);

        count_scen = length(HSB_linearIndices)*length(DT_linearIndices)*length(DY_linearIndices);

        if count_scen<Nens
            corr_thd = corr_thd-0.01;
        else
            count = 1;
            for DT_id = 1:length(DT_linearIndices)
                for DY_id = 1:length(DY_linearIndices)
                    for HSB_id = 1:length(HSB_linearIndices)
                        % CA_pred_idxs(tran, count) = CA_pred_idx(ca_id);
                        % CA_window_idxs(tran, count) = CA_window_idx(ca_id);
                        % CE_pred_idxs(tran, count) = CE_pred_idx(ce_id);
                        % CE_window_idxs(tran, count) = CE_window_idx(ce_id);
                        % PHI_pred_idxs(tran, count) = PHI_pred_idx(phi_id);
                        % PHI_window_idxs(tran, count) = PHI_window_idx(phi_id);
                        % 
                        if count <= 200
                            for var_idx=1:length(targ_vars)
                                targ_var = targ_vars{var_idx};
                                pred_var = eval(['pred_vars{' targ_var '_pred_idx(' targ_var '_id)}']);
                                window = eval([targ_var '_window_idx(' targ_var '_id)']);
                                predictor = eval(['squeeze(' pred_var '(window, tran, ismember(w.t, t_cali)))']);
                                targ = eval([targ_var '(tran, ismember(t_output, t_cali))']);
                                targ = targ';
                                [brob,stats] = robustfit(predictor(~isnan(predictor)), targ(~isnan(predictor)));

                                % Apply the fitted model to the full time series
                                model = brob(1)+brob(2)*eval(['squeeze(' pred_var '(window, tran, ismember(w.t, t1)))']);
                                eval([targ_var '_cali(tran, count, :) = model;'])
                            end
                        end

                        count = count+1;

                    end
                end
            end
        end

    end
end


%%
%=========================================================================
% Pred targs
%=========================================================================
% corr_thd = 0.7; % Only combinations of corr above the threshold will be evaluated
% 
% for var_idx=1:length(targ_vars)
%     targ_var = targ_vars{var_idx};
%     var_corr = eval([targ_var '_corrs']);
% 
%     %corr_thd = prctile(abs(var_corr), 75, [1,2]); % Use 99 percentile as the threshold
%     corr_thd = prctile(abs(var_corr), 100-Nens/(size(var_corr, 1)*size(var_corr, 2))*100, [1,2]); % Use 99 percentile as the threshold
%     var_filter = abs(var_corr)>=corr_thd;
%     linearIndices = find(var_filter);
%     [pred_idxs, window_idxs, tran_idxs] = ind2sub(size(var_corr), linearIndices);
%     num_scens = sum(tran_idxs==1); % Number of variable combinations for each transect after filtering
% 
%     model_var = NaN(num_tran, num_scens, num_step); % Slope of linear regression
%     for tran = 1:num_tran
%         pred_idx = pred_idxs(tran_idxs==tran);
%         window_idx = window_idxs(tran_idxs==tran);
%         for scen = 1:num_scens
%             pred_var = pred_vars{pred_idx(scen)};
%             window = window_idx(scen);
%             predictor = eval(['squeeze(' pred_var '(window, tran, ismember(w.t, t_cali)))']);
%             targ = eval([targ_var '(tran, ismember(t_output, t_cali))']);
%             targ = targ';
% 
%             % Fit linear regression
%             [brob,stats] = robustfit(predictor(~isnan(predictor)), targ(~isnan(predictor)));
% 
%             % Apply the fitted model to the full time series
%             model = brob(1)+brob(2)*eval(['squeeze(' pred_var '(window, tran, ismember(w.t, t1)))']);
%             model_var(tran, scen, :) = model;
% 
%             % figure
%             % scatter(x,y,'filled')
%             % hold on
%             % plot(x,brob(1)+brob(2)*x,'g')
%             % hold off
%             % xlabel('x')
%             % ylabel('y')
%             % legend('Data','Outlier','Ordinary Least Squares','Robust Regression')
%             % grid on
% 
%         end
%     end
% 
%     % Generate samples from model_var
%     % ensemble_var = NaN(num_tran, Nens, num_step);
%     % 
%     % for i = 1:num_tran
%     %     for j = 1:num_step
%     %         % For each transect and time step, sample from the num_sample domain
%     %         ensemble_var(i, :, j) = randsample(model_var(i, :, j), Nens, true);
%     %     end
%     % end
% 
%     %eval([targ_var '_cali = ensemble_var;'])
%     eval([targ_var '_cali = model_var;'])
%     %eval([targ_var '_cali = squeeze(mean(model_var, 2, "omitnan"));']) % Average over different scens
% end


% Plot results
for id_tr=1:num_tran

    f1=figure(1); 
    set(f1,'PaperPositionMode','Auto','Position',[100 0 800 800],'color','w', 'Visible','off');

    ALPHA=0.25;

    tmore=60;

    xp0=0.1;
    yp0=0.7;
    width=0.805;
    height=0.2;
    ysep=0.05;

    for var_idx=1:length(targ_vars)
        targ_var = targ_vars{var_idx};
        var_corr = eval([targ_var '_corrs']);
        var_cali = eval([targ_var '_cali']);

        [M, I] = max(squeeze(abs(var_corr(:, :, id_tr))), [], 'all');
        [I_row, I_col] = ind2sub(size(squeeze(abs(var_corr(:, :, id_tr)))), I); 
        wave_var = pred_vars{I_row};
        wave_labels = {'$\bar{H_{s,b}}$', '$std(H_{s,b})$', '$\bar{H_{s,b}^{2}}$', '$std(H_{s,b}^{2})$', '$\bar{\Omega}$', '$std(\Omega)$',  };
        wave_label = wave_labels{I_row};


        roll_window = roll_windows(I_col);
        predictor = squeeze(eval([wave_var '(I_col, id_tr, ismember(w.t, t1))']))';
        target = eval([targ_var '(id_tr, :)']);
        [brob,stats] = robustfit(predictor(ismember(t1, t_cali)), target(ismember(t_output, t_cali)));
        model = brob(1)+brob(2)*predictor;

        nplot=var_idx; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;

        han3=plot(t1, squeeze(var_cali(id_tr, :, :)),'-r'); 
        han2=plot(t1, model,'-r'); 
        han1=plot(t_output, target,'-k');


        set(han3, 'Color', [1 0 0 0.01]);
        set(han1,'linewidth',2);
        han3=xline(tforecast, '--k'); set(han3,'linewidth',2);

        ylabel(targ_var);

        yyaxis right
        han1_1= plot(t1, predictor,'-b'); 
        ylabel(wave_label, 'Interpreter', 'latex')
        set(gca, 'YColor', 'k');
        if var_idx == 1
            legend([han1, han2, han1_1], {'Target', 'Model', 'Wave'}, 'Location', 'northwest');
        end
        datetick('x','keeplimits','keepticks');
        title(sprintf('%s (Corr: %.2f)', targ_var, M))
    end
    saveas(figure(f1),[OUTPUT_DIR,filesep,'figures\CoSMoS-COAST_corrs',num2str(ID(id_tr),'%0.5d'),'.png'],'png');
    close all;
end

