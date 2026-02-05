%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save CoSMoS-COAST figures (in parallel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('ids2plot','var') % if this 'ids2plot' variable exists, then plot those transects
    ids=1:Ntr;
else                        % else plot all of them
    [~,ids]=ismember(ids2plot,ID);
    ids(ids==0)=[]; % remote entries where the MOP ID doesn't exist
end

% create parallel pool
poolsize=8;
parpool(poolsize);

parfor id_tr=(ids)'

    if strcmp(transects(id_tr).model_type,'full model') || strcmp(transects(id_tr).model_type,'cross-shore only') || strcmp(transects(id_tr).model_type,'rate only')
        
        fprintf('working on transect ID # %d (%d of %d) ... \n',ID(id_tr),find(id_tr==ids),length(ids)); % display progress
        
        f1=1;
        f2=2;
        f3=3;
        f4=4;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Figure 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%a
        hanf1=figure(f1); set(hanf1,'PaperPositionMode','Auto','Position',[50 100 920 1000],'color','w','Visible','off');
        
        ALPHA=0.25;
        
        tmore=60;
        
        xp0=0.15;
        yp0=0.845;
        width=0.82;
        height=0.12;
        ysep=0.014;
        
        [idw]=ismember(t1(n),t);
        
        nplot=1; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;
        han1=plot(t,Hs(id_tr,:,1),'-b',t(idw),Hs(id_tr,idw,1),'ro'); set(han1,'MarkerFaceColor','r','MarkerSize',10);  set(gca,'XTick',linspace(t(1),t(min(n+tmore,length(t))),12),'Xlim',[t(1) t(min(n+tmore,length(t)))]); datetick('x','keeplimits','keepticks');
        
        vmin=0;
        vmax=1.1*max(Hs(id_tr,:,1));
        han1=plot([tforecast tforecast],[vmin vmax],'--k'); set(han1,'linewidth',2);
        han1=plot([tforecast2 tforecast2],[vmin vmax],'--k'); set(han1,'linewidth',2);
        axis([t0 t1(min(n+tmore,length(t1))) vmin vmax]);
        
        ylabel({'wave','height [m]'});
        set(gca,'FontSize',14,'XTickLabel',[]);
        ax = gca; ax.YRuler.Exponent = 0;
        title(sprintf('CoSMoS-COAST: model components (transect ID: %d)',ID(id_tr)),'FontSize',14);
        
        nplot=2; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;
        
        nn=~isnan(Y_obs(id_tr,:)) & SAT(id_tr,:)==1;
        han1=errorbar(t_obs(nn),Y_obs(id_tr,nn)-Y00(id_tr),2*Y_rms(id_tr,nn),'b'); set(han1,'linestyle','none','color',[0.4 0.6 1]);
        han1=plot(t_obs(nn),Y_obs(id_tr,nn)-Y00(id_tr),'bo'); set(han1,'MarkerSize',4,'MarkerFaceColor',[0.4 0.6 1]); datetick('x','keeplimits','keepticks')
        
        nn=~isnan(Y_obs(id_tr,:)) & SAT(id_tr,:)==0;
        han1=errorbar(t_obs(nn),Y_obs(id_tr,nn)-Y00(id_tr),2*Y_rms(id_tr,nn),'m'); set(han1,'linestyle','none','color',[0.278 0 0.7]);
        han1=plot(t_obs(nn),Y_obs(id_tr,nn)-Y00(id_tr),'mo'); set(han1,'MarkerSize',4,'MarkerFaceColor',[0.58 0.30 1],'MarkerEdgeColor',[0.278 0 0.7]); datetick('x','keeplimits','keepticks')
        
        nn=~isnan(Y_obs(id_tr,:));
        
        han1=fill([t_output(~isnan(t_output)) fliplr(t_output(~isnan(t_output)))],[YCI1(id_tr,~isnan(t_output))-Y00(id_tr) fliplr(YCI2(id_tr,~isnan(t_output))-Y00(id_tr))],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
        plot(t_output,YY(id_tr,:)-Y00(id_tr),'-r'); set(gca,'XTick',linspace(t1(1),t1(min(n+tmore,length(t1))),12),'Xlim',[t1(1) t1(min(n+tmore,length(t1)))]);
        
        han1=plot([t1(1) t1(min(n+tmore,length(t1)))],[Ymin(id_tr)-Y00(id_tr) Ymin(id_tr)-Y00(id_tr)],'k--'); set(han1,'LineWidth',2);
        
        vmin=min(YCI1(id_tr,~isnan(t_output))-Y00(id_tr)); if ~isempty(Y_obs(id_tr,nn)); vmin=min(vmin,min(Y_obs(id_tr,nn)-Y00(id_tr)-2*Y_rms(id_tr,nn))); end; if vmin<0; vmin=1.05*vmin; else; vmin=0.95*vmin; end
        vmax=max(YCI2(id_tr,~isnan(t_output))-Y00(id_tr)); if ~isempty(Y_obs(id_tr,nn)); vmax=max(vmax,max(Y_obs(id_tr,nn)-Y00(id_tr)+2*Y_rms(id_tr,nn))); end; if vmax<0; vmax=0.95*vmax; else; vmax=1.05*vmax; end
        
        han1=plot([tforecast tforecast],[vmin vmax],'--k'); set(han1,'linewidth',2);
        han1=plot([tforecast2 tforecast2],[vmin vmax],'--k'); set(han1,'linewidth',2);
        axis([t0 t1(min(n+tmore,length(t1))) vmin vmax]);
        
        ylabel({'total','shoreline','position [m]'});
        set(gca,'FontSize',14,'XTickLabel',[]);
        ax = gca; ax.YRuler.Exponent = 0;
        
        nplot=3; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;
        han1=fill([t_output(~isnan(t_output)) fliplr(t_output(~isnan(t_output)))],[YSTCI1(id_tr,~isnan(t_output)) fliplr(YSTCI2(id_tr,~isnan(t_output)))],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
        plot(t_output,YST(id_tr,:),'-r'); set(gca,'XTick',linspace(t1(1),t1(min(n+tmore,length(t1))),12),'Xlim',[t1(1) t1(min(n+tmore,length(t1)))]); datetick('x','keeplimits','keepticks');
        
        vmin=min(YSTCI1(id_tr,~isnan(t_output))); if vmin<0; vmin=1.2*vmin; else; vmin=0.9*vmin; end
        vmax=max(YSTCI2(id_tr,~isnan(t_output))); if vmax<0; vmax=0.9*vmax; else; vmax=1.2*vmax; end
        han1=plot([tforecast tforecast],[vmin vmax],'--k'); set(han1,'linewidth',2);
        han1=plot([tforecast2 tforecast2],[vmin vmax],'--k'); set(han1,'linewidth',2);
        axis([t0 t1(min(n+tmore,length(t1))) vmin vmax]);
        
        axis([t0 t1(min(n+tmore,length(t1))) vmin vmax]);
        ylabel({'short-term','position [m]'});
        set(gca,'FontSize',14,'XTickLabel',[]);
        ax = gca; ax.YRuler.Exponent = 0;
        
        nplot=4; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;
        han1=fill([t_output(~isnan(t_output)) fliplr(t_output(~isnan(t_output)))],[YBRUCI1(id_tr,~isnan(t_output)) fliplr(YBRUCI2(id_tr,~isnan(t_output)))],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
        plot(t_output,YBRU(id_tr,:),'-r'); set(gca,'XTick',linspace(t1(1),t1(min(n+tmore,length(t1))),12),'Xlim',[t1(1) t1(min(n+tmore,length(t1)))]); datetick('x','keeplimits','keepticks');
        
        vmin=min(YBRUCI1(id_tr,~isnan(t_output))); if vmin<0; vmin=1.2*vmin; else; vmin=0.9*vmin; end
        vmax=max(YBRUCI2(id_tr,~isnan(t_output))); if vmax<0; vmax=0.9*vmax; else; vmax=1.2*vmax; end
        han1=plot([tforecast tforecast],[vmin vmax],'--k'); set(han1,'linewidth',2);
        han1=plot([tforecast2 tforecast2],[vmin vmax],'--k'); set(han1,'linewidth',2);
        axis([t0 t1(min(n+tmore,length(t1))) vmin vmax]);
        
        ylabel({'Bruunian','recession [m]'});
        set(gca,'FontSize',14,'XTickLabel',[]);
        ax = gca; ax.YRuler.Exponent = 0;
        
        nplot=5; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;
        han1=fill([t_output(~isnan(t_output)) fliplr(t_output(~isnan(t_output)))],[YVLTCI1(id_tr,~isnan(t_output)) fliplr(YVLTCI2(id_tr,~isnan(t_output)))],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
        plot(t_output,YVLT(id_tr,:),'-r'); set(gca,'XTick',linspace(t1(1),t1(min(n+tmore,length(t1))),12),'Xlim',[t1(1) t1(min(n+tmore,length(t1)))]); datetick('x','keeplimits','keepticks');
        
        vmin=min(YVLTCI1(id_tr,~isnan(t_output))); if vmin<0; vmin=1.2*vmin; else; vmin=0.9*vmin; end
        vmax=max(YVLTCI2(id_tr,~isnan(t_output))); if vmax<0; vmax=0.9*vmax; else; vmax=1.2*vmax; end
        han1=plot([tforecast tforecast],[vmin vmax],'--k'); set(han1,'linewidth',2);
        han1=plot([tforecast2 tforecast2],[vmin vmax],'--k'); set(han1,'linewidth',2);
        axis([t0 t1(min(n+tmore,length(t1))) vmin vmax]);
        
        ylabel({'long-term','rate term [m]'});
        set(gca,'FontSize',14,'XTickLabel',[]);
        ax = gca; ax.YRuler.Exponent = 0;
        
        nplot=6; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;
        
        han1=fill([t_output(~isnan(t_output)) fliplr(t_output(~isnan(t_output)))],[YLSTCI1(id_tr,~isnan(t_output)) fliplr(YLSTCI2(id_tr,~isnan(t_output)))],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
        plot(t_output,YLST(id_tr,:),'-r'); set(gca,'XTick',linspace(t1(1),t1(min(n+tmore,length(t1))),12),'Xlim',[t1(1) t1(min(n+tmore,length(t1)))]); datetick('x','keeplimits','keepticks');
        
        YLST_steps=YLST(id_tr,:)-YLST_only(id_tr,:);
        plot(t_output,YLST_steps,'-m'); set(gca,'XTick',linspace(t1(1),t1(min(n+tmore,length(t1))),12),'Xlim',[t1(1) t1(min(n+tmore,length(t1)))]); datetick('x','keeplimits','keepticks');
        
        vmin=min(YLSTCI1(id_tr,~isnan(t_output))); vmin=min(vmin,min(YLST_steps)); if vmin<0; vmin=1.2*vmin; else; vmin=0.9*vmin; end
        vmax=max(YLSTCI2(id_tr,~isnan(t_output))); vmax=max(vmax,max(YLST_steps)); if vmax<0; vmax=0.9*vmax; else; vmax=1.2*vmax; end
        han1=plot([tforecast tforecast],[vmin vmax],'--k'); set(han1,'linewidth',2);
        han1=plot([tforecast2 tforecast2],[vmin vmax],'--k'); set(han1,'linewidth',2);
        axis([t0 t1(min(n+tmore,length(t1))) vmin vmax]);
        
        ylabel({'longshore','position [m]'});
        set(gca,'FontSize',14,'XTickLabel',[]);
        ax = gca; ax.YRuler.Exponent = 0;
        
        nplot=7; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;
        
        han1=fill([t_output(~isnan(t_output)) fliplr(t_output(~isnan(t_output)))],[YLST_only_CI1(id_tr,~isnan(t_output)) fliplr(YLST_only_CI2(id_tr,~isnan(t_output)))],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
        
        plot(t_output,YLST_only(id_tr,:),'-r',t_output,-YLST_steps,'-m'); set(gca,'XTick',linspace(t1(1),t1(min(n+tmore,length(t1))),12),'Xlim',[t1(1) t1(min(n+tmore,length(t1)))]); datetick('x','keeplimits','keepticks');
        
        vmin=min(YLST_only_CI1(id_tr,~isnan(t_output))); if vmin<0; vmin=1.2*vmin; else; vmin=0.9*vmin; end
        vmax=max(YLST_only_CI2(id_tr,~isnan(t_output))); if vmax<0; vmax=0.9*vmax; else; vmax=1.2*vmax; end
        han1=plot([tforecast tforecast],[vmin vmax],'--k'); set(han1,'linewidth',2);
        han1=plot([tforecast2 tforecast2],[vmin vmax],'--k'); set(han1,'linewidth',2);
        axis([t0 t1(min(n+tmore,length(t1))) vmin vmax]);
        
        ylabel({'longshore','position','(no assim.) [m]'});
        set(gca,'FontSize',14,'XTickLabel',[]);
        ax = gca; ax.YRuler.Exponent = 0;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Figure 2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%a
        hanf2=figure(f2); set(hanf2,'PaperPositionMode','Auto','Position',[985 100 920 1000],'color','w','Visible','off');
        
        ALPHA=0.25;
        
        tmore=60;
        
        xp0=0.165;
        yp0=0.865;
        width=0.805;
        height=0.102;
        ysep=0.0125;
        
        nplot=1; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;
        han1=plot(t,Hs(id_tr,:,1),'-b',t(idw),Hs(id_tr,idw,1),'ro'); set(han1,'MarkerFaceColor','r','MarkerSize',10);  set(gca,'XTick',linspace(t(1),t(min(n+tmore,length(t))),12),'Xlim',[t(1) t(min(n+tmore,length(t)))]); datetick('x','keeplimits','keepticks');

        vmin=0;
        vmax=1.1*max(Hs(id_tr,:,1));
        han1=plot([tforecast tforecast],[vmin vmax],'--k'); set(han1,'linewidth',2);
        han1=plot([tforecast2 tforecast2],[vmin vmax],'--k'); set(han1,'linewidth',2);
        
        axis([t0 tforecast2+2*365 vmin vmax]);
        
        ylabel({'wave','height [m]'});
        set(gca,'FontSize',14,'XTickLabel',[]);
        ax = gca; ax.YRuler.Exponent = 0;
        title(sprintf('CoSMoS-COAST: model parameters (transect ID: %d)',ID(id_tr)),'FontSize',14);
        
        nplot=2; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;
        
        nn=~isnan(Y_obs(id_tr,:)) & SAT(id_tr,:)==1;
        han1=errorbar(t_obs(nn),Y_obs(id_tr,nn)-Y00(id_tr),2*Y_rms(id_tr,nn),'b'); set(han1,'linestyle','none','color',[0.4 0.6 1]);
        han1=plot(t_obs(nn),Y_obs(id_tr,nn)-Y00(id_tr),'bo'); set(han1,'MarkerSize',4,'MarkerFaceColor',[0.4 0.6 1]); datetick('x','keeplimits','keepticks')
        
        nn=~isnan(Y_obs(id_tr,:)) & SAT(id_tr,:)==0;
        han1=errorbar(t_obs(nn),Y_obs(id_tr,nn)-Y00(id_tr),2*Y_rms(id_tr,nn),'m'); set(han1,'linestyle','none','color',[0.278 0 0.7]);
        han1=plot(t_obs(nn),Y_obs(id_tr,nn)-Y00(id_tr),'mo'); set(han1,'MarkerSize',4,'MarkerFaceColor',[0.58 0.30 1],'MarkerEdgeColor',[0.278 0 0.7]); datetick('x','keeplimits','keepticks')
        
        nn=~isnan(Y_obs(id_tr,:));
        
        han1=fill([t_output(~isnan(t_output)) fliplr(t_output(~isnan(t_output)))],[YCI1(id_tr,~isnan(t_output))-Y00(id_tr) fliplr(YCI2(id_tr,~isnan(t_output))-Y00(id_tr))],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
        plot(t_output,YY(id_tr,:)-Y00(id_tr),'-r'); set(gca,'XTick',linspace(t1(1),t1(min(n+tmore,length(t1))),12),'Xlim',[t1(1) t1(min(n+tmore,length(t1)))]);
        
        han1=plot([t1(1) t1(min(n+tmore,length(t1)))],[Ymin(id_tr)-Y00(id_tr) Ymin(id_tr)-Y00(id_tr)],'k--'); set(han1,'LineWidth',2);
        
        vmin=min(YCI1(id_tr,~isnan(t_output))-Y00(id_tr)); if ~isempty(Y_obs(id_tr,nn)); vmin=min(vmin,min(Y_obs(id_tr,nn)-Y00(id_tr)-2*Y_rms(id_tr,nn))); end; if vmin<0; vmin=1.05*vmin; else; vmin=0.95*vmin; end
        vmax=max(YCI2(id_tr,~isnan(t_output))-Y00(id_tr)); if ~isempty(Y_obs(id_tr,nn)); vmax=max(vmax,max(Y_obs(id_tr,nn)-Y00(id_tr)+2*Y_rms(id_tr,nn))); end; if vmax<0; vmax=0.95*vmax; else; vmax=1.05*vmax; end
        han1=plot([tforecast tforecast],[vmin vmax],'--k'); set(han1,'linewidth',2);
        han1=plot([tforecast2 tforecast2],[vmin vmax],'--k'); set(han1,'linewidth',2);
        
        axis([t0 tforecast2+2*365 vmin vmax]);
        
        ylabel({'total','shoreline','position [m]'});
        set(gca,'FontSize',14,'XTickLabel',[]);
        ax = gca; ax.YRuler.Exponent = 0;
        
        nplot=3; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;
        han1=fill([t_output(~isnan(t_output)) fliplr(t_output(~isnan(t_output)))],[DTCI1(id_tr,~isnan(t_output)) fliplr(DTCI2(id_tr,~isnan(t_output)))],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
        plot(t_output,DTOUT(id_tr,:),'-r'); set(gca,'XTick',linspace(t1(1),t1(min(n+tmore,length(t1))),12),'Xlim',[t1(1) t1(min(n+tmore,length(t1)))]); datetick('x','keeplimits','keepticks');
        
        vmin=min(DTCI1(id_tr,~isnan(t_output))); if vmin<0; vmin=1.2*vmin; else; vmin=0.9*vmin; end
        vmax=max(DTCI2(id_tr,~isnan(t_output))); if vmax<0; vmax=0.9*vmax; else; vmax=1.2*vmax; end
        
        if isnan(vmin); vmin=lb_DT; end
        if isnan(vmax); vmax=ub_DT; end
        
        han1=plot([tforecast tforecast],[vmin vmax],'--k'); set(han1,'linewidth',2);
        han1=plot([tforecast2 tforecast2],[vmin vmax],'--k'); set(han1,'linewidth',2);
        
        axis([t0 tforecast2+2*365 vmin vmax]);
        
        ylabel({'time-scale','parameter','DT [s]'});
        set(gca,'FontSize',14,'XTickLabel',[]);
        ax = gca; ax.YRuler.Exponent = 0;
        
        nplot=4; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;
        
        han1=fill([t_output(~isnan(t_output)) fliplr(t_output(~isnan(t_output)))],[DYCI1(id_tr,~isnan(t_output)) fliplr(DYCI2(id_tr,~isnan(t_output)))],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
        plot(t_output,DYOUT(id_tr,:),'-r'); set(gca,'XTick',linspace(t1(1),t1(min(n+tmore,length(t1))),12),'Xlim',[t1(1) t1(min(n+tmore,length(t1)))]); datetick('x','keeplimits','keepticks');
        
        vmin=min(DYCI1(id_tr,~isnan(t_output))); if vmin<0; vmin=1.2*vmin; else; vmin=0.9*vmin; end
        vmax=max(DYCI2(id_tr,~isnan(t_output))); if vmax<0; vmax=0.9*vmax; else; vmax=1.2*vmax; end
        
        if isnan(vmin); vmin=lb_DY; end
        if isnan(vmax); vmax=ub_DY; end
        
        han1=plot([tforecast tforecast],[vmin vmax],'--k'); set(han1,'linewidth',2);
        han1=plot([tforecast2 tforecast2],[vmin vmax],'--k'); set(han1,'linewidth',2);
        
        axis([t0 tforecast2+2*365 vmin vmax]);
        
        ylabel({'excursion','parameter','DY [m]'});
        set(gca,'FontSize',14,'XTickLabel',[]);
        ax = gca; ax.YRuler.Exponent = 0;
        
        nplot=5; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;
        
        han1=plot([t1(1) t1(min(n+tmore,length(t1)))],[HSB0(id_tr) HSB0(id_tr)],'--b'); set(han1,'LineWidth',2);
        
        han1=fill([t_output(~isnan(t_output)) fliplr(t_output(~isnan(t_output)))],[HSBCI1(id_tr,~isnan(t_output)) fliplr(HSBCI2(id_tr,~isnan(t_output)))],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
        plot(t_output,HSBOUT(id_tr,:),'-r'); set(gca,'XTick',linspace(t1(1),t1(min(n+tmore,length(t1))),12),'Xlim',[t1(1) t1(min(n+tmore,length(t1)))]); datetick('x','keeplimits','keepticks');
        
        vmin=min(HSBCI1(id_tr,~isnan(t_output))); vmin=min(vmin,HSB0(id_tr)); if vmin<0; vmin=1.2*vmin; else; vmin=0.9*vmin; end
        vmax=max(HSBCI2(id_tr,~isnan(t_output))); vmax=max(vmax,HSB0(id_tr)); if vmax<0; vmax=0.9*vmax; else; vmax=1.2*vmax; end
        
        if isnan(vmin); vmin=lb_HSB; end
        if isnan(vmax); vmax=ub_HSB; end
        
        han1=plot([tforecast tforecast],[vmin vmax],'--k'); set(han1,'linewidth',2);
        han1=plot([tforecast2 tforecast2],[vmin vmax],'--k'); set(han1,'linewidth',2);
        
        axis([t0 tforecast2+2*365 vmin vmax]);
        
        ylabel({'wave height','parameter','(H_s)_b [m]'});
        set(gca,'FontSize',14,'XTickLabel',[]);
        ax = gca; ax.YRuler.Exponent = 0;
        
        nplot=6; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;
        han1=fill([t_output(~isnan(t_output)) fliplr(t_output(~isnan(t_output)))],[KCI1(id_tr,~isnan(t_output)) fliplr(KCI2(id_tr,~isnan(t_output)))],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
        plot(t_output,KOUT(id_tr,:),'-r'); set(gca,'XTick',linspace(t1(1),t1(min(n+tmore,length(t1))),12),'Xlim',[t1(1) t1(min(n+tmore,length(t1)))]); datetick('x','keeplimits','keepticks');
        
        vmin=min(KCI1(id_tr,~isnan(t_output))); if vmin<0; vmin=1.2*vmin; else; vmin=0.9*vmin; end
        vmax=max(KCI2(id_tr,~isnan(t_output))); if vmax<0; vmax=0.9*vmax; else; vmax=1.2*vmax; end
        
        if isnan(vmin); vmin=lb_K; end
        if isnan(vmax); vmax=ub_K; end
        
        han1=plot([tforecast tforecast],[vmin vmax],'--k'); set(han1,'linewidth',2);
        han1=plot([tforecast2 tforecast2],[vmin vmax],'--k'); set(han1,'linewidth',2);
        
        axis([t0 tforecast2+2*365 vmin vmax]);
        
        ylabel({'longshore','transport','parameter','K [-]'});
        set(gca,'FontSize',14,'XTickLabel',[]);
        ax = gca; ax.YRuler.Exponent = 0;
        
        nplot=7; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;
        
        han1=plot([t1(1) t1(min(n+tmore,length(t1)))],[LTER(id_tr) LTER(id_tr)],'--b'); set(han1,'LineWidth',2);
        
        han1=fill([t_output(~isnan(t_output)) fliplr(t_output(~isnan(t_output)))],[VLTCI1(id_tr,~isnan(t_output)) fliplr(VLTCI2(id_tr,~isnan(t_output)))],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
        plot(t_output,VLTOUT(id_tr,:),'-r'); set(gca,'XTick',linspace(t1(1),t1(min(n+tmore,length(t1))),12),'Xlim',[t1(1) t1(min(n+tmore,length(t1)))]); datetick('x','keeplimits','keepticks');
        
        vmin=min(VLTCI1(id_tr,~isnan(t_output))); vmin=min(vmin,LTER(id_tr)); if vmin<0; vmin=1.2*vmin; else; vmin=0.9*vmin; end
        vmax=max(VLTCI2(id_tr,~isnan(t_output))); vmax=max(vmax,LTER(id_tr)); if vmax<0; vmax=0.9*vmax; else; vmax=1.2*vmax; end
        
        if isnan(vmin); vmin=lb_vlt; end
        if isnan(vmax); vmax=ub_vlt; end
        
        han1=plot([tforecast tforecast],[vmin vmax],'--k'); set(han1,'linewidth',2);
        han1=plot([tforecast2 tforecast2],[vmin vmax],'--k'); set(han1,'linewidth',2);
        
        axis([t0 tforecast2+2*365 vmin vmax]);
        
        ylabel({'long-term rate','parameter','v_{lt} [m/yr]'});
        set(gca,'FontSize',14,'XTickLabel',[]);
        ax = gca; ax.YRuler.Exponent = 0;
        
        nplot=8; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;
        
        han1=plot([t1(1) t1(min(n+tmore,length(t1)))],[sigma00 sigma00],'--b'); set(han1,'LineWidth',2);
        
        han1=fill([t_output(~isnan(t_output)) fliplr(t_output(~isnan(t_output)))],[SIGMACI1(id_tr,~isnan(t_output)) fliplr(SIGMACI2(id_tr,~isnan(t_output)))],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
        plot(t_output,SIGMAOUT(id_tr,:),'-r'); set(gca,'XTick',linspace(t1(1),t1(min(n+tmore,length(t1))),12),'Xlim',[t1(1) t1(min(n+tmore,length(t1)))]); datetick('x','keeplimits','keepticks');
        
        vmin=min(SIGMACI1(id_tr,~isnan(t_output))); if vmin<0; vmin=1.2*vmin; else; vmin=0.9*vmin; end
        vmax=max(SIGMACI2(id_tr,~isnan(t_output))); if vmax<0; vmax=0.9*vmax; else; vmax=1.2*vmax; end
        
        if isnan(vmin); vmin=lb_sigma; end
        if isnan(vmax); vmax=ub_sigma; end
        
        han1=plot([tforecast tforecast],[vmin vmax],'--k'); set(han1,'linewidth',2);
        han1=plot([tforecast2 tforecast2],[vmin vmax],'--k'); set(han1,'linewidth',2);
        
        axis([t0 tforecast2+2*365 vmin vmax]);
        
        ylabel({'noise','parameter','\sigma [var.]'});
        set(gca,'FontSize',14);
        ax = gca; ax.YRuler.Exponent = 0;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Figure 3
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%a
        hanf3=figure(f3); set(hanf3,'PaperPositionMode','Auto','Position',[50 100 920 1000],'color','w','Visible','off');
        
        ALPHA=0.25;
        
        tmore=60;
        
        xp0=0.16;
        yp0=0.84;
        width=0.81;
        height=0.120;
        ysep=0.01;
        
        NTR=ID;
        
        nn=find(t_obs<=t1(min(n+1,length(t1))),1,'last');
        nnn=find(t_output<=t1(min(n+1,length(t1))),1,'last');
        
        if PLOT_SECTIONS==1
            if strcmp(transects(id_tr).model_type,'full model')
                id=find(id_full_model_end>=id_tr,1,'first');
                id1=id_full_model_start(id);
                id2=id_full_model_end(id);
            elseif strcmp(transects(id_tr).model_type,'cross-shore only')
                id=find(id_cross_shore_only_end>=id_tr,1,'first');
                id1=max(id_cross_shore_only_start(id)-1,1);
                id2=min(id_cross_shore_only_end(id)+1,Ntr);
            elseif  strcmp(transects(id_tr).model_type,'rate only')
                id=find(id_rate_only_end>=id_tr,1,'first');
                id1=max(id_rate_only_start(id)-1,1);
                id2=min(id_rate_only_end(id)+1,Ntr);
            end
            xmin=ID(id1);
            xmax=ID(id2);
        else
            xmin=min(NTR);
            xmax=max(NTR);
        end
        
        nplot=1; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;
        
        if ~isempty(nn) && DATA_ASSIMILATION
            han1=errorbar(NTR,Y_obs(:,nn)-Y00,2*Y_rms(:,nn),'b'); set(han1,'linestyle','none','color',[0.4 0.6 1]);
            han1=plot(NTR,Y_obs(:,nn)-Y00,'bo'); set(han1,'MarkerSize',4,'MarkerFaceColor',[0.4 0.6 1]);
        end
        
        VAR=YY(:,nnn)-Y00;
        VARCI1=YCI1(:,nnn)-Y00;
        VARCI2=YCI2(:,nnn)-Y00;
        VAROBS=Y_obs(:,nn)-Y00;
        VARMIN=Ymin-Y00;
        
        VAR(bool_cliff_only | bool_no_prediction)=NaN;
        VARCI1(bool_cliff_only | bool_no_prediction)=NaN;
        VARCI2(bool_cliff_only | bool_no_prediction)=NaN;
        
        id=find(~isnan(VARCI1));
        idx=find(diff(id)~=1);
        if ~isempty(idx)
            A=[idx(1);diff(idx);numel(id)-idx(end)];
        else
            A=length(id);
        end
        C1=mat2cell(NTR(id),A,1);
        C2=mat2cell(VAR(id),A,1);
        C3=mat2cell(VARCI1(id),A,1);
        C4=mat2cell(VARCI2(id),A,1);
        
        for ii=1:length(C1)
            han1=fill([C1{ii}; flipud(C1{ii})],[C3{ii}; flipud(C4{ii})],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
            plot(C1{ii},C2{ii},'-r');
        end
        
        vmin=min(VARCI1(id1:id2)); if vmin<0; vmin=1.2*vmin; else; vmin=0.9*vmin; end
        vmax=max(VARCI2(id1:id2)); if vmax<0; vmax=0.9*vmax; else; vmax=1.2*vmax; end
        
        vmin=min(vmin,min(VAROBS(id1:id2))-25);
        vmax=max(vmax,max(VAROBS(id1:id2))+25);
        
        han1=plot([ID(id_tr) ID(id_tr)],[vmin vmax],'--r'); set(han1,'linewidth',2);
        
        han1=plot(NTR,VARMIN,'k--'); set(han1,'LineWidth',2);
        
        title(sprintf('CoSMoS-COAST: alongshore model parameters (transect ID: %d)',ID(id_tr)),'FontSize',14);
        
        ylabel({'total','shoreline','position [m]'});
        set(gca,'FontSize',14,'XTickLabel',[],'Xlim',[xmin xmax],'Ylim',[vmin vmax]);
        ax = gca; ax.YRuler.Exponent = 0;
        
        nplot=2; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;
        
        VAR=DTOUT(:,nnn);
        VARCI1=DTCI1(:,nnn);
        VARCI2=DTCI2(:,nnn);
        
        VAR(bool_cliff_only | bool_no_prediction | bool_rate_only)=NaN;
        VARCI1(bool_cliff_only | bool_no_prediction | bool_rate_only)=NaN;
        VARCI2(bool_cliff_only | bool_no_prediction | bool_rate_only)=NaN;
        
        id=find(~isnan(VARCI1));
        idx=find(diff(id)~=1);
        if ~isempty(idx)
            A=[idx(1);diff(idx);numel(id)-idx(end)];
        else
            A=length(id);
        end
        C1=mat2cell(NTR(id),A,1);
        C2=mat2cell(VAR(id),A,1);
        C3=mat2cell(VARCI1(id),A,1);
        C4=mat2cell(VARCI2(id),A,1);
        
        for ii=1:length(C1)
            han1=fill([C1{ii}; flipud(C1{ii})],[C3{ii}; flipud(C4{ii})],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
            plot(C1{ii},C2{ii},'-r');
        end
        
        vmin=min(VARCI1(id1:id2)); if vmin<0; vmin=1.2*vmin; else; vmin=0.9*vmin; end
        vmax=max(VARCI2(id1:id2)); if vmax<0; vmax=0.9*vmax; else; vmax=1.2*vmax; end
        
        if isnan(vmin); vmin=0; end
        if isnan(vmax); vmax=ub_DT; end

        han1=plot([ID(id_tr) ID(id_tr)],[vmin vmax],'--r'); set(han1,'linewidth',2);
        
        ylabel({'time-scale','parameter','DT [s]'});
        set(gca,'FontSize',14,'XTickLabel',[],'Xlim',[xmin xmax],'Ylim',[vmin vmax]);
        ax = gca; ax.YRuler.Exponent = 0;
        
        nplot=3; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;
        
        VAR=DYOUT(:,nnn);
        VARCI1=DYCI1(:,nnn);
        VARCI2=DYCI2(:,nnn);
        
        VAR(bool_cliff_only | bool_no_prediction | bool_rate_only)=NaN;
        VARCI1(bool_cliff_only | bool_no_prediction | bool_rate_only)=NaN;
        VARCI2(bool_cliff_only | bool_no_prediction | bool_rate_only)=NaN;
        
        id=find(~isnan(VARCI1));
        idx=find(diff(id)~=1);
        if ~isempty(idx)
            A=[idx(1);diff(idx);numel(id)-idx(end)];
        else
            A=length(id);
        end
        C1=mat2cell(NTR(id),A,1);
        C2=mat2cell(VAR(id),A,1);
        C3=mat2cell(VARCI1(id),A,1);
        C4=mat2cell(VARCI2(id),A,1);
        
        for ii=1:length(C1)
            han1=fill([C1{ii}; flipud(C1{ii})],[C3{ii}; flipud(C4{ii})],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
            plot(C1{ii},C2{ii},'-r');
        end
        
        vmin=min(VARCI1(id1:id2)); if vmin<0; vmin=1.2*vmin; else; vmin=0.9*vmin; end
        vmax=max(VARCI2(id1:id2)); if vmax<0; vmax=0.9*vmax; else; vmax=1.2*vmax; end
        
        if isnan(vmin); vmin=0; end
        if isnan(vmax); vmax=ub_DY; end
        
        han1=plot([ID(id_tr) ID(id_tr)],[vmin vmax],'--r'); set(han1,'linewidth',2);
        
        ylabel({'excursion','parameter','DY [m]'});
        set(gca,'FontSize',14,'XTickLabel',[],'Xlim',[xmin xmax],'Ylim',[vmin vmax]);
        ax = gca; ax.YRuler.Exponent = 0;
        
        nplot=4; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;
        
        VAR=HSBOUT(:,nnn);
        VARCI1=HSBCI1(:,nnn);
        VARCI2=HSBCI2(:,nnn);
        
        VAR(bool_cliff_only | bool_no_prediction | bool_rate_only)=NaN;
        VARCI1(bool_cliff_only | bool_no_prediction | bool_rate_only)=NaN;
        VARCI2(bool_cliff_only | bool_no_prediction | bool_rate_only)=NaN;
        
        id=find(~isnan(VARCI1));
        idx=find(diff(id)~=1);
        if ~isempty(idx)
            A=[idx(1);diff(idx);numel(id)-idx(end)];
        else
            A=length(id);
        end
        C1=mat2cell(NTR(id),A,1);
        C2=mat2cell(VAR(id),A,1);
        C3=mat2cell(VARCI1(id),A,1);
        C4=mat2cell(VARCI2(id),A,1);
        
        for ii=1:length(C1)
            han1=fill([C1{ii}; flipud(C1{ii})],[C3{ii}; flipud(C4{ii})],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
            plot(C1{ii},C2{ii},'-r');
        end
        
        vmin=min(VARCI1(id1:id2)); if vmin<0; vmin=1.2*vmin; else; vmin=0.9*vmin; end
        vmax=max(VARCI2(id1:id2)); if vmax<0; vmax=0.9*vmax; else; vmax=1.2*vmax; end
        
        if isnan(vmin); vmin=lb_HSB; end
        if isnan(vmax); vmax=ub_HSB; end
        
        han1=plot(NTR,HSB0,'-b'); set(han1,'LineWidth',2);
        
        han1=plot([ID(id_tr) ID(id_tr)],[vmin vmax],'--r'); set(han1,'linewidth',2);
        
        ylabel({'wave height','parameter','(H_s)_b [m]'});
        set(gca,'FontSize',14,'XTickLabel',[],'Xlim',[xmin xmax],'Ylim',[vmin vmax]);
        ax = gca; ax.YRuler.Exponent = 0;
        
        nplot=5; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;
        
        VAR=KOUT(:,nnn);
        VARCI1=KCI1(:,nnn);
        VARCI2=KCI2(:,nnn);
        
        VAR(bool_cliff_only | bool_no_prediction | bool_rate_only | bool_cross_shore_only)=NaN;
        VARCI1(bool_cliff_only | bool_no_prediction | bool_rate_only | bool_cross_shore_only)=NaN;
        VARCI2(bool_cliff_only | bool_no_prediction | bool_rate_only | bool_cross_shore_only)=NaN;
        
        id=find(~isnan(VARCI1));
        idx=find(diff(id)~=1);
        if ~isempty(idx)
            A=[idx(1);diff(idx);numel(id)-idx(end)];
        else
            A=length(id);
        end
        C1=mat2cell(NTR(id),A,1);
        C2=mat2cell(VAR(id),A,1);
        C3=mat2cell(VARCI1(id),A,1);
        C4=mat2cell(VARCI2(id),A,1);
        
        for ii=1:length(C1)
            han1=fill([C1{ii}; flipud(C1{ii})],[C3{ii}; flipud(C4{ii})],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
            plot(C1{ii},C2{ii},'-r');
        end
        
        vmin=min(VARCI1(id1:id2)); if vmin<0; vmin=1.2*vmin; else; vmin=0.9*vmin; end
        vmax=max(VARCI2(id1:id2)); if vmax<0; vmax=0.9*vmax; else; vmax=1.2*vmax; end
        
        if isnan(vmin); vmin=0; end
        if isnan(vmax); vmax=ub_K; end
        
        han1=plot([ID(id_tr) ID(id_tr)],[vmin vmax],'--r'); set(han1,'linewidth',2);
        
        ylabel({'longshore','transport','parameter','K [-]'});
        set(gca,'FontSize',14,'XTickLabel',[],'Xlim',[xmin xmax],'Ylim',[vmin vmax]);
        ax = gca; ax.YRuler.Exponent = 0;
        
        nplot=6; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;
        
        VAR=VLTOUT(:,nnn);
        VARCI1=VLTCI1(:,nnn);
        VARCI2=VLTCI2(:,nnn);
        
        VAR(bool_cliff_only | bool_no_prediction)=NaN;
        VARCI1(bool_cliff_only | bool_no_prediction)=NaN;
        VARCI2(bool_cliff_only | bool_no_prediction)=NaN;
        
        id=find(~isnan(VARCI1));
        idx=find(diff(id)~=1);
        if ~isempty(idx)
            A=[idx(1);diff(idx);numel(id)-idx(end)];
        else
            A=length(id);
        end
        C1=mat2cell(NTR(id),A,1);
        C2=mat2cell(VAR(id),A,1);
        C3=mat2cell(VARCI1(id),A,1);
        C4=mat2cell(VARCI2(id),A,1);
        
        for ii=1:length(C1)
            han1=fill([C1{ii}; flipud(C1{ii})],[C3{ii}; flipud(C4{ii})],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
            plot(C1{ii},C2{ii},'-r');
        end
        
        vmin=min(VARCI1(id1:id2)); if vmin<0; vmin=1.2*vmin; else; vmin=0.9*vmin; end
        vmax=max(VARCI2(id1:id2)); if vmax<0; vmax=0.9*vmax; else; vmax=1.2*vmax; end
        
        vmin=min(vmin,min(LTER(id1:id2))-0.1);
        vmax=max(vmax,max(LTER(id1:id2))+0.1);
        
        han1=plot(NTR,LTER,'-b'); set(han1,'LineWidth',2);
        
        han1=plot([ID(id_tr) ID(id_tr)],[vmin vmax],'--r'); set(han1,'linewidth',2);
        
        ylabel({'long-term rate','parameter','v_{lt} [m/yr]'});
        set(gca,'FontSize',14,'XTickLabel',[],'Xlim',[xmin xmax],'Ylim',[vmin vmax]);
        ax = gca; ax.YRuler.Exponent = 0;
        
        nplot=7; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;
        
        VAR=SIGMAOUT(:,nnn);
        VARCI1=SIGMACI1(:,nnn);
        VARCI2=SIGMACI2(:,nnn);
        
        VAR(bool_cliff_only | bool_no_prediction)=NaN;
        VARCI1(bool_cliff_only | bool_no_prediction)=NaN;
        VARCI2(bool_cliff_only | bool_no_prediction)=NaN;
        
        id=find(~isnan(VARCI1));
        idx=find(diff(id)~=1);
        if ~isempty(idx)
            A=[idx(1);diff(idx);numel(id)-idx(end)];
        else
            A=length(id);
        end
        C1=mat2cell(NTR(id),A,1);
        C2=mat2cell(VAR(id),A,1);
        C3=mat2cell(VARCI1(id),A,1);
        C4=mat2cell(VARCI2(id),A,1);
        
        han1=plot([NTR(1) NTR(end)],[sigma00 sigma00],'--b'); set(han1,'LineWidth',2);
        
        for ii=1:length(C1)
            han1=fill([C1{ii}; flipud(C1{ii})],[C3{ii}; flipud(C4{ii})],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
            plot(C1{ii},C2{ii},'-r');
        end
        
        vmin=min(VARCI1(id1:id2)); if vmin<0; vmin=1.2*vmin; else; vmin=0.9*vmin; end
        vmax=max(VARCI2(id1:id2)); if vmax<0; vmax=0.9*vmax; else; vmax=1.2*vmax; end
        
        if isnan(vmin); vmin=0; end
        if isnan(vmax); vmax=ub_sigma; end
        
        han1=plot([ID(id_tr) ID(id_tr)],[vmin vmax],'--r'); set(han1,'linewidth',2);
        
        xlabel('Transect #');
        ylabel({'noise','parameter','\sigma [var.]'});
        set(gca,'FontSize',14,'Xlim',[xmin xmax],'Ylim',[vmin vmax]);
        ax = gca; ax.YRuler.Exponent = 0;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Figure 4
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%a
        hanf4=figure(f4); set(hanf4,'PaperPositionMode','Auto','Position',[100 100 920 700],'color','w','Visible','off');
        
        ALPHA=0.2;
        
        xp0=0.08;
        yp0=0.08;
        
        width=0.9;
        height=0.38;
        
        ysep=0.12;
        
        subplot('Position',[xp0 yp0+height+ysep width height]); hold on; box on;
        
        han1=plot(t,Hs(id_tr,:,1),'-b');

        vmin=0;
        vmax=1.1*max(Hs(id_tr,:,1));
        axis([t0 t1(min(n+tmore,length(t1))) vmin vmax]);
        
        set(gca,'Xlim',[t0 tforecast2+1*365.25],'Ylim',[vmin vmax]);
        set(gca,'xtick',datenum(1995:2:2021,1,1)); datetick('x','keeplimits','keepticks');
        set(gca,'Xlim',[t0 tforecast2+1*365.25],'Ylim',[vmin vmax]);
        set(gca,'FontSize',14);
        
        title(sprintf('CoSMoS-COAST: model calibration/validation (transect ID: %d)',ID(id_tr)),'FontSize',14);
        
        ylabel({'Wave height [m]'});
        
        subplot('Position',[xp0 yp0 width height]); hold on; box on;
        
        nn=~isnan(Y_obs(id_tr,:)) & SAT(id_tr,:)==1;
        han1=errorbar(t_obs(nn),Y_obs(id_tr,nn)-Y00(id_tr),2*Y_rms(id_tr,nn),'b'); set(han1,'linestyle','none','color',[0.4 0.6 1]);
        han1=plot(t_obs(nn),Y_obs(id_tr,nn)-Y00(id_tr),'bo'); set(han1,'MarkerSize',4,'MarkerFaceColor',[0.4 0.6 1]); datetick('x','keeplimits','keepticks')
        
        nn=~isnan(Y_obs(id_tr,:)) & SAT(id_tr,:)==0;
        han1=errorbar(t_obs(nn),Y_obs(id_tr,nn)-Y00(id_tr),2*Y_rms(id_tr,nn),'m'); set(han1,'linestyle','none','color',[0.278 0 0.7]);
        han1=plot(t_obs(nn),Y_obs(id_tr,nn)-Y00(id_tr),'mo'); set(han1,'MarkerSize',4,'MarkerFaceColor',[0.58 0.30 1],'MarkerEdgeColor',[0.278 0 0.7]); datetick('x','keeplimits','keepticks')
        
        nn=~isnan(Y_obs(id_tr,:));
        
        han1=fill([t_output(~isnan(t_output)) fliplr(t_output(~isnan(t_output)))],[YCI1(id_tr,~isnan(t_output))-Y00(id_tr) fliplr(YCI2(id_tr,~isnan(t_output))-Y00(id_tr))],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
        plot(t_output,YY(id_tr,:)-Y00(id_tr),'-r'); set(gca,'XTick',linspace(t1(1),t1(min(n+tmore,length(t1))),12),'Xlim',[t1(1) t1(min(n+tmore,length(t1)))]);
        
        han1=plot([t0 tforecast2+1*365.25],[Ymin(id_tr)-Y00(id_tr) Ymin(id_tr)-Y00(id_tr)],'k--'); set(han1,'LineWidth',2);
        
        vmin=min(YCI1(id_tr,~isnan(t_output))-Y00(id_tr)); if ~isempty(Y_obs(id_tr,nn)); vmin=min(vmin,min(Y_obs(id_tr,nn)-Y00(id_tr)-2*Y_rms(id_tr,nn))); end; if vmin<0; vmin=1.05*vmin; else; vmin=0.95*vmin; end
        vmax=max(YCI2(id_tr,~isnan(t_output))-Y00(id_tr)); if ~isempty(Y_obs(id_tr,nn)); vmax=max(vmax,max(Y_obs(id_tr,nn)-Y00(id_tr)+2*Y_rms(id_tr,nn))); end; if vmax<0; vmax=0.95*vmax; else; vmax=1.05*vmax; end
        
        % han1=plot([tforecast tforecast],[vmin vmax],'--k'); set(han1,'linewidth',2);
        % han1=plot([tforecast2 tforecast2],[vmin vmax],'--k'); set(han1,'linewidth',2);
        axis([t0 t1(min(n+tmore,length(t1))) vmin vmax]);
        
        set(gca,'Xlim',[t0 tforecast2+1*365.25],'Ylim',[vmin vmax]);
        set(gca,'xtick',datenum(1999:2:2021,1,1)); datetick('x','keeplimits','keepticks'); axis tight;
        set(gca,'Xlim',[t0 tforecast2+1*365.25],'Ylim',[vmin vmax]);
        set(gca,'FontSize',14);
        xlabel('time'); ylabel('Shoreline position [m]');
        
        han1=legend('observed C.I.','observed','95% CI','ensemble mean');
        set(han1,'Position',[0.0909193847350033 0.395361616706461 0.185869561330132 0.127857139280864]);
        
        annotation(gcf,'line',[0.746444746376811 0.746444746376811],[0.96047619047619 0.069165017807677],'LineStyle','--','LineWidth',2);
        annotation(gcf,'line',[0.946512681159419 0.946512681159419],[0.96047619047619 0.069165017807677],'LineStyle','--','LineWidth',2);
        
        annotation(gcf,'doublearrow',[0.724660326086956 0.767368659420289],[0.506848568790397 0.506925207756233]);
        annotation(gcf,'doublearrow',[0.925203804347827 0.96791213768116],[0.506935628545047 0.507012267510882]);
        
        han1=annotation(gcf,'textbox','Position',[0.625448369565217 0.473173196148265 0.07 0.0700000000000002],'String',{'Hindcast','(Calibration)'},'FitBoxToText','off','LineStyle','none','FontSize',14,'HorizontalAlignment','center');
        han1=annotation(gcf,'textbox','Position',[0.820298913043479 0.474601767576835 0.07 0.0700000000000001],'String',{'Hindcast','(Validation)'},'FitBoxToText','off','LineStyle','none','FontSize',14,'HorizontalAlignment','center');
        %han1=annotation(gcf,'textbox','Position',[0.926354166666666 0.461428923185154 0.07 0.0700000000000001],'String',{},'FitBoxToText','off','LineStyle','none','FontSize',12,'HorizontalAlignment','center','Rotation',90);
        annotation('textarrow',[0.891485507246376 0.858786231884057],[0.553111739745404 0.610749646393211],'String','Forecast','HeadStyle','none','LineStyle','none', 'TextRotation',90,'FontSize',12);

        text(0.5*(tforecast+tforecast2)-2.5*365,vmin+0.2*(vmax-vmin),['RMSE (sat) = ',num2str(RMSE_sat(id_tr),3),' m'],'FontSize',14);
        text(0.5*(tforecast+tforecast2)-2.5*365,vmin+0.075*(vmax-vmin),['within CI = ',num2str(PCT_sat(id_tr),3),'%'],'FontSize',14);
        
        % save to files
        %print(figure(f1),[OUTPUT_DIR,filesep,'figures\CoSMoS-COAST_results_components_tr_ID_',num2str(ID(id_tr),'%0.5d'),'.png'],'-dpng','-r200');
        %print(figure(f2),[OUTPUT_DIR,filesep,'figures\CoSMoS-COAST_results_parameters_tr_ID_',num2str(ID(id_tr),'%0.5d'),'.png'],'-dpng','-r200');
        %print(figure(f3),[OUTPUT_DIR,filesep,'figures\CoSMoS-COAST_results_alongshore_tr_ID_',num2str(ID(id_tr),'%0.5d'),'.png'],'-dpng','-r200');
        %print(figure(f4),[OUTPUT_DIR,filesep,'figures\CoSMoS-COAST_calibration_validation_tr_ID_',num2str(ID(id_tr),'%0.5d'),'.png'],'-dpng','-r200');
        
        saveas(figure(f1),[OUTPUT_DIR,filesep,'figures\CoSMoS-COAST_results_components_tr_ID_',num2str(ID(id_tr),'%0.5d'),'.png'],'png');
        saveas(figure(f2),[OUTPUT_DIR,filesep,'figures\CoSMoS-COAST_results_parameters_tr_ID_',num2str(ID(id_tr),'%0.5d'),'.png'],'png');
        saveas(figure(f3),[OUTPUT_DIR,filesep,'figures\CoSMoS-COAST_results_alongshore_tr_ID_',num2str(ID(id_tr),'%0.5d'),'.png'],'png');
        saveas(figure(f4),[OUTPUT_DIR,filesep,'figures\CoSMoS-COAST_calibration_validation_tr_ID_',num2str(ID(id_tr),'%0.5d'),'.png'],'png');
        
        close all;
        
    end

end

delete(gcp('nocreate'))