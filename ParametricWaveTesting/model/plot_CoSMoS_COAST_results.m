%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT the simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); clf;

ALPHA=0.25;

tmore=60;

xp0=0.15;
yp0=0.845;
width=0.82;
height=0.12;
ysep=0.014;

[idw]=ismember(t1(n),t);

nplot=1; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;
han1=plot(t,Hs(WAVEID(id_tr),:,1),'-b',t(idw),Hs(WAVEID(id_tr),idw,1),'ro'); set(han1,'MarkerFaceColor','r','MarkerSize',10);  set(gca,'XTick',linspace(t(1),t(min(n+tmore,length(t))),12),'Xlim',[t(1) t(min(n+tmore,length(t)))]); datetick('x','keeplimits','keepticks');

vmin=0;
vmax=1.3*max(Hs(id_tr,1:min(n+tmore,length(t))));
han1=plot([tforecast tforecast],[vmin vmax],'--k'); set(han1,'linewidth',2);
han1=plot([tforecast2 tforecast2],[vmin vmax],'--k'); set(han1,'linewidth',2);
axis([t0 t1(min(n+tmore,length(t1))) vmin vmax]);

ylabel({'wave','height [m]'});
set(gca,'FontSize',14,'XTickLabel',[]);
ax = gca; ax.YRuler.Exponent = 0;
title(sprintf('CoSMoS-COAST %s : model components (transect ID: %d)',Model_name,ID(id_tr)),'FontSize',14);

nplot=2; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;

nn=~isnan(Y_obs(id_tr,:)) & SAT(id_tr,:)==1;
han1=errorbar(t_obs(nn),Y_obs(id_tr,nn)-Y00(id_tr),2*Y_rms(id_tr,nn),'b'); set(han1,'linestyle','none','color',[0.4 0.6 1]);
han1=plot(t_obs(nn),Y_obs(id_tr,nn)-Y00(id_tr),'bo'); set(han1,'MarkerSize',4,'MarkerFaceColor',[0.4 0.6 1]); datetick('x','keeplimits','keepticks')

nn=~isnan(Y_obs(id_tr,:)) & SAT(id_tr,:)==0;
han1=errorbar(t_obs(nn),Y_obs(id_tr,nn)-Y00(id_tr),2*Y_rms(id_tr,nn),'m'); set(han1,'linestyle','none','color',[0.278 0 0.7]);
han1=plot(t_obs(nn),Y_obs(id_tr,nn)-Y00(id_tr),'mo'); set(han1,'MarkerSize',4,'MarkerFaceColor',[0.58 0.30 1],'MarkerEdgeColor',[0.278 0 0.7]); datetick('x','keeplimits','keepticks')

nn=~isnan(Y_obs(id_tr,:));

han1=fill([t_output(~isnan(t_output)) fliplr(t_output(~isnan(t_output)))],[YCI1(id_tr,~isnan(t_output))-Y00(id_tr) fliplr(YCI2(id_tr,~isnan(t_output))-Y00(id_tr))],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
plot(t_output,YY(id_tr,:)-Y00(id_tr),'-r'); set(gca,'XTick',linspace(t(1),t(min(n+tmore,length(t))),12),'Xlim',[t(1) t(min(n+tmore,length(t)))]);

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
plot(t_output,YST(id_tr,:),'-r'); set(gca,'XTick',linspace(t(1),t(min(n+tmore,length(t))),12),'Xlim',[t(1) t(min(n+tmore,length(t)))]); datetick('x','keeplimits','keepticks');

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
plot(t_output,YBRU(id_tr,:),'-r'); set(gca,'XTick',linspace(t(1),t(min(n+tmore,length(t))),12),'Xlim',[t(1) t(min(n+tmore,length(t)))]); datetick('x','keeplimits','keepticks');

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
plot(t_output,YVLT(id_tr,:),'-r'); set(gca,'XTick',linspace(t(1),t(min(n+tmore,length(t))),12),'Xlim',[t(1) t(min(n+tmore,length(t)))]); datetick('x','keeplimits','keepticks');

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
plot(t_output,YLST(id_tr,:),'-r'); set(gca,'XTick',linspace(t(1),t(min(n+tmore,length(t))),12),'Xlim',[t(1) t(min(n+tmore,length(t)))]); datetick('x','keeplimits','keepticks');

YLST_steps=YLST(id_tr,:)-YLST_only(id_tr,:);
plot(t_output,YLST_steps,'-m'); set(gca,'XTick',linspace(t(1),t(min(n+tmore,length(t))),12),'Xlim',[t(1) t(min(n+tmore,length(t)))]); datetick('x','keeplimits','keepticks');

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

plot(t_output,YLST_only(id_tr,:),'-r',t_output,-YLST_steps,'-m'); set(gca,'XTick',linspace(t(1),t(min(n+tmore,length(t))),12),'Xlim',[t(1) t(min(n+tmore,length(t)))]); datetick('x','keeplimits','keepticks');

vmin=min(YLST_only_CI1(id_tr,~isnan(t_output))); if vmin<0; vmin=1.2*vmin; else; vmin=0.9*vmin; end
vmax=max(YLST_only_CI2(id_tr,~isnan(t_output))); if vmax<0; vmax=0.9*vmax; else; vmax=1.2*vmax; end
han1=plot([tforecast tforecast],[vmin vmax],'--k'); set(han1,'linewidth',2);
han1=plot([tforecast2 tforecast2],[vmin vmax],'--k'); set(han1,'linewidth',2);
axis([t0 t1(min(n+tmore,length(t1))) vmin vmax]);

ylabel({'longshore','position','(no assim.) [m]'});
set(gca,'FontSize',14,'XTickLabel',[]);
ax = gca; ax.YRuler.Exponent = 0;

drawnow

if MAKE_ANIMATION % make gif animation
    
    % grab frame
    frame = getframe(gcf);
    im=frame2im(frame);
    
    [imind,cm]=rgb2ind(im,256);
    
    if count == 1
        imwrite(imind,cm,filename1,'gif', 'Loopcount',inf,'delaytime',1/fps);
    else
        imwrite(imind,cm,filename1,'gif','WriteMode','append','delaytime',1/fps);
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2); clf;

ALPHA=0.25;

tmore=60;

xp0=0.165;
yp0=0.865;
width=0.805;
height=0.102;
ysep=0.0125;

nplot=1; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;
han1=plot(t,Hs(WAVEID(id_tr),:,1),'-b',t(idw),Hs(WAVEID(id_tr),idw,1),'ro'); set(han1,'MarkerFaceColor','r','MarkerSize',10);  set(gca,'XTick',linspace(t(1),t(min(n+tmore,length(t))),12),'Xlim',[t(1) t(min(n+tmore,length(t)))]); datetick('x','keeplimits','keepticks');

vmin=0;
vmax=1.3*max(Hs(id_tr,1:min(n+tmore,length(t))));
han1=plot([tforecast tforecast],[vmin vmax],'--k'); set(han1,'linewidth',2);
han1=plot([tforecast2 tforecast2],[vmin vmax],'--k'); set(han1,'linewidth',2);

if MAKE_ANIMATION
    axis([t0 t1(min(n+tmore,length(t1))) vmin vmax]);
else
    axis([t0 tforecast2+2*365 vmin vmax]);
end

ylabel({'wave','height [m]'});
set(gca,'FontSize',14,'XTickLabel',[]);
ax = gca; ax.YRuler.Exponent = 0;
title(sprintf('CoSMoS-COAST %s: model parameters (transect ID: %d)',Model_name,ID(id_tr)),'FontSize',14);

nplot=2; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;

nn=~isnan(Y_obs(id_tr,:)) & SAT(id_tr,:)==1;
han1=errorbar(t_obs(nn),Y_obs(id_tr,nn)-Y00(id_tr),2*Y_rms(id_tr,nn),'b'); set(han1,'linestyle','none','color',[0.4 0.6 1]);
han1=plot(t_obs(nn),Y_obs(id_tr,nn)-Y00(id_tr),'bo'); set(han1,'MarkerSize',4,'MarkerFaceColor',[0.4 0.6 1]); datetick('x','keeplimits','keepticks')

nn=~isnan(Y_obs(id_tr,:)) & SAT(id_tr,:)==0;
han1=errorbar(t_obs(nn),Y_obs(id_tr,nn)-Y00(id_tr),2*Y_rms(id_tr,nn),'m'); set(han1,'linestyle','none','color',[0.278 0 0.7]);
han1=plot(t_obs(nn),Y_obs(id_tr,nn)-Y00(id_tr),'mo'); set(han1,'MarkerSize',4,'MarkerFaceColor',[0.58 0.30 1],'MarkerEdgeColor',[0.278 0 0.7]); datetick('x','keeplimits','keepticks')

nn=~isnan(Y_obs(id_tr,:));

han1=fill([t_output(~isnan(t_output)) fliplr(t_output(~isnan(t_output)))],[YCI1(id_tr,~isnan(t_output))-Y00(id_tr) fliplr(YCI2(id_tr,~isnan(t_output))-Y00(id_tr))],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
plot(t_output,YY(id_tr,:)-Y00(id_tr),'-r'); set(gca,'XTick',linspace(t(1),t(min(n+tmore,length(t))),12),'Xlim',[t(1) t(min(n+tmore,length(t)))]);

han1=plot([t(1) t(min(n+tmore,length(t)))],[Ymin(id_tr)-Y00(id_tr) Ymin(id_tr)-Y00(id_tr)],'k--'); set(han1,'LineWidth',2);

vmin=min(YCI1(id_tr,~isnan(t_output))-Y00(id_tr)); if ~isempty(Y_obs(id_tr,nn)); vmin=min(vmin,min(Y_obs(id_tr,nn)-Y00(id_tr)-2*Y_rms(id_tr,nn))); end; if vmin<0; vmin=1.05*vmin; else; vmin=0.95*vmin; end
vmax=max(YCI2(id_tr,~isnan(t_output))-Y00(id_tr)); if ~isempty(Y_obs(id_tr,nn)); vmax=max(vmax,max(Y_obs(id_tr,nn)-Y00(id_tr)+2*Y_rms(id_tr,nn))); end; if vmax<0; vmax=0.95*vmax; else; vmax=1.05*vmax; end
han1=plot([tforecast tforecast],[vmin vmax],'--k'); set(han1,'linewidth',2);
han1=plot([tforecast2 tforecast2],[vmin vmax],'--k'); set(han1,'linewidth',2);

if MAKE_ANIMATION
    axis([t0 t1(min(n+tmore,length(t1))) vmin vmax]);
else
    axis([t0 tforecast2+2*365 vmin vmax]);
end

ylabel({'total','shoreline','position [m]'});
set(gca,'FontSize',14,'XTickLabel',[]);
ax = gca; ax.YRuler.Exponent = 0;

nplot=3; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;
han1=fill([t_output(~isnan(t_output)) fliplr(t_output(~isnan(t_output)))],[DTCI1(id_tr,~isnan(t_output)) fliplr(DTCI2(id_tr,~isnan(t_output)))],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
plot(t_output,DTOUT(id_tr,:),'-r'); set(gca,'XTick',linspace(t(1),t(min(n+tmore,length(t))),12),'Xlim',[t(1) t(min(n+tmore,length(t)))]); datetick('x','keeplimits','keepticks');

vmin=min(DTCI1(id_tr,~isnan(t_output))); if vmin<0; vmin=1.2*vmin; else; vmin=0.9*vmin; end
vmax=max(DTCI2(id_tr,~isnan(t_output))); if vmax<0; vmax=0.9*vmax; else; vmax=1.2*vmax; end

if isnan(vmin); vmin=lb_DT; end
if isnan(vmax); vmax=ub_DT; end

han1=plot([tforecast tforecast],[vmin vmax],'--k'); set(han1,'linewidth',2);
han1=plot([tforecast2 tforecast2],[vmin vmax],'--k'); set(han1,'linewidth',2);

if MAKE_ANIMATION
    axis([t0 t1(min(n+tmore,length(t1))) vmin vmax]);
else
    axis([t0 tforecast2+2*365 vmin vmax]);
end

ylabel({'time-scale','parameter','DT [s]'});
set(gca,'FontSize',14,'XTickLabel',[]);
ax = gca; ax.YRuler.Exponent = 0;

nplot=4; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;

han1=fill([t_output(~isnan(t_output)) fliplr(t_output(~isnan(t_output)))],[DYCI1(id_tr,~isnan(t_output)) fliplr(DYCI2(id_tr,~isnan(t_output)))],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
plot(t_output,DYOUT(id_tr,:),'-r'); set(gca,'XTick',linspace(t(1),t(min(n+tmore,length(t))),12),'Xlim',[t(1) t(min(n+tmore,length(t)))]); datetick('x','keeplimits','keepticks');

vmin=min(DYCI1(id_tr,~isnan(t_output))); if vmin<0; vmin=1.2*vmin; else; vmin=0.9*vmin; end
vmax=max(DYCI2(id_tr,~isnan(t_output))); if vmax<0; vmax=0.9*vmax; else; vmax=1.2*vmax; end

if isnan(vmin); vmin=lb_DY; end
if isnan(vmax); vmax=ub_DY; end

han1=plot([tforecast tforecast],[vmin vmax],'--k'); set(han1,'linewidth',2);
han1=plot([tforecast2 tforecast2],[vmin vmax],'--k'); set(han1,'linewidth',2);

if MAKE_ANIMATION
    axis([t0 t1(min(n+tmore,length(t1))) vmin vmax]);
else
    axis([t0 tforecast2+2*365 vmin vmax]);
end

ylabel({'excursion','parameter','DY [m]'});
set(gca,'FontSize',14,'XTickLabel',[]);
ax = gca; ax.YRuler.Exponent = 0;

nplot=5; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;

han1=plot([t(1) t(min(n+tmore,length(t)))],[HSB0(id_tr) HSB0(id_tr)],'--b'); set(han1,'LineWidth',2);

han1=fill([t_output(~isnan(t_output)) fliplr(t_output(~isnan(t_output)))],[HSBCI1(id_tr,~isnan(t_output)) fliplr(HSBCI2(id_tr,~isnan(t_output)))],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
plot(t_output,HSBOUT(id_tr,:),'-r'); set(gca,'XTick',linspace(t(1),t(min(n+tmore,length(t))),12),'Xlim',[t(1) t(min(n+tmore,length(t)))]); datetick('x','keeplimits','keepticks');

vmin=min(HSBCI1(id_tr,~isnan(t_output))); vmin=min(vmin,HSB0(id_tr)); if vmin<0; vmin=1.2*vmin; else; vmin=0.9*vmin; end
vmax=max(HSBCI2(id_tr,~isnan(t_output))); vmax=max(vmax,HSB0(id_tr)); if vmax<0; vmax=0.9*vmax; else; vmax=1.2*vmax; end

if isnan(vmin); vmin=lb_HSB; end
if isnan(vmax); vmax=ub_HSB; end

han1=plot([tforecast tforecast],[vmin vmax],'--k'); set(han1,'linewidth',2);
han1=plot([tforecast2 tforecast2],[vmin vmax],'--k'); set(han1,'linewidth',2);

if MAKE_ANIMATION
    axis([t0 t1(min(n+tmore,length(t1))) vmin vmax]);
else
    axis([t0 tforecast2+2*365 vmin vmax]);
end

ylabel({'wave height','parameter','(H_s)_b [m]'});
set(gca,'FontSize',14,'XTickLabel',[]);
ax = gca; ax.YRuler.Exponent = 0;

nplot=6; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;
han1=fill([t_output(~isnan(t_output)) fliplr(t_output(~isnan(t_output)))],[KCI1(id_tr,~isnan(t_output)) fliplr(KCI2(id_tr,~isnan(t_output)))],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
plot(t_output,KOUT(id_tr,:),'-r'); set(gca,'XTick',linspace(t(1),t(min(n+tmore,length(t))),12),'Xlim',[t(1) t(min(n+tmore,length(t)))]); datetick('x','keeplimits','keepticks');

vmin=min(KCI1(id_tr,~isnan(t_output))); if vmin<0; vmin=1.2*vmin; else; vmin=0.9*vmin; end
vmax=max(KCI2(id_tr,~isnan(t_output))); if vmax<0; vmax=0.9*vmax; else; vmax=1.2*vmax; end

if isnan(vmin); vmin=lb_K; end
if isnan(vmax); vmax=ub_K; end

han1=plot([tforecast tforecast],[vmin vmax],'--k'); set(han1,'linewidth',2);
han1=plot([tforecast2 tforecast2],[vmin vmax],'--k'); set(han1,'linewidth',2);

if MAKE_ANIMATION
    axis([t0 t1(min(n+tmore,length(t1))) vmin vmax]);
else
    axis([t0 tforecast2+2*365 vmin vmax]);
end

ylabel({'longshore','transport','parameter','K [-]'});
set(gca,'FontSize',14,'XTickLabel',[]);
ax = gca; ax.YRuler.Exponent = 0;

nplot=7; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;

han1=plot([t(1) t(min(n+tmore,length(t)))],[LTER(id_tr) LTER(id_tr)],'--b'); set(han1,'LineWidth',2);

han1=fill([t_output(~isnan(t_output)) fliplr(t_output(~isnan(t_output)))],[VLTCI1(id_tr,~isnan(t_output)) fliplr(VLTCI2(id_tr,~isnan(t_output)))],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
plot(t_output,VLTOUT(id_tr,:),'-r'); set(gca,'XTick',linspace(t(1),t(min(n+tmore,length(t))),12),'Xlim',[t(1) t(min(n+tmore,length(t)))]); datetick('x','keeplimits','keepticks');

vmin=min(VLTCI1(id_tr,~isnan(t_output))); vmin=min(vmin,LTER(id_tr)); if vmin<0; vmin=1.2*vmin; else; vmin=0.9*vmin; end
vmax=max(VLTCI2(id_tr,~isnan(t_output))); vmax=max(vmax,LTER(id_tr)); if vmax<0; vmax=0.9*vmax; else; vmax=1.2*vmax; end

if isnan(vmin); vmin=lb_vlt; end
if isnan(vmax); vmax=ub_vlt; end

han1=plot([tforecast tforecast],[vmin vmax],'--k'); set(han1,'linewidth',2);
han1=plot([tforecast2 tforecast2],[vmin vmax],'--k'); set(han1,'linewidth',2);

if MAKE_ANIMATION
    axis([t0 t1(min(n+tmore,length(t1))) vmin vmax]);
else
    axis([t0 tforecast2+2*365 vmin vmax]);
end

ylabel({'long-term rate','parameter','v_{lt} [m/yr]'});
set(gca,'FontSize',14,'XTickLabel',[]);
ax = gca; ax.YRuler.Exponent = 0;

nplot=8; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;

han1=plot([t(1) t(min(n+tmore,length(t)))],[sigma00 sigma00],'--b'); set(han1,'LineWidth',2);

han1=fill([t_output(~isnan(t_output)) fliplr(t_output(~isnan(t_output)))],[SIGMACI1(id_tr,~isnan(t_output)) fliplr(SIGMACI2(id_tr,~isnan(t_output)))],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
plot(t_output,SIGMAOUT(id_tr,:),'-r'); set(gca,'XTick',linspace(t(1),t(min(n+tmore,length(t))),12),'Xlim',[t(1) t(min(n+tmore,length(t)))]); datetick('x','keeplimits','keepticks');

vmin=min(SIGMACI1(id_tr,~isnan(t_output))); if vmin<0; vmin=1.2*vmin; else; vmin=0.9*vmin; end
vmax=max(SIGMACI2(id_tr,~isnan(t_output))); if vmax<0; vmax=0.9*vmax; else; vmax=1.2*vmax; end

if isnan(vmin); vmin=lb_sigma; end
if isnan(vmax); vmax=ub_sigma; end

han1=plot([tforecast tforecast],[vmin vmax],'--k'); set(han1,'linewidth',2);
han1=plot([tforecast2 tforecast2],[vmin vmax],'--k'); set(han1,'linewidth',2);

if MAKE_ANIMATION
    axis([t0 t1(min(n+tmore,length(t1))) vmin vmax]);
else
    axis([t0 tforecast2+2*365 vmin vmax]);
end

ylabel({'noise','parameter','\sigma [var.]'});
set(gca,'FontSize',14);
ax = gca; ax.YRuler.Exponent = 0;

drawnow

if MAKE_ANIMATION % make gif animation
    
    % grab frame
    frame = getframe(gcf);
    im=frame2im(frame);
    
    [imind,cm]=rgb2ind(im,256);
    
    if count == 1
        imwrite(imind,cm,filename2,'gif', 'Loopcount',inf,'delaytime',1/fps);
    else
        imwrite(imind,cm,filename2,'gif','WriteMode','append','delaytime',1/fps);
    end
    
    count=count+1;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3); clf;

ALPHA=0.25;

tmore=60;

xp0=0.16;
yp0=0.84;
width=0.81;
height=0.120;
ysep=0.01;

NTR=ids((1:Ntr))';

%Y00=nanmean(Y0,2);

nn=find(t_obs<=t1(min(n+1,length(t))),1,'last');
nnn=find(t_output<=t1(min(n+1,length(t))),1,'last');

nplot=1; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;

if ~isempty(nn) && DATA_ASSIMILATION
    han1=errorbar(NTR,Y_obs(:,nn)-Y00,2*Y_rms(:,nn),'b'); set(han1,'linestyle','none','color',[0.4 0.6 1]);
    han1=plot(NTR,Y_obs(:,nn)-Y00,'bo'); set(han1,'MarkerSize',4,'MarkerFaceColor',[0.4 0.6 1]);
end

id_nan=isnan(YCI1(:,nnn));
idsec1=find(id_nan(1:end-1) & ~id_nan(2:end))+1; % if a nan, then becomes not a nan
idsec2=find(~id_nan(1:end-1) & id_nan(2:end));   % if not a nan, then becomes a nan
if ~id_nan(1)                % if the first entry is not a nan
    idsec1=cat(1,1,idsec1);  % account for it
end
if ~id_nan(end)            % if the last entry is not a nan
    idsec2=cat(1,idsec2,Ntr); % account for it
end

Nsections=length(idsec1);

for iiii=1:Nsections;
    id_not_nan=(idsec1(iiii):idsec2(iiii));
    han1=fill([NTR(id_not_nan); flipud(NTR(id_not_nan))],[YCI1(id_not_nan,nnn)-Y00(id_not_nan); flipud(YCI2(id_not_nan,nnn)-Y00(id_not_nan))],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
end

plot(NTR,YY(:,nnn)-Y00,'-r');

han1=plot(NTR,Ymin-Y00,'k--'); set(han1,'LineWidth',2);

if PLOT_SECTIONS==1
    if strcmp(transects(id_tr).model_type,'full model')
        id=find(id_full_model_end>=id_tr,1,'first');
        xmin=id_full_model_start(id);
        xmax=id_full_model_end(id);
    elseif strcmp(transects(id_tr).model_type,'cross-shore only')
        id=find(id_cross_shore_only_end>=id_tr,1,'first');
        xmin=id_cross_shore_only_start(id)-1;
        xmax=id_cross_shore_only_end(id)+1;
    elseif  strcmp(transects(id_tr).model_type,'rate only')
        id=find(id_rate_only_end>=id_tr,1,'first');
        xmin=id_rate_only_start(id)-1;
        xmax=id_rate_only_end(id)+1;
    end
else
    xmin=min(NTR);
    xmax=max(NTR);
end

vmin=min(YCI1(id_not_nan,nnn)-Y00(id_not_nan)); if vmin<0; vmin=1.2*vmin; else; vmin=0.9*vmin; end
vmax=max(YCI2(id_not_nan,nnn)-Y00(id_not_nan)); if vmax<0; vmax=0.9*vmax; else; vmax=1.2*vmax; end

han1=plot([ids(id_tr) ids(id_tr)],[vmin vmax],'--r'); set(han1,'linewidth',2);

title(sprintf('CoSMoS-COAST %s: alongshore model parameters (transect ID: %d)',Model_name,ID(id_tr)),'FontSize',14);

ylabel({'total','shoreline','position [m]'});
set(gca,'FontSize',14,'XTickLabel',[],'Xlim',[xmin xmax],'Ylim',[vmin vmax]);
ax = gca; ax.YRuler.Exponent = 0;

nplot=2; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;

id_nan=isnan(DTOUT(:,nnn));
idsec1=find(id_nan(1:end-1) & ~id_nan(2:end))+1; % if a nan, then becomes not a nan
idsec2=find(~id_nan(1:end-1) & id_nan(2:end));   % if not a nan, then becomes a nan
if ~id_nan(1)              % if the first entry is not a nan
    idsec1=cat(1,1,idsec1);  % account for it
end
if ~id_nan(end)            % if the last entry is not a nan
    idsec2=cat(1,idsec2,Ntr); % account for it
end

Nsections=length(idsec1);

vmin=nanmean(DTOUT(:,nnn));
vmax=vmin;

for iiii=1:Nsections
    id_not_nan=(idsec1(iiii):idsec2(iiii));
    han1=fill([NTR(id_not_nan); flipud(NTR(id_not_nan))],[DTCI1(id_not_nan,nnn); flipud(DTCI2(id_not_nan,nnn))],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
    
    vmin=min(vmin,min(DTCI1(id_not_nan,nnn)));
    vmax=max(vmax,max(DTCI2(id_not_nan,nnn)));
end

plot(NTR,DTOUT(:,nnn),'-r');

han1=plot([ids(id_tr) ids(id_tr)],[vmin vmax],'--r'); set(han1,'linewidth',2);

ylabel({'time-scale','parameter','DT [s]'});
set(gca,'FontSize',14,'XTickLabel',[],'Xlim',[xmin xmax],'Ylim',[vmin vmax]);
ax = gca; ax.YRuler.Exponent = 0;

nplot=3; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;

id_nan=isnan(DYOUT(:,nnn));
idsec1=find(id_nan(1:end-1) & ~id_nan(2:end))+1; % if a nan, then becomes not a nan
idsec2=find(~id_nan(1:end-1) & id_nan(2:end));   % if not a nan, then becomes a nan
if ~id_nan(1)              % if the first entry is not a nan
    idsec1=cat(1,1,idsec1);  % account for it
end
if ~id_nan(end)            % if the last entry is not a nan
    idsec2=cat(1,idsec2,Ntr); % account for it
end

Nsections=length(idsec1);

vmin=nanmean(DYOUT(:,nnn));
vmax=vmin;

for iiii=1:Nsections
    id_not_nan=(idsec1(iiii):idsec2(iiii));
    han1=fill([NTR(id_not_nan); flipud(NTR(id_not_nan))],[DYCI1(id_not_nan,nnn); flipud(DYCI2(id_not_nan,nnn))],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
    
    vmin=min(vmin,min(DYCI1(id_not_nan,nnn)));
    vmax=max(vmax,max(DYCI2(id_not_nan,nnn)));
end

plot(NTR,DYOUT(:,nnn),'-r');

han1=plot([ids(id_tr) ids(id_tr)],[vmin vmax],'--r'); set(han1,'linewidth',2);

ylabel({'excursion','parameter','DY [m]'});
set(gca,'FontSize',14,'XTickLabel',[],'Xlim',[xmin xmax],'Ylim',[vmin vmax]);
ax = gca; ax.YRuler.Exponent = 0;

nplot=4; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;

id_nan=isnan(HSBOUT(:,nnn));
idsec1=find(id_nan(1:end-1) & ~id_nan(2:end))+1; % if a nan, then becomes not a nan
idsec2=find(~id_nan(1:end-1) & id_nan(2:end));   % if not a nan, then becomes a nan
if ~id_nan(1)              % if the first entry is not a nan
    idsec1=cat(1,1,idsec1);  % account for it
end
if ~id_nan(end)            % if the last entry is not a nan
    idsec2=cat(1,idsec2,Ntr); % account for it
end

han1=plot(NTR,HSB0,'-b'); set(han1,'LineWidth',2);

Nsections=length(idsec1);

for iiii=1:Nsections
    id_not_nan=(idsec1(iiii):idsec2(iiii));
    han1=fill([NTR(id_not_nan); flipud(NTR(id_not_nan))],[HSBCI1(id_not_nan,nnn); flipud(HSBCI2(id_not_nan,nnn))],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
end

plot(NTR,HSBOUT(:,nnn),'-r');

vmin=min(HSBCI1(id_not_nan,nnn)); vmin=min(vmin,min(HSB0)); if vmin<0; vmin=1.2*vmin; else; vmin=0.9*vmin; end
vmax=max(HSBCI2(id_not_nan,nnn)); vmax=max(vmax,max(HSB0)); if vmax<0; vmax=0.9*vmax; else; vmax=1.2*vmax; end

han1=plot([ids(id_tr) ids(id_tr)],[vmin vmax],'--r'); set(han1,'linewidth',2);

ylabel({'wave height','parameter','(H_s)_b [m]'});
set(gca,'FontSize',14,'XTickLabel',[],'Xlim',[xmin xmax],'Ylim',[vmin vmax]);
ax = gca; ax.YRuler.Exponent = 0;

nplot=5; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;

id_nan=isnan(KOUT(:,nnn));
idsec1=find(id_nan(1:end-1) & ~id_nan(2:end))+1; % if a nan, then becomes not a nan
idsec2=find(~id_nan(1:end-1) & id_nan(2:end));   % if not a nan, then becomes a nan
if ~id_nan(1)              % if the first entry is not a nan
    idsec1=cat(1,1,idsec1);  % account for it
end
if ~id_nan(end)            % if the last entry is not a nan
    idsec2=cat(1,idsec2,Ntr); % account for it
end

Nsections=length(idsec1);

vmin=nanmean(KOUT(:,nnn));
vmax=vmin;

for iiii=1:Nsections
    id_not_nan=(idsec1(iiii):idsec2(iiii));
    han1=fill([NTR(id_not_nan); flipud(NTR(id_not_nan))],[KCI1(id_not_nan,nnn); flipud(KCI2(id_not_nan,nnn))],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
    
    vmin=min(vmin,min(KCI1(id_not_nan,nnn)));
    vmax=max(vmax,max(KCI2(id_not_nan,nnn)));
end

plot(NTR,KOUT(:,nnn),'-r');

han1=plot([ids(id_tr) ids(id_tr)],[vmin vmax],'--r'); set(han1,'linewidth',2);

ylabel({'longshore','transport','parameter','K [-]'});
set(gca,'FontSize',14,'XTickLabel',[],'Xlim',[xmin xmax],'Ylim',[vmin vmax]);
ax = gca; ax.YRuler.Exponent = 0;

nplot=6; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;

id_nan=isnan(VLTOUT(:,nnn));
idsec1=find(id_nan(1:end-1) & ~id_nan(2:end))+1; % if a nan, then becomes not a nan
idsec2=find(~id_nan(1:end-1) & id_nan(2:end));   % if not a nan, then becomes a nan
if ~id_nan(1)              % if the first entry is not a nan
    idsec1=cat(1,1,idsec1);  % account for it
end
if ~id_nan(end)            % if the last entry is not a nan
    idsec2=cat(1,idsec2,Ntr); % account for it
end

Nsections=length(idsec1);

han1=plot(NTR,LTER,'-b'); set(han1,'LineWidth',2);

for iiii=1:Nsections
    id_not_nan=(idsec1(iiii):idsec2(iiii));
    han1=fill([NTR(id_not_nan); flipud(NTR(id_not_nan))],[VLTCI1(id_not_nan,nnn); flipud(VLTCI2(id_not_nan,nnn))],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
end

plot(NTR,VLTOUT(:,nnn),'-r');

vmin=min(VLTCI1(id_not_nan,nnn)); vmin=min(vmin,min(LTER)); if vmin<0; vmin=1.2*vmin; else; vmin=0.9*vmin; end
vmax=max(VLTCI2(id_not_nan,nnn)); vmax=max(vmax,max(LTER)); if vmax<0; vmax=0.9*vmax; else; vmax=1.2*vmax; end

han1=plot([ids(id_tr) ids(id_tr)],[vmin vmax],'--r'); set(han1,'linewidth',2);

ylabel({'long-term rate','parameter','v_{lt} [m/yr]'});
set(gca,'FontSize',14,'XTickLabel',[],'Xlim',[xmin xmax],'Ylim',[vmin vmax]);
ax = gca; ax.YRuler.Exponent = 0;

nplot=7; subplot('position',[xp0 yp0-(height+ysep)*(nplot-1) width height]); hold on; box on;

id_nan=isnan(SIGMAOUT(:,nnn));
idsec1=find(id_nan(1:end-1) & ~id_nan(2:end))+1; % if a nan, then becomes not a nan
idsec2=find(~id_nan(1:end-1) & id_nan(2:end));   % if not a nan, then becomes a nan
if ~id_nan(1)              % if the first entry is not a nan
    idsec1=cat(1,1,idsec1);  % account for it
end
if ~id_nan(end)            % if the last entry is not a nan
    idsec2=cat(1,idsec2,Ntr); % account for it
end

Nsections=length(idsec1);

han1=plot([min(NTR) max(NTR)],[sigma00 sigma00],'-b'); set(han1,'LineWidth',2);

for iiii=1:Nsections
    id_not_nan=(idsec1(iiii):idsec2(iiii));
    han1=fill([NTR(id_not_nan); flipud(NTR(id_not_nan))],[SIGMACI1(id_not_nan,nnn); flipud(SIGMACI2(id_not_nan,nnn))],'r');  set(han1,'FaceAlpha',ALPHA,'EdgeColor','none');
end

plot(NTR,SIGMAOUT(:,nnn),'-r');

vmin=min(SIGMACI1(id_not_nan,nnn)); if vmin<0; vmin=1.2*vmin; else; vmin=0.9*vmin; end
vmax=max(SIGMACI2(id_not_nan,nnn)); if vmax<0; vmax=0.9*vmax; else; vmax=1.2*vmax; end

han1=plot([ids(id_tr) ids(id_tr)],[vmin vmax],'--r'); set(han1,'linewidth',2);

xlabel('Transect #');
ylabel({'noise','parameter','\sigma [var.]'});
set(gca,'FontSize',14,'Xlim',[xmin xmax],'Ylim',[vmin vmax]);
ax = gca; ax.YRuler.Exponent = 0;

drawnow

if MAKE_ANIMATION % make gif animation
    
    % grab frame
    frame2 = getframe(gcf);
    im2=frame2im(frame2);
    
    [imind2,cm2]=rgb2ind(im2,256);
    
    if count2 == 1
        imwrite(imind2,cm2,filename3,'gif', 'Loopcount',inf,'delaytime',1/fps2);
    else
        imwrite(imind2,cm2,filename3,'gif','WriteMode','append','delaytime',1/fps2);
    end
    
    count2=count2+1;
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4); clf;

ALPHA=0.2;

xp0=0.08;
yp0=0.08;

width=0.9;
height=0.38;

ysep=0.12;

subplot('Position',[xp0 yp0+height+ysep width height]); hold on; box on;

han1=plot(t,Hs(WAVEID(id_tr),:,1),'-b');

vmin=0;
vmax=1.1*max(Hs(id_tr,:));
axis([t0 t1(min(n+tmore,length(t1))) vmin vmax]);

set(gca,'Xlim',[t0 tforecast2+1*365.25],'Ylim',[vmin vmax]);
set(gca,'xtick',datenum(1995:2:2021,1,1)); datetick('x','keeplimits','keepticks');
set(gca,'Xlim',[t0 tforecast2+1*365.25],'Ylim',[vmin vmax]);
set(gca,'FontSize',14);

title(sprintf('CoSMoS-COAST %s: model calibration/validation (transect ID: %d)',Model_name,ID(id_tr)),'FontSize',14);

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
plot(t_output,YY(id_tr,:)-Y00(id_tr),'-r'); set(gca,'XTick',linspace(t(1),t(min(n+tmore,length(t))),12),'Xlim',[t(1) t(min(n+tmore,length(t)))]);

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

% RMS error during forecast period
tobs=t_obs(nn);
obs=Y_obs(id_tr,nn)-Y0(id_tr);

id=(tobs<tforecast | tobs>tforecast2);
tobs(id)=[];
obs(id)=[];

tmod=t_output(~isnan(t_output));
model=YY(id_tr,~isnan(t_output))-Y0(id_tr);
obs_CI1=obs-2*y_rms_sat;
obs_CI2=obs+2*y_rms_sat;

if length(model)>=3
    mod_obs=interp1(tmod,model,tobs);
else
    mod_obs=NaN*obs;
end

y_error=mod_obs-obs;
RMSE(i)=sqrt(nanmean((y_error).^2));
within_CI=(mod_obs>=obs_CI1 & mod_obs<=obs_CI2);
PCT(i)=(sum(within_CI)./length(mod_obs))*100;

text(0.5*(tforecast+tforecast2)-2.5*365,vmin+0.2*(vmax-vmin),['RMSE = ',num2str(RMSE(i),3),' m'],'FontSize',14);
text(0.5*(tforecast+tforecast2)-2.5*365,vmin+0.075*(vmax-vmin),['within CI = ',num2str(PCT(i),3),'%'],'FontSize',14);

drawnow
