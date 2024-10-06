t_output(count_output)=t1(n+1)+increment;

Yout=Y; % Yout is the main shoreline results variable

if HOLD_THE_LINE   
    % clip to minimum (non-erodible) shoreline
    Yout=max(Yout,Ymin);
end

YY(:,count_output)=quantile(Yout,0.5,2);
YCI1(:,count_output)=quantile(Yout,0.025,2)-0.0001;  % 95% confidence bands
YCI2(:,count_output)=quantile(Yout,0.975,2)+0.0001;

YST(:,count_output)=quantile(Yst,0.5,2);
YSTCI1(:,count_output)=quantile(Yst,0.025,2)-0.0001;  % 95% confidence bands
YSTCI2(:,count_output)=quantile(Yst,0.975,2)+0.0001;

YLST(:,count_output)=quantile(Ylst,0.5,2);
YLSTCI1(:,count_output)=quantile(Ylst,0.025,2)-0.0001;  % 95% confidence bands
YLSTCI2(:,count_output)=quantile(Ylst,0.975,2)+0.0001;

QCUMOUT(:,count_output)=quantile(Qcum,0.5,2);
QCUMOUTCI1(:,count_output)=quantile(Qcum,0.025,2)-1e-6; % 95% confidence bands
QCUMOUTCI2(:,count_output)=quantile(Qcum,0.975,2)+1e-6;

% YLST_only(:,count_output)=quantile(Ylst_only,0.5  ,2);
% YLST_only_CI1(:,count_output)=quantile(Ylst_only,0.025,2)-0.0001;  % 95% confidence bands
% YLST_only_CI2(:,count_output)=quantile(Ylst_only,0.975,2)+0.0001;

YBRU(:,count_output)=quantile(Ybru,0.5,2);
YBRUCI1(:,count_output)=quantile(Ybru,0.025,2)-0.0001;  % 95% confidence bands
YBRUCI2(:,count_output)=quantile(Ybru,0.975,2)+0.0001;

YVLT(:,count_output)=quantile(Yvlt,0.5,2);
YVLTCI1(:,count_output)=quantile(Yvlt,0.025,2)-0.0001;  % 95% confidence bands
YVLTCI2(:,count_output)=quantile(Yvlt,0.975,2)+0.0001;

DTOUT(:,count_output)=quantile(DT,0.5,2);
DTCI1(:,count_output)=quantile(DT,0.025,2)-1e-6;  % 95% confidence bands
DTCI2(:,count_output)=quantile(DT,0.975,2)+1e-6;

DYOUT(:,count_output)=quantile(DY,0.5,2);
DYCI1(:,count_output)=quantile(DY,0.025,2)-1e-6;  % 95% confidence bands
DYCI2(:,count_output)=quantile(DY,0.975,2)+1e-6;

HSBOUT(:,count_output)=quantile(HSB,0.5,2);
HSBCI1(:,count_output)=quantile(HSB,0.025,2)-1e-6;  % 95% confidence bands
HSBCI2(:,count_output)=quantile(HSB,0.975,2)+1e-6;

COUT(:,count_output)=quantile(c,0.5,2);
CCI1(:,count_output)=quantile(c,0.025,2)-1e-6;  % 95% confidence bands
CCI2(:,count_output)=quantile(c,0.975,2)+1e-6;

KOUT(:,count_output)=quantile(K,0.5,2);
KCI1(:,count_output)=quantile(K,0.025,2)-1e-6;  % 95% confidence bands
KCI2(:,count_output)=quantile(K,0.975,2)+1e-6;

VLTOUT(:,count_output)=quantile(vlt,0.5,2);
VLTCI1(:,count_output)=quantile(vlt,0.025,2)-1e-6;  % 95% confidence bands
VLTCI2(:,count_output)=quantile(vlt,0.975,2)+1e-6;

SIGMAOUT(:,count_output)=quantile(sigma,0.5,2);
SIGMACI1(:,count_output)=quantile(sigma,0.025,2)-1e-6;  % 95% confidence bands
SIGMACI2(:,count_output)=quantile(sigma,0.975,2)+1e-6;

count_output=count_output+1;