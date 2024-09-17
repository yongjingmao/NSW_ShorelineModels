t_output(count_output)=t1(n+1)+increment;

Yout=Y; % Yout is the main shoreline results variable

if HOLD_THE_LINE   
    % clip to minimum (non-erodible) shoreline
    Yout=max(Yout,Ymin);
end

phi_scale = 10.^phi;
dummyC = sqrt(2229)*24*3600;
ca_scale = ca/dummyC;
ce_scale = ce/dummyC;
dummy_K = 1/16/(2.65-1)*sqrt(9.81/0.78)*24*3600;
K_scale = K/dummy_K;


% phi_scale = phi;
% ca_scale = ca;
% ce_scale = ce;


YY(:,count_output)=quantile(Yout,0.5,2);
YCI1(:,count_output)=quantile(Yout,0.025,2)-0.0001;  % 95% confidence bands
YCI2(:,count_output)=quantile(Yout,0.975,2)+0.0001;

YST(:,count_output)=quantile(Yst,0.5,2);
YSTCI1(:,count_output)=quantile(Yst,0.025,2)-0.0001;  % 95% confidence bands
YSTCI2(:,count_output)=quantile(Yst,0.975,2)+0.0001;

YLST(:,count_output)=quantile(Ylst,0.5,2);
YLSTCI1(:,count_output)=quantile(Ylst,0.025,2)-0.0001;  % 95% confidence bands
YLSTCI2(:,count_output)=quantile(Ylst,0.975,2)+0.0001;

YLST_only(:,count_output)=quantile(Ylst_only,0.5  ,2);
YLST_only_CI1(:,count_output)=quantile(Ylst_only,0.025,2)-0.0001;  % 95% confidence bands
YLST_only_CI2(:,count_output)=quantile(Ylst_only,0.975,2)+0.0001;

YBRU(:,count_output)=quantile(Ybru,0.5,2);
YBRUCI1(:,count_output)=quantile(Ybru,0.025,2)-0.0001;  % 95% confidence bands
YBRUCI2(:,count_output)=quantile(Ybru,0.975,2)+0.0001;

YVLT(:,count_output)=quantile(Yvlt,0.5,2);
YVLTCI1(:,count_output)=quantile(Yvlt,0.025,2)-0.0001;  % 95% confidence bands
YVLTCI2(:,count_output)=quantile(Yvlt,0.975,2)+0.0001;

CAOUT(:,count_output)=quantile(ca_scale,0.5,2);
CACI1(:,count_output)=quantile(ca_scale,0.025,2);  % 95% confidence bands
CACI2(:,count_output)=quantile(ca_scale,0.975,2);

CEOUT(:,count_output)=quantile(ce_scale,0.5,2);
CECI1(:,count_output)=quantile(ce_scale,0.025,2);  % 95% confidence bands
CECI2(:,count_output)=quantile(ce_scale,0.975,2);

PHIOUT(:,count_output)=quantile(phi_scale,0.5,2);
PHICI1(:,count_output)=quantile(phi_scale,0.025,2)-1e-6;  % 95% confidence bands
PHICI2(:,count_output)=quantile(phi_scale,0.975,2)+1e-6;

COUT(:,count_output)=quantile(c,0.5,2);
CCI1(:,count_output)=quantile(c,0.025,2)-1e-6;  % 95% confidence bands
CCI2(:,count_output)=quantile(c,0.975,2)+1e-6;

KOUT(:,count_output)=quantile(K_scale,0.5,2);
KCI1(:,count_output)=quantile(K_scale,0.025,2)-1e-6;  % 95% confidence bands
KCI2(:,count_output)=quantile(K_scale,0.975,2)+1e-6;

VLTOUT(:,count_output)=quantile(vlt,0.5,2);
VLTCI1(:,count_output)=quantile(vlt,0.025,2)-1e-6;  % 95% confidence bands
VLTCI2(:,count_output)=quantile(vlt,0.975,2)+1e-6;

SIGMAOUT(:,count_output)=quantile(sigma,0.5,2);
SIGMACI1(:,count_output)=quantile(sigma,0.025,2)-1e-6;  % 95% confidence bands
SIGMACI2(:,count_output)=quantile(sigma,0.975,2)+1e-6;

count_output=count_output+1;