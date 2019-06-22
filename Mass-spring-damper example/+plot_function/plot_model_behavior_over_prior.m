function plot_model_behavior_over_prior(envir_data,t_data,z_measure)

b = envir_data.b;
y0 = [envir_data.z0; envir_data.v0];
c2 = y0(1);
Hk = envir_data.H_k;
Ht = envir_data.H_t;
x_true = @(k,t) exp(-0.5*b*t).*(2*(y0(2)+0.5*b*c2)/sqrt(4*k-b^2)*sin(0.5*sqrt(4*k-b^2)*t)+c2*cos(0.5*sqrt(4*k-b^2)*t));
x_mi = @(k,t) y0(2)/sqrt(k)*sin(sqrt(k)*t)+c2*cos(sqrt(k)*t);
ymeasure = z_measure{1};
eps = envir_data.eps(1);
eps = eps*ones(length(t_data),1);
nSample = 1e2;
nS = 1e4;
kSample = linspace(Hk(1),Hk(2),nSample);
tsample = linspace(Ht(1),Ht(2),nS);
z_true = zeros(nSample,nS);
z_inadequate = zeros(nSample,nS);
for i = 1:nSample
   z_true(i,:) = x_true(kSample(i),tsample);
   z_inadequate(i,:) = x_mi(kSample(i),tsample);
end
yshift = 0.005;
y_shift = x_mi(0.20001,tsample)+yshift;
f = figure('Name','Model behavior over prior uncertainty of k','Units','inches');
lp = [0.08 0.55];
bp = 0.12;
ww = 0.43;
hh = 0.85;
f.Position = [2,1,8,4.5];
subplot(1,2,1);
hold on
for i = 1:nSample
   plot(tsample,z_true(i,:),'c','LineWidth',2);
end
errorbar(t_data,ymeasure,eps,'k.','LineWidth',2,'MarkerSize',1e-5);
hold off
box on
set(gca,'Position',[lp(1) bp ww hh],'FontSize',12,'LineWidth',1.5,...
   'YLim',[-1.7,2],'YTick',-2:2,'XTick',0:3);
text(-0.1242,0.3383,'Displacement','Units','normalized','FontSize',15,'Rotation',90);
text(0.4860,-0.0963,'\itt','Units','normalized','FontSize',15);
subplot(1,2,2);
hold on
for i = 1:nSample
   plot(tsample,z_inadequate(i,:),'c','LineWidth',2);
end
plot(tsample,y_shift,'r--','LineWidth',3);
errorbar(t_data,ymeasure,eps,'k.','LineWidth',2,'MarkerSize',1e-5);
hold off
set(gca,'Position',[lp(2) bp ww hh],'FontSize',12,'LineWidth',1.5,'YTickLabel',[],...
   'YLim',[-1.7,2],'YTick',-2:2,'YTickLabel',[],'XTick',0:3);
text(0.4860,-0.0963,'\itt','Units','normalized','FontSize',15);
box on
plotID = tsample >= 2.3 & tsample <=2.8;
f = figure('Name','Inset 2','Units','inches');
f.Position = [2,1,2,2];
hold on
for i = 1:nSample
   plot(tsample(plotID),z_true(i,plotID),'c','LineWidth',2);
end
errorbar(t_data(end-1:end),ymeasure(end-1:end),eps(end-1:end),'k.','LineWidth',3.5,'MarkerSize',1e-5);
hold off
box on
set(gca,'Position',[0.03 0.03 0.94 0.94],'FontSize',12,'LineWidth',1.5,...
   'XTick',[],'YTick',[],'XLim',[2.3,2.8]);
f = figure('Name','Inset 1','Units','inches');
f.Position = [2,1,2,2];
hold on
for i = 1:nSample
   plot(tsample(plotID),z_inadequate(i,plotID),'c','LineWidth',2);
end
plot(tsample(plotID),y_shift(plotID),'r--','LineWidth',3);
errorbar(t_data(end-1:end),ymeasure(end-1:end),eps(end-1:end),'k.','LineWidth',3.5,'MarkerSize',1e-5);
hold off
box on
set(gca,'Position',[0.03 0.03 0.94 0.94],'FontSize',12,'LineWidth',1.5,...
   'XTick',[],'YTick',[],'XLim',[2.3,2.8]);