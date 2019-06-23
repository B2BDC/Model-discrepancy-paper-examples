function plot_Confounding(confounding_result)

figure('Units','inches','Position',[2 1 8 6]);
subplot(2,2,1);
X = confounding_result.xmesh;
Y = confounding_result.ymesh;
Z = confounding_result.zmesh;
[x,y,z] = supplementary_function.createPatch3D(X,Y,Z);
lims = [min(X(:)), max(X(:)); min(Y(:)), max(Y(:)); min(Z(:)) max(Z(:))];
xtick = 0.2:0.05:0.3;
ytick = [0.01 0.03];
ztick = -0.09:0.03:0;


patch(x,y,z,10*z,'FaceAlpha',0.4,'FaceColor','interp','EdgeColor','interp');
view(3)
box off
grid on
cl = colorbar;
cl.Ticks = 10*ztick;
cl.TickLabels = ztick;
cl.Units = 'normalized';
cl.Position = [0.4570 0.5850 0.0195 0.3750];
text(0.5335,-0.0405,'\itk','FontSize',15,'Units','normalized');
text(-0.1097,0.0504,'{\itc}_0','FontSize',15,'Units','normalized');
text(-0.2063,0.6496,'{\itc}_1','FontSize',15,'Units','normalized');
set(gca,'LineWidth',1.5,'FontSize',12,'XLim',lims(1,:),'YLim',lims(2,:),...
   'ZLim',lims(3,:),'XTick',xtick,'YTick',ytick,'ZTick',ztick,...
   'Position',[0.09,0.55,0.35,0.42],'View',[-17.229028831511030,23.343466376044585]);
% k vs c0
subplot(2,2,2);
X = confounding_result.xdata;
Y = confounding_result.yVSx;
[x,y] = supplementary_function.createPatch2D(X,Y);
patch(x,y,zeros(size(y)),'FaceAlpha',1,'FaceColor','r','EdgeColor','red');
box on
text(0.4863,-0.1415,'\itk','FontSize',15,'Units','normalized');
text(-0.1199,0.5281,'{\itc}_0','FontSize',15,'Units','normalized');
set(gca,'LineWidth',1.5,'FontSize',12,'XLim',lims(1,:),'YLim',lims(2,:),...
   'XTick',xtick,'YTick',ytick,'Position',[0.58,0.57,0.38,0.41]);
% k vs c1
subplot(2,2,3);
X = confounding_result.xdata;
Y = confounding_result.zVSx;
[x,y] = supplementary_function.createPatch2D(X,Y);
patch(x,y,zeros(size(y)),'FaceAlpha',1,'FaceColor','r','EdgeColor','red');
box on
text(0.4863,-0.1415,'\itk','FontSize',15,'Units','normalized');
text(-0.1747,0.5407,'{\itc}_1','FontSize',15,'Units','normalized');
set(gca,'LineWidth',1.5,'FontSize',12,'XLim',lims(1,:),'YLim',lims(3,:),...
   'XTick',xtick,'YTick',ztick,'Position',[0.1,0.08,0.38,0.41]);
% c0 vs c1
subplot(2,2,4);
X = confounding_result.ydata;
Y = confounding_result.zVSy;
[x,y] = supplementary_function.createPatch2D(X,Y);
patch(x,y,zeros(size(y)),'FaceAlpha',1,'FaceColor','r','EdgeColor','red');
box on
text(0.4863,-0.1203,'{\itc}_0','FontSize',15,'Units','normalized');
text(-0.1747,0.5407,'{\itc}_1','FontSize',15,'Units','normalized');
set(gca,'LineWidth',1.5,'FontSize',12,'XLim',lims(2,:),'YLim',lims(3,:),...
   'XTick',ytick,'YTick',ztick,'Position',[0.58,0.08,0.38,0.41]);