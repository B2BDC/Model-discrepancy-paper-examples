function plot_True_and_inadequate_models(envir_data)

k = envir_data.k;
b = envir_data.b;
lambda = envir_data.lambda;
c1 = envir_data.c1;
c2 = envir_data.c2;
z_true = @(t) exp(-0.5*b*t).*(c1*sin(lambda*t)+c2*cos(lambda*t));
z_inadequate = @(t) (c1*sin(sqrt(k)*t)+c2*cos(sqrt(k)*t));
delta = @(t) z_true(t)-z_inadequate(t);
f = figure('Name','z* and z and \delta*');
set(f,'Units','inches','Position',[2,1,8,4.5]);
xx = 0:0.01:4;
yyaxis left
hold on
plot(xx,z_true(xx),'r-','LineWidth',2);
plot(xx,z_inadequate(xx),'b-','LineWidth',2);
hold off
ylabel('Displacement','Position',[-0.2914 0.4948]);
text(2.7362,1.3994,'{\itz}^*({\itk}^*, {\itt})','FontSize',15,'Color','r')
text(2.1954,1.8540,'{\itz}({\itk}^*, {\itt})','FontSize',15,'Color','b')
set(gca,'FontSize',15,'LineWidth',1.5,'XLim',[0 4],'XTick',1:4,...
   'Position',[0.1 0.12 0.8 0.84],'YLim',[-2 3],'YTick',-2:2.5:3)
box on;
yyaxis right
plot(xx,delta(xx),'k-','LineWidth',2);
text(2.7134,-0.1455,'\delta^*({\itt})','FontSize',15)
lx = xlabel('\itt');
lx.Position = [2.00000190734863, -2.34435260640688];
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
ylabel('Model discrepancy','Position',[4.2428 -0.1000]);
set(gca,'FontSize',15,'LineWidth',1.5,'XLim',[0 4],'XTick',1:4,'YLim',[-0.25 0.05],'YTick',-0.25:0.15:0.05)
box on;