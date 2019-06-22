function plot_B2Bprediction_with_true_model(zPred,tpred,zTrue)

xLab = {'0.10','0.05'};
xx = [1, 2];
xl = [0.5 2.5];
yl = {0.19:0.02:0.23, 1.81:0.01:1.83, 2.13:0.02:2.17};
f = figure('Name','Prediction with M*');
set(f,'Units','inches','Position',[4 2 8 3.2]);
bp = repmat(0.15,1,3);
lp = [0.09 0.41 0.72];
hh = 0.8;
ww = 0.25;
for i = 1:3
   yy = zPred(:,:,i);
   subplot(1,3,i);
   hold on;
   plot(xl,repmat(zTrue(i),1,2),'r--','LineWidth',2);
   errorbar(xx,mean(yy,2),0.5*diff(yy,[],2),'b.','LineWidth',2,'MarkerSize',1e-5);
   hold off
   if i == 1
      text(-0.2854,0.2707,'Predicted  {\itz}^*','Units','normalized','FontSize',15,'Rotation',90);
   end
   text(0.4698, -0.0967,'\epsilon','Units','normalized','FontSize',15);
   set(gca,'FontSize',12,'LineWidth',1.5,'XLim',xl,'XTick',xx,'XTickLabel',xLab,...
      'YLim',[yl{i}(1) yl{i}(end)],'YTick',yl{i},'Position',[lp(i),bp(i),ww,hh]);
   box on;
   text(0.65,0.9395,['{\itt} = ' num2str(tpred(i),'%2.1f')],'FontSize',15,'Units','normalized');
end