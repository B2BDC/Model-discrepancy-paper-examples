function plot_prior_constrain_on_delta(pred_with_prior_constraint,pred_no_prior,envir_data,ntest)

zPred = pred_with_prior_constraint.prediction;
alpha_opt = pred_with_prior_constraint.minimal_alpha;
alpha = pred_with_prior_constraint.maximal_alpha;
tpred = envir_data.t_pred;
zTrue = envir_data.z_pred_true;
zInf = permute(pred_no_prior(1,:,:),[3,2,1]);
xl = {[0:0.004:0.02 0.024 0.028 0.032],[0:0.01:0.05 0.06 0.07 0.08],[0:0.02:0.1 0.12 0.14 0.16]};
yl = {0.13:0.07:0.27, 1.4:0.5:2.4, 0.6:1.6:3.8};
f = figure('Name','Prediction with prior constraint on \delta');
set(f,'Units','inches','Position',[4 2 8 6]);
lp = 0.12;
bp = [0.7 0.4 0.1];
ww = 0.84;
hh = 0.25;
for i = 1:3
   yy = zPred(:,:,i);
   xx = linspace(alpha_opt(i),alpha(i),ntest+1);
   xx = xx(2:end);
   subplot(3,1,i);
   hold on;
   plot([xl{i}(1) xl{i}(end)],repmat(zTrue(i),1,2),'r--','LineWidth',2);
   errorbar(xx,mean(yy'),0.5*diff(yy,[],2),'b.',...
      'LineWidth',2,'MarkerSize',1e-5);
   errorbar(xl{i}(end-1),mean(zInf(i,:)),0.5*diff(zInf(i,:)),'b.',...
      'LineWidth',2,'MarkerSize',1e-5);
   hold off
   if i == 3
      text(0.4922, -0.2611,'\alpha','FontSize',15,'Units','normalized');
   end
   if i == 2
      text(-0.0969, -0.2306,'Predicted  {\itz}({\itk}, {\itt})  +  \delta({\itt})','FontSize',15,'Units','normalized','Rotation',90);
   end
   set(gca,'FontSize',12,'LineWidth',1.5,'XLim',[xl{i}(1), xl{i}(end)],'XTick',xl{i}(1:end-1),...
      'XTickLabel',[num2cell(xl{i}(1:end-3)), {'...','\infty'}],...
      'YLim',[yl{i}(1) yl{i}(end)],'YTick',yl{i},'Position',[lp,bp(i),ww,hh]);
   box on;
   text(0.0114,0.8611,['{\itt} = ' num2str(tpred(i),'%2.1f')],'Units','normalized','FontSize',15);
end
