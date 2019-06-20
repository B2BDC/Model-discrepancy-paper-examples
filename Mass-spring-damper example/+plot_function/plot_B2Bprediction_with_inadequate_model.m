function plot_B2Bprediction_with_inadequate_model(zPred,tPred,zTrue)

eplab = {'0.10','0.05'};
ww = 0.38;
hh = 0.27;
lp = [0.15 0.15+0.44];
bp = [0.95-hh, 0.95-2*hh-0.02, 0.1];
f = figure('Name','Prediction with inadequate model');
set(f,'Units','inches','Position',[2 0.2 8 11]);
n_poly = length(zPred);
n_pred = length(tPred);
xx = 1:n_poly;
xl = [0.5 n_poly+0.5];
yl = {0.18:0.06:0.3; 1.4:0.5:2.4; 0.6:1.5:3.6;};
cc = 1;
for i = 1:n_pred
   yPred = zeros(n_poly,2,2);
   for j = 1:n_poly
      yPred(j,:,1) = zPred{j}(1,:,i);
      yPred(j,:,2) = zPred{j}(2,:,i);
   end
   for j = 1:2
      yy = yPred(:,:,j);
      infflag = isinf(yy(:,1));
      subplot(3,2,cc);
      hold on;
      plot(xl,repmat(zTrue(i),1,2),'r--','LineWidth',2);
      errorbar(xx(~infflag),mean(yy(~infflag,:)'),0.5*diff(yy(~infflag,:),[],2),'b.',...
         'LineWidth',2,'MarkerSize',1e-5);
      plot(xx(infflag),repmat(zTrue(i),1,sum(infflag)),'xk','LineWidth',2,'MarkerSize',11);
      hold off
      set(gca,'FontSize',12,'LineWidth',1.5,'XLim',xl,'XTick',[],...
         'YLim',[yl{i}(1) yl{i}(end)],'YTick',yl{i},'Position',[lp(j),bp(i),ww,hh]);
      box on;
      if i == 1
         text(0.3732,1.0694,['\epsilon = ' eplab{j}],'FontSize',15,'Units','normalized');
      elseif i == 3
         text(0.4829,-0.1640,'{\itn}','FontSize',15,'Units','normalized');
         set(gca,'XTick',xx,'XTickLabel',xx-1);
      end
      if j == 1
         text(-0.1782,0.3562,['at {\itt} = ' num2str(tPred(i),'%.1f')],'FontSize',15,'Units','normalized','Rotation',90);
         if i == 2
            text(-0.3072,0.1283,'Predicted  {\itz}({\itk}, {\itt}) + \delta({\itt})','FontSize',15,'Units','normalized','Rotation',90);
         end
      end
      cc = cc+1;
   end
end