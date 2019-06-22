function plot_B2Bposterior_delta(post_delta,envir_data)

pmin = post_delta.min;
pmax = post_delta.max;
tRange = post_delta.tRange;
count = 1;
while isempty(pmin{count})
   count = count+1;
end
ndata = length(pmin{count});
eps = envir_data.eps;
n_eps = length(eps);
k = envir_data.k;
b = envir_data.b;
lambda = 0.5*sqrt(4*k-b^2);
c1 = envir_data.c1;
c2 = envir_data.c2;
x1 = @(t) exp(-0.5*b*t).*(c1*sin(lambda*t)+c2*cos(lambda*t));
x2 = @(t) (envir_data.v0/sqrt(envir_data.k)*sin(sqrt(k)*t)+c2*cos(sqrt(k)*t));
delta = @(t) x1(t)-x2(t);
xx = linspace(tRange(1),tRange(2),ndata);
yy = delta(xx);
yl = {-0.7:0.5:0.3,-2:2:2};
f = figure('Name','Posterior \delta(t) over t \in [0, 4]');
set(f,'Units','inches','Position',[4 2 8 5]);
lp = [0.1 0.57];
bp = [0.55 0.1];
ww = 0.4;
hh = 0.4;
cc = 1;
pId = [1 3 2 4];
for i = 1:n_eps
   for j = 1:2
      subplot(2,2,pId(cc));
      plot(xx,yy,'r-','LineWidth',2);
      hold on
      plot(xx,pmin{i,j+2},'b-','LineWidth',2);
      plot(xx,pmax{i,j+2},'b-','LineWidth',2);
      plot([3 3],[-10 10],'k--','LineWidth',2);
      hold off
      set(gca,'FontSize',12,'LineWidth',1.5,'XLim',[0 4],'XTick',0:4,...
         'YLim',[yl{j}(1) yl{j}(end)],'YTick',yl{j},'Position',[lp(i),bp(j),ww,hh]);
      text(0.0371,0.0938,['{\itn} = ' num2str(j+2)],'Units','normalized','FontSize',15);
      box on;
      if j == 1
         set(gca,'XTickLabel',[])
         text(0.3795,1.0682,['\epsilon = ' num2str(eps(i),'%.2f')],'Units','normalized','FontSize',15);
      end
      if j == 2
         text(0.4909, -0.1761,'\itt','FontSize',15,'Units','normalized');
      end
      if cc == 1
         text(-0.1878,-0.7344,'Model discrepancy function \delta','FontSize',15,'Units','normalized','Rotation',90)
      end
      cc = cc+1;
   end
end