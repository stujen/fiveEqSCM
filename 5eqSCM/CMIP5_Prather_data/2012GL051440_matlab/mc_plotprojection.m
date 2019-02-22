function mc_plotprojection( species, N )

param=readparams(species);

fname85 = sprintf( 'output/output_%s_%s_%d.mat', species, 'RCP85', N );
fname60 = sprintf( 'output/output_%s_%s_%d.mat', species, 'RCP60', N );
fname45 = sprintf( 'output/output_%s_%s_%d.mat', species, 'RCP45', N );
fname26 = sprintf( 'output/output_%s_%s_%d.mat', species, 'RCP26', N );

if (~exist(fname85,'file'))
    fprintf(1,'First run monte carlo for %s, %s, %d iterations\n',species,'RCP85',N)
    return
end
if (~exist(fname60,'file'))
    fprintf(1,'First run monte carlo for %s, %s, %d iterations\n',species,'RCP60',N)
    return
end
if (~exist(fname45,'file'))
    fprintf(1,'First run monte carlo for %s, %s, %d iterations\n',species,'RCP45',N)
    return
end
if (~exist(fname26,'file'))
    fprintf(1,'First run monte carlo for %s, %s, %d iterations\n',species,'RCP26',N)
    return
end

load(fname85,'Ercp')
load(fname85,'cF')
load(fname85,'tt')
E85=Ercp;
C85=cF;
t=tt;
load(fname60,'Ercp')
load(fname60,'cF')
E60=Ercp;
C60=cF;
load(fname45,'Ercp')
load(fname45,'cF')
E45=Ercp;
C45=cF;
load(fname26,'Ercp')
load(fname26,'cF')
E26=Ercp;
C26=cF;

if (mean(E85) < 1)
    unitA='ppt';
    unitE='Gg/a';
    E85 = E85*1e3;
    E60 = E60*1e3;
    E45 = E45*1e3;
    E26 = E26*1e3;
    C85 = C85*1e3;
    C60 = C60*1e3;
    C45 = C45*1e3;
    C26 = C26*1e3;
else
    unitA='ppb';
    unitE='Tg/a';

end

%FIGURE : SIMPLE EMISSIONS
figure(1); clf;
plot(t,mean(E26),'g-*')
hold on
plot(t,mean(E45),'b-o')
plot(t,mean(E60),'r-s')
plot(t,mean(E85),'k-+')
plot(t,param.Ercp26(t),'g')
plot(t,param.Ercp45(t),'b')
plot(t,param.Ercp60(t),'r')
plot(t,param.Ercp85(t),'k')
hold off
title('RCP emissions')
legend('RCP2.6&','RCP4.5&','RCP6.0&','RCP8.5&',...
    'RCP2.6','RCP4.5','RCP6.0','RCP8.5','Location','NorthWest')
ylabel(sprintf('anthro %s emissions, %s',species,unitE))

saveas(gcf,sprintf('figs/%s_emiss_mc%i',lower(species),N),'png')


%FIGURE : SIMPLE ABUNDANCE
figure(2); clf;
plot(t,mean(C26),'g-*')
hold on
plot(t,mean(C45),'b-o')
plot(t,mean(C60),'r-s')
plot(t,mean(C85),'k-+')
plot(t,median(C26),'g--');
plot(t,median(C45),'b--');
plot(t,median(C60),'r--');
plot(t,median(C85),'k--');
plot(t,param.Crcp26(t),'g')
plot(t,param.Crcp45(t),'b')
plot(t,param.Crcp60(t),'r')
plot(t,param.Crcp85(t),'k')

hold off
title('RCP projections')
legend('RCP2.6&','RCP4.5&','RCP6.0&','RCP8.5&',...
    'RCP2.6','RCP4.5','RCP6.0','RCP8.5','Location','NorthWest')
ylabel(sprintf('%s abundance, %s',species,unitA))

saveas(gcf,sprintf('figs/%s_conc_mc%i',lower(species),N),'png')


%FIGURE : EMISSIONS WITH UNCERTAINTY
figure(3); clf;
subplot(2,1,1)
h1=fill([t,fliplr(t)],[mean(E26)+std(E26),fliplr(mean(E26)-std(E26))],'g');
hold on
h2=fill([t,fliplr(t)],[mean(E45)+std(E45),fliplr(mean(E45)-std(E45))],'b');
h3=fill([t,fliplr(t)],[mean(E60)+std(E60),fliplr(mean(E60)-std(E60))],'r');
h4=fill([t,fliplr(t)],[mean(E85)+std(E85),fliplr(mean(E85)-std(E85))],'k');
h5=plot(t,mean(E26),'g-*');
h6=plot(t,mean(E45),'b-o');
h7=plot(t,mean(E60),'r-s');
h8=plot(t,mean(E85),'k-+');
h9=plot(t,param.Ercp26(t),'g','LineWidth',2);
h10=plot(t,param.Ercp45(t),'b','LineWidth',2);
h11=plot(t,param.Ercp60(t),'r','LineWidth',2);
h12=plot(t,param.Ercp85(t),'k','LineWidth',2);
hold off
set(h1,'EdgeColor','none')
set(h2,'EdgeColor','none')
set(h3,'EdgeColor','none')
set(h4,'EdgeColor','none')
%alpha(0.3)
title('RCP emissions')
legend([h5 h6 h7 h8 h9 h10 h11 h12],...
    'RCP2.6&','RCP4.5&','RCP6.0&','RCP8.5&',...
    'RCP2.6','RCP4.5','RCP6.0','RCP8.5','Location','NorthWest')
ylabel(sprintf('anthro %s emissions, %s',species,unitE))

%FIGURE : ABUNDANCES WITH UNCERTAINTY
subplot(2,1,2)
h1=fill([t,fliplr(t)],[mean(C26)+std(C26),fliplr(mean(C26)-std(C26))],'g');
hold on
h2=fill([t,fliplr(t)],[mean(C45)+std(C45),fliplr(mean(C45)-std(C45))],'b');
h3=fill([t,fliplr(t)],[mean(C60)+std(C60),fliplr(mean(C60)-std(C60))],'r');
h4=fill([t,fliplr(t)],[mean(C85)+std(C85),fliplr(mean(C85)-std(C85))],'k');
h5=plot(t,mean(C26),'g-*');
h6=plot(t,mean(C45),'b-o');
h7=plot(t,mean(C60),'r-s');
h8=plot(t,mean(C85),'k-+');
h9 =plot(t,param.Crcp26(t),'g','LineWidth',2);
h10=plot(t,param.Crcp45(t),'b','LineWidth',2);
h11=plot(t,param.Crcp60(t),'r','LineWidth',2);
h12=plot(t,param.Crcp85(t),'k','LineWidth',2);
hold off
set(h1,'EdgeColor','none')
set(h2,'EdgeColor','none')
set(h3,'EdgeColor','none')
set(h4,'EdgeColor','none')
%alpha(0.3)
title('RCP projections')
legend([h5 h6 h7 h8 h9 h10 h11 h12],...
    'RCP2.6&','RCP4.5&','RCP6.0&','RCP8.5&',...
    'RCP2.6','RCP4.5','RCP6.0','RCP8.5','Location','NorthWest')
ylabel(sprintf('%s abundance, %s',species,unitA))

% save as eps
margin = .25; ppos = [ margin, margin, 6,10];
set(gcf, 'PaperUnits','inches','PaperPosition', ppos);
print('-depsc2',sprintf('figs/%s_mc%i',lower(species),N));

% add transparency, then save as png
alpha(0.3)
subplot(2,1,1)
alpha(0.3)
saveas(gcf,sprintf('figs/%s_mc%i',lower(species),N),'png')


%FIGURE : ABUNDANCES PDF in 2100
%figure(1); clf;
%[f26,x26]=ksdensity(C26(:,end));
%[f45,x45]=ksdensity(C45(:,end));
%[f60,x60]=ksdensity(C60(:,end));
%[f85,x85]=ksdensity(C85(:,end));

%pdf
%plot(x26,f26,'g')
%hold on
%plot(x45,f45,'b')
%plot(x60,f60,'r')
%plot(x85,f85,'k')

% mean
%[x,i]=min(abs(x26-mean(C26(:,end)))); plot(x26(i),f26(i),'g.','MarkerSize',20)
%[x,i]=min(abs(x45-mean(C45(:,end)))); plot(x45(i),f45(i),'b.','MarkerSize',20)
%[x,i]=min(abs(x60-mean(C60(:,end)))); plot(x60(i),f60(i),'r.','MarkerSize',20)
%[x,i]=min(abs(x85-mean(C85(:,end)))); plot(x85(i),f85(i),'k.','MarkerSize',20)
%
% median
%[x,i]=min(abs(x26-median(C26(:,end)))); plot(x26(i),f26(i),'go')
%[x,i]=min(abs(x45-median(C45(:,end)))); plot(x45(i),f45(i),'bo')
%[x,i]=min(abs(x60-median(C60(:,end)))); plot(x60(i),f60(i),'ro')
%[x,i]=min(abs(x85-median(C85(:,end)))); plot(x85(i),f85(i),'ko')
%
%hold off
%title('Abundance y2100')
%ylabel('PDF')
%xlabel(sprintf('%s, %s',species,unitA))
%legend('RCP2.6&','RCP4.5&','RCP6.0&','RCP8.5&',...
%    'o=median','.=mean','Location','NorthEast')
%
%% save as eps, png
%margin = .25; ppos = [ margin, margin, 6,4];
%set(gcf, 'PaperUnits','inches','PaperPosition', ppos);
%print('-depsc2',sprintf('figs/%s_pdf_y2100_mc%i',lower(species),N));
%saveas(gcf,sprintf('figs/%s_pdf_y2100_mc%i',lower(species),N),'png')

end
