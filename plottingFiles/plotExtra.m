fig1 = figure(1);
plot(R_f_array,final_lysis_time,'sk','linewidth',1.5);
set(fig1,  'units', 'centimeters','pos',[5, 5, 8.44, 5.908]);
set(gca,'fontsize',10)
xlhand = get(gca,'xlabel');
ylhand = get(gca,'ylabel');
set(xlhand,'string','Radius (nm)','fontsize',10);
set(ylhand,'string','Time till 95% Lysis (s)','fontsize',10);
%print(strcat(fname,'minAP'),'-depsc','-tiff','-r600')

fig2 = figure(2);
plot(t,lysis_curve_store,'linewidth',1.5)
set(fig2,  'units', 'centimeters','pos',[5, 5, 8.44, 5.908]);
set(gca,'fontsize',10)
xlhand = get(gca,'xlabel');
ylhand = get(gca,'ylabel');
set(xlhand,'string','Time (s)','fontsize',10);
set(ylhand,'string','Extent of Lysis','fontsize',10);

fig3 = figure(3);
plot(R_f_array,N_CS_array,'sk','linewidth',1.5)
set(fig3,  'units', 'centimeters','pos',[5, 5, 8.44, 5.908]);
set(gca,'fontsize',10)
xlhand = get(gca,'xlabel');
ylhand = get(gca,'ylabel');
set(xlhand,'string','Radius (nm)','fontsize',10);
set(ylhand,'string','Number of Cross-Sections','fontsize',10);

fig4 = figure(4);
plot(t,BS_store,'linewidth',1.5)
set(fig4,  'units', 'centimeters','pos',[5, 5, 8.44, 5.908]);
set(gca,'fontsize',10)
xlhand = get(gca,'xlabel');
ylhand = get(gca,'ylabel');
set(xlhand,'string','Time (s)','fontsize',10);
set(ylhand,'string','Binding Site Concentration (\muM)','fontsize',10);
    
    
    