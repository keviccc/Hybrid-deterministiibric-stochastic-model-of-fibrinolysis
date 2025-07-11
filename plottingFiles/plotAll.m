% Plotting m file

t = linspace(startTime, finalTime, nt);

lysis = 1 - ploti.extentLysis;

% ========================================================================
% EVALUATING LYSIS
% ========================================================================
% Tasklist
% plot binding site change
% plot cross-section change

lysisLine = ones(1,nt)*N_0*(1-epsilon_0)*lysisLevel/N_CS;

fig1 = figure(1);
width = 12.5;
height = 6;
set(fig1, 'units','centimeters','pos',[0 0 width height])

subplot(1,2,1)
plot(t, ploti.crossSection,t,lysisLine,'-k')
xlabel('Time (s)');
ylabel('Number of Binding Sites');
title('A');

subplot(1,2,2)
plot(t,ploti.extentLysis,'-k',t,ploti.ntot/(n_0*(1-epsilon_0)),'--k')
xlabel('Time (s)');
ylabel('Fraction');
title('B');

% Saving file
saveas(fig1,'hybridTransLysisPlot.fig');

% ========================================================================
% EVALUATING CONCENTRATIONS
% ========================================================================

fig2 = figure(2);
set(fig2, 'units','centimeters','pos',[0 0 width height*2])

subplot(2,2,1)
plot(t,ploti.epsilon,'-k')
xlabel('Time (s)');
ylabel('Voidage');
title('A')

subplot(2,2,2)
plot(t,ploti.CtPA,'-k',t,ploti.ntPA,'--k')
xlabel('Time (s)');
ylabel('Concentration (\muM)');
title('B')

subplot(2,2,3)
plot(t,ploti.CPLG/50,'-k',t,ploti.nPLG,'--k')
xlabel('Time (s)');
ylabel('Concentration (\muM)');
title('C')

subplot(2,2,4)
plot(t,ploti.CPLS/20,'-k',t,ploti.nPLS,'--k',t,ploti.LPLS/10,'-.k')
xlabel('Time (s)');
ylabel('Concentration (\muM)');
title('D')

% Saving file
saveas(fig2,'hybridTransConcPlot.fig');

% ========================================================================
% PLOTTING LYSIS DISTRIBUTION FOR DIFFERENT EXTENTS OF LYSIS
% ========================================================================
fig3 = figure(3);
set(fig2, 'units','centimeters','pos',[0 0 width height*2])

xmax = N_CS/8;

subplot(2,3,1)
bar(ploti.FDPDist(:,find(lysis< 0.99,1)));
axis([0 xmax 0 20])
%xlabel('Length (Units)');
%ylabel('Number of Fibres');
title('1% Lysis')

subplot(2,3,2)
bar(ploti.FDPDist(:,find(lysis< 0.95,1)));
axis([0 xmax 0 20])
%xlabel('Length (Units)');
%ylabel('Number of Fibres');
title('5% Lysis')

subplot(2,3,3)
bar(ploti.FDPDist(:,find(lysis< 0.90,1)));
axis([0 xmax 0 20])
%xlabel('Length (Units)');
%ylabel('Number of Fibres');
title('10% Lysis')

subplot(2,3,4)
bar(ploti.FDPDist(:,find(lysis< 0.85,1)));
axis([0 xmax 0 20])
%xlabel('Length (Units)');
%ylabel('Number of Fibres');
title('15% Lysis')

subplot(2,3,5)
bar(ploti.FDPDist(:,find(lysis< 0.80,1)));
axis([0 xmax 0 20])
%xlabel('Length (Units)');
%ylabel('Number of Fibres');
title('20% Lysis')

subplot(2,3,6)
bar(ploti.FDPDist(:,find(lysis< 0.50,1)));
axis([0 xmax 0 20])
%xlabel('Length (Units)');
%ylabel('Number of Fibres');
title('50% Lysis')

% Saving file
saveas(fig3,'hybridTransFDPBar.fig');