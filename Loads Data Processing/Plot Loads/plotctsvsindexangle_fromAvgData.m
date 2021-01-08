clear CT_data CTerr CP_data CPerr CTlo CTloerr CPlo CPloerr CTup CTuperr CPup CPuperr ctcp ctcperr col col_uni
load('colors.mat')
AvgData_corr = AvgData;

RPM_des = 1200;
col_des = 12;
diffcol_des = 0;
seperate = false;
upcolor = colors{5};
locolor = colors{3};
totcolor = colors{1};

% 4-BLADED ROTOR
% CT_90 = 0.0830;
% CP_90 = 0.0084;

%% GET CONSTANTS
diffcols = MeanData.diffcols;
phis = MeanData.phis;
col = MeanData.meancols;
RPMs = MeanData.RPMs;
phis_uni = unique(phis);
err = [AvgData_corr.err_cts_outer{:}]';

%% CHECK FOR OUTLIERS
% figure(12)
% loc = (col == col_des);
% errorbar(phis(loc),[AvgData_corr.avg_cts_total{loc}],[AvgData_corr.err_cts_total{loc}],'o')
% % errorbar(phis(loc),[AvgData_corr.avg_ctcp{loc}],[AvgData_corr.err_ctcp{loc}],'o')
% hold on
% loc = (col == col_des) &([AvgData_corr.err_cts_total{:}]'<.05);
% errorbar(phis(loc),[AvgData_corr.avg_cts_total{loc}],[AvgData_corr.err_cts_total{loc}],'k.')
% hold off

%% AVERAGE REPEATING DATA
CT_data=[];
CP_data=[];
CTlo=[];
CPlo=[];
CTup=[];
CPup=[];
CTloerr=[];
CPloerr=[];
CTuperr=[];
CPuperr=[];
CTerr=[];
CPerr=[];
ctcperr=[];
phis_plot=[];

for i = 1:length(phis_uni)
    loc =(RPMs> RPM_des*.98)&(RPMs < RPM_des*1.02);
    loc = (phis_uni(i) == phis)&loc & (col == col_des) & (diffcols == diffcol_des);
    loc=loc&(err<0.01);
    
    if sum(loc)>0
        phis_plot = [phis_plot,phis_uni(i)];
        CT_data = [CT_data, mean([AvgData_corr.avg_cts_total{loc}])];
        CP_data = [CP_data, mean([AvgData_corr.avg_cps_total{loc}])];
        CTlo = [CTlo, mean([AvgData_corr.avg_cts_inner{loc}])];
        CPlo = [CPlo, mean([AvgData_corr.avg_cps_inner{loc}])];
        CTup = [CTup, mean([AvgData_corr.avg_cts_outer{loc}])];
        CPup = [CPup, mean([AvgData_corr.avg_cps_outer{loc}])];
        if sum(loc)<4
            CTerr = [CTerr, sumsquares([AvgData_corr.err_cts_total{loc}])];
            CPerr = [CPerr, sumsquares([AvgData_corr.err_cps_total{loc}])];
            CTloerr = [CTloerr, sumsquares([AvgData_corr.err_cts_inner{loc}])];
            CPloerr = [CPloerr, sumsquares([AvgData_corr.err_cps_inner{loc}])];
            CTuperr = [CTuperr, sumsquares([AvgData_corr.err_cts_outer{loc}])];
            CPuperr = [CPuperr, sumsquares([AvgData_corr.err_cps_outer{loc}])];
        else
            Nrevs = mean([SortedData.nrevs{loc}]);
            CTerr = [CTerr,1.96* std([AvgData_corr.avg_cts_total{loc}])/sum(loc)];
            CPerr = [CPerr,1.96* std([AvgData_corr.avg_cps_total{loc}])/sum(loc)];
            CTloerr = [CTloerr,1.96* std([AvgData_corr.avg_cts_inner{loc}])/sum(loc)];
            CPloerr = [CPloerr,1.96* std([AvgData_corr.avg_cps_inner{loc}])/sum(loc)];
            CTuperr = [CTuperr,1.96* std([AvgData_corr.avg_cts_outer{loc}])/sum(loc)];
            CPuperr = [CPuperr,1.96* std([AvgData_corr.avg_cps_outer{loc}])/sum(loc)];
        end
        
%         loc = (col == col_des) & (phis_uni(i) == phis) & ([AvgData_corr.err_ctcp{:}]'<0.25);
        ctcperr = [ctcperr,sumsquares([AvgData_corr.err_ctcp{loc}])];
    end
end
% phis_plot = [phis_plot, -phis_plot];
% CT_data = [CT_data, CT_data];
% CP_data = [CP_data, CP_data];
% CTlo = [CTlo, CTlo];
% CPlo = [CPlo, CPlo];
% CTup = [CTup, CTup];
% CPup = [CPup, CPup];
% CTerr = [CTerr, CTerr];
% CPerr = [CPerr, CPerr];
% CTloerr = [CTloerr, CTloerr];
% CPloerr = [CPloerr, CPloerr];
% CTuperr = [CTuperr, CTuperr];
% CPuperr = [CPuperr, CPuperr];
% 
% [phis_plot,I] = sort(phis_plot);
% CT_data = CT_data(I);
% CP_data = CP_data(I);
% CTlo = CTlo(I);
% CPlo = CPlo(I);
% CTup = CTup(I);
% CPup = CPup(I);
% CTerr = CTerr(I);
% CPerr = CPerr(I);
% CTloerr = CTloerr(I);
% CPloerr = CPloerr(I);
% CTuperr = CTuperr(I);
% CPuperr = CPuperr(I);


% add -90 deg case
if (false)
    if sum(phis_uni==90)>0
        phis_uni(end+1) = -90;
        loc = (phis_uni == 90);
        CT_data(end+1) = CT_data(loc);
        CTerr(end+1) = CTerr(loc);
        CP_data(end+1) = CP_data(loc);
        CPerr(end+1) = CPerr(loc);
        
        CTlo(end+1) = CTlo(loc);
        CTloerr(end+1) = CTloerr(loc);
        CPlo(end+1) = CPlo(loc);
        CPloerr(end+1) = CPloerr(loc);
        
        CTup(end+1) = CTup(loc);
        CTuperr(end+1) = CTuperr(loc);
        CPup(end+1) = CPup(loc);
        CPuperr(end+1) = CPuperr(loc);
        
        ctcperr(end+1) = ctcperr(loc);
    end
end

%%
if seperate
    % LOWER
    figure(1)
    subplot(2,1,2)
    hold on
    errorbar(phis_plot,CTlo,CTloerr, 'o','color',locolor,'MarkerEdgeColor',locolor,'MarkerFaceColor',locolor,'LineWidth', 1)
    % plot([-95,95],[CT_90,CT_90], '--','color',[0 0 0]+0.7, 'linewidth',1.2)
    ylim([0.05, 0.17])
    xlim([-95 95])
    xticks([-90:15:90])
    yticks([0.05:0.02:0.17])
    ylabel('C_T/ \sigma')
    xlabel('Index Angle, deg')
    grid on
    grid minor
    set(gca, 'fontsize',14)
    
    % UPPER
    subplot(2,1,1)
    hold on
    errorbar(phis_plot,CTup,CTuperr,'o','color',upcolor,'MarkerEdgeColor',upcolor,'MarkerFaceColor',upcolor,'LineWidth', 1)
    % plot([-95,95],[CT_90,CT_90], '--','color',[0 0 0]+0.7, 'linewidth',1.2)
    ylim([0.05, 0.17])
    xlim([-95 95])
    xticks([-90:15:90])
    yticks([0.05:0.02:0.17])
    xlabel('')
    ylabel('C_T/ \sigma')
    grid on
    set(gca, 'fontsize',14)
    legend('CFD','VVPM','Exp','location',[.85 .88 .1 .1])
    
    % TOTAL
    figure(2)
    hold on
    errorbar(phis_plot,CT_data,CTerr, 'o','color',totcolor,'MarkerEdgeColor',totcolor,'MarkerFaceColor',totcolor,'LineWidth', 1)
    % plot([-95,95],[CT_90,CT_90], '--','color',[0 0 0]+0.7, 'linewidth',1.2)
    xlabel('Index Angle, deg')
    ylabel('C_T/ \sigma')
    set(gca,'FontSize',18)
    grid on
    hold on
    ylim([0.05, 0.17])
    xlim([-95 95])
    xticks([-90:15:90])
    yticks([0.05:0.02:0.17])
    legend('CFD','VVPM','Exp','location',[.85 .86 .1 .1])
    
    % CT/CP
    figure(3)
    hold on
    errorbar(phis_plot,CT_data./CP_data,ctcperr, 'o','color',totcolor,'MarkerEdgeColor',totcolor,'MarkerFaceColor',totcolor,'LineWidth', 1)
    % plot([-95,95],[CT_90./CP_90,CT_90./CP_90], '--','color',[0 0 0]+0.7, 'linewidth',1.2)
    xlabel('Index Angle, deg')
    ylabel('C_T/ C_P')
    set(gca,'FontSize',18)
    grid on
    hold on
    ylim([5 12])
    xlim([-95 95])
    xticks([-90:15:90])
    % yticks([0.05:0.02:0.17])
    legend('CFD','VVPM','Exp','location',[.88 .88 .1 .1])
    %%
else
    % LOWER
    figure(1)
    hold on
    errorbar(phis_plot,CTlo,CTloerr, 's','color',locolor,'MarkerEdgeColor',locolor,'MarkerFaceColor',locolor,'LineWidth', 1)
    % plot([-95,95],[CT_90,CT_90], '--','color',[0 0 0]+0.7, 'linewidth',1.2)
    
    % UPPER
    hold on
    errorbar(phis_plot,CTup,CTuperr,'^','color',upcolor,'MarkerEdgeColor',upcolor,'MarkerFaceColor',upcolor,'LineWidth', 1)
    % plot([-95,95],[CT_90,CT_90], '--','color',[0 0 0]+0.7, 'linewidth',1.2)
    
    % TOTAL
    hold on
    errorbar(phis_plot,CT_data,CTerr, 'o','color',totcolor,'MarkerEdgeColor',totcolor,'MarkerFaceColor',totcolor,'LineWidth', 1)
    % plot([-95,95],[CT_90,CT_90], '--','color',[0 0 0]+0.7, 'linewidth',1.2)
    xlabel('Index Angle, deg')
    ylabel('C_T/ \sigma')
    set(gca,'FontSize',18)
    grid on
    hold on
    xlim([-95 95])
    xticks([-90:15:90])
    ylim([0.04, 0.14])
    %     yticks([0.05:0.02:0.17])
    legend('Lower', 'Upper','Total', 'location', 'northwest')
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % LOWER
    figure(2)
    hold on
    errorbar(phis_plot,CPlo,CPloerr, 's','color',locolor,'MarkerEdgeColor',locolor,'MarkerFaceColor',locolor,'LineWidth', 1)
    % plot([-95,95],[CT_90,CT_90], '--','color',[0 0 0]+0.7, 'linewidth',1.2)
    
    % UPPER
    hold on
    errorbar(phis_plot,CPup,CPuperr,'^','color',upcolor,'MarkerEdgeColor',upcolor,'MarkerFaceColor',upcolor,'LineWidth', 1)
    % plot([-95,95],[CT_90,CT_90], '--','color',[0 0 0]+0.7, 'linewidth',1.2)
    
    % TOTAL
    hold on
    errorbar(phis_plot,CP_data,CPerr, 'o','color',totcolor,'MarkerEdgeColor',totcolor,'MarkerFaceColor',totcolor,'LineWidth', 1)
    % plot([-95,95],[CT_90,CT_90], '--','color',[0 0 0]+0.7, 'linewidth',1.2)
    xlabel('Index Angle, deg')
    ylabel('C_P/ \sigma')
    set(gca,'FontSize',18)
    grid on
    hold on
    xlim([-95 95])
    xticks([-90:15:90])
    %     ylim([0.05, 0.17])
    %     yticks([0.05:0.02:0.17])
    legend('Lower', 'Upper','Total', 'location', 'northwest')
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % LOWER
    figure(3)
    hold on
    plot(phis_plot,CTlo./CPlo, '-s','color',locolor,'MarkerEdgeColor',locolor,'MarkerFaceColor',locolor,'LineWidth', 1)
    % plot([-95,95],[CT_90,CT_90], '--','color',[0 0 0]+0.7, 'linewidth',1.2)
    
    % UPPER
    hold on
    plot(phis_plot,CTup./CPup,'-^','color',upcolor,'MarkerEdgeColor',upcolor,'MarkerFaceColor',upcolor,'LineWidth', 1)
    % plot([-95,95],[CT_90,CT_90], '--','color',[0 0 0]+0.7, 'linewidth',1.2)
    
    % TOTAL
    hold on
    plot(phis_plot,CT_data./CP_data, '-o','color',totcolor,'MarkerEdgeColor',totcolor,'MarkerFaceColor',totcolor,'LineWidth', 1)
    % plot([-95,95],[CT_90,CT_90], '--','color',[0 0 0]+0.7, 'linewidth',1.2)
    xlabel('Index Angle, deg')
    ylabel('C_T/C_P')
    set(gca,'FontSize',18)
    grid on
    hold on
    xlim([-95 95])
    xticks([-90:15:90])
    ylim([5, 15])
    yticks([5:2:15])
    legend('Lower', 'Upper','Total', 'location', 'northwest')
end

%%
function x = sumsquares(y)
for i=1: length(y)
    x(i) = (y(i))^2;
end
x = sqrt(sum(x));
end