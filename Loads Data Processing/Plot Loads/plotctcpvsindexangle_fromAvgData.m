clear CT_data CTerr CP_data CPerr CTlo CTloerr CPlo CPloerr CTup CTuperr CPup CPuperr ctcp ctcperr col col_uni phis_uni
load('colors.mat')
AvgData_corr = AvgData;
for i = 1:length(AvgData_corr.avg_cps_inner)
    AvgData_corr.avg_cts_inner{i} = -AvgData_corr.avg_cts_inner{i};
    AvgData_corr.avg_cts_total{i} = (AvgData_corr.avg_cts_inner{i}+ AvgData_corr.avg_cts_outer{i})./2;
    
    AvgData_corr.avg_cps_inner{i} = -AvgData_corr.avg_cps_inner{i};
    AvgData_corr.avg_cps_total{i} = (AvgData_corr.avg_cps_inner{i}+ AvgData_corr.avg_cps_outer{i})./2;
end

RPM_des = 1250;
col_des = 8;
diffcol_des = 0;
seperate = false;
upcolor = colors{2};
locolor = colors{5};
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
phis_uni = phis_uni(phis_uni~=2);
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
    %     loc=loc&(err<0.001);
    
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


%% add -90 deg case
if (true)
    if sum(phis_uni==90)>0
        loc = (phis_uni == 90)';
        phis_uni(end+1) = -90;
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
[phis_plot,loc] = sort(phis_uni);
CT_data = CT_data(loc);
CTerr = CTerr(loc);
CP_data = CP_data(loc);
CPerr = CPerr(loc);
CTlo = CTlo(loc);
CTloerr = CTloerr(loc);
CPlo = CPlo(loc);
CPloerr = CPloerr(loc);
CTup = CTup(loc);
CTuperr = CTuperr(loc);
CPup = CPup(loc);
CPuperr = CPuperr(loc);
ctcperr = ctcperr(loc);



%%
% LOWER
figure(3)
hold on
CTCPerrlo = sqrt( (CTloerr./CPlo).^2 + (CPloerr .* CTlo./CPlo.^2).^2);
errorbar(phis_plot,CTlo./CPlo,CTCPerrlo, '-s','color',locolor,'MarkerEdgeColor',locolor,'MarkerFaceColor',locolor,'LineWidth', 1)
% plot([-95,95],[CT_90,CT_90], '--','color',[0 0 0]+0.7, 'linewidth',1.2)

% UPPER
hold on
CTCPerrup = sqrt( (CTuperr./CPup).^2 + (CPuperr .* CTup./CPup.^2).^2);
errorbar(phis_plot,CTup./CPup,CTCPerrup,'-^','color',upcolor,'MarkerEdgeColor',upcolor,'MarkerFaceColor',upcolor,'LineWidth', 1)
% plot([-95,95],[CT_90,CT_90], '--','color',[0 0 0]+0.7, 'linewidth',1.2)

% TOTAL
hold on
CTCPerr = sqrt( (CTerr./CP_data).^2 + (CPerr .* CT_data./CP_data.^2).^2);
errorbar(phis_plot,CT_data./CP_data,CTCPerr, '-o','color',totcolor,'MarkerEdgeColor',totcolor,'MarkerFaceColor',totcolor,'LineWidth', 1)
% plot([-95,95],[CT_90,CT_90], '--','color',[0 0 0]+0.7, 'linewidth',1.2)
xlabel('Index Angle, deg')
ylabel('C_T/C_P')
set(gca,'FontSize',18)
grid on
hold on
xlim([-90 90])
xticks([-90:15:90])
ylim([5, 15])
yticks([5:2:15])
legend('Lower', 'Upper','Total', 'location', 'northwest')


%%
function x = sumsquares(y)
for i=1: length(y)
    x(i) = (y(i))^2;
end
x = sqrt(sum(x));
end