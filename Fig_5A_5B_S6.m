% fix the random seed for reproducibility
rndseed = 1; rng(rndseed);

specific_file = true; % set to true to only use one data file
specific_file_name = "dose_response_gefit_osim_analysis";

load(specific_file_name);

specific_drugs = true; % set to true to pick a list of drugs
normalize_by_dmso_sta = true; % normalize with dmso / sta or max and min
specific_plates = true; % use only specific plates  

figure(1)
tiledlayout(3, 2);
figure(2)
tiledlayout(3, 2);
figure(3)
tiledlayout(2, 3);
figure(4)
tiledlayout(2, 3);
figure(5)
figure(6)
figure(7)
figure(8)

% create variables storing results
classno = zeros(8,18);
class2no = zeros(18,1);
labels = strings(18,1);

growth_scores = zeros(2,6,8);
viability_scores = zeros(2,6,8);
EC50s_alive = zeros(2,6);
EC50s_growth = zeros(2,6);
EC50s_ctg = zeros(2,6);
CCTRs = zeros(2,6);
growth_scores_gefit = zeros(8,6);
growth_scores_osim = zeros(8,6);
viability_scores_gefit = zeros(8,6);
viability_scores_osim = zeros(8,6);
viability_score_vec = zeros(size(drug_numbers,1),1);

% constraints for Hill curve fitting
lower_lims = [0, 0, 0, 1];
upper_lims = [Inf, 5, 1, 1];

count = 0;
count_drug = 0;

% colors for phenotypic categories
cmap = [
255/255, 128/255, 114/255;
255/255,211/255,44/255;
163/255,121/255,217/255;
%1 1 1;
207/255,235/255,241/255;
100/255,149/255,237/255;
0 0 200/255];

opts = optimset('Display','off');

for j = [4,6]
    count_drug = count_drug + 1;
    count_organoid = 0;

    % load CTG data
    load("ctg_data.mat");
    if j == 4
        ctg_data = gefitinib_ctg;
    else
        ctg_data = osimertinib_ctg;
    end

    for i = 1:6
        count_organoid = count_organoid+1;
        count = count+1;
    
        figure(count_drug);
        nexttile;
        a = find(drug_numbers == j & organoid_values == i); % find drug - organoid pair
        specific_drugs_list = [1,3,j]; % drugs to display
        specific_plates_list_single_file = [plate_numbers(a,1)]; % plates to use
        % generate GV plot for drug - organoid pair
        [perc_flokkar, growth_score, doses_part] = plot_script_drug_response_function(specific_file, specific_file_name, specific_drugs, specific_drugs_list, ...
        normalize_by_dmso_sta, specific_plates, specific_plates_list_single_file, [], plate_numbers, dose_values);
        title(append(organoids(i)," ",drugs(j)));

        % make sure doses are in ascending order
        if ~issorted(doses_part(3:end))
            error('Doses not in ascending order')
            quit
        end

        % compute and store viability scores for drug - organoid pair;
        vec = (perc_flokkar(3:end,1)-perc_flokkar(2,1))/(perc_flokkar(1,1)-perc_flokkar(2,1));
        viability_scores(count_drug,count_organoid,:) = vec;
        viability_score_vec(a) = vec;
        % store growth scores for drug - organoid pair;
        growth_scores(count_drug,count_organoid,:) = growth_score(3:end);
    
        % generate DRC for viability score
        opt = Inf;
        hill_fit = @(b,x)  b(4)+(b(3)-b(4)).*x.^b(2)./(b(1).^b(2)+x.^b(2));
        b0 = [1,1,1,1];                           
        for k = 1:1000
            [B_temp, val_temp] = lsqcurvefit(hill_fit, b0, doses_part(3:end), vec, lower_lims, upper_lims, opts);
            if val_temp < opt
                opt = val_temp;
                B = B_temp;
            end
        end

        % compute IC50 dose
        vec_check = [10^(-4):10^(-5):10^3];
        for ell = 1:size(vec_check,2)
            x = vec_check(ell);
            if hill_fit(B,x) <= 0.5
                break;
            end
        end
        EC50s_alive(count_drug,count_organoid) = x;

        % plot DRC for viability score
        figure(2+count_drug);
        nexttile;
        scatter(doses_part(3:end),vec,30,'','red',"^");
        set(gca,'ylim',[-0.5 1.5])
        set(gca,'XScale','log');
        hold on
        fplot(@(x) hill_fit(B,x),'Color','red','LineWidth',2.5)
        set(gca,'color', [0.93 0.93 0.93]);
        
        yline(0,':','Color','black','LineWidth',1);
        yline(0.5,':','Color','black','LineWidth',1);
        yline(1,':','Color','black','LineWidth',1);
        xline(x,':','Color','red','LineWidth',2);
    
        % generate DRC for growth score
        opt = Inf;
        hill_fit = @(b,x)  b(4)+(b(3)-b(4)).*x.^b(2)./(b(1).^b(2)+x.^b(2));
        b0 = [1,1,1,1];
        for k = 1:1000
            [B_temp, val_temp] = lsqcurvefit(hill_fit, b0, [doses_part(3:end)], [growth_score(3:end)], lower_lims, upper_lims, opts);
            if val_temp < opt
                opt = val_temp;
                C = B_temp;
            end
        end

        % compute IC50 dose
        vec_check = [10^(-4):10^(-5):10^3];
        for ell = 1:size(vec_check,2)
            x = vec_check(ell);
            if hill_fit(C,x) <= 0.5
                break;
            end
        end
        EC50s_growth(count_drug,count_organoid) = x;

        % plot DRC for growth score
        scatter(doses_part(3:end),growth_score(3:end),50,'','blue',"square");
        fplot(@(x) hill_fit(C,x),'Color','blue','LineWidth',2.5);
        xline(x,':','Color','blue','LineWidth',2);
        set(gca,'xlim',[0.002 70]);
        xticks([10^(-2),10^(-1),10^0,10^1]);

        % calculate CCTRs
        CCTRs(count_drug,count_organoid) = log10(EC50s_alive(count_drug,count_organoid)) - log10(EC50s_growth(count_drug,count_organoid));

        % generate DRC for CTG data
        ctg_data_here = ctg_data(8*(count_organoid-1)+1:8*(count_organoid-1)+8,:);
        doses_here = [doses_part(3:end), doses_part(3:end), doses_part(3:end)];
        empty = find(ctg_data_here == 0);
        ctg_data_here(empty) = [];
        doses_here(empty) = [];

        opt = Inf;
        hill_fit = @(b,x)  b(4)+(b(3)-b(4)).*x.^b(2)./(b(1).^b(2)+x.^b(2));
        b0 = [1,1,1,1];
        for k = 1:1000
            [B_temp, val_temp] = lsqcurvefit(hill_fit, b0, doses_here, ctg_data_here, lower_lims, upper_lims, opts);
            if val_temp < opt
                opt = val_temp;
                D = B_temp;
            end
        end

        % compute IC50 dose
        vec_check = [10^(-4):10^(-5):10^3];
        for ell = 1:size(vec_check,2)
            x = vec_check(ell);
            if hill_fit(D,x) <= 0.5
                break;
            end
        end
        EC50s_ctg(count_drug,count_organoid) = x;

        % plot DRC for CTG data
        means_ctg = zeros(8,1);
        for k = 1:8
            ind = find(doses_here == doses_part(2+k));
            means_ctg(k) = mean(ctg_data_here(ind));
        end
        scatter(doses_part(3:end),means_ctg,30,[0,111,60]/255,"filled");
        fplot(@(x) hill_fit(D,x),'Color',[0,111,60]/255,'LineWidth',2.5);
        xline(x,':','Color',[0,111,60]/255,'LineWidth',2);

        set(gca,'FontSize',12)
        xlabel('\bf Dose (μM)')
        ylabel('\bf Normalized RLU and GV score')
        title(append(organoids(i)," ",drugs(j)),'FontSize',18);

       % compute phenotypic categories across doses
       livelive = perc_flokkar(:, 1);
       livedead = perc_flokkar(:, 2);
       deaddead = perc_flokkar(:, 3);
       cytotoxic_bool = deaddead >= 0.7;
       late_cytotoxic_bool = (livedead >= 0.8 - 0.866667.*livelive) & (livedead >= 0.15);
       cytostatic_bool = livelive >= 0.75 & livedead <= 0.15;
       cytotoxic_cytostatic_bool = ~(cytotoxic_bool | late_cytotoxic_bool | cytostatic_bool);

       bool_vec = [cytotoxic_bool, late_cytotoxic_bool, cytostatic_bool, cytotoxic_cytostatic_bool];
       for k=1:8
            if bool_vec(k+2,:) == [true, false, false, false]
                classno(k+2,count) = 1;
          elseif bool_vec(k+2,:) == [false, true, false, false]
             classno(k+2,count) = 3; 
          elseif bool_vec(k+2,:) == [false, false, true, false]
             if growth_score(k+2) >= 0.7
                classno(k+2,count) = 4;
             elseif growth_score(k+2) >= 0.3
                 classno(k+2,count) = 5;
             else
                 classno(k+2,count) = 6;
             end
          elseif bool_vec(k+2,:) == [false, false, false, true]
             classno(k+2,count) = 2;
            end
       end

       labels(count) = append(organoids(i)," ",drugs(j));
    end

    % plot phenotypic categories
    figure(4+count_drug);
    imagesc(classno(3:end,6*(count_drug-1)+1:6*(count_drug-1)+6)','CDataMapping','direct');
    colormap(cmap);
    xticklabels(doses_part(3:10))
    yticks([1:18]);
    yticklabels(labels(6*(count_drug-1)+1:6*(count_drug-1)+6));
    xlabel('Doses');
    title('Phenotypic categories across doses');
    colorbar('Ticks',[1:1:7],'TickLabels',["Cytotoxic" "Cytost.+cytot." "Late cytot." "Cytost. low to no effect" "Cytost. medium" "Cytost. strong"]);
    xline([0.5:1:8.5]);
    yline([0.5:1:16.5]);

    % create table with growth and viability scores
    fprintf(append('Growth scores for ',drugs(j),':'));
    table(squeeze(growth_scores(count_drug,1,:)),squeeze(growth_scores(count_drug,2,:)),squeeze(growth_scores(count_drug,3,:)),squeeze(growth_scores(count_drug,4,:)),squeeze(growth_scores(count_drug,5,:)),squeeze(growth_scores(count_drug,6,:)),'VariableNames',organoids)
    if count_drug == 1
        growth_scores_gefit = squeeze(growth_scores(count_drug,:,:))';
    else
        growth_scores_osim = squeeze(growth_scores(count_drug,:,:))';
    end
    fprintf(append('Viability scores for ',drugs(j),':'));
    table(squeeze(viability_scores(count_drug,1,:)),squeeze(viability_scores(count_drug,2,:)),squeeze(viability_scores(count_drug,3,:)),squeeze(viability_scores(count_drug,4,:)),squeeze(viability_scores(count_drug,5,:)),squeeze(viability_scores(count_drug,6,:)),'VariableNames',organoids)
    if count_drug == 1
        growth_scores_gefit = squeeze(growth_scores(count_drug,:,:))';
    else
        growth_scores_osim = squeeze(growth_scores(count_drug,:,:))';
    end
end

% plot CCTRs
figure(7);
b = bar(organoids,CCTRs(1,:),0.6);
set(gca,'ylim',[0 1.5])
b.FaceColor = 'flat';
b.CData = [0.2 0.2 0.2; 0.5 0.5 0.5; 0.4 0.4 0.4; 0.8 0.8 0.8; 0.4 0.4 0.4; 0.7 0.7 0.7];
ylabel({'CCTR';'log_{10}(IC_{50} viability) − log_{10}(IC_{50} growth)'})
set(gca,'FontSize',14)
title('Gefitinib','FontSize',20)
figure(8);
b = bar(organoids,CCTRs(2,:),0.6);
set(gca,'ylim',[0 1.5])
b.FaceColor = 'flat';
b.CData = [0.2 0.2 0.2; 0.5 0.5 0.5; 0.4 0.4 0.4; 0.8 0.8 0.8; 0.4 0.4 0.4; 0.7 0.7 0.7];
ylabel({'CCTR';'log_{10}(IC_{50} viability) − log_{10}(IC_{50} growth)'})
set(gca,'FontSize',14)
title('Osimertinib','FontSize',20)

% calculate confidence intervals
limits = zeros(size(drug_numbers,1),2);
for i = 1:size(drug_numbers,1)
    if drug_numbers(i) > 3
        a = i;
        b = find(drug_numbers == 1 & plate_numbers == plate_numbers(a));
        c = find(drug_numbers == 3 & plate_numbers == plate_numbers(a));
        idx = [a,b,c];
        growth_score_vec = zeros(10000,1);
        for j = 1:size(growth_score_vec,1)
            sample = randi(no_organoids(a),no_organoids(a),1);
            params_here = params;
            params_here(1:no_organoids(a),:,a) = params(sample,:,a);
            params_2_here = params_2;
            params_2_here(1:no_organoids(a),:,a) = params_2(sample,:,a);
            error_vec_bic_here = error_vec_bic;
            error_vec_bic_here(1:no_organoids(a),a)  = error_vec_bic(sample,a);
            error_vec_2_bic_here = error_vec_2_bic;
            error_vec_2_bic_here(1:no_organoids(a),a) = error_vec_2_bic(sample,a);
            [growth_score, growth_score_ab, growth_score_a] = calculate_growth_score_dose_response(no_organoids(idx), params_here(:,:,idx), params_2_here(:,:,idx), error_vec_bic_here(:,idx), error_vec_2_bic_here(:,idx), normalize_by_dmso_sta, drug_numbers(idx), plate_numbers(idx));
            growth_score_vec(j) = growth_score(1);
        end
        limits(i,:) = [prctile(growth_score_vec,2.5), prctile(growth_score_vec,97.5)];
    end
end

% store all growth and viability scores
growth_score_vec = zeros(size(drug_numbers,1),1);
for i = 1:size(drug_numbers,1)
    if drug_numbers(i) > 3
        a = i;
        b = find(drug_numbers == 1 & plate_numbers == plate_numbers(a));
        c = find(drug_numbers == 3 & plate_numbers == plate_numbers(a));
        idx = [a,b,c];
        [growth_score_here, growth_score_ab_here, growth_score_a_here] = calculate_growth_score_dose_response(no_organoids(idx), params(:,:,idx), params_2(:,:,idx), error_vec_bic(:,idx), error_vec_2_bic(:,idx), normalize_by_dmso_sta, drug_numbers(idx), plate_numbers(idx));
        growth_score_vec(i) = growth_score_here(1);
    end
end

summary_table = [drugs(drug_numbers)', organoids(organoid_values), dose_values, viability_score_vec, growth_score_vec, limits, no_organoids'];
summary_table = summary_table(drug_numbers > 3,:);