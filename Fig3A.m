% fix the random seed for reproducibility
rndseed = 1; rng(rndseed);

combine = true; % set to true to combine wells

specific_file_name = "up_high_dose_analysis.mat"; 

normalize_by_dmso_sta = true; % normalize with dmso / sta or max and min
specific_drugs = false; % set to true to pick a list of drugs
specific_drugs_list = []; % drugs to display

figure;
data_files = specific_file_name;

file_name = data_files;
load(file_name);
perc_flokkar = perc_flokkar_fully; % only use fully tracked organoids

% rename the variable for consistancy
plate_numbers = plate_number;
clear plate_number;

% group organoids by wells
if combine
  [no_organoids, perc_flokkar, params, params_2, error_vec_bic, error_vec_2_bic, plate_numbers, drug_numbers] = combine_wells(no_organoids, perc_flokkar, params, params_2, error_vec_bic, error_vec_2_bic, plate_numbers, drug_numbers);
end 

% remove groups that have 0 organoids
if sum(no_organoids == 0) > 0
  fprintf('The groups with indices ');
  fprintf('%d ', find(no_organoids == 0));
  fprintf('are zero, they will be removed.\n');
end
filter = find(no_organoids ~= 0);
[no_organoids, perc_flokkar, params, params_2, error_vec_bic, error_vec_2_bic, plate_numbers, drug_numbers] = select_groups(filter, no_organoids, perc_flokkar, params, params_2, error_vec_bic, error_vec_2_bic, plate_numbers, drug_numbers);

% calculate the growth scores
if normalize_by_dmso_sta
  [growth_score, ~, ~] = calculate_growth_score(no_organoids, params, params_2, error_vec_bic, error_vec_2_bic, true, drug_numbers, plate_numbers);
else
  [growth_score, ~, ~] = calculate_growth_score(no_organoids, params, params_2, error_vec_bic, error_vec_2_bic);
end

% select chosen drugs
if specific_drugs
  filter = find(ismember(drug_numbers, specific_drugs_list));
  [no_organoids, perc_flokkar, params, params_2, error_vec_bic, error_vec_2_bic, plate_numbers, drug_numbers] = select_groups(filter, no_organoids, perc_flokkar, params, params_2, error_vec_bic, error_vec_2_bic, plate_numbers, drug_numbers);
  growth_score = growth_score(filter);
end

load("drug_list.mat");
% generate labels for groups (they are different depending on whether we combine wells)
if combine
    labels = strings(1,size(drug_numbers,1));
    for ii = 1:size(labels,2)
        labels(ii) = drug_list(drug_numbers(ii));
    end
else
    labels = cellstr(num2str(drug_numbers));
    labels = append(labels, '/');
    labels = append(labels, cellstr(num2str(plate_numbers)));
end

% normalize perc_flokkar
perc_flokkar = perc_flokkar./sum(perc_flokkar')'; 

% plot the data
view_vec = [22.1, 18.9];
plot_data_with_controls(growth_score, perc_flokkar, drug_numbers, view_vec);

% add a legend
h0 = plot(NaN, NaN, 'o', 'MarkerSize',50, 'MarkerFaceColor', [0.9 0.9 0.9], 'MarkerEdgeColor', 'black', 'DisplayName', '166 anti-cancer compounds');
h5 = plot(NaN, NaN, "^", 'MarkerSize',100, 'MarkerFaceColor', [0 0 200/255], 'MarkerEdgeColor', 'black', 'DisplayName', 'DMSO');
h6 = plot(NaN, NaN, "v", 'MarkerSize',100, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'black', 'DisplayName', 'STA');
h7 = plot(NaN, NaN, "pentagram", 'MarkerSize',150, 'MarkerFaceColor', [36/255 217/255 238/255], 'MarkerEdgeColor', 'black', 'DisplayName', '5FU');
h8 = plot(NaN, NaN, "diamond", 'MarkerSize',100,'MarkerFaceColor', [255 165 0]/255, 'MarkerEdgeColor', 'black', 'DisplayName', 'SN38');

lgd = legend([h0, h5, h7, h6, h8],'Location','NorthOutside','Orientation', 'horizontal');
set(gca,'fontsize', 19);

% put control groups in focus
ax = gca;
ax.Children = ax.Children([10:14,1:9]);
ax.SortMethod = 'childorder';

% calculate confidence intervals
limits = zeros(size(drug_numbers,1),2);
for i = 1:size(drug_numbers,1)
    if drug_numbers(i) > 4
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
            [growth_score_here, growth_score_ab_here, growth_score_a_here] = calculate_growth_score(no_organoids(idx), params_here(:,:,idx), params_2_here(:,:,idx), error_vec_bic_here(:,idx), error_vec_2_bic_here(:,idx), normalize_by_dmso_sta, drug_numbers(idx), plate_numbers(idx));
            growth_score_vec(j) = growth_score_here(1);
        end
        limits(i,:) = [prctile(growth_score_vec,2.5), prctile(growth_score_vec,97.5)];
    end
end

summary_table = [drug_numbers, perc_flokkar, growth_score, limits, no_organoids'];
summary_table = summary_table(drug_numbers > 4,:);