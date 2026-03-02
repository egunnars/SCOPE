combine = true; % set to true to combine wells

specific_file_name = "up_high_dose_analysis.mat"; 

normalize_by_dmso_sta = true; % normalize with dmso / sta or max and min
specific_drugs = false; % set to true to pick a list of drugs
specific_drugs_list = [114, 111, 74, 169, 119, 11]; % drugs to display

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
    labels = append(labels, cellstr(num2str(plate_numbers)))
end

% normalize perc_flokkar
perc_flokkar = perc_flokkar./sum(perc_flokkar')'; 

% classify into phenotypic categories
livelive = perc_flokkar(:, 1);
livedead = perc_flokkar(:, 2);
deaddead = perc_flokkar(:, 3);

cytotoxic_bool = deaddead >= 0.7;
late_cytotoxic_bool = (livedead >= 0.8 - 0.866667.*livelive) & (livedead >= 0.15);
cytostatic_bool = livelive >= 0.75 & livedead <= 0.15;
cytotoxic_cytostatic_bool = ~(cytotoxic_bool | late_cytotoxic_bool | cytostatic_bool);

ind_cytotoxic = find(cytotoxic_bool == 1);
ind_cytotoxic_cytostatic = find(cytotoxic_cytostatic_bool == 1);
ind_late_cytotoxic = find(late_cytotoxic_bool == 1);
ind_cytostatic = find(cytostatic_bool == 1);

% plot the data
view_vec = [22.1, 18.9];
plot_data(growth_score, perc_flokkar, labels, drug_numbers, no_organoids, ind_late_cytotoxic, ind_cytotoxic_cytostatic, ind_cytotoxic, [], view_vec, specific_drugs);

% add a legend
h1 = plot(NaN, NaN, 'o', 'MarkerSize',50, 'MarkerFaceColor', [255/255, 128/255, 114/255], 'MarkerEdgeColor', 'black', 'DisplayName', 'cytotoxic');
h2 = plot(NaN, NaN, 'o', 'MarkerSize',50, 'MarkerFaceColor', [255,211,44]/255, 'MarkerEdgeColor', 'black', 'DisplayName', 'cytotoxic + cytostatic');
h3 = plot(NaN, NaN, 'o', 'MarkerSize',50, 'MarkerFaceColor', [163,121,217]/255, 'MarkerEdgeColor', 'black', 'DisplayName', 'late cytotoxic');
h4 = plot(NaN, NaN, 'o', 'MarkerSize',50, 'MarkerFaceColor', [100/255, 149/255, 237/255], 'MarkerEdgeColor', 'black', 'DisplayName', 'cytostatic');

set(gca,'fontsize', 19);
lgd = legend([h1, h2, h3, h4],'Location','NorthOutside','Orientation', 'horizontal');
