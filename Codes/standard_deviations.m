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

%Put control groups in focus
ax = gca;
ax.Children = ax.Children([10:14,1:9]);
ax.SortMethod = 'childorder';


limits = zeros(size(drug_numbers,1),2);
for i=5:170
    i
    a = find(drug_numbers == i);
    b = find(drug_numbers == 1 & plate_numbers == plate_numbers(a));
    c = find(drug_numbers == 3 & plate_numbers == plate_numbers(a));
    idx = [a,b,c];
    growth_score_vec = zeros(1000,1);
    for j=1:1000
        sample = randi(no_organoids(a),no_organoids(a),1);
        [growth_score, growth_score_ab, growth_score_a] = calculate_growth_score(no_organoids(idx), params(sample,:,idx), params_2(sample,:,idx), error_vec_bic(sample,idx), error_vec_2_bic(sample,idx), normalize_by_dmso_sta , drug_numbers(idx), plate_numbers(idx));
        growth_score_vec(j) = growth_score(1);
    end
    limits(i,:) = [prctile(growth_score_vec,5), prctile(growth_score_vec,95)];
end
%}

%{
load(file_name);
plate_numbers = plate_number;
triplicates = zeros(size(drug_numbers,1),3);
for i=5:170
    i
    d = find(drug_numbers == i);
    for k=1:3
        a = d(k);
        b = find(drug_numbers == 1 & plate_numbers == plate_numbers(a));
        c = find(drug_numbers == 3 & plate_numbers == plate_numbers(a));
        idx = [a;b;c];
        [growth_score, growth_score_ab, growth_score_a] = calculate_growth_score_stdev_1(no_organoids(idx), params(:,:,idx), params_2(:,:,idx), error_vec_bic(:,idx), error_vec_2_bic(:,idx), normalize_by_dmso_sta , drug_numbers(idx), plate_numbers(idx));
        triplicates(i,k) = growth_score(1);
    end
end
%}
