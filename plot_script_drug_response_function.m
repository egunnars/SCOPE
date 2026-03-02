function [perc_flokkar, growth_score, doses_part] = plot_script_drug_response_function(specific_file, specific_file_name, specific_drugs, specific_drugs_list, ...
    normalize_by_dmso_sta, specific_plates, specific_plates_list_single_file, specific_plates_list_all_files, plate_numbers, dose_values)

data_files = specific_file_name;

% iterate through files
for file_no = 1:numel(data_files)
   %nexttile;
   file_name = data_files(file_no);
   load(file_name);

   perc_flokkar = perc_flokkar_fully; % only use fully tracked organoids

   % select chosen plates
   if specific_plates
      if specific_file
         filter = find(ismember(plate_numbers, specific_plates_list_single_file));
         [no_organoids, perc_flokkar, params, params_2, error_vec_bic, error_vec_2_bic, plate_numbers, drug_numbers, dose_values] = select_groups_dose_response(filter, no_organoids, perc_flokkar, params, params_2, error_vec_bic, error_vec_2_bic, plate_numbers, drug_numbers, dose_values);
      else 
         filter = find(ismember(plate_numbers, specific_plates_list_all_files(file_no, :)));
         [no_organoids, perc_flokkar, params, params_2, error_vec_bic, error_vec_2_bic, plate_numbers, drug_numbers, dose_values] = select_groups_dose_response(filter, no_organoids, perc_flokkar, params, params_2, error_vec_bic, error_vec_2_bic, plate_numbers, drug_numbers, dose_values);
      end
   end

   % remove groups that have 0 organoids
   if sum(no_organoids == 0) > 0
      fprintf('The groups with indices ');
      fprintf('%d ', find(no_organoids == 0));
      fprintf('are zero, they will be removed.\n');
   end
   filter = find(no_organoids ~= 0);

   % calculate the growth scores
   if normalize_by_dmso_sta
      [growth_score, ~, ~] = calculate_growth_score_dose_response(no_organoids, params, params_2, error_vec_bic, error_vec_2_bic, true, drug_numbers, plate_numbers);
   else
      [growth_score, ~, ~] = calculate_growth_score_dose_response(no_organoids, params, params_2, error_vec_bic, error_vec_2_bic);
   end

   % select chosen drugs
   if specific_drugs
      filter = find(ismember(drug_numbers, specific_drugs_list));
      [no_organoids, perc_flokkar, params, params_2, error_vec_bic, error_vec_2_bic, plate_numbers, drug_numbers, dose_values] = select_groups_dose_response(filter, no_organoids, perc_flokkar, params, params_2, error_vec_bic, error_vec_2_bic, plate_numbers, drug_numbers, dose_values);
      growth_score = growth_score(filter);
   end

   % normalize perc_flokkar
   perc_flokkar = perc_flokkar./sum(perc_flokkar')'; 

   %plot the data
   plot_data_dose_response(growth_score, perc_flokkar, drug_numbers, no_organoids);
   view(15,10)

   doses_part = dose_values;
end

end
