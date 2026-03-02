% A function that takes an array if indicies, filter, and chooses only those groups
function [no_organoids_new, perc_flokkar_new, params_new, params_2_new, error_vec_bic_new, error_vec_2_bic_new, plate_numbers_new, drug_numbers_new, dose_values_new] = select_groups_dose_response(filter, no_organoids_old, perc_flokkar_old, params_old, params_2_old, error_vec_bic_old, error_vec_2_bic_old, plate_numbers_old, drug_numbers_old, dose_values_old)
   no_organoids_new = no_organoids_old(filter);
   perc_flokkar_new = perc_flokkar_old(filter, :);
   params_new = params_old(:, :, filter);
   params_2_new = params_2_old(:, :, filter);
   error_vec_bic_new = error_vec_bic_old(:, filter);
   error_vec_2_bic_new = error_vec_2_bic_old(:, filter);
   plate_numbers_new = plate_numbers_old(filter);
   drug_numbers_new = drug_numbers_old(filter);
   dose_values_new = dose_values_old(filter);
end
