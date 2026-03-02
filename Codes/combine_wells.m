% A function to combine the wells such that if two wells on the same plate have the same drug they are combined.
% NOTE: this function expects perc_flokkar to be a count and not a percentage. The normalization must not be done before this function is called.
function [no_organoids_new, perc_flokkar_new, params_new, params_2_new, error_vec_bic_new, error_vec_2_bic_new, plate_numbers_new, drug_numbers_new] = combine_wells(no_organoids_old, perc_flokkar_old, params_old, params_2_old, error_vec_bic_old, error_vec_2_bic_old, plate_numbers_old, drug_numbers_old)

   % security check to make sure sizes fit
   if ~isequal(size(no_organoids_old, 2), size(params_old, 3), size(params_2_old, 3), size(error_vec_bic_old, 2), size(error_vec_2_bic_old, 2), size(plate_numbers_old, 1), size(drug_numbers_old, 1))
      error('combine_wells: sizes of inputs do not match');
      return
   end

   % calculate the number of organoids of each group
   no_organoids_new = [];
   for plate_no = 1:max(plate_numbers_old)
      for drug_no = 1:max(drug_numbers_old)
         well_idxs = find(plate_numbers_old == plate_no & drug_numbers_old == drug_no);
         if numel(well_idxs) == 0
            continue
         end
         no_organoids_new = [no_organoids_new, sum(no_organoids_old(well_idxs))];
      end
   end

   % find the maximum number of organoids in order to format the data correctly
   max_no_org = max(no_organoids_new);
   no_groups = size(no_organoids_new, 2);

   % create arrays to store the data
   params_new = zeros(max_no_org, 2, no_groups);
   params_2_new = zeros(max_no_org, 3, no_groups);
   error_vec_bic_new = zeros(max_no_org, no_groups);
   error_vec_2_bic_new = zeros(max_no_org, no_groups);

   perc_flokkar_new = [];
   plate_numbers_new = [];
   drug_numbers_new = [];

   % go through all plate, drug pairs and fill the arrays
   new_idx = 0;
   for plate_no = 1:max(plate_numbers_old)
      for drug_no = 1:max(drug_numbers_old)
         well_idxs = find((plate_numbers_old == plate_no & drug_numbers_old == drug_no));
         if numel(well_idxs) == 0
            continue;
         end
         new_idx = new_idx + 1;

         % go through each well in group
         no_organoids_curr = 0;
         for i = 1:numel(well_idxs)
            old_idx = well_idxs(i);
            orgs = no_organoids_old(old_idx);

            % append parameter data and error vec data 
            left_idxs = no_organoids_curr + 1 : no_organoids_curr + orgs;
            right_idxs = 1 : no_organoids_old(old_idx);
            params_new(left_idxs, :, new_idx) = params_old(right_idxs, :, old_idx);
            params_2_new(left_idxs, :, new_idx) = params_2_old(right_idxs, :, old_idx);

            error_vec_bic_new(left_idxs, new_idx) = error_vec_bic_old(right_idxs, old_idx);
            error_vec_2_bic_new(left_idxs, new_idx) = error_vec_2_bic_old(right_idxs, old_idx);

            no_organoids_curr = no_organoids_curr + orgs;
         end

         plate_numbers_new = [plate_numbers_new; plate_no];
         drug_numbers_new = [drug_numbers_new; drug_no];
         perc_flokkar_new = [perc_flokkar_new; sum(perc_flokkar_old(well_idxs, :), 1)];
      end
   end

   % security check to make sure sizes fit
   if ~isequal(size(no_organoids_new, 2), size(params_new, 3), size(params_2_new, 3), size(error_vec_bic_new, 2), size(error_vec_2_bic_new, 2), size(plate_numbers_new, 1), size(drug_numbers_new, 1))
      error('combine_wells: sizes of outputs do not match');
      return
   end

end
