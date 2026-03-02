% Calculates the growth score provided the number of organoids in each well, the parameters of the model, and the bic errors of the model fits. 
% The parameters norm_DMSO_STA, drug_numbers and plate_numbers are optional and only if we want normalization. 
function [growth_score, growth_score_ab, growth_score_a] = calculate_growth_score_dose_response(no_organoids, params, params_2, error_vec_bic, error_vec_2_bic, norm_DMSO_STA, drug_numbers, plate_numbers)

   % make sure the sizes of the arrays fit
   if ~isequal(size(no_organoids, 2), size(params, 3), size(params_2, 3), size(error_vec_bic, 2), size(error_vec_2_bic, 2))
      error('calculate_growth_score: sizes of input arrays do not match')
      return
   end

   % make sure all data points have organoids
   if ismember(0, no_organoids)
      error('calculate_growth_score: some point has no organoids')
      return
   end

   no_points = size(no_organoids, 2);

   % logical arrays to check which model to use for each organoid
   fit_by_gomp = error_vec_2_bic < error_vec_bic;
   fit_by_exp = error_vec_bic < error_vec_2_bic;

   % arrays to store the center values of the parameters for each data point
   gomp_ab = [];
   exp_a = [];

   % number of organoids of each type
   no_gomp = [];
   no_exp = [];

   for i = 1:no_points
      no_gomp = [no_gomp ; sum(fit_by_gomp(:, i))];
      no_exp = [no_exp ; sum(fit_by_exp(:, i))]; 

      gomp_ab = [gomp_ab ; interquartile_mean(params_2(fit_by_gomp(:, i), 1, i)./params_2(fit_by_gomp(:, i), 2, i))];
      exp_a = [exp_a ; interquartile_mean(params(fit_by_exp(:, i), 1, i))];
   end

   % make sure the number of organoids makes sense
   if no_gomp + no_exp ~= no_organoids'
      error('calculate_growth_score: no_gom, no_ex and no_organoids do not match')
      return
   end

   % normalization using DMSO and STA
   if exist('norm_DMSO_STA', 'var') && norm_DMSO_STA 
      growth_score_ab = zeros(no_points, 1);
      growth_score_a = zeros(no_points, 1);
      for i = 1:max(plate_numbers)
         p_no = i;
         dmso_ab = mean(gomp_ab(drug_numbers == 1 & plate_numbers == p_no));
         dmso_a = mean(exp_a(drug_numbers == 1 & plate_numbers == p_no));
         sta_ab = mean(gomp_ab(drug_numbers == 3 & plate_numbers == p_no));
         sta_a = mean(exp_a(drug_numbers == 3 & plate_numbers == p_no));
         p_idxs = plate_numbers == p_no;
         growth_score_ab(p_idxs) = (gomp_ab(p_idxs) - sta_ab)./(dmso_ab - sta_ab);
         growth_score_a(p_idxs) = (exp_a(p_idxs) - sta_a)./(dmso_a - sta_a);
      end
   else 
      % otherwise normalize using highest and lowest scores
      growth_score_ab = (gomp_ab-min(gomp_ab))./(max(gomp_ab)-min(gomp_ab));
      growth_score_a = (exp_a-min(exp_a))./(max(exp_a)-min(exp_a));
   end

   % handle the possibility of only having one type of organoid 
   growth_score_ab(no_gomp == 0) = 0;    
   growth_score_a(no_exp == 0) = 0; % (these scores will be multiplied by zero later)

   % combine the growth scores via weighted average
   growth_score = (no_gomp .* growth_score_ab + no_exp .* growth_score_a) ./ no_organoids';

   % normalize the growth_score again:
   if ~exist('norm_DMSO_STA', 'var') || ~norm_DMSO_STA 
      growth_score = (growth_score-min(growth_score))./(max(growth_score)-min(growth_score));
   end

   % helper function to calculate the interquartile mean
   function iq_mean = interquartile_mean(data)
       data = data(:)';
       sorted_data = sort(data);
       Q1 = prctile(sorted_data, 25);
       Q3 = prctile(sorted_data, 75);
       interquartile_data = sorted_data(sorted_data >= Q1 & sorted_data <= Q3);
       iq_mean = mean(interquartile_data);
   end

end
