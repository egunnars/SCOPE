% A function to plot the data. The parameters ind_green, ind_red, ... are optional
function plot_data_dose_response(growth_score, perc_flokkar, drug_numbers, no_organoids, ind_green, ind_yellow, ind_red, ind_magenta)
   % generate gray color scale for doses
   colors = zeros(size(perc_flokkar, 1), 3);
   c = gray(9);
   c = flip(c);
   for i = 3:size(perc_flokkar, 1)
      colors(i,:) = c(i-2,:);
   end

   if ~exist('ind_green','var') 
      ind_green = find(drug_numbers == 1); % DMSO
   end
   for i = 1:size(ind_green, 1)
      colors(ind_green(i), :) = [0 0 1]; % blue
   end
    
   if ~exist('ind_red','var')
      ind_red = find(drug_numbers == 3); % STA
   end
   for i = 1:size(ind_red, 1)
      colors(ind_red(i), :) = [1 0 0]; % red
   end
   
   % make GV plot
   z = [perc_flokkar(:, 1), perc_flokkar(:, 2), growth_score];
   for k = 3:size(perc_flokkar,1)
       scatter3(z(k, 1), z(k, 2), z(k, 3), 150+(k-3)*40, colors(k,:), "filled", 'MarkerEdgeColor', [0 0 0]);
       hold on
   end
   % add DMSO and STA
   z = [perc_flokkar(ind_green, 1), perc_flokkar(ind_green, 2), growth_score(ind_green)];
   scatter3(z(:, 1), z(:, 2), z(:, 3), 250, colors(ind_green,:), "filled","^", 'MarkerEdgeColor', [0 0 0], 'LineWidth',1.2);
   z = [perc_flokkar(ind_red, 1), perc_flokkar(ind_red, 2), growth_score(ind_red)];
   scatter3(z(:, 1), z(:, 2), z(:, 3), 250, colors(ind_red,:), "filled","v", 'MarkerEdgeColor', [0 0 0], 'LineWidth',1.2);
    
   % axis labels and scales
   xlabel('Live-Live');
   ylabel('Live-Dead');
   zlabel('Growth score');
   xlim([0 1]);
   ylim([0 1]);
   zlim([-0.5,1.25]);

   % add boundaries
   set(gca,'color',[0.95 0.95 0.95]);

   face_alpha_value = 0.3;
   edge_alpha_value = 1;
   min_height = 0;
   max_height = 1;
   colors = [0.9,0.9,0.9];
   fill3([0,0,0,0],[0,1,1,0],[min_height,min_height,max_height,max_height],colors,'FaceAlpha',face_alpha_value,'EdgeAlpha',edge_alpha_value);
   fill3([0,1,1,0],[1,0,0,1],[min_height,min_height,max_height,max_height],colors,'FaceAlpha',face_alpha_value,'EdgeAlpha',edge_alpha_value);
   fill3([0,1,0],[0,0,1],1*[1,1,1],colors,'FaceAlpha',0.2,'EdgeAlpha',0.7);
   fill3([0,1,0],[0,0,1],[0,0,0],colors,'FaceAlpha',0.2,'EdgeAlpha',0.7);

   set(gca,'fontsize', 19);
end
