% A function to plot the data. The parameters ind_green, ind_red, ... are optional
function plot_data(growth_score, perc_flokkar, labels, drug_numbers, no_organoids, ind_green, ind_yellow, ind_red, ind_magenta, view_vec, specific_drugs)
   % set the colors of the phenotypic categories
   colors = zeros(size(perc_flokkar, 1), 3);
    
   for i = 1:size(perc_flokkar, 1)
       colors(i, :) = [100/255, 149/255, 255/255];
   end
    
   for i = 1:size(ind_green, 1)
      colors(ind_green(i), :) = [163,121,217]/255;
   end
    
   for i = 1:size(ind_yellow, 1)
      colors(ind_yellow(i), :) = [255,211,44]/255; 
   end
    
   for i = 1:size(ind_red, 1)
      colors(ind_red(i), :) = [255/255, 128/255, 114/255];
   end
    
   for i = 1:size(ind_magenta, 1)
      colors(ind_magenta(i), :) = [1 0 1];
   end
   
   % make 3D plot
   a = find(drug_numbers > 4);
   z = [perc_flokkar(a, 1), perc_flokkar(a, 2), growth_score(a)];
   scatter3(z(:, 1), z(:, 2), z(:, 3), 120, colors(a,:), "filled", 'MarkerEdgeColor', [0 0 0], 'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',1);
   hold on

   if specific_drugs
        % draw labels next to each dot
        z = [perc_flokkar(a, 1), perc_flokkar(a, 2), growth_score(a)];
        x = z(:, 1); y = z(:, 2); z = z(:, 3);
        dx = 0.01; dy = 0.01; dz = 0.01; % displacement so the text does not overlay the data points
        text(x+dx, y+dy, z+dz, labels, 'fontsize', 18);
        set(gca,'fontsize', 16);
   end
    
   % axis labels and scales
   xlabel('Live-Live');
   ylabel('Live-Dead');
   zlabel('Growth score');
   xlim([0 1]);
   ylim([0 1]);
   zlim([-0.1 1.2]);
   set(gca,'color',[0.9 0.9 0.9]);

   % color background according to phenotypic categories
   color_cytostatic = [100/255, 149/255, 237/255];
   color_cytotoxic = [255/255, 128/255, 114/255];
   color_cytos_cytotox = [255,211,44]/255;
   color_late_cytotoxic = [163,121,217]/255;

   face_alpha_value = 0.3;
   edge_alpha_value = 1;
   min_height = -0.1;
   max_height = 1.2;
   fill3([0,0.3,0],[0,0,0.3],min_height*ones(1,3),color_cytotoxic,'FaceAlpha',face_alpha_value,'EdgeAlpha',edge_alpha_value);
   fill3([0,0,0,0],[0,0.3,0.3,0],[min_height,min_height,max_height,max_height],color_cytotoxic,'FaceAlpha',face_alpha_value,'EdgeAlpha',edge_alpha_value);
   fill3([0,0.3,0.75,0.75,0],[0.3,0,0,0.15,0.8],min_height*ones(1,5),color_cytos_cytotox,'FaceAlpha',face_alpha_value,'EdgeAlpha',edge_alpha_value);
   fill3([0,0,0,0],[0.3,0.8,0.8,0.3],[min_height,min_height,max_height,max_height],color_cytos_cytotox,'FaceAlpha',face_alpha_value,'EdgeAlpha',edge_alpha_value);
   fill3([0.75,1,0.85,0.75],[0,0,0.15,0.15],min_height*ones(1,4),color_cytostatic,'FaceAlpha',face_alpha_value,'EdgeAlpha',edge_alpha_value);
   fill3([1,1,0.85,0.85],[0,0,0.15,0.15],[min_height,max_height,max_height,min_height],color_cytostatic,'FaceAlpha',face_alpha_value,'EdgeAlpha',edge_alpha_value);
   fill3([0,0.75,0.85,0],[0.8,0.15,0.15,1],min_height*ones(1,4),color_late_cytotoxic,'FaceAlpha',face_alpha_value,'EdgeAlpha',edge_alpha_value);
   fill3([0,0,0,0],[0.8,1,1,0.8],[min_height,min_height,max_height,max_height],color_late_cytotoxic,'FaceAlpha',face_alpha_value,'EdgeAlpha',edge_alpha_value);
   fill3([0,0.85,0.85,0],[1,0.15,0.15,1],[min_height,min_height,max_height,max_height],color_late_cytotoxic,'FaceAlpha',face_alpha_value,'EdgeAlpha',edge_alpha_value);
   
   view(view_vec)
   set(gcf, 'Position',  [100, 100, 1000, 600]);
end
