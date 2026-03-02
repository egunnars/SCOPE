% A function to plot the data. The parameters ind_green, ind_red, ... are optional
function plot_data_with_controls(growth_score, perc_flokkar, drug_numbers, view_vec)
   % assign all drugs blank color
   colors = zeros(size(perc_flokkar, 1), 3);

   for i = 1:size(perc_flokkar, 1)
       colors(i, :) = [1 1 1];
   end
   
   % draw plot
   a = find(drug_numbers >= 5);
   z = [perc_flokkar(a, 1), perc_flokkar(a, 2), growth_score(a)];
   scatter3(z(:, 1), z(:, 2), z(:, 3), 130, colors(a,:), "filled", 'LineWidth', 1.1, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceAlpha',0,'MarkerEdgeAlpha',1);
   hold on

    % add DMSO dot
    ind_green = find(drug_numbers == 1);
    for i = 1:size(ind_green, 1)
    colors(ind_green(i), :) = [0 0 200/255]; % green
    end
    z = [perc_flokkar(ind_green, 1), perc_flokkar(ind_green, 2), growth_score(ind_green)];
    z = mean(z);
    h = scatter3(z(:, 1), z(:, 2), z(:, 3), 350, colors(ind_green(1),:), "filled","^", 'MarkerEdgeColor', [0 0 0], 'LineWidth',1.6);
    uistack(h,'top');
    
    % add 5FU dot
    ind_yellow = find(drug_numbers == 2);
    for i = 1:size(ind_yellow, 1)
    colors(ind_yellow(i), :) = [36/255 217/255 238/255]; % yellow 
    end
    z = [perc_flokkar(ind_yellow, 1), perc_flokkar(ind_yellow, 2), growth_score(ind_yellow)];
    z = mean(z);
    h = scatter3(z(:, 1), z(:, 2), z(:, 3), 550, colors(ind_yellow(1),:), "filled", "pentagram", 'MarkerEdgeColor', [0 0 0], 'LineWidth',1.6);
    uistack(h,'top');
    
    % add STA dot
    ind_red = find(drug_numbers == 3); 
    for i = 1:size(ind_red, 1)
    colors(ind_red(i), :) = [1 0 0]; % red
    end
    z = [perc_flokkar(ind_red, 1), perc_flokkar(ind_red, 2), growth_score(ind_red)];
    z = mean(z);
    h = scatter3(z(:, 1), z(:, 2), z(:, 3), 350, colors(ind_red(1),:), "filled","v", 'MarkerEdgeColor', [0 0 0], 'LineWidth',1.6);
    uistack(h,'top');
    
    % add SN38 dot
    ind_magenta = find(drug_numbers == 4);
    for i = 1:size(ind_magenta, 1)
    colors(ind_magenta(i), :) = [255 165 0]/255; %[246/255 194/255 23/255]; % magenta
    end
    z = [perc_flokkar(ind_magenta, 1), perc_flokkar(ind_magenta, 2), growth_score(ind_magenta)];
    z = mean(z);
    h = scatter3(z(:, 1), z(:, 2), z(:, 3), 400, colors(ind_magenta(1),:), "filled", "diamond", 'MarkerEdgeColor', [0 0 0], 'LineWidth',1.6);
    uistack(h,'top');

   % axis labels and scales
   xlabel('Live-Live');
   ylabel('Live-Dead');
   zlabel('Growth score');
   xlim([0 1]);
   ylim([0 1]);
   zlim([-0.1 1.2]);
   set(gca,'color',[0.95 0.95 0.95]);

   % draw boundary
   face_alpha_value = 0.3;
   edge_alpha_value = 1;
   min_height = 0;
   max_height = 1;
   colors = [0.9 0.9 0.9];
   fill3([0,0,0,0],[0,1,1,0],[min_height,min_height,max_height,max_height],colors,'FaceAlpha',face_alpha_value,'EdgeAlpha',edge_alpha_value);
   fill3([0,1,1,0],[1,0,0,1],[min_height,min_height,max_height,max_height],colors,'FaceAlpha',face_alpha_value,'EdgeAlpha',edge_alpha_value);
   fill3([0,1,0],[0,0,1],1*[1,1,1],colors,'FaceAlpha',0.2,'EdgeAlpha',0.7);
   fill3([0,1,0],[0,0,1],[0,0,0],colors,'FaceAlpha',0.2,'EdgeAlpha',0.7);

   view(view_vec);
   set(gcf, 'Position',  [100, 100, 1000, 600]);
end