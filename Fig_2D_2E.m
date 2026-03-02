%Fix the random seed for reproducibility
rndseed = 1; rng(rndseed);

doses = [0.00001,0.0001,0.001,0.005,0.01,0.1,1,2.55];

colors = [orderedcolors("gem");0.1 0.3 0.2];
figure;
tiledlayout(1,2);

for ell=1:2

    nexttile(ell);
    if ell == 1
        %Parameters for UA dose response curve
        IC50_vi = 0.0102;
        hill_coef_vi = 0.6386;
        Emin_vi = 0.0336;
        IC50_gr = 0.0056;
        hill_coef_gr = 1.9483;
        Emin_gr = 0.0232;
    else
        %Parameters for V8 dose response curve
        IC50_vi = 0.3156;
        hill_coef_vi = 1.0723;
        Emin_vi = 0.0139;
        IC50_gr = 0.0006;
        hill_coef_gr = 2.8032;
        Emin_gr = 0.0658;
    end
    
    %Generate dose response curves
    percent_inhibition_vi = hill_gr_inhibition(doses, IC50_vi, hill_coef_vi, Emin_vi);
    percent_inhibition_gr = hill_gr_inhibition(doses, IC50_gr, hill_coef_gr, Emin_gr);
    size_end = zeros(1,size(doses,2));
    
    DT = 100; %time point of drug dose
    
    N = 10^4; %number of subclones
    T = 1000; %number of time points
    r = 0.03; % mean exp growth rate
    timevec = linspace(0,20,T);
    count = zeros(size(doses,2),1);
    X = 10^4*ones(N,T,size(doses,2));
    death_time = ceil(900*rand(1,N)); %time at which each subclone dies
    r = r + (rand(N,1)-0.5)*0.01; 
    
    for j=1:size(doses,2) %drug
        rd = r.*(percent_inhibition_gr(j))/100; %all clones growth rate reduced by drug dose
           
        %percent of clones killed by drug dose
        for i=1:DT-1 %time before drug applied
            X(:,i, j) = X(:,1,j).*exp(r*timevec(i));
        end
        count(j) = 0;
        for k=1:N
            tmp = rand;
            if rand*100 < 100-percent_inhibition_vi(j)
                %killing of subclones
                for i=DT:DT+death_time(k)-1
                    X(k,i, j) = X(k,DT-1,j).*exp(rd(k).*(timevec(i) - timevec(DT-1)));
                end
                X(k,DT+death_time(k), j) = 0;
                count(j) = count(j)+1;
            else
                for i=DT:DT+death_time(k)
                    X(k,i, j) = X(k,DT-1,j).*exp(rd(k).*(timevec(i) - timevec(DT-1)));
                end
            end
            for i=DT+death_time(k)+1:T %time after drug applied
                X(k,i, j) = X(k,DT+death_time(k),j).*exp(rd(k).*(timevec(i) - timevec(DT+death_time(k))));
            end
        end
    
    end
    
    % make plots
    for i=1:size(doses,2)
        plot(timevec, sum(X(:,:,i),1), 'color', colors(i,:), 'linewidth', 3);
        hold on;
    end
    
    for j=1:size(doses,2)
        size_end(j) = sum(X(:,end,j));
        if ell == 1
            p = plot(timevec(end),size_end(j),'square','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',colors(j,:), 'MarkerSize',16);
            uistack(p, 'top');
        else
            p = plot(timevec(end),size_end(j),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',colors(j,:),'MarkerSize',14);
            uistack(p, 'top');
        end
    end
    
    xlabel('Time (days)');
    ylabel('Predicted tumor size')
    lgd = legend('0.00001','0.0001','0.001','0.005','0.01','0.1','1','2.55','Orientation','horizontal');
    lgd.NumColumns = 2;
    
    set(gca,'Fontsize',16)
    if ell == 1
        title('UA organoid','FontSize',20);
        else
        title('V8 organoid','FontSize',20);
    end
    ylim([10^6 2*10^8]);
end

function inhibition = hill_gr_inhibition(d, IC50, hill_coef, Emin, Emax)
    if nargin < 4
        hill_coef = 1;
    end
    if nargin < 5
        Emax = 100;
    end

    inhibition = Emax + (Emin - Emax) .* (d .^ hill_coef) ./ (IC50 ^ hill_coef + d .^ hill_coef);
end