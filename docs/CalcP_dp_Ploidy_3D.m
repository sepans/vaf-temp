function [ww, wf, pred_outs] = CalcP_dp_Ploidy_3D (p0, freq, d, ploidy, df_ci, dp_ci)

ddf = 0.005;
ddp = 0.01;

if (nargin < 4) 
    error ('purity, VAF, total depth, and ploidy must be entered. Optional: VAF CI, purity CI.');
elseif (nargin < 5)
    df_ci = 0.01;
    dp_ci = 0.05;
elseif (nargin < 6)
    dp_ci = 0.05;
end    

Colors = {'c', 'k', 'r', 'g',  'b', [1 0.8 0.2], [0.5 0.8 0.5], [0.2 0.3 0.5], ...
    [0.5 0.5 0.5], [0.5 0.3 0.3], [0.5 0.8 1], [1 0.5 1], ...
    [0.8 0.1 0.8], [0.5 0.1 0.8], [0 0.1 .5], [0.5 0.1 0]};

types = {};
types{1} = 'Somatic LOH CNmut = 1';
if (ploidy > 1)
    l = length(types);
    for i=1:ploidy;
        types{l+i} = sprintf ('Somatic, CNmut = %i', i);
        if (i > 1); types{l+i} = sprintf ('Somatic LOH, CNmut = %i', i); end
    end
end
l = length(types);
types{l+1} = 'Germline LOH CNmut = 1';
if (ploidy > 1)
    l = length(types);
    for i=1:ploidy;
        types{l+i} = sprintf ('Germline, CNmut = %i', i);
        if (i > 1); types{l+i} = sprintf ('Germline LOH, CNmut = %i', i); end
    end
end

[~, df] = binofit(round(d*freq), d, df_ci);

fs = df(1) : ddf : df(2); 
fs = sort([fs(1:end-1), df(2)]);
dp = sort(p0*(1-dp_ci) : ddp : p0*(1+dp_ci));

aics = zeros (length(fs), length(dp), size(types, 2));

for j=1:length(dp)
    p = dp(j);
    for k=1:length(fs)
        f = fs(k);
        
        aic = [];
        aic(1) = 2 - 2 * log (binopdf (round(d*f), d, (p)/(2*(1-p)+1*p))); %somatic LOH
        if (ploidy > 1)
            l = length(aic);
            for i=1:ploidy
                aic(l+i) = 2 - 2 * log (binopdf (round(d*f), d, (i*p)/(2*(1-p)+ploidy*p))); %somatic LOH high CN
            end
        end
        l = length(aic);
        aic(l+1) = 2 - 2 * log (binopdf (round(d*f), d, (1-p+p)/(2*(1-p)+1*p))); %germline LOH high CN;        
        if (ploidy > 1)
            l = length(aic);
            for i=1:ploidy
                aic(l+i) = 2 - 2 * log (binopdf (round(d*f), d, (1-p+i*p)/(2*(1-p)+ploidy*p))); %germline LOH high CN;
            end
        end

        w = zeros(1, size(types, 2));
        if (j == 1); ws = zeros(length(fs), length(dp), size(types, 2)); end;
        not_nan = ~isnan(aic);

        D = sum(exp(-0.5*(aic(not_nan)-min(aic(not_nan)))));
        w(not_nan) = exp(-0.5*(aic(not_nan)-min(aic(not_nan)))) / D;

        ws (k, j, 1:size(types, 2)) = w;
        aics (k, j, 1:size(types, 2)) = aic;
        clear aic w w_sorted_i;
    end
end

not_nan = ~isnan(aics);
D = sum(exp(-0.5*(aics(not_nan)-min(aics(not_nan)))));
ww = exp(-0.5*(aics-min(aics(not_nan)))) / D;

sum_ww = zeros (1, size(types, 2));
for i=1:size(types, 2)
    w = ww(:,:,i); 
    sum_ww (i) = sum(w(~isnan(w)));
    pred_outs{i} = sprintf ('%s\t%2.2e', types{i}, sum_ww (i));
    leg{i} = sprintf ('%s, w = %2.2f', types{i}, sum_ww (i));
end

aics_p = zeros (length(fs), size(types, 2));
wf = zeros (length(fs), length(dp), size(types, 2));
for j=1:length(dp)
    for i=1:size(types, 2)
        aics_p(1:length(fs), i) = aics(:,j,i); 
    end
    
    for k=1:length(fs)
        aic_f = aics_p(k,:);
        not_nan = ~isnan(aic_f);
        D = sum(exp(-0.5*(aic_f(not_nan)-min(aic_f(not_nan)))));
        wf(k, j, 1:size(types, 2)) = exp(-0.5*(aic_f-min(aic_f(not_nan)))) / D;
    end
end
    
[~, i] = sort (sum_ww, 'descend');
fprintf ('%i\t%s\t%2.2e\t%s\t%2.2e\n', round(p0*100), types{i(1)}, sum_ww(i(1)), types{i(2)}, sum_ww(i(2)));

pred_out = sprintf ('%s (%2.2f)', types{i(1)}, sum_ww(i(1)));
if (max(sum_ww) < 0.99) 
    pred_out = sprintf ('%s, %s (%2.2f)', pred_out, types{i(2)}, sum_ww(i(2)));
end

figure;
for k=1:size(types, 2)
    up = min(mean(wf(:,:,k), 2)+std(wf(:,:,k), 0, 2), 1)-mean(wf(:,:,k), 2);
    down = -max(mean(wf(:,:,k), 2)-std(wf(:,:,k), 0, 2), 0)+mean(wf(:,:,k), 2);
    ps = shadedErrorBar (fs, mean(wf(:,:,k), 2), [up down], '.', Colors{k}); 
    if (size(ps.patch,1))
        ps.patch.FaceColor = Colors{k}; ps.edge(1).Color = Colors{k}; ps.edge(2).Color = Colors{k};
    end
    hold on;
    pp(k) = plot (fs, mean(wf(:,:,k), 2), 'o-', 'Color', Colors{k}, 'MarkerSize', 5, 'LineWidth', 2);  
    set(pp(k),'MarkerEdgeColor','k','MarkerFaceColor', Colors{k});
    hold on; 
end

title (sprintf ('VAF: %2.2f - Depth: %i - Purity %2.2f (%i%% CI)\n%s', freq, d, p0, round(dp_ci*100), pred_out));

xlabel (sprintf ('VAF (%0.0f%% Confidence Interval)', 100-100*df_ci));
ylabel ('W');

title (sprintf ('VAF: %2.2f - Depth: %i - Purity %2.2f (%i%% CI)\n%s', freq, d, p0, round(dp_ci*100), pred_out));

hold on;
plot ([freq freq], [0 1], '.-', 'Color', [.7 .7 .7], 'LineWidth', 2);
xlim ([fs(1)-ddf fs(end)+ddf]);
ylim ([-.01 1.01]);
legend (pp, leg, 'location', 'northeastoutside');
