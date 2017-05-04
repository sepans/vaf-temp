function [ws] = CalcP_dp_Ploidy_3D_OnlyTheLoop (p0, freq, d, ploidy, df_ci, dp_ci)

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
l = length(types);

ddf = 0.005;
ddp = 0.01;

fs = sort(freq*(1-df_ci) : ddf : freq*(1+df_ci));
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
