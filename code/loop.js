//function [ws] = CalcP_dp_Ploidy_3D_OnlyTheLoop (p0, freq, d, ploidy, df_ci, dp_ci)

function CalcP_dp_Ploidy_3D(p0, freq, d, ploidy, df_ci, dp_ci) {


// types = {};
// types{1} = 'Somatic LOH CNmut = 1';
// if (ploidy > 1)
//     l = length(types);
//     for i=1:ploidy;
//         types{l+i} = sprintf ('Somatic, CNmut = %i', i);
//         if (i > 1); types{l+i} = sprintf ('Somatic LOH, CNmut = %i', i); end
//     end
// end
// l = length(types);
// types{l+1} = 'Germline LOH CNmut = 1';
// if (ploidy > 1)
//     l = length(types);
//     for i=1:ploidy;
//         types{l+i} = sprintf ('Germline, CNmut = %i', i);
//         if (i > 1); types{l+i} = sprintf ('Germline LOH, CNmut = %i', i); end
//     end
// end
// l = length(types);

// ddf = 0.005;
// ddp = 0.01;

// fs = sort(freq*(1-df_ci) : ddf : freq*(1+df_ci));
// dp = sort(p0*(1-dp_ci) : ddp : p0*(1+dp_ci));

// aics = zeros (length(fs), length(dp), size(types, 2));