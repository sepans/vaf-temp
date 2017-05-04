
function CalcP_dp_Ploidy_3D(p0, freq, d, ploidy, df_ci, dp_ci) {

	var ddf = 0.005
	var ddp = 0.01

	if (arguments.length < 4) {
	    var err = 'purity, VAF, total depth, and ploidy must be entered. Optional: VAF CI, purity CI.'
	    console.error(err)
		throw(err)	    
	 }
	else if (arguments.length < 5) {
	    df_ci = 0.01;
	    dp_ci = 0.05;

	}
	else if (arguments.length < 6) {
	    dp_ci = 0.05;
	}

// Colors = {'c', 'k', 'r', 'g',  'b', [1 0.8 0.2], [0.5 0.8 0.5], [0.2 0.3 0.5], ...
//     [0.5 0.5 0.5], [0.5 0.3 0.3], [0.5 0.8 1], [1 0.5 1], ...
//     [0.8 0.1 0.8], [0.5 0.1 0.8], [0 0.1 .5], [0.5 0.1 0]};

// types = {};
// types{1} = 'Somatic LOH CNmut = 1';
// if (ploidy > 1)
//     l = length(types);
//     for i=1:ploidy;
//         types{l+i} = sprintf ('Somatic, CNmut = %i', i);
//         if (i > 1); types{l+i} = sprintf ('Somatic LOH, CNmut = %i', i); end
//     end
// end

	types = [];
	types[1] = 'Somatic LOH CNmut = 1';
	if (ploidy > 1) {
    	var l = types.length;
    	for(var i=1; i<ploidy; i++) {
        	types[l+i] = 'Somatic, CNmut = ' + i
        	if (i > 1) {
        		types[l+i] = 'Somatic LOH, CNmut = ' + i	
        	} 
        }
    }

// [~, df] = binofit(round(d*freq), d, df_ci);

// fs = df(1) : ddf : df(2); 
// fs = sort([fs(1:end-1), df(2)]);
// dp = sort(p0*(1-dp_ci) : ddp : p0*(1+dp_ci));
    



}

CalcP_dp_Ploidy_3D(1, 2, 3 ,4)