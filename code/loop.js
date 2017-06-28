//function [ws] = CalcP_dp_Ploidy_3D_OnlyTheLoop (p0, freq, d, ploidy, df_ci, dp_ci)

//hack to make code run both with node and browser
if(typeof require!=='undefined') {
    var jStat = require('jstat').jStat
    var BigNumber = require('bignumber.js')
    var math = require('mathjs');
}

function CalcP_dp_Ploidy_3D(p0, freq, d, ploidy, df_ci, dp_ci) {

    console.log('called with args', arguments)
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

	var types = [];
	types.push('Somatic LOH CNmut = 1') 
	if (ploidy > 1) {
        for(var i=0; i<ploidy; i++) {
        	types.push('Somatic, CNmut = ' + (i + 1))
        }
    }

	types.push('Germline LOH CNmut = 1') 
	if (ploidy > 1) {
        for(var i=0; i<ploidy; i++) {
        	types.push('Germline, CNmut = ' + (i + 1))
        }
    }

    var l = types.length
    console.log(types)

    var ddf = 0.005
	var ddp = 0.01



    //[~, df] = binofit(round(d*freq), d, df_ci);
    // fs = df(1) : ddf : df(2); 
    // fs = sort([fs(1:end-1), df(2)]);

    var binoret = binofit(d*freq, d, (1 - df_ci) * 100)
    console.log(d*freq, d, df_ci, binoret)

    var ddp2 = [p0*(1-dp_ci), p0*(1+dp_ci)]

	//fs = [freq*(1-df_ci), ddf, freq*(1+df_ci)].sort()
	//dp = [p0*(1-dp_ci) : ddp : p0*(1+dp_ci)].sort()
	var fs = sequence(binoret[0], binoret[1], ddf)
    //var fs = sequence(freq*(1-df_ci), freq*(1+df_ci), ddf)
    var dp = sequence(p0*(1-dp_ci), p0*(1+dp_ci), ddp)

    // var fs = [freq*(1-df_ci)]//sequence(freq*(1-df_ci), freq*(1+df_ci), ddf)
    // var dp = [p0*(1-dp_ci)]//sequence(p0*(1-dp_ci), p0*(1+dp_ci), ddp)

	console.log('fs ', fs, fs.length)

	console.log('dp ', dp, dp.length)

    var rectangle = [[binoret[0], binoret[1]], [ddp2[0], ddp2[1]] ]


    //new chart 
    /*

        for j=1:length(dp);
            p = dp (j);
            ff(1) = (p)/(2*(1-p)+1*p); %somatic LOH
            if (ploidy > 1)
                l = size(ff);
                for i=1:ploidy
                    ff(l+i) = (i*p)/(2*(1-p)+ploidy*p); %somatic LOH high CN
                end
            end
            l = length(ff);
            ff(l+1) = (1-p+p)/(2*(1-p)+1*p); %germline LOH high CN;        
            if (ploidy > 1)
                l = length(ff);
                for i=1:ploidy
                    ff(l+i) = (1-p+i*p)/(2*(1-p)+ploidy*p); %germline LOH high CN;
                end
            end
            all_freq(j,1:size(types, 2)) = ff;
            clear ff;
        end

        figure;

        for i=1:size(types, 2)
            pp(i) = plot (all_freq(:, i), dp, '-', 'Color', c_type{i}, 'MarkerSize', 5, 'LineWidth', 2);  
            hold on;
        end


    */
    var dp2 = sequence(0, 1, 0.01)

    all_freq = []
    for(var i=0; i<types.length; i++) {
        all_freq.push([])
    }
    console.log(all_freq)
    for(var j = 0; j < dp2.length; j++) {
        var ff = []
        p = dp2[j]
        ff.push( (p)/(2*(1-p)+1*p) )
        if (ploidy > 1) {
            for(var i=0; i<ploidy; i++) {
                ff.push((i*p)/(2*(1-p)+ploidy*p))
            }
        }
        ff.push( (1-p+p)/(2*(1-p)+1*p) ) // %germline LOH high CN; 
        if (ploidy > 1) {
            for(var i=0; i<ploidy; i++) {
                ff.push((1-p+i*p)/(2*(1-p)+ploidy*p))
            }
        }
        ff.forEach(function(fff, i) {
            all_freq[i].push(fff)
        })
        //all_freq.push(ff)

    }   

	//aics = zeros (length(fs), length(dp), size(types, 2));

	aics = zeros(fs.length, dp.length, types.length)


	//console.log(zeros(3,2,2))


	/*
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
	*/

    var ws = zeros(fs.length, dp.length, types.length)
    var aics = zeros(fs.length, dp.length, types.length)

	for(var j = 0; j < dp.length; j++) {
		var p = dp[j]
		for(var k= 0; k < fs.length; k++ ) {
			f = fs[k]

			var aic = []
			//aic[0] = 2 - 2 * Math.log( binopdf( Math.round(d*f), d, (p)/(2*(1-p)+1*p))) 

			if (ploidy > 0) {
				for(var i = 1; i <= ploidy; i++) {
					aic.push(2 - 2 * Math.log( binopdf( Math.round(d*f), d, (i * p)/(2*(1-p)+ploidy*p))) )
				}
			}

        //aic(l+1) = 2 - 2 * log (binopdf (round(d*f), d, (1-p+p)/(2*(1-p)+1*p))); %germline LOH high CN;        
			//aic.push(2 - 2 * Math.log( binopdf( Math.round(d*f), d, (1-p + p)/(2*(1-p)+1*p))) )

			if (ploidy > 0) {
				for(var i = 1; i <= ploidy; i++) {
					aic.push(2 - 2 * Math.log( binopdf( Math.round(d*f), d, (1-p + i*p)/(2*(1-p)+ploidy*p))) )
				}
			}
            //console.log('aic', aic)

			// w = zeros2D(1, types.length)

			// if(j===0) {
			// 	ws = zeros(fs.length, dp.length, types.length)
			// }

			//var not_nan = notNanArr(aic)
			var notNanAic = aic.filter(function(d) { return !isNaN(d)})

			var aicMin = getMinOfArray(notNanAic)

			var D = aic.reduce(function(acc, curr)  {
				return acc + Math.exp( -0.5 * (curr - aicMin))
			}, 0)

            //console.log('D', D)

			var w = aic.map(function(curr) {
				return Math.exp( -0.5 * (curr - aicMin))/ D
			})
            // rounding needed to match results with matlab
            w = w.map(function(curr) {
                return Math.round(curr * 1000) / 1000
            })
            // console.log(j, k)
            // console.log('w', w)
            // console.log('aic', aic)

			ws[k][j] = w

            //console.log('ws', ws)
			
			aics[k][j] = aic

            //console.log('aics', aics)

			// console.log('notNanAic', j, k)
			// console.log(aic)

			

		}
	}

    /*

        not_nan = ~isnan(aics);

        // min is a 3 dim min returning 1 number
        // sum is also 3 dim sum
        D = sum(exp(-0.5*(aics(not_nan)-min(aics(not_nan)))));
        ww = exp(-0.5*(aics-min(aics(not_nan)))) / D;
        disp(aics)
    */

    //var notNanAics = aics.filter(function(dd) { return dd.filter(function(d) { return !isNaN(d) })}) //needs 3 dims
    //console.log('notNanAics', notNanAics)
    console.log(aics.length, aics[0].length, aics[0][0].length)

    aicsMin = getMinOf3DMatrix(aics)
    console.log(aicsMin)

    var D = reduceAll3D(aics, function(acc, curr) {
        return acc + (isNaN(curr) ? 0 : Math.exp( -0.5 * (curr - aicsMin)))
    }, 0)  

    console.log('D', D)      

    var ww = map3D(aics, function(curr) {
        return Math.exp( - 0.5 * (curr - aicsMin)) / D 
    })

    console.log(ww.length, ww[0].length, ww[0][0].length)
    //console.log(ww)
    // var D = aic.reduce(function(acc, curr)  {
    //     return acc + Math.exp( -0.5 * (curr - aicMin))
    // }, 0)


    /*

        sum_ww = zeros (1, size(types, 2));
        for i=1:size(types, 2)
            w = ww(:,:,i);
            // sum_ww is the same size as type vector
            sum_ww (i) = sum(w(~isnan(w)));
            pred_outs{i} = sprintf ('%s\t%2.2e', types{i}, sum_ww (i));
            leg{i} = sprintf ('%s, w = %2.2f', types{i}, sum_ww (i));
        end
    */
    
    var sum_ww = sumOfFirst2Dims(ww)

    var pred_outs = sum_ww.map((d, i) => `${types[i]}   ${d.toFixed(2)}`)
    var leg = sum_ww.map((d, i) => `${types[i]}, w = ${d.toFixed(2)}`)

    // console.log(sum_ww)
    // console.log(pred_outs)
    // console.log(leg)

    /*
    // NEW LOOP (above code with loop instead of matrix operations)
    aics_p = zeros (length(fs), size(types, 2));
    wf = zeros (length(fs), length(dp), size(types, 2));
    for j=1:length(dp)
        for i=1:size(types, 2)
            for k=1:length(fs)
                aics_p(k, i) = aics(k,j,i);
            end
        end

        for k=1:length(fs)
            for i=1:size(types, 2)
                aic_f(i) = aics_p(k,i);
            end
            not_nan = ~isnan(aic_f);
            D = sum(exp(-0.5*(aic_f(not_nan)-min(aic_f(not_nan)))));
            wf(k, j, 1:size(types, 2)) = exp(-0.5*(aic_f-min(aic_f(not_nan)))) / D;
        end
    end
*/
    var aics_p = zeros2D(fs.length, types.length)
    var wf = zeros(fs.length, dp.length, types.length)
    var aic_f = []
    //console.log('wf', wf)
    for(var j = 0; j<dp.length; j++) {
        
        for(var i = 0; i< types.length; i++) {
            for(var k = 0; k < fs.length; k++) {
                aics_p[k][i] = aics[k][j][i]
            }
        }
       // console.log('aics_p', aics_p)
        for(var k = 0; k < fs.length; k++) {
            for(var i = 0; i< types.length; i++) {
                aic_f[i] = aics_p[k][i]
            }

            var D = aic_f.reduce(function(acc, curr) {
                return acc + (isNaN(curr) ? 0 : Math.exp( -0.5 * (curr - aicsMin)))
            }, 0) 

            //console.log('k D', k, D)
            wf_element = aic_f.map(function(curr) {
                return Math.exp( - 0.5 * (curr - aicsMin)) / D 
            })
            //console.log('wf_element', wf_element)
            wf[k][j] = wf_element

            

        }
    }
    //console.log('wf after', wf)

    console.log('wf', wf.length, wf[0].length, wf[0][0].length)

    /*
    [~, i] = sort (sum_ww, 'descend');
    fprintf ('%i\t%s\t%2.2e\t%s\t%2.2e\n', round(p0*100), types{i(1)}, sum_ww(i(1)), types{i(2)}, sum_ww(i(2)));

    pred_out = sprintf ('%s (%2.2f)', types{i(1)}, sum_ww(i(1)));
    if (max(sum_ww) < 0.99) 
        pred_out = sprintf ('%s, %s (%2.2f)', pred_out, types{i(2)}, sum_ww(i(2)));
    end
*/

/*
for k=1:size(types, 2)
    up = min(mean(wf(:,:,k), 2)+std(wf(:,:,k), 0, 2), 1)-mean(wf(:,:,k), 2);
    down = -max(mean(wf(:,:,k), 2)-std(wf(:,:,k), 0, 2), 0)+mean(wf(:,:,k), 2);
    %disp('up')
    %disp(up)
    %disp('down')
    %disp(down)
    ps = shadedErrorBar (fs, mean(wf(:,:,k), 2), [up down], '.', Colors{k}); 
    if (size(ps.patch,1))
        ps.patch.FaceColor = Colors{k}; ps.edge(1).Color = Colors{k}; ps.edge(2).Color = Colors{k};
    end
    hold on;
    pp(k) = plot (fs, mean(wf(:,:,k), 2), 'o-', 'Color', Colors{k}, 'MarkerSize', 5, 'LineWidth', 2);  
    set(pp(k),'MarkerEdgeColor','k','MarkerFaceColor', Colors{k});
    hold on; 
end
*/
    //var means = means = wf.map(dim1 => dim1.map(dim2 => { console.log(dim2); return d3.mean(dim2)}))
    var plotData = []
    for(var k = 0; k< types.length; k++) {
        //wfmat = math.matrix(wf)
        //newMat = wfmat.subset([0, 16], [0, 8], k)
        var kSlice = slice3DDim3(wf, k)
        console.log(kSlice)
        var means2ndDim = kSlice.map(function(dim1) { return d3.mean(dim1)})
        var std2ndDim = kSlice.map(function(dim1) { return d3.deviation(dim1)})
        var up = means2ndDim.map(function(mean, i) {
            return Math.min(means2ndDim[i] + std2ndDim[i], 1) //- mean
        })
        var down = means2ndDim.map(function(mean, i) {
            return Math.max(means2ndDim[i] - std2ndDim[i], 0) //+ mean
        })
        plotData.push({
            up: up,
            down: down,
            line: means2ndDim,
            lable: leg[k]
        })
        //console.log(wf)
        //var data = sumOfFirst2Dims(wf).map(function(d) { return d / types.length})
        //console.log('data', data)


    }
    console.log(plotData)


    return [fs, plotData, all_freq, rectangle]


}// what is this?


//sumOfFirst2Dims([[[1, 2, 4],[2, 3, 5]],[[4, 0, 2],[-1, 4, 3]], [[1, 0, 1],[1, 1, 1]], [[0, 0, 2],[1, 1, 3]]])
function sumOfFirst2Dims(mat) {
    var dim1 = mat.length
        dim2 = mat[0].length
        dim3 = mat[0][0].length

    var sum = new Array(dim3)
    for(var k=0; k< dim3; k++) {
        sum[k] = 0
        for(var i=0; i< dim1; i++) {
            for(var j=0; j< dim2; j++) {
                sum[k] = sum[k] + mat[i][j][k]
            }
        }

    }
    return sum;  

}

var mat = [[[1, 2, 4],[2, 3, 5]],[[4, 0, 2],[-1, 4, 3]], [[1, 0, 1],[1, 1, 1]], [[0, 0, 2],[1, 1, 3]]]
var aa = slice3DDim3(mat, 0)
console.log(mat)
console.log(aa)
function slice3DDim3(mat, dim3slice) {
    var ret = []
    mat.forEach(function(dim1, i) {
        ret[i] = []
        dim1.forEach(function(dim2, j) {
            dim2.forEach(function(el, k) {
                if(k===dim3slice) {
                    ret[i][j] = el
                }
            })
        })
    })
    return ret;

}

function size3D(mat) {
    return [mat.length, mat[0].length, mat[0][0].length]
}

 

//console.log(notNanArr([-2, -1, 0, 1, 2].map(d => 0/d)))
function notNanArr(arr) {
	return arr.map(function(d) { return !isNaN(d)})
}

function sequence(from, to, step) {
	var length = Math.abs(to - from) / step + 1
	var ret = []
	for(var i = 0; i< length; i++) {
		ret.push(from + i * step)
	}
	return ret
}

//console.log(zeros(2,3,2))
function zeros(xSize, ySize, zSize) {

	var ret = new Array(xSize);
	
	for(var i = 0; i< xSize; i++)  {
		ret[i] = new Array(ySize)
		for(var j = 0; j< ySize; j++) {
			ret[i][j] = new Array(zSize)
			for(var k = 0; k < zSize; k++) {
				ret[i][j][k] = 0
			}
		}

	}

	return ret
}

//console.log(zeros2D(2,3))
function zeros2D(xSize , ySize) {

	var ret = new Array(xSize);
	
	for(var i = 0; i< xSize; i++)  {
		ret[i] = new Array(ySize)
		for(var j = 0; j< ySize; j++) {
			ret[i][j] = 0
		}

	}

	return ret

}

function binopdf(k, n, p) {
    //use BigNumber to avoid underflowing (e.g. (1-p)**(n-k) == 0 where 1-p is small and n-k is large)
    var pow1 = new BigNumber((1-p).toString()).toPower(n - k)
    return (p === 0 || p === 1) ?
      ((n * p) === k ? 1 : 0) :
      pow1.times(Math.pow(p, k).toString()).times(jStat.combination(n, k).toString()).toNumber();

}



//console.log(binopdf(8, 12, 95))
function binofit(vx, vN, p) {
	var ret = []

	var tails = calcTails(p)
	var vTL=tails[0];
	var vTU=tails[1];// var vCL=95

	var vP = vx/vN
    if(vx==0) { 
    	ret[0] = v
    } else { 
    	var v=vP/2;
    	vsL=0;
    	vsH=vP;
    	var p=vTL/100
        while((vsH-vsL)>1e-5) {
        	if(BinP(vN,v,vx,vN)>p) {
        	 	vsH=v;
        	 	v=(vsL+v)/2 
        	} else {
        		vsL=v;
        		v=(v+vsH)/2
        	}
        }
        ret[0] = v
    }
    if(vx==vN) {
    	 ret[1] = 1
    } else { 
    	var v=(1+vP)/2;
    	vsL=vP;
    	vsH=1; 
    	var p=vTU/100
        while((vsH-vsL)>1e-5) { 
        	if(BinP(vN,v,0,vx)<p) { 
        		vsH=v; 
        		v=(vsL+v)/2 
        	} else {
        	 vsL=v; 
        	 v=(v+vsH)/2 
        	} 
        }
        ret[1] = v
    }
    return ret
 
}

function calcTails(vCL) {
    //vCL = eval(form.CL.value)
    vTU = (100-vCL)/2
    vTL = vTU
    return[vTU, vTL]
}

function BinP(N,p,x1,x2) {
    var q=p/(1-p);
    var k=0;
    var v = 1;
    var s=0;
    var tot=0
    while(k<=N) {
        tot=tot+v
        if(k>=x1 & k<=x2) {
        	s=s+v
        }
        if(tot>1e30) {
        	s=s/1e30;
        	tot=tot/1e30; 
        	v=v/1e30
        }
        k=k+1;
        v=v*q*(N+1-k)/k
    }
    return s/tot
}

function getMinOfArray(numArray) {
  return Math.min.apply(null, numArray);
}

function getMinOf3DMatrix(mat) {
   var flatten = [].concat.apply([], [].concat.apply([], mat))
   return getMinOfArray(flatten)
}

function reduceAll3D(mat, fn, init) {
   var flatten = [].concat.apply([], [].concat.apply([], mat))
   return flatten.reduce(fn, init ? init : 0)
}

function map3D(mat, fn) {
    return mat.map(function(dim0) {
        return dim0.map(function(dim1) {
            return dim1.map(function(el) {
                return fn(el)
            })
        })
    })

}



//getMinOf3DMatrix([[[1, 2],[2, 3]],[[4, 0],[-1, 4]]])
//console.log(map3D([[[1, 2],[2, 3]],[[4, 0],[-1, 4]]], function(d) { return d+1}))

//console.log(reduceAll3D([[[1, 2],[2, 3]],[[4, 0],[-1, 4]]], function(acc, curr) { return acc + curr}, 0))

//CalcP_dp_Ploidy_3D(0.6, 0.31, 1000, 2, 0.01, 0.05)