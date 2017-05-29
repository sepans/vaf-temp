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

	var types = [];
	types.push('Somatic LOH CNmut = 1') 
	if (ploidy > 1) {
        for(var i=0; i<ploidy; i++) {
        	types.push('Somatic LOH, CNmut = ' + (i + 1))
        }
    }

	types.push('Germline LOH CNmut = 1') 
	if (ploidy > 1) {
        for(var i=0; i<ploidy; i++) {
        	types.push('Germline, CNmut = = ' + (i + 1))
        }
    }

    var l = types.length
    console.log(types)

    var ddf = 0.005
	var ddp = 0.01

	//fs = [freq*(1-df_ci), ddf, freq*(1+df_ci)].sort()
	//dp = [p0*(1-dp_ci) : ddp : p0*(1+dp_ci)].sort()
	var fs = sequence(freq*(1-df_ci), freq*(1+df_ci), ddf)
	var dp = sequence(p0*(1-dp_ci), p0*(1+dp_ci), ddp)

	console.log(fs.length, fs[0], fs[1], fs[fs.length-1])

	console.log(dp.length, dp[0], dp[1], dp[dp.length-1])

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

	for(var j = 0; j < dp.length; j++) {
		var p = dp[j]
		for(var k= 0; k < fs.length; k++ ) {
			f = fs[k]

			var aic = []
			aic[0] = 2 - 2 * Math.log( binopdf( Math.round(d*f), d, (p)/(2*(1-p)+1*p))) 
			if (ploidy > 1) {
				for(var i = 1; i < ploidy; i++) {
					aic.push(2 - 2 * Math.log( binopdf( Math.round(d*f), d, (i * p)/(2*(1-p)+ploidy*p))) )
				}
			}

			aic.push(2 - 2 * Math.log( binopdf( Math.round(d*f), d, (1-p + p)/(2*(1-p)+1*p))) )

			if (ploidy > 1) {
				for(var i = 1; i < ploidy; i++) {
					aic.push(2 - 2 * Math.log( binopdf( Math.round(d*f), d, (1-p + i*p)/(2*(1-p)+ploidy*p))) )
				}
			}

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

			var w = aic.map(function(curr) {
				return Math.exp( -0.5 * (curr - aicMin))/ D
			})

			ws[k][j][types.length] = w
			
			aics[k][j][types.length] = aic

			// console.log('notNanAic', j, k)
			// console.log(aic)

			

		}
	}
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



//console.log(binopdf(8, 12, 95))
function binopdf(vx, vN, p) {
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




CalcP_dp_Ploidy_3D(1, 1, 3, 5, 6, 7)