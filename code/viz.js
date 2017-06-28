var variables = ['Purity', 'VAF', 'Depth', 'ploidy', 'df_ci', 'dp_ci']
var defaults = [0.6, 0.5, 1000, 2, 0.01, 0.05]
//var defaults = [0.6, 0.5, 1000, 2, 0.01, 0.05]

drawInputs()

var calculating = d3.select('#calculating')

var myWorker;

if(window.Worker) {

	myWorker = new Worker('calculateWorker.js')

	myWorker.onmessage = function(e) {
	  
	  console.log('Message received from worker', e.data);

	  drawCharts(e.data)

	}
}

var margin = {top: 20, right: 180, bottom: 40, left: 45},
    w = 500,
    h = 400,
    width = w - margin.left - margin.right,
    height = h - margin.top - margin.bottom


function drawCharts(allPlotData) {

		var fs = allPlotData[0]
			plotData = allPlotData[1]
			all_freq = allPlotData[2]
			rect = allPlotData[3]

		calculating.style('display', 'none')

		d3.select('#chart1 g').remove()
		d3.select('#chart2').selectAll('g').remove()

		console.log('lables', plotData.map(d => d.lable))

		drawChart(fs, plotData)

		drawChart2(all_freq, rect, plotData.map(d => d.lable))

}


function init() {
	//var allPlotData = [[0.2728157043457031,0.27781570434570313,0.28281570434570313,0.28781570434570314,0.29281570434570314,0.29781570434570315,0.30281570434570315,0.30781570434570316,0.3128157043457031,0.3178157043457031,0.3228157043457031,0.3278157043457031,0.3328157043457031,0.3378157043457031,0.34281570434570313,0.34781570434570314,0.35281570434570314],[{"up":[4.932558694326207e-16,6.2693874893369864e-15,7.968472357116566e-14,1.0127958370024586e-12,1.2872573557710332e-11,1.6360826375135302e-10,2.0794150908031073e-9,2.6428525356087214e-8,3.3589212498564793e-7,0.000004268914130640371,0.000054246665123777895,0.0006881795533953284,0.008550144419692492,0.08457152221995344,0.2803310748139693,0.3704916249970779,0.43859057624096603],"down":[0.9999999999999998,0.9999999999999978,0.9999999999999716,0.9999999999996383,0.9999999999953992,0.9999999999414709,0.9999999992553591,0.9999999905255105,0.9999998794406599,0.9999984658028411,0.9999804768984198,0.9997519406215629,0.9969114816340691,0.9692008838773263,0.8897394935923141,0.8065485711199569,0.7265042965899404],"line":[1.7575637498574978e-16,2.235398410227832e-15,2.843288333186359e-14,3.616686809491535e-13,4.600732931235348e-12,5.852905854046107e-11,7.446409509933829e-10,9.474489556394089e-9,1.2055934004638456e-7,0.0000015341971588910085,0.000019523101580197086,0.0002480593784371129,0.003088518365930901,0.030799116122673737,0.11026050640768588,0.1934514288800431,0.27349570341005963],"lable":"Somatic LOH CNmut = 1, w = 0.00"},{"up":[2.220446049250313e-16,2.220446049250313e-15,2.842170943040401e-14,3.61710661422876e-13,4.600764214046649e-12,5.852907047909639e-11,7.446409044931102e-10,9.474489526262175e-9,1.2055934006660607e-7,0.0000015341971588700432,0.00001952310158026549,0.00024805937843708215,0.0030885183659309368,0.03079911612267372,0.11026050640768592,0.19345142888004319,0.2734957034100596],"down":[2.220446049250313e-16,2.220446049250313e-15,2.842170943040401e-14,3.61710661422876e-13,4.600764214046649e-12,5.852907047909639e-11,7.446409044931102e-10,9.474489526262175e-9,1.2055934006660607e-7,0.0000015341971588700432,0.00001952310158026549,0.00024805937843708215,0.0030885183659309368,0.03079911612267372,0.11026050640768592,0.19345142888004319,0.2734957034100596],"line":[0.9999999999999998,0.9999999999999978,0.9999999999999716,0.9999999999996383,0.9999999999953992,0.9999999999414709,0.9999999992553591,0.9999999905255105,0.9999998794406599,0.9999984658028411,0.9999804768984197,0.9997519406215629,0.9969114816340691,0.9692008838773263,0.8897394935923141,0.8065485711199568,0.7265042965899404],"lable":"Somatic LOH, CNmut = 1, w = 1.00"},{"up":[1.5160665314820943e-79,6.166747084912657e-77,2.5083839378770797e-74,1.0203093866076319e-71,4.150206914910459e-69,1.688136707661042e-66,6.866658849253146e-64,2.793079498846213e-61,1.1361110071888485e-58,4.621189446808537e-56,1.8794470084661803e-53,7.631089986876954e-51,3.034681570214056e-48,9.614014425645878e-46,1.0261009454517183e-43,4.0223204366542207e-42,1.3102794509221379e-40],"down":[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],"line":[5.360142080382253e-80,2.180290887121158e-77,8.868549557987236e-75,3.607370722679929e-72,1.4673339919446157e-69,5.968527588423503e-67,2.427758358773598e-64,9.875150261196784e-62,4.016813320475739e-59,1.6338612004231323e-56,6.644954660147487e-54,2.698044581362309e-51,1.0729426433643694e-48,3.399158245456407e-46,3.6281862059378403e-44,1.4229152900988593e-42,4.637018208043723e-41],"lable":"Somatic LOH, CNmut = 2, w = 0.00"},{"up":[0,0,1.4150788185096167e-160,9.566054426979348e-157,6.466905198044845e-153,4.371798544172186e-149,2.9554511511916395e-145,1.9979628206249343e-141,1.3506743988141192e-137,9.130813534096082e-134,6.171796410500561e-130,4.164797482706296e-126,2.7526196351981494e-122,1.4493157302558581e-118,2.5707910283572594e-115,1.6747093441527367e-112,9.066281355726479e-110],"down":[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],"line":[1.094778540129342e-168,7.4009896044010235e-165,5.00326277680404e-161,3.382336657674529e-157,2.2865481564209075e-153,1.5457664344499621e-149,1.0449785865587027e-145,7.064328438261889e-142,4.7756668571278223e-138,3.2284398936954434e-134,2.1822007617885536e-130,1.4725733146029413e-126,9.732619480989676e-123,5.124532777622947e-119,9.091784617349588e-116,5.929847254097844e-113,3.212876902978368e-110],"lable":"Germline LOH CNmut = 1, w = 0.00"},{"up":[2.315979537959975e-45,9.93075468121381e-44,4.2112205738078276e-42,1.7415804561391529e-40,6.796671560467261e-39,2.350732382639958e-37,1.013473455657963e-35,1.2255338140091524e-33,1.4670571311385903e-31,1.584685077984498e-29,1.6349063780616288e-27,1.650906798468138e-25,1.6191496324200515e-23,1.2695175356719932e-21,4.569137956323071e-20,1.7432919873085957e-18,5.665342599589032e-17],"down":[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],"line":[1.379830409206092e-45,6.571363299127691e-44,3.2202403290144175e-42,1.6430782556591453e-40,8.864159257043992e-39,5.14256903049536e-37,3.252716734204074e-35,2.2537145752528848e-33,1.6993147678680992e-31,1.3739844193635515e-29,1.171211822528942e-27,1.036062905905292e-25,9.288023064631295e-24,7.2948374713217295e-22,3.318560550829306e-20,1.24151764298647e-18,3.932322925296916e-17],"lable":"Germline, CNmut = = 1, w = 0.00"},{"up":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"down":[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],"line":[6.698038225002441e-261,4.3192617006901255e-256,2.785296383306884e-251,1.796111576995268e-246,1.1582310654755027e-241,7.468907934355221e-237,4.816360666720559e-232,3.105852797697506e-227,2.002822146777812e-222,1.2915149622508197e-217,8.327215640132111e-213,5.360192242284173e-208,3.379333198968487e-203,1.69725502415241e-198,2.871846018872556e-194,1.78484619122592e-190,9.217992441551726e-187],"lable":"Germline, CNmut = = 2, w = 0.00"}]]
	//var allPlotData = CalcP_dp_Ploidy_3D(0.6, 0.31, 1000, 2, 0.01, 0.05)
	var allPlotData;

	if(myWorker) {
		myWorker.postMessage([defaults])
	}
	else {
		allPlotData = CalcP_dp_Ploidy_3D.apply(this, defaults)
		drawCharts(allPlotData)
	}

	 

	//console.log(JSON.stringify(allPlotData))

}

init()

function drawChart2(all_freq, rect, lables) {

	var svg = d3.select("#chart2"),
	    g = svg.append("g").attr("transform", "translate(" + margin.left + "," + margin.top + ")");

	svg.attr('width', w)
		.attr('height', h)

	var x = d3.scaleLinear().range([0, width]),//d3.scaleTime().range([0, width]),
	    y = d3.scaleLinear().range([height, 0]),
	    z = d3.scaleOrdinal(d3.schemeCategory10);
//	    z2 = d3.scaleOrdinal(d3.schemeCategory10);

	var line = d3.line()
	    .curve(d3.curveCardinal)
	    .y(function(d, i) { return y(i* 0.01); })
	    .x(function(d) { return x(d); });

	x.domain([0, 1]);

	  //y.domain([d3.min(data, d => d3.min(d.up)), d3.max(data, d => d3.max(d.down))]);
	y.domain([0, 1])

	 z.domain(all_freq.length);

	  g.append("g")
	      .attr("class", "axis axis--x")
	      .attr("transform", "translate(0," + height + ")")
	      .call(d3.axisBottom(x));

	  g.append("g")
	      .attr("class", "axis axis--y")
	      .call(d3.axisLeft(y))
	    .append("text")
	      .attr("transform", "rotate(-90)")
	      .attr("y", 6)
	      .attr("dy", "0.71em")
	      .attr("fill", "#000")
	      //.text("Temperature, ºF");
	  g.append('rect')
	  	.attr('class', 'rectangle')
	  	.attr('x', x(rect[0][0]))
	  	.attr('y', y(rect[1][1]))
	  	.attr('width',  x(rect[1][0]) - x(rect[0][0]))
	  	.attr('height', y(rect[0][1]) - y(rect[1][1]))

	  console.log('rect', rect)

	  g.append('line')
	  	.attr('class', 'cross')
	  	.attr('x1', x((rect[0][0] + rect[1][0])/2))
	  	.attr('y1', y((rect[0][1] + rect[1][1])/2) - 5)
	  	.attr('x2', x((rect[0][0] + rect[1][0])/2))
	  	.attr('y2', y((rect[0][1] + rect[1][1])/2) + 5)

	  g.append('line')
	  	.attr('class', 'cross')
	  	.attr('x1', x((rect[0][0] + rect[1][0])/2) - 5)
	  	.attr('y1', y((rect[0][1] + rect[1][1])/2))
	  	.attr('x2', x((rect[0][0] + rect[1][0])/2) + 5)
	  	.attr('y2', y((rect[0][1] + rect[1][1])/2))


	  var city = g.selectAll(".city")
	    .data(all_freq)
	    .enter().append("g")
	      .attr("class", "city");

	  city.append("path")
	      .attr("class", function(d, i) { return "line"})
	      .attr("d", function(d, i) { return line(d); })
	      .style("stroke", function(d, i) { return z(i); });	 

	  city.append("text")
	      .data(lables)
	      .attr("transform", function(d, i) { return "translate("+(width + 10)+"," + (i*15) + ")"; })
	      .attr("x", 3)
	      .attr("dy", "0.35em")
	      .style("font", "10px sans-serif")
	      .style("fill", function(d, i) { return z(i); })
	      .text(function(d) { return d; });

	    svg.append("text")
	    	.attr('class', 'axis-lable axis-lable-x')
	    	.attr("transform",
	            "translate(" + (width/2 + margin.left) + " ," + 
	                           h + ")")
	        .style("text-anchor", "middle")
	        .text('x axis');


	  // text label for the y axis
		svg.append("text")
	    	.attr('class', 'axis-lable axis-lable-y')
		    .attr("transform", "translate(10,"+height/2+") rotate(-90) ")
		    // .attr("y", 0 - margin.left)
		    // .attr("x",0 - (height / 2))
		    // .attr("dy", "1em")
		    .style("text-anchor", "middle")
		    .text("y axis"); 		    

}


function drawChart(fs, data) {

	var svg = d3.select("#chart1"),
	    g = svg.append("g").attr("transform", "translate(" + margin.left + "," + margin.top + ")");

	svg.attr('width', w)
		.attr('height', h)

	data.map(d => {
		return d.updown = d.up.map((dd, i) => { return {up: d.up[i], down: d.down[i]}}) 
	})
	
	//remove zeros and ones
	data = data.filter(d => {
		var sum = Math.round(d3.sum(d.line) * 1000) / 1000
		var zeroOrOne = sum===0 || sum===d.line.length
		console.log(d.line, sum ,zeroOrOne)
		return !zeroOrOne
	})

	console.log('data',data)
	console.log('DATA',data[0].line)

	var x = d3.scaleLinear().range([0, width]),//d3.scaleTime().range([0, width]),
	    y = d3.scaleLinear().range([height, 0]),
	    z = d3.scaleOrdinal(d3.schemeCategory10);
//	    z2 = d3.scaleOrdinal(d3.schemeCategory10);

	var dataLenght = data[0].line.length

	var line = d3.line()
	    .curve(d3.curveCardinal)
	    .x(function(d, i) { return x(fs[i]); })
	    .y(function(d) { return y(d); });


	var upArea = d3.area()
		    .curve(d3.curveBasis)
		    .x0(function(d, i) { return x(fs[i]); })
		    .x1(function(d, i) { return x(fs[i]); })
	        .y0( function(d) { return y(d.up) } )
          	.y1(  function(d) { return y(d.down) } );	    

	// var downArea = d3.area()
	// 	    .curve(d3.curveBasis)
	// 	    .x0(function(d, i) { console.log('ar', d, i); return x(fs[i]); })
	// 	    .x1(function(d, i) { return x(fs[i]); })
	//         .y0( function(d) { return y(y.domain()[1]) } )
 //          	.y1(  function(d) { return y(d.down) } );	    

	    /*
	  var cities = Object.keys(data[0]).slice(1).map(function(id) {
	    return {
	      id: id,
	      values: data.map(function(d) {
	        return {date: d.date, temperature: d[id]};
	      })
	    };
	  });
		*/

	  //x.domain(d3.extent(data, function(d) { return d.date; }));
	  x.domain(d3.extent(fs));

	  //y.domain([d3.min(data, d => d3.min(d.up)), d3.max(data, d => d3.max(d.down))]);
	  y.domain([0, 1])

	  z.domain(data.length);

	  g.append("g")
	      .attr("class", "axis axis--x")
	      .attr("transform", "translate(0," + height + ")")
	      .call(d3.axisBottom(x));

	  g.append("g")
	      .attr("class", "axis axis--y")
	      .call(d3.axisLeft(y))
	    .append("text")
	      .attr("transform", "rotate(-90)")
	      .attr("y", 6)
	      .attr("dy", "0.71em")
	      .attr("fill", "#000")
	      //.text("Temperature, ºF");

	  var city = g.selectAll(".city")
	    .data(data)
	    .enter().append("g")
	      .attr("class", "city");

	  city.append("path")
	      .attr("class", "line")
	      .attr("d", function(d, i) { return line(d.line); })
	      .style("stroke", function(d, i) { return z(i); });

	   city.append('path')
	   		.attr('class', ' area up')   
		     .attr("d", function(d, i) { return upArea(d.updown) })
	      	.style("fill", function(d, i) { return z(i); })
	      	.style("stroke", function(d, i) { return z(i); })
	      	.style("stroke-opacity", 0.2)
	      	.style("fill-opacity", 0.1);

	  city.selectAll('.circle')
	  	.data(d => d.line)
	  	.enter().append("circle")
	  	.attr('class', 'circle')
	  	.attr('r', function(d) { return 3})
	    .style("stroke", function(d, i, a) { 
	    	return z(data.indexOf(this.parentNode.__data__))
	    })
	    .attr('cx', function(d, i) { return x(fs[i]); })
	    .attr('cy', function(d) { return y(d); });

	  city.append("text")
	      .datum(function(d) { return d})
	      .attr("transform", function(d, i) { return "translate("+(width + 10)+"," + (i*15) + ")"; })
	      .attr("x", 3)
	      .attr("dy", "0.35em")
	      .style("font", "10px sans-serif")
	      .style("fill", function(d, i) { return z(i); })
	      .text(function(d) { return d.lable; });

	    svg.append("text")
	    	.attr('class', 'axis-lable axis-lable-x')
	    	.attr("transform",
	            "translate(" + (width/2 + margin.left) + " ," + 
	                           h + ")")
	        .style("text-anchor", "middle")
	        .text('VAF');


	  // text label for the y axis
		svg.append("text")
	    	.attr('class', 'axis-lable axis-lable-y')
		    .attr("transform", "translate(10,"+height/2+") rotate(-90) ")
		    // .attr("y", 0 - margin.left)
		    // .attr("x",0 - (height / 2))
		    // .attr("dy", "1em")
		    .style("text-anchor", "middle")
		    .text("W"); 	      
	      

}

function drawInputs() {
	var variableEl = d3.select('#variables')
	variableItems = variableEl.selectAll('.variable')
		.data(variables)
		.enter()
		.append('div')
		.classed('variable', true)

	variableItems.append('label')
		.text(d => d+':')

	variableItems.append('input')
		.attr('id', d => toInputId(d))
		.attr('value', (d,i) => defaults[i])

	variableEl.append('button')
		.text( 'draw')
		.on('click', function(e) {
			console.log('calculating opacity 1 on')
			calculating.style('display', 'block')

			var values = variables.map(d => +document.getElementById(toInputId(d)).value)
			console.log('values', values)
			if(myWorker) {
				myWorker.postMessage([values])
			}
			else {
				var allPlotData = CalcP_dp_Ploidy_3D.apply(this, values)

				drawCharts(allPlotData)

			}
		})

}

function toInputId(d) {
	console.log(d, 'input-'+d.toLowerCase().split(' ').join('-'))
	return 'input-'+d.toLowerCase().split(' ').join('-')
}
