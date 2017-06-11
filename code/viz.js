var plotData = CalcP_dp_Ploidy_3D(0.6, 0.31, 1000, 2, 0.01, 0.05)

var mockdata = [
	{date: 20111001, 'New York': 63.8 , 'San Francisco': 62.7 , 'Austin': 72.2},
	{date: 20111002, 'New York': 58.8 , 'San Francisco': 60.7 , 'Austin': 71.2},
	{date: 20111003, 'New York': 53.8 , 'San Francisco': 52.7 , 'Austin': 65.2},
	{date: 20111004, 'New York': 55.8 , 'San Francisco': 62.7 , 'Austin': 45.2},
	{date: 20111005, 'New York': 57.8 , 'San Francisco': 62.7 , 'Austin': 74.2},
]

drawChart(plotData)

function drawChart(data) {

	console.log(data)
	var svg = d3.select("svg"),
	    margin = {top: 20, right: 80, bottom: 30, left: 50},
	    width = svg.attr("width") - margin.left - margin.right,
	    height = svg.attr("height") - margin.top - margin.bottom,
	    g = svg.append("g").attr("transform", "translate(" + margin.left + "," + margin.top + ")");

	var parseTime = d3.timeParse("%Y%m%d");

	data.map(d => {
		return d.updown = d.up.concat( d.down)
	})

	var x = d3.scaleLinear().range([0, width]),//d3.scaleTime().range([0, width]),
	    y = d3.scaleLinear().range([height, 0]),
	    z = d3.scaleOrdinal(d3.schemeCategory10);

	var dataLenght = data[0].line.length

	var line = d3.line()
	    .curve(d3.curveBasis)
	    .x(function(d, i) { return x(i); })
	    .y(function(d) { return y(d); });

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
	  x.domain([0, dataLenght]);

	  y.domain([d3.min(data, d => d3.min(d.up)), d3.max(data, d => d3.max(d.down))]);

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
	      .attr("d", function(d, i) { console.log(d, i); return line(d.line); })
	      .style("stroke", function(d, i) { return z(i); });

	  city.append("path")
	      .attr("class", "line")
	      .attr("d", function(d, i) { console.log(d, i); return line(d.up); })
	      .style('stroke-dasharray', '3, 3')
	      .style("stroke", function(d, i) { return z(i); });

	  city.append("path")
	      .attr("class", "line")
	      .attr("d", function(d, i) { console.log(d, i); return line(d.down); })
	      .style('stroke-dasharray', '3, 3')
	      .style("stroke", function(d, i) { return z(i); });

	  city.append("text")
	      .datum(function(d) { console.log(d); return d})
	      .attr("transform", function(d, i) { console.log(d, i); return "translate(" + x(d.line.length-3) + "," + (i*15) + ")"; })
	      .attr("x", 3)
	      .attr("dy", "0.35em")
	      .style("font", "10px sans-serif")
	      .style("fill", function(d, i) { return z(i); })
	      .text(function(d) { console.log(d); return d.lable; });
	      

}


function type(d, _, columns) {
  d.date = parseTime(d.date);
  for (var i = 1, n = columns.length, c; i < n; ++i) d[c = columns[i]] = +d[c];
  return d;
}