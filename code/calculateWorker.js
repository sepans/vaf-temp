importScripts('loop.js', 'lib/bignumber.min.js', 'lib/math.min.js', 'lib/jstat.min.js', 'lib/d3.min.js'); 

onmessage = function(e) {
	console.log('Message received from main script', e.data);
	var defaults = e.data[0]
	
  	var workerResult = CalcP_dp_Ploidy_3D.apply(this, defaults)
  	console.log('Posting message back to main script');
  	postMessage(workerResult);
}