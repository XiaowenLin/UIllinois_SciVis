/**
 * get input from UI
 * 
 * */

//-------------------------------------------------------
// Global variables

var x_extent=[-1.0,1.0];
var y_extent=[-1.0,1.0];

//------------------------------------------------------
//MAIN
function main() {
	render();
}

//--Function: render-------------------------------------
//Main drawing function

function render(canvas){
	var L = parseFloat(document.getElementById("L").value);
	var canvas = document.getElementById('example');
	
	// Get the rendering context for 2DCG <- (2)
	var ctx = canvas.getContext('2d');

	// Draw the scalar data using an image representation
	var imgData=ctx.getImageData(0, 0, canvas.width, canvas.height);
	
	
	// Choose the scalar function
	var scalar_func = gaussian;
	//if (document.getElementById("show_hedgehog").checked)
	//  scalar_func = sin2D;
	//Determine the data range...useful for the color mapping
	//1. map to data range
	//2. use get value
	//3. get data range, mn and mx
	var mn = scalar_func(pixel2pt(canvas.width, canvas.height, x_extent, y_extent, 0, 0));
	var mx = mn;
	/* iterate thru pixels */
	for (var y=0; y<canvas.height; y++)
		for (var x=0; x<canvas.width; x++)
	  	{
			var pt = pixel2pt(canvas.width, canvas.height, x_extent, y_extent, x, y);
	  		var fval = scalar_func(pt);
	  		
	  		if (fval < mn)
	  			mn=fval;
	  		if (fval>mx)
	  			mx=fval;
	  	}
	
	var noise_texture = get_noise_texture(canvas, Math.random); 
	var texture_after_lic = get_lic(noise_texture, L, lic_at_pt, gaussian, gaussian_grad);
	var color_func = greyscale_map;
	for (var y=0;y<canvas.height;y++)
		for (var x=0;x<canvas.width;x++)
	  	{
	  		var fval = texture_after_lic[y][x];
	  		var color = color_func(fval,mn,mx);
	  		i = (y*canvas.width + x) * 4;
	  		imgData.data[i]=color[0];
	  		imgData.data[i+1]= color[1];
	  		imgData.data[i+2]= color[2];
	  		imgData.data[i+3]= color[3];
	     }

	ctx.putImageData(imgData,0,0);

  // Draw the grid if necessary
	if (document.getElementById("show_hedgehog").checked) {
		var grid_size = parseFloat(document.getElementById("grid_size").value);
		var k = parseFloat(document.getElementById("k").value);
		var sample_func = sample_uniform_grid;
		if (document.getElementById("sampling_r").checked)
			sample_func = sample_random;
		render_hedgehog(ctx, grid_size, k, sample_func, gaussian_grad, canvas);
	}
}

//--------------------------------------------------------
// Map a point in pixel coordinates to the 2D function domain
function pixel2pt(width,height,x_extent,y_extent, p_x,p_y){
	var pt = [0,0];
	xlen=x_extent[1]-x_extent[0]
	ylen=y_extent[1]-y_extent[0]
	pt[0]=(p_x/width)*xlen + x_extent[0];
	pt[1]=(p_y/height)*ylen + y_extent[0];
	return pt;
	}

//--------------------------------------------------------
// Map a point in domain coordinates to pixel coordinates
function pt2pixel(width,height,x_extent,y_extent, p_x,p_y){
	var pt = [0,0];

	var xlen = (p_x-x_extent[0])/(x_extent[1]-x_extent[0]);
  var ylen = (p_y-y_extent[0])/(y_extent[1]-y_extent[0]);

	pt[0]=Math.round(xlen*width);
	pt[1]=Math.round(ylen*height);
	return pt;
	}

//contours string to array
function string2arr(string){
	var arr = string.split(',');
	for (var i = 0; i < arr.length; i++)
		arr[i] = arr[i].trim();
	return arr;
}

function linearInterp(p_a, p_b, val_a_ind, val_b_ind, five_val, borders, contour){
	var t = (five_val[val_a_ind] - contour) / 1.0 / (five_val[val_a_ind] - five_val[val_b_ind]);
	var x = borders[p_a[0]] + t * (borders[p_b[0]] - borders[p_a[0]]);
	var y = borders[p_a[1]] + t * (borders[p_b[1]] - borders[p_a[1]]);
	return([x, y]);
}

function drawALine(ctx, width, col, p0, p1){
	ctx.beginPath();              
	ctx.lineWidth = width;
	ctx.strokeStyle = col;  // Green path
	ctx.moveTo(p0[0], p0[1]);
	ctx.lineTo(p1[0], p1[1]);
	ctx.stroke();  // Draw it
}

// 2-positive code to 2 codes 
function splitCode(bts){
	var codes = ['', ''];
	var inds = [];
	for (var i = 0; i < bts.length; i++)
		if (bts[i] == 1){
			inds.push(i);
			console.log(i);
		};
	for (var i = 0; i < inds.length; i++)
		for (var j = 0; j < 5; j++)
		if (inds[i] == j){
			codes[i] = codes[i] + '1';
		}else{
			codes[i] = codes[i] + '0';
		};
	return codes;
}

function marchingSqaure(bts, five_val, borders, canvas, contour, wid){
	//mark corners of rectangular as: lt: 0, rt: 1, lb: 2, rb: 3, cent: 4
	//mark edges of rectangular as: left:0, right:1, top:2, bot:3
	var mapCodeToEdge_One = {
		'1000': [[0,2], [1,2], [0,2], [0,3], 0, 1, 0, 2],
		'0100': [[0,2], [1,2], [1,2], [1,3], 0, 1, 1, 3],
		'0010': [[0,2], [0,3], [0,3], [1,3], 0, 2, 2, 3],
		'0001': [[1,2], [1,3], [0,3], [1,3], 1, 3, 2, 3]
	};
	var mapCodeToEdge_Three = {
		'0111': [[0,2], [1,2], [0,2], [0,3], 0, 1, 0, 2],
		'1011': [[0,2], [1,2], [1,2], [1,3], 0, 1, 1, 3],
		'1101': [[0,2], [0,3], [0,3], [1,3], 0, 2, 2, 3],
		'1110': [[1,2], [1,3], [0,3], [1,3], 1, 3, 2, 3]
	};
	var mapCodeToEdge_Two_SameSide = {
		'1100': [[0,2], [0,3], [1,2], [1,3], 0, 2, 1, 3],
		'0101': [[1,2], [0,2], [0,3], [1,3], 1, 0, 2, 3],
		'0011': [[0,2], [0,3], [1,2], [1,3], 0, 2, 1, 3],
		'1010': [[0,2], [1,2], [0,3], [1,3], 0, 1, 2, 3]
	};
//	var mapEdgeToPoints_one ={
//		'1010': [[0,0], []]
//			
//	}
	var ctx = canvas.getContext('2d');
	var code = bts.reduce(function(a, b){return(''+a+b)});
	code = code.substr(0,4);
	var sm = 0;
	for (var i = 0; i < 4; i++)
		sm = sm + bts[i];
	//is sm is 0, pass
	if (sm == 1){
		//find the two points
		var edges = mapCodeToEdge_One[code];
		var p0 = linearInterp(edges[0], edges[1], edges[4], edges[5], five_val, borders, contour);
		var p1 = linearInterp(edges[2], edges[3], edges[6], edges[7], five_val, borders, contour);
		//draw
		drawALine(ctx, wid, "green", p0, p1);
	}else if (sm == 2) {
		//if two points are on the same side
		if (typeof(mapCodeToEdge_Two_SameSide[code]) != "undefined"){
			var edges = mapCodeToEdge_Two_SameSide[code];
			var p0 = linearInterp(edges[0], edges[1], edges[4], edges[5], five_val, borders, contour);
			var p1 = linearInterp(edges[2], edges[3], edges[6], edges[7], five_val, borders, contour);
			drawALine(ctx, wid, "green", p0, p1);
		}else{
			//if two points are diagonal
			// if center is green
			if (bts[4] == 1){
				var r_bts = bts.map(function(x){return((x == 1) ? 0 : 1)});
				var codes = splitCode(r_bts);
				for (var i = 0; i < codes.length; i++){
					code = codes[i].substr(0,4);
					var edges = mapCodeToEdge_One[code];
					var p0 = linearInterp(edges[0], edges[1], edges[4], edges[5], five_val, borders, contour);
					var p1 = linearInterp(edges[2], edges[3], edges[6], edges[7], five_val, borders, contour);
					drawALine(ctx, wid, "green", p0, p1);
				}
			}else{
				//if center is white
				var codes = splitCode(bts);
				for (var i = 0; i < codes.length; i++){
					code = codes[i].substr(0,4);
					var edges = mapCodeToEdge_One[code];
					var p0 = linearInterp(edges[0], edges[1], edges[4], edges[5], five_val, borders, contour);
					var p1 = linearInterp(edges[2], edges[3], edges[6], edges[7], five_val, borders, contour);
					drawALine(ctx, wid, "green", p0, p1);
				}
			}		
		}
	}else if (sm == 3) {
		//find the points
		code = code.substr(0,4);
		var edges = mapCodeToEdge_Three[code];
		var p0 = linearInterp(edges[0], edges[1], edges[4], edges[5], five_val, borders, contour);
		var p1 = linearInterp(edges[2], edges[3], edges[6], edges[7], five_val, borders, contour);
		//draw
		drawALine(ctx, wid, "green", p0, p1);
	}
}

function cellDrawContour(five_val, borders, canvas, contour){
	var bts = [];
	for (var i = 0; i < 5; i++)
		bts.push( (five_val[i] >= contour) + 0);
	marchingSqaure(bts, five_val, borders, canvas, contour, "4");
}

function drawContour(contour, canvas, res, scalar_func){
	//get pixel locs
	var ctx = canvas.getContext('2d');
	loc=[0,0];
	delta_x = canvas.width/res;
	delta_y = canvas.height/res;
	for (var i=0;i<res;i++)
		for (var j=0;j<res;j++)
		{
			//get pixel locs
			var left = loc[0] + i * delta_x;//min x
			var right = left + delta_x;//max x
			var top = loc[1] + j * delta_y;
			var bot = top + delta_y;
			//get four points, get mapped locs and get fval
			var lt_fval = scalar_func(pixel2pt(canvas.width,canvas.height,x_extent,y_extent,left,top));
			var rt_fval = scalar_func(pixel2pt(canvas.width,canvas.height,x_extent,y_extent,right,top));
			var lb_fval = scalar_func(pixel2pt(canvas.width,canvas.height,x_extent,y_extent,left,bot));
			var rb_fval = scalar_func(pixel2pt(canvas.width,canvas.height,x_extent,y_extent,right,bot));
			var cent = scalar_func(pixel2pt(canvas.width,canvas.height,x_extent,y_extent,(left + right)/2.0,(top + bot) / 2.0));
			var five_val = [lt_fval, rt_fval, lb_fval, rb_fval, cent];
			var borders = [left, right, top, bot];
			//draw contour for each cell
			cellDrawContour(five_val, borders, canvas, contour);
		}
};

var get_noise_texture = function(canvas, noise_func) {
	var noise_arr = new Array(canvas.height);
	for (var i = 0; i < canvas.height; i++) {
	  noise_arr[i] = new Array(canvas.width);
	}
	for (var y = 0; y < canvas.height; y++) {
		for (var x = 0; x < canvas.width; x++) {
			var pt = pixel2pt(canvas.width, canvas.height, x_extent, y_extent, x, y)
	  		noise_arr[y][x] = noise_func(pt);
	    }
	}
	return noise_arr;
};

/**
 * weight_func and field_func should be in pt space
 * get all cells on trace
 * calculate weights on trace
 * get LIC value
 * i, j: 0 - 400, in space of px
 * L, in space of px
 * weight_func, field_func in the space of pt
 * */
var lic_at_pt = function(i, j, noise_texture, L, weight_func, field_grad_func){
	var l = 0;
	var cur_i = i;
	var cur_i_ = i;
	var cur_j = j;
	var cur_j_ = j;
	var i_bnd = noise_texture.length;
	var j_bnd = noise_texture[0].length;
	var pt = pixel2pt(j_bnd, i_bnd, x_extent, y_extent, cur_j, cur_i);
	var orig_pt = pt;
	var noises = [];
	var weights = [];
	while(l <= L && cur_i < i_bnd && cur_i >= 0 && cur_j < j_bnd && cur_j >=0) {
		noises.push(noise_texture[cur_i][cur_j]);
		weights.push(weight_func(pt_diff(pt, orig_pt)));
		if (pt[0] == 0 && pt[1] == 0)
			break;
		var vec = field_grad_func(pt);
		cur_j_ = cur_j_ + vec[0];
		cur_i_ = cur_i_ + vec[1];
		cur_j = Math.round(cur_j_);
		cur_i = Math.round(cur_i_);
		if (cur_j == j_bnd || cur_i == i_bnd || cur_j == 0 || cur_i == 0)
			break;
		l = Math.sqrt(Math.pow(cur_j_ - j, 2) + Math.pow(cur_i_ - i, 2));
		pt = pixel2pt(j_bnd, i_bnd, x_extent, y_extent, cur_j, cur_i);
	}
	l = 0;
	cur_i = i;
	cur_j = j;
	cur_i_ = i;
	cur_j_ = j;
	pt = pixel2pt(j_bnd, i_bnd, x_extent, y_extent, cur_j_, cur_i_);
	while(l < L && cur_i > 0 && cur_i < i_bnd && cur_j > 0 && cur_j < j_bnd) {
		if (pt[0] == 0 && pt[1] == 0)
			break;
		var vec = field_grad_func(pt);
//		var vec_px = pt2pixel(j_bnd, i_bnd, x_extent, y_extent, vec[0], vec[1]);
//		vec_px = normalize(vec_px);
		cur_j_ = cur_j_ - vec[0];
		cur_i_ = cur_i_ - vec[1];
		cur_j = Math.round(cur_j_);
		cur_i = Math.round(cur_i_);
		if (cur_j == j_bnd || cur_i == i_bnd || cur_j == 0 || cur_i == 0)
			break;
		l = Math.sqrt(Math.pow(cur_j_ - j, 2) + Math.pow(cur_i_ - i, 2));
		pt = pixel2pt(j_bnd, i_bnd, x_extent, y_extent, cur_j, cur_i);
		noises.push(noise_texture[cur_i][cur_j]);
		weights.push(weight_func(pt_diff(pt, orig_pt)));
	}
	var w_sum = 0;
	var w_noise = 0;
	for (var i = 0; i < weights.length; i++) {
		w_sum = w_sum + weights[i];
		w_noise = w_noise + weights[i] * noises[i]; 
	}
	return w_noise/w_sum;
};

var get_lic = function(noise_texture, L, lic_func, weight_func, field_grad_func){
	var lic_arr = new Array(noise_texture.length);
	for (var i = 0; i < noise_texture[0].length; i++) {
	  lic_arr[i] = new Array(noise_texture.length);
	}
	for (var i = 0; i < noise_texture.length; i++) {
		for (var j = 0; j < noise_texture[0].length; j++) {
			lic_arr[i][j] = lic_func(i, j, noise_texture, L, weight_func, field_grad_func);
	    }
	}
	return lic_arr;
};

var render_hedgehog = function(ctx, grid_size, k, sample_func, field_grad_func, canvas) {
	/* take sample */
	var samples = sample_func(grid_size, canvas);
	for (var i = 0; i < samples.length; i++) {
		/** 
		 * draw line (x, x + kv(x))
		 * */
		var wid = "4";
		var col = "red";
		var pt = pixel2pt(canvas.width, canvas.height, x_extent, y_extent, samples[i][0], samples[i][1]);
		var vec = field_grad_func(pt);
		var kv = pt_times(vec, k);
		drawALine(ctx, wid, col, samples[i], pt_add(samples[i], kv));
		
	}
	
}

var sample_uniform_grid = function(grid_size, canvas) {
	var ctx = canvas.getContext('2d');
	var loc=[0,0];
	var delta_x = canvas.width/(grid_size - 1);
	var delta_y = canvas.height/(grid_size - 1);
	var samples = [];
	for (var i = 0; i <= grid_size; i++)
		for (var j = 0; j <= grid_size; j++)
		{
			var left = loc[0] + i * delta_x;//min x
			var top = loc[1] + j * delta_y;
			var cur = [left, top];
			samples.push(cur);
		}
	return samples;
}

var sample_random = function(grid_size, canvas) {
	var n = grid_size * grid_size;
	var samples = [];
	for (var i = 0; i < n; i++){
		var px = new Array(2);
		px[0] = Math.random() * canvas.width;
		px[1] = Math.random() * canvas.height;
		samples.push(px);
	}
	return samples;
}

function rainbow_colormap(fval,fmin,fmax){
	var dx=0.8;
	var fval_nrm = (fval-fmin)/(fmax-fmin);
	var g = (6.0-2.0*dx)*fval_nrm +dx;
	var R = Math.max(0.0,(3.0-Math.abs(g-4.0)-Math.abs(g-5.0))/2.0 )*255;
	var G = Math.max(0.0,(4.0-Math.abs(g-2.0)-Math.abs(g-4.0))/2.0 )*255;
	var B = Math.max(0.0,(3.0-Math.abs(g-1.0)-Math.abs(g-2.0))/2.0 )*255;
	color = [Math.round(R),Math.round(G),Math.round(B),255];
	return color;
}

//function gray_colormap(fval,fmin,fmax){
//	var dx=0.8;
//	var fval_nrm = (fval-fmin)/(fmax-fmin);
//	var g = (6.0-2.0*dx)*fval_nrm +dx;
//	var R = Math.max(0.0,(3.0-Math.abs(g-4.0)-Math.abs(g-5.0))/2.0 )*255;
//	var G = Math.max(0.0,(3.0-Math.abs(g-4.0)-Math.abs(g-5.0))/2.0 )*255;
//	var B = Math.max(0.0,(3.0-Math.abs(g-4.0)-Math.abs(g-5.0))/2.0 )*255;
//	color = [Math.round(R),Math.round(G),Math.round(B),255];
//	return color;
//}

function greyscale_map(fval,fmin,fmax){
	  var c=255*((fval-fmin)/(fmax-fmin));
	  var color = [Math.round(c),Math.round(c),Math.round(c),255];
		return color;
	}

/**
 * math functions
 * */
var normalize = function(vec){
	var sq_sum = Math.sqrt(vec[0] * vec[0] + vec[1] * vec[1]);
	var ret = new Array(2);
	ret[0] = vec[0] / sq_sum;
	ret[1] = vec[1] / sq_sum;
	return ret;
}

var gaussian = function(pt){
	return Math.exp(-(pt[0]*pt[0]+pt[1]*pt[1]));
};

var gaussian_grad = function(pt) {
	var sq_sum = pt[0]*pt[0]+pt[1]*pt[1];
	if (sq_sum == 0){
		return [0, 0];
	}
	var ret = new Array(2);
	ret[0] = -Math.exp(-(sq_sum)) * pt[0] * Math.sqrt(sq_sum);
	ret[1] = -Math.exp(-(sq_sum)) * pt[1] * Math.sqrt(sq_sum);
	var norm = Math.sqrt(ret[0]*ret[0] + ret[1]*ret[1]);
	ret[0] = ret[0] / norm;
	ret[1] = ret[1] / norm;
	return ret;
};

var pt_diff = function(a, b){
	var ret = new Array(2);
	ret[0] = a[0] - b[0];
	ret[1] = a[1] - b[1];
	return ret;
}

var pt_add = function(a, b){
	var ret = new Array(2);
	ret[0] = a[0] + b[0];
	ret[1] = a[1] + b[1];
	return ret;
}

var pt_times = function(a, scal){
	var ret = new Array(2);
	ret[0] = scal * a[0];
	ret[1] = scal * a[1];
	return ret;
}

