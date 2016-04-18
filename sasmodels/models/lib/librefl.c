double err_mod_func(double n_sub, double ind, double nu);
double linearfunc(double n_sub, double ind, double nu);
double power_r(double n_sub, double ind, double nu);
double power_l(double n_sub, double ind, double nu);
double exp_r(double n_sub, double ind, double nu);
double exp_l(double n_sub, double ind, double nu);
double intersldfunc(int fun_type, double n_sub, double i, double nu, double sld_l, double sld_r);
double interfunc(int fun_type, double n_sub, double i, double sld_l, double sld_r);

// normalized and modified erf
//   |
// 1 +                __  - - - -
//   |             _
//	 |            _
//   |        __
// 0 + - - -
//   |-------------+------------+--
//   0           center       n_sub    --->
//                                     ind
//
// n_sub = total no. of bins(or sublayers)
// ind = x position: 0 to max
// nu = max x to integration
double err_mod_func(double n_sub, double ind, double nu)
{
  double center, func;
  if (nu == 0.0)
		nu = 1.0e-14;
	if (n_sub == 0.0)
		n_sub = 1.0;


	//ind = (n_sub-1.0)/2.0-1.0 +ind;
	center = n_sub/2.0;
	// transform it so that min(ind) = 0
	ind -= center;
	// normalize by max limit
	ind /= center;
	// divide by sqrt(2) to get Gaussian func
	nu /= sqrt(2.0);
	ind *= nu;
	// re-scale and normalize it so that max(erf)=1, min(erf)=0
	func = erf(ind)/erf(nu)/2.0;
	// shift it by +0.5 in y-direction so that min(erf) = 0
	func += 0.5;

	return func;
}
double linearfunc(double n_sub, double ind, double nu)
{
  double bin_size, func;
	if (n_sub == 0.0)
		n_sub = 1.0;

	bin_size = 1.0/n_sub;  //size of each sub-layer
	// rescale
	ind *= bin_size;
	func = ind;

	return func;
}
// use the right hand side from the center of power func
double power_r(double n_sub, double ind, double nu)
{
  double bin_size,func;
	if (nu == 0.0)
		nu = 1.0e-14;
	if (n_sub == 0.0)
		n_sub = 1.0;

	bin_size = 1.0/n_sub;  //size of each sub-layer
	// rescale
	ind *= bin_size;
	func = pow(ind, nu);

	return func;
}
// use the left hand side from the center of power func
double power_l(double n_sub, double ind, double nu)
{
  double bin_size, func;
	if (nu == 0.0)
		nu = 1.0e-14;
	if (n_sub == 0.0)
		n_sub = 1.0;

	bin_size = 1.0/n_sub;  //size of each sub-layer
	// rescale
	ind *= bin_size;
	func = 1.0-pow((1.0-ind),nu);

	return func;
}
// use 1-exp func from x=0 to x=1
double exp_r(double n_sub, double ind, double nu)
{
  double bin_size, func;
	if (nu == 0.0)
		nu = 1.0e-14;
	if (n_sub == 0.0)
		n_sub = 1.0;

	bin_size = 1.0/n_sub;  //size of each sub-layer
	// rescale
	ind *= bin_size;
	// modify func so that func(0) =0 and func(max)=1
	func = 1.0-exp(-nu*ind);
	// normalize by its max
	func /= (1.0-exp(-nu));

	return func;
}

// use the left hand side mirror image of exp func
double exp_l(double n_sub, double ind, double nu)
{
  double bin_size, func;
	if (nu == 0.0)
		nu = 1.0e-14;
	if (n_sub == 0.0)
		n_sub = 1.0;

	bin_size = 1.0/n_sub;  //size of each sub-layer
	// rescale
	ind *= bin_size;
	// modify func
	func = exp(-nu*(1.0-ind))-exp(-nu);
	// normalize by its max
	func /= (1.0-exp(-nu));

	return func;
}

// To select function called
// At nu = 0 (singular point), call line function
double intersldfunc(int fun_type, double n_sub, double i, double nu, double sld_l, double sld_r)
{
	double sld_i, func;
	// this condition protects an error from the singular point
	if (nu == 0.0){
		nu = 1.0e-13;
	}
	// select func
	switch(fun_type){
		case 1 :
			func = power_r(n_sub, i, nu);
			break;
		case 2 :
			func = power_l(n_sub, i, nu);
			break;
		case 3 :
			func = exp_r(n_sub, i, nu);
			break;
		case 4 :
			func = exp_l(n_sub, i, nu);
			break;
		case 5 :
			func = linearfunc(n_sub, i, nu);
			break;
		default:
			func = err_mod_func(n_sub, i, nu);
			break;
	}
	// compute sld
	if (sld_r>sld_l){
		sld_i = (sld_r-sld_l)*func+sld_l; //sld_cal(sld[i],sld[i+1],n_sub,dz,thick);
	}
	else if (sld_r<sld_l){
		func = 1.0-func;
		sld_i = (sld_l-sld_r)*func+sld_r; //sld_cal(sld[i],sld[i+1],n_sub,dz,thick);
	}
	else{
		sld_i = sld_r;
	}
	return sld_i;
}


// used by refl.c
double interfunc(int fun_type, double n_sub, double i, double sld_l, double sld_r)
{
	double sld_i, func;
	switch(fun_type){
		case 0 :
			func = err_mod_func(n_sub, i, 2.5);
			break;
		default:
			func = linearfunc(n_sub, i, 1.0);
			break;
	}
	if (sld_r>sld_l){
		sld_i = (sld_r-sld_l)*func+sld_l; //sld_cal(sld[i],sld[i+1],n_sub,dz,thick);
	}
	else if (sld_r<sld_l){
		func = 1.0-func;
		sld_i = (sld_l-sld_r)*func+sld_r; //sld_cal(sld[i],sld[i+1],n_sub,dz,thick);
	}
	else{
		sld_i = sld_r;
	}
	return sld_i;
}
