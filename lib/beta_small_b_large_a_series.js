'use strict';

// MODULES //

var abs = require( 'math-abs' );
var factorial = require( 'math-factorial' );
var gammainc = require( 'math-gammainc');
var ln = require( 'math-ln' );
var log1p = require( 'math-log1p' );
var pow = require( 'math-power' );


// FUNCTIONS //

var full_igamma_prefix = require( './full_igamma_prefix.js' );
var regularised_gamma_prefix = require( './regularised_gamma_prefix.js' );
var tgamma_delta_ratio = require( './tgamma_delta_ratio.js');


// CONSTANTS //

var MIN_VALUE = require( 'const-min-float64' );
var EPSILON = 2.220446049250313e-16;


// BGRAT ROUTINE //

/**
* FUNCTION: beta_small_b_large_a_series( a, b, x, y, s0, mult, normalised )
*	This is DiDonato and Morris's BGRAT routine, see Eq's 9 through 9.6. It calculates I_x(a,b) for the case `a > b`.
*
* @param {Number} a - function parameter
* @param {Number} b - function parameter
* @param {Number} x - function input
* @param {Number} y - equal to `1 - x`
* @param {Number} s0 - initial value
* @param {Number} mult - multiplication term of series
* @param {Boolean} normalised - boolean indicating whether to evaluate the regularized or non-regularized incomplete beta function
* @returns {Number} function value
*/
function beta_small_b_large_a_series( a, b, x, y, s0, mult, normalised ) {
	var bm1;
	var lx, lx2;
	var prefix;
	var h;
	var j;
	var m;
	var n;
	var p;
	var r;
	var t;
	var u;
	var sum;
	var lxp;
	var tnp1;
	var tmp1;
	var t4;
	var b2n;
	var mbn;

	bm1 = b - 1;
	t = a + bm1 / 2;
	if (y < 0.35) {
		lx = log1p(-y);
	} else {
		lx = ln(x);
	}
	u = -t * lx;
	// And from from 9.2:
	h = regularised_gamma_prefix( b, u );
	if( h <= MIN_VALUE ) {
		return s0;
	}
	if ( normalised ) {
		prefix = h / tgamma_delta_ratio( a, b );
		prefix /= pow(t, b);
	} else {
		prefix = full_igamma_prefix( b, u ) / pow( t, b );
	}
	prefix *= mult;

	/*
		Now we need the quantity Pn, unfortunatately this is computed
		recursively, and requires a full history of all the previous values
		so no choice but to declare a big table and hope it's big enough...
	*/
	p = new Array( 30 );
	p[ 0 ] = 1;  // see 9.3.

	// Now an initial value for J, see 9.6:
	// Call: gammainc( x, s, regularized, upper )
	j = gammainc( u, b, true, true );
	j /= h;

	// Now we can start to pull things together and evaluate the sum in Eq 9:
	sum = s0 + prefix * j;  // Value at N = 0
	// Some variables we'll need:
	tnp1 = 1; // 2*N+1
	lx2 = lx / 2;
	lx2 *= lx2;
	lxp = 1;
	t4 = 4 * t * t;
	b2n = b;

	for( n = 1; n < p.length; ++n ) {
		// Begin by evaluating the next Pn from Eq 9.4:
		tnp1 += 2;
		p[n] = 0;
		mbn = b - n;
		tmp1 = 3;
		for( m = 1; m < n; ++m ) {
			mbn = m * b - n;
			p[n] += mbn * p[n-m] / factorial(tmp1);
			tmp1 += 2;
		}
		p[n] /= n;
		p[n] += bm1 / factorial(tnp1);
		// Now we want Jn from Jn-1 using Eq 9.6:
		j = ( b2n * (b2n + 1) * j + (u + b2n + 1) * lxp ) / t4;
		lxp *= lx2;
		b2n += 2;
		// Pull it together with Eq 9:
		r = prefix * p[n] * j;
		sum += r;
		if ( r > 1 ) {
			if( abs(r) < abs(EPSILON * sum) ) {
				break;
			}
		} else {
			if( abs(r / EPSILON) < abs(sum) ) {
				break;
			}
		}
	}
	return sum;
} // end FUNCTION beta_small_b_large_a_series()


// EXPORTS //

module.exports = beta_small_b_large_a_series;
