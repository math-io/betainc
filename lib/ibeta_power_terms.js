'use strict';

// MODULES //

var abs = require( 'math-abs' );
var exp = require( 'math-exp' );
var ln = require( 'math-ln' );
var log1p = require( 'math-log1p' );
var pow = require( 'math-power' );


// FUNCTIONS //

var lower_gamma_series = require( './lower_gamma_series.js' );
var upper_gamma_fraction = require( './upper_gamma_fraction.js' );


// CONSTANTS //

var EPSILON = 2.220446049250313e-16;
var LOG_MAX_VALUE = 709.0;
var LOG_MIN_VALUE = -708.0;


// IBETA POWER TERMS //

/**
* FUNCTION: ibeta_power_terms( a, b, x, y, normalised )
*	Compute the leading power terms in the incomplete beta function:
*		(x^a)(y^b)/Beta(a,b) when normalised, and
*		(x^a)(y^b) otherwise.
*
* @param {Number} a - function parameter
* @param {Number} b - function parameter
* @param {Number} x - function input
* @param {Number} y - equal to `1 - x`
* @param {Boolean} normalised - boolean indicating whether to evaluate the regularized or non-regularized incomplete beta function
* @returns {Number} leading power terms
*/
function ibeta_power_terms( a, b, x, y, normalised ) {
	var result;
	var c;
	var la, lb, lc;
	var b1, b2, e1, lb1, lb2;
	var sa, sb, sc;
	var p1, p2, p3;

	if ( !normalised ) {
		return pow( x, a ) * pow( y, b );
	}
	result = 0;
	c = a + b;
	// integration limits for the gamma functions:
	la = a + 5;
	lb = b + 5;
	lc = a + b + 5;
	// gamma function partials:
	sa = lower_gamma_series( a, la ) / a;
	sa += upper_gamma_fraction( a, la, EPSILON );
	sb = lower_gamma_series( b, lb ) / b;
	sb += upper_gamma_fraction( b, lb, EPSILON );
	sc = lower_gamma_series( c, lc ) / c;
	sc += upper_gamma_fraction( c, lc, EPSILON );
	// gamma function powers combined with incomplete beta powers:
	b1 = ( x * lc ) / la;
	b2 = ( y * lc ) / lb;
	e1 = lc - la - lb;
	lb1 = a * ln( b1 );
	lb2 = b * ln( b2 );

	if((lb1 >= LOG_MAX_VALUE ) || (lb1 <= LOG_MIN_VALUE ) || (lb2 >= LOG_MAX_VALUE ) ||
		(lb2 <= LOG_MIN_VALUE) || (e1 >= LOG_MAX_VALUE ) || (e1 <= LOG_MIN_VALUE )
	) {
		result = exp( lb1 + lb2 - e1 );
	} else {
		if( (abs(b1 - 1) * a < 10) && (a > 1) ) {
			p1 = exp( a * log1p( (x * b - y * la) / la ) );
		} else {
			p1 = pow( b1, a );
		}
		if( ( abs(b2 - 1) * b < 10 ) && ( b > 1 ) ) {
			p2 = exp( b * log1p( (y * a - x * lb) / lb ) );
		} else {
			p2 = pow( b2, b );
		}
		p3 = exp( e1 );
		result = p1 * p2 / p3;
	}
	// and combine with the remaining gamma function components:
	result /= sa * sb / sc;
	return result;
} // end FUNCTION ibeta_power_terms()


// EXPORTS //

module.exports = ibeta_power_terms;
