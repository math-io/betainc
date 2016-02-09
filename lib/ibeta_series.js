'use strict';

// MODULES //

var exp = require( 'math-exp' );
var ln = require( 'math-ln' );
var log1p = require( 'math-log1p' );
var pow = require( 'math-power' );


// FUNCTIONS //

var lower_gamma_series = require( './lower_gamma_series.js' );
var upper_gamma_fraction = require( './upper_gamma_fraction.js' );
var sum_series = require( './sum_series.js' );
var ibeta_series_t = require( './ibeta_series_t.js' );


// CONSTANTS //

var EPSILON = 2.220446049250313e-16;
var LOG_MAX_VALUE = 709.0;
var LOG_MIN_VALUE = -708.0;
var MIN_VALUE = require( 'const-min-float64' );


// IBETA SERIES //

/**
* FUNCTION: ibeta_series( a, b, x, s0, normalised )
*	Creates a function to evaluate a series expansion of the incomplete beta function.
*
* @param {Number} a - function parameter
* @param {Number} b - function parameter
* @param {Number} x - function input
* @param {Number} s0 - initial value
* @param {Boolean} normalised - boolean indicating whether to evaluate the regularized or non-regularized incomplete beta function
* @returns {Number} sum of terms of lower incomplete beta expansion
*/
function ibeta_series( a, b, x, s0, normalised ) {
	var result;
	var c;
	var la, lb, lc;
	var sa, sb, sc;
	var b1, b2;
	var e1;
	var lb1, lb2;
	var p;
	var s;
	var max_iter;

	if ( normalised ) {
		c = a + b;
		// figure out integration limits for the gamma function:
		la = a + 5;
		lb = b + 5;
		lc = a + b + 5;

		// calculate the gamma parts:
		sa = lower_gamma_series( a, la ) / a;
		sa += upper_gamma_fraction( a, la, EPSILON );
		sb = lower_gamma_series( b, lb ) / b;
		sb += upper_gamma_fraction( b, lb, EPSILON );
		sc = lower_gamma_series( c, lc ) / c;
		sc += upper_gamma_fraction( c, lc, EPSILON );

		// and their combined power-terms:
		b1 = (x * lc) / la;
		b2 = lc/lb;
		e1 = lc - la - lb;
		lb1 = a * ln(b1);
		lb2 = b * ln(b2);

	if( ( lb1 >= LOG_MAX_VALUE ) || ( lb1 <= LOG_MIN_VALUE ) || ( lb2 >= LOG_MAX_VALUE ) ||
		( lb2 <= LOG_MIN_VALUE ) || ( e1 >= LOG_MAX_VALUE ) || ( e1 <= LOG_MIN_VALUE ) ) {
			p = lb1 + lb2 - e1;
			result = exp( p );
		} else {
			result = pow(b1, a);
			if ( a * b < lb * 10 ) {
				result *= exp( b * log1p( a / lb ) );
			} else {
				result *= pow( b2, b );
			}
			result /= exp( e1 );
		}
		// And combine the results:
		result /= sa * sb / sc;
	} else {
		// Non-normalised, just compute the power:
		result = pow( x, a );
	}
	if ( result < MIN_VALUE ) {
		return s0; // Safeguard: series can't cope with denorms.
	}
	s = ibeta_series_t( a, b, x, result );
	max_iter = 100; //CHANGE IT!!
	result = sum_series( s, EPSILON, max_iter, s0 );
	return result;
} // end FUNCTION ibeta_series()


// EXPORTS //

module.exports = ibeta_series;
