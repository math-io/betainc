'use strict';

// FUNCTIONS //

var continued_fraction_b = require( './continued_fraction_b.js' );
var ibeta_fraction2_t = require( './ibeta_fraction2_t' );
var ibeta_power_terms = require( './ibeta_power_terms.js' );


// CONSTANTS //

var EPSILON = 2.220446049250313e-16;


// IBETA FRACTION 2 //

/**
* FUNCTION: ibeta_fraction2( a, b, x, y, normalised )
*	Evaluates the incomplete beta function via the continued fraction representation.
*
* @param {Number} a - function parameter
* @param {Number} b - function parameter
* @param {Number} x - input value
* @param {Number} y - equal to `1 - x`
* @param {Boolean} normalised - boolean indicating whether to evaluate the regularized or non-regularized incomplete beta function
* @returns {Number} function value
*/
function ibeta_fraction2( a, b, x, y, normalised ) {
	var f;
	var result;
	var fract;

	result = ibeta_power_terms( a, b, x, y, normalised );
	if ( result === 0 ) {
		return result;
	}
	f = ibeta_fraction2_t( a, b, x, y );
	fract = continued_fraction_b( f, EPSILON );
	return result / fract;
} // end FUNCTION ibeta_fraction2()


// EXPORTS //

module.exports = ibeta_fraction2;
