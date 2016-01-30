'use strict';

// MODULES //

var ibeta_power_terms = require( './ibeta_power_terms.js' );


// IBETA A STEP //

/**
* FUNCTION: ibeta_a_step( a, b, x, y, k, normalised, p_derivative )
*	Computes the difference between ibeta(a,b,x) and ibeta(a+k,b,x).
*
* @param {Number} a - function parameter
* @param {Number} b - function parameter
* @param {Number} x - function input
* @param {Number} y - equal to `1 - x`
* @param {Boolean} normalised - boolean indicating whether to evaluate the regularized or non-regularized incomplete beta function
* @param {Object} p_derivative - object to hold the derivative of `p`.
* @returns {Number} difference of ibeta(a,b,x) and ibeta(a+k,b,x)
*/
function ibeta_a_step( a, b, x, y, k, normalised, p_derivative ) {
	var prefix;
	var sum;
	var term;
	var i;

	prefix = ibeta_power_terms( a, b, x, y, normalised );
	if ( p_derivative ) {
		p_derivative.value = prefix;
	}
	prefix /= a;
	if ( prefix === 0 ) {
		return prefix;
	}
	sum = 1;
	term = 1;
	// series summation from 0 to k-1:
	for( i = 0; i < k-1; ++i ) {
		term *= ( a + b + i ) * x / ( a + i + 1 );
		sum += term;
	}
	prefix *= sum;
	return prefix;
} // end FUNCTION ibeta_a_step()


// EXPORTS //

module.exports = ibeta_a_step;
