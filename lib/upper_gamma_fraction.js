'use strict';

// MODULES //

var continued_fraction = require( 'math-continued-fraction' );


// FUNCTIONS //

var upper_incomplete_gamma_fract = require( './upper_incomplete_gamma_fract' );


// UPPER GAMMA FRACTION //

/**
* FUNCTION: lower_incomplete_gamma_series( a, z, eps )
*	Evaluate the lower incomplete gamma integral via a series expansion and divide by gamma(z) to normalise.
*
* @param {Number} a - function parameter
* @param {Number} z - function parameter
* @param {Number} eps - tolerance for continued fraction
* @returns {Number} function value
*/
function upper_gamma_fraction( a, z, eps ) {
	var f = upper_incomplete_gamma_fract( a, z );
	return 1 / ( z - a + 1 + continued_fraction( f, { 'tolerance': eps } ) );
} // end FUNCTION upper_gamma_fraction()


// EXPORTS //

module.exports = upper_gamma_fraction;
