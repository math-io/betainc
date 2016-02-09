'use strict';

/**
* FUNCTION: lower_incomplete_gamma_series( a1, z1 )
*	Creates a function to evaluate a series expansion of the lower incomplete gamma function.
*
* @param {Number} a1 - function parameter
* @param {Number} z1 - function parameter
* @returns {Function} series function
*/
function upper_incomplete_gamma_fract( a1, z1 ) {
	var z = z1 - a1 + 1,
		a = a1,
		k = 0;
	return function() {
		++k;
		z += 2;
		return {
			'a': k * (a - k),
			'b': z
		};
	};
} // end FUNCTION upper_incomplete_gamma_fract()


// EXPORTS //

module.exports = upper_incomplete_gamma_fract;
