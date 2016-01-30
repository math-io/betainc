'use strict';

/**
* FUNCTION: ibeta_fraction2_t( a, b, x, y )
*	Creates a term function used in evaluating the continued fraction expansion for the incomplete beta function.
*
* @param {Number} a - function parameter
* @param {Number} b - function parameter
* @param {Number} x - input value
* @param {Number} y - equal to `1 - x`
* @returns {Function} function returning the `a` and `b` terms of the continued fraction expansion
*/
function ibeta_fraction2_t( a, b, x, y ) {
	var m = 0;
	return function() {
		var aN;
		var bN;
		var denom;

		aN = (a + m - 1) * (a + b + m - 1) * m * (b - m) * x * x;
		denom = (a + 2 * m - 1);
		aN /= denom * denom;

		bN = m;
		bN += (m * (b - m) * x) / (a + 2*m - 1);
		bN += ( (a + m) * ( a * y - b * x + 1 + m *(2 - x) ) ) / (a + 2*m + 1);

		++m;
		return {
			'a': aN,
			'b': bN
		};
	};
} // end FUNCTION ibeta_fraction2_t()


// EXPORTS //

module.exports = ibeta_fraction2_t;
