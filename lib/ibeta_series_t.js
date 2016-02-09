'use strict';

// IBETA SERIES //

/**
* FUNCTION: ibeta_series_t( a_, b_, x_, mult  )
*	Creates a function to evaluate a series expansion of the incomplete beta function.
*
* @param {Number} a_ - function parameter
* @param {Number} b_ - function parameter
* @param {Number} x_ - function input
* @param {Number} mult - multiplication term of series
* @returns {Function} series function
*/
function ibeta_series_t( a_, b_, x_, mult ) {
	var result = mult,
		x = x_,
		apn = a_,
		poch = 1 - b_,
		n = 1;
	return function() {
		var r = result / apn;
		apn += 1;
		result *= poch * x / n;
		++n;
		poch += 1;
		return r;
	};
} // end FUNCTION ibeta_series_t()


// EXPORTS //

module.exports = ibeta_series_t;
