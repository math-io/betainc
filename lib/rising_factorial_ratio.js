'use strict';

// RISING FACTORIAL RATIO //

/**
* FUNCTION: rising_factorial_ratio( a, b, k )
*	Calculates
*		(a)(a+1)(a+2)...(a+k-1)
*		_______________________
*		(b)(b+1)(b+2)...(b+k-1)
*	This function is only called with small k, for large k
*	it is grossly inefficient. This function is only needed for the non-regular incomplete beta,
*	it computes the delta in: beta(a,b,x) = prefix + delta * beta(a+k,b,x)
*
* @param {Number} a - input value
* @param {Number} b - input value
* @param {Number} k - input value
* @returns {Number} ratio value
*/
function rising_factorial_ratio( a, b, k ) {
	var result;
	var i;

	if ( k === 0 ) {
		return 1;
	}
	result = 1;
	for ( i = 0; i < k; ++i ) {
		result *= ( a + i ) / ( b + i );
	}
	return result;
} // end FUNCTION rising_factorial_ratio()


// EXPORTS //

module.exports = rising_factorial_ratio;
