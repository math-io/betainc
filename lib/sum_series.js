'use strict';

// MODULES //

var abs = require( 'math-abs' );


// SUM SERIES //

/**
* FUNCTION: sum_series( func, factor, max_terms, init_value )
*	Sums the element of the series given by the supplied function.
*
* @param {Function} func - series function
* @param {Number} factor - further terms are only added as long as factor*result is smaller than the next term
* @param {Number} max_terms - maximum number of terms to be added
* @param {Number} init_value - initial value of the resulting sum
* @returns {Number} sum of all series terms
*/
function sum_series( func, factor, max_terms, init_value ) {
	var counter;
	var result;
	var next_term;

	counter = max_terms;
	result = init_value;
	do {
		next_term = func();
		result += next_term;
	}
	while( ( abs(factor * result) < abs(next_term) ) && --counter );

	// Set max_terms to the actual number of terms of the series evaluated:
	max_terms = max_terms - counter;
	return result;
} // end FUNCTION sum_series()


// EXPORTS //

module.exports = sum_series;
