'use strict';

// MODULES //

var abs = require( 'math-abs' );


// CONSTANTS //

var TINY = require( 'compute-const-smallest-float32' ).VALUE;


// CONTINUED FRACTION B //

/**
* FUNCTION: continued_fraction_b( g, factor )
*	Evaluates
*	  b0 +	   a1
*	  ---------------
*	  b1 +	 a2
*		   ----------
*		   b2 +   a3
*				-----
*				b3 + ...
*
* @param {Function} g - function giving terms of continued fraction expansion
* @param {Number} factor - further terms are only added as long as factor*result is smaller than the next term
* @returns {Number} evaluated expansion
*/
function continued_fraction_b( g, factor ) {
	var v = g();
	var f;
	var C;
	var D;
	var delta;
	f = v.b;
	if( f === 0 ) {
		f = TINY;
	}
	C = f;
	D = 0;

	var max_iter = 1000;
	do {
		v = g();
		D = v.b + v.a * D;
		if ( D === 0 ) {
			D = TINY;
		}
		C = v.b + v.a / C;
		if ( C === 0 ) {
			C = TINY;
		}
		D = 1 / D;
		delta = C * D;
		f = f * delta;
	} while( ( abs( delta - 1 ) > factor ) && --max_iter );
	return f;
} // end FUNCTION continued_fraction_b()


// EXPORTS //

module.exports = continued_fraction_b;
