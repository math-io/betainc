'use strict';

// FUNCTIONS //

var ibeta_imp = require( './ibeta_imp.js' );


// BETAINC //

/**
* FUNCTION: betainc( x, a, b, regularized, upper )
*	Evaluates the incomplete beta function.
*
* @param {Number} x - function parameter
* @param {Number} a - function parameter
* @param {Number} b - function parameter
* @param {Boolean} [regularized=true] - boolean indicating if the function should evaluate the regularized or non-regularized incomplete beta function
* @param {Boolean} [upper=false] - boolean indicating if the function should return the upper tail of the incomplete beta function
* @returns {Number} function value
*/
function betainc( x, a, b, regularized, upper ) {
	if ( regularized !== false ) {
		return upper ? ibeta_imp( a, b, x, undefined, true, true ) : ibeta_imp( a, b, x, false, true );
	}
	return upper ? ibeta_imp( a, b, x, undefined, true, false ) : ibeta_imp( a, b, x, false, false );
} // end FUNCTION betainc()


// EXPORTS //

module.exports = betainc;
