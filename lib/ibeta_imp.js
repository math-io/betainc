'use strict';

// MODULES //

var asin = Math.asin;
var beta = require( 'math-beta' );
var exp = require( 'math-exp' );
var expm1 = require( 'math-expm1' );
var floor = require( 'math-floor' );
var log1p = require( 'math-log1p' );
var max = Math.max;
var min = Math.min;
var pow = require( 'math-power' );
var sqrt = require( 'math-sqrt' );


// FUNCTIONS //

var beta_small_b_large_a_series = require( './beta_small_b_large_a_series.js' );
var binomial_ccdf = require( './binomial_ccdf.js' );
var ibeta_a_step = require( './ibeta_a_step.js' );
var ibeta_series = require( './ibeta_series.js' );
var ibeta_fraction2 = require( './ibeta_fraction2.js');
var rising_factorial_ratio = require( './rising_factorial_ratio.js' );


// CONSTANTS //

var MAX_INT32 = require( 'compute-const-max-int32' );
var PI = require( 'const-pi' );
var HALF_PI = PI / 2;


// IBETA IMPLEMENTATION //

/**
* FUNCTION: ibeta_imp( a, b, x, invert, normalised )
*	Evaluates the incomplete beta function. This function divides up the
*	input range and selects the right implementation method for each domain
*
* @param {Number} a - function parameter
* @param {Number} b - function parameter
* @param {Number} x - function input
* @param {Boolean} invert - boolean indicating if the function should return the upper tail of the incomplete beta function instead
* @param {Boolean} normalised - boolean indicating if the function should evaluate the regularized incomplete beta function
* @returns {Number} function value
*/
function ibeta_imp( a, b, x, invert, normalised ) {
	// jshint maxstatements: 400
	var fract;
	var y;
	var p;
	var tmp;
	var prefix;
	var lambda;
	var k;
	var n;
	var bbar;
	y = 1 - x;
	if ( (x < 0) || (x > 1) ) {
		return NaN;
	}
	if ( normalised ) {
		if (a < 0) { return NaN; }
		if (b < 0) { return NaN; }
		// extend to a few very special cases:
		if ( a === 0) {
			if (b === 0) {
				return NaN;
			}
			if (b > 0) {
				return invert ? 0 : 1;
			}
		} else if (b === 0) {
			if (a > 0) {
				return invert ? 1 : 0;
			}
		}
	} else {
		if ( a <= 0 || b <= 0 ) {
			return NaN;
		}
	 }
	if ( x === 0 ) {
		return (invert ? (normalised ? 1 : beta(a, b) ) : 0 );
	}
	if ( x === 1 ) {
		return (invert === 0 ? (normalised ? 1 : beta(a, b)) : 0);
	}
	if( (a === 0.5) && (b === 0.5) ) {
		// We have an arcsine distribution:
		p = invert ? asin( sqrt(y) ) / HALF_PI : asin( sqrt(x) ) / HALF_PI;
		if ( !normalised ) {
			p *= PI;
		}
		return p;
	}
	if ( a === 1 ) {
		tmp = b;
		b = a;
		a = tmp;

		tmp = y;
		y = x;
		x = tmp;

		invert = !invert;
	}
	if ( b === 1 ) {
		// Special case see: http://functions.wolfram.com/GammaBetaErf/BetaRegularized/03/01/01/
		if ( a === 1 ) {
			return invert ? y : x;
		}
		if ( y < 0.5 ) {
			p = invert ? -expm1(a * log1p(-y)) : exp(a * log1p(-y));
		} else {
			p = invert ? -( pow( x, a ) - 1 ) : pow( x, a );
	 	}
		if( !normalised ) {
			p /= a;
		}
		return p;
	}
	if ( min(a, b) <= 1) {
		if(x > 0.5) {
			tmp = b;
			b = a;
			a = tmp;

			tmp = y;
			y = x;
			x = tmp;

			invert = !invert;
		}
		if ( max(a, b) <= 1 ) {
			// Both a,b < 1:
			if( (a >= min( 0.2, b ) ) || ( pow(x, a) <= 0.9 ) ) {
				if ( !invert ) {
					fract = ibeta_series(a, b, x, 0, normalised );
				} else {
					fract = -(normalised ? 1 : beta( a, b ) );
					invert = false;
					fract = -ibeta_series(a, b, x, fract, normalised );
				}
			} else {
				tmp = b;
				b = a;
				a = tmp;

				tmp = y;
				y = x;
				x = tmp;

				invert = !invert;
				if ( y >= 0.3 ) {
					if ( !invert ) {
						fract = ibeta_series(a, b, x, 0, normalised );
					} else {
						fract = -( normalised ? 1 : beta( a, b ) );
						invert = false;
						fract = -ibeta_series(a, b, x, fract, normalised );
					}
				} else {
					// Sidestep on a, and then use the series representation:
					if ( !normalised ) {
						prefix = rising_factorial_ratio( a + b, a, 20 );
					} else {
						prefix = 1;
					}
					fract = ibeta_a_step( a, b, x, y, 20, normalised );
					if ( !invert ) {
						fract = beta_small_b_large_a_series( a + 20, b, x, y, fract, prefix, normalised );
					} else {
						fract -= ( normalised ? 1 : beta( a, b ) );
						invert = false;
						fract = -beta_small_b_large_a_series( a + 20, b, x, y, fract, prefix, normalised );
					}
				}
			}
		} else {
			// One of a, b < 1 only:
			if( (b <= 1) || ( (x < 0.1) && ( pow(b * x, a) <= 0.7 ) ) ) {
				if ( !invert ) {
					fract = ibeta_series( a, b, x, 0, normalised );
				} else {
					fract = -( normalised ? 1 : beta( a, b ) );
					invert = false;
					fract = -ibeta_series( a, b, x, fract, normalised );
				}
			} else {
				tmp = b;
				b = a;
				a = tmp;

				tmp = y;
				y = x;
				x = tmp;
				invert = !invert;

				if ( y >= 0.3 ) {
					if (!invert) {
						fract = ibeta_series(a, b, x, 0, normalised );
					} else {
						fract = -(normalised ? 1 : beta( a, b ));
						invert = false;
						fract = -ibeta_series(a, b, x, fract, normalised );
					}
				} else if ( a >= 15 ) {
					if(!invert) {
						fract = beta_small_b_large_a_series( a, b, x, y, 0, 1, normalised );
					} else {
						fract = -(normalised ? 1 : beta( a, b ));
						invert = false;
						fract = -beta_small_b_large_a_series( a, b, x, y, fract, 1, normalised );
					}
				} else {
					// Sidestep to improve errors:
					if ( !normalised ) {
						prefix = rising_factorial_ratio( a + b, a, 20 );
					} else {
						prefix = 1;
					}
					fract = ibeta_a_step( a, b, x, y, 20, normalised );
					if ( !invert ) {
						fract = beta_small_b_large_a_series( a + 20, b, x, y, fract, prefix, normalised );
					} else {
						fract -= ( normalised ? 1 : beta( a, b ) );
						invert = false;
						fract = -beta_small_b_large_a_series( a + 20, b, x, y, fract, prefix, normalised );
					}
				}
			}
		}
	} else {
		// Both a,b >= 1:
		if(a < b) {
			lambda = a - (a + b) * x;
		} else {
			lambda = (a + b) * y - b;
		}
		if ( lambda < 0 ) {
			tmp = b;
			b = a;
			a = tmp;

			tmp = y;
			y = x;
			x = tmp;
			invert = !invert;
		}
		if (b < 40) {
			if ( (floor(a) === a) && (floor(b) === b) && (a < MAX_INT32 - 100) ) {
				// relate to the binomial distribution and use a finite sum:
				k = a - 1;
				n = b + k;
				fract = binomial_ccdf(n, k, x, y);
				if (!normalised) {
					fract *= beta( a, b );
				}
			} else if (b * x <= 0.7) {
				if( !invert ) {
					fract = ibeta_series( a, b, x, 0, normalised );
				} else {
					fract = -( normalised ? 1 : beta( a, b ) );
					invert = false;
					fract = -ibeta_series( a, b, x, fract, normalised );
				}
			} else if ( a > 15 ) {
				// sidestep so we can use the series representation:
				n = floor(b);
				if (n === b) {
					--n;
				}
				bbar = b - n;
				if ( !normalised ) {
					prefix = rising_factorial_ratio( a + bbar, bbar, n );
				} else {
					prefix = 1;
				}
				fract = ibeta_a_step( bbar, a, y, x, n, normalised );
				fract = beta_small_b_large_a_series( a, bbar, x, y, fract, 1, normalised );
				fract /= prefix;
			} else if( normalised ) {
				n = floor(b);
				bbar = b - n;
				if ( bbar <= 0 ) {
					--n;
					bbar += 1;
				}
				fract = ibeta_a_step( bbar, a, y, x, n, normalised );
				fract += ibeta_a_step( a, bbar, x, y, 20, normalised );
				if (invert) {
					fract -= 1;
				}
				fract = beta_small_b_large_a_series( a + 20,  bbar, x, y, fract, 1, normalised );
				if(invert) {
					fract = -fract;
					invert = false;
				}
			} else {
				fract = ibeta_fraction2( a, b, x, y, normalised );
			}
		} else {
			fract = ibeta_fraction2( a, b, x, y , normalised );
		}
	}
	return invert ? ( normalised ? 1 : beta( a, b ) ) - fract : fract;
} // end FUNCTION ibeta_imp()


// EXPORTS //

module.exports = ibeta_imp;
