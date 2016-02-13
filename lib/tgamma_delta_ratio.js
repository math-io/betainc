'use strict';

// MODULES //

var abs = require( 'math-abs' );
var exp = require( 'math-exp' );
var factorial = require( 'factorial' );
var floor = require( 'math-floor' );
var log1p = require( 'math-log1p' );
var pow = require( 'math-power' );


// FUNCTIONS //

var lower_gamma_series = require( './lower_gamma_series.js' );
var upper_gamma_fraction = require( './upper_gamma_fraction.js' );


// CONSTANTS //

var EPSILON = 2.220446049250313e-16;
var MAX_FACTORIAL = 170;


// TGAMMA DELTA RATIO IMPLEMENTATION //

function tgamma_delta_ratio_imp_lanczos( z, delta ) {
	var prefix;
	var zd;
	var sum;
	var sumd;
	//
	// The upper gamma fraction is *very* slow for z < 6, actually it's very
	// slow to converge everywhere but recursing until z > 6 gets rid of the
	// worst of it's behaviour.
	//
	prefix = 1;
	zd = z + delta;
	while((zd < 6) && (z < 6)) {
		prefix /= z;
		prefix *= zd;
		z += 1;
		zd += 1;
	}
	if ( delta < 10 ) {
		prefix *= exp( -z * log1p( delta / z ) );
	} else {
		prefix *= pow( z / zd, z );
	}
	prefix *= pow( Math.E / zd, delta);
	sum = lower_gamma_series( z, z ) / z;
	sum += upper_gamma_fraction(z, z, EPSILON );
	sumd = lower_gamma_series( zd, zd ) / zd;
	sumd += upper_gamma_fraction(zd, zd, EPSILON );
	sum /= sumd;
	return sum * prefix;
} // end FUNCTION tgamma_delta_ratio_imp_lanczos()


function tgamma_delta_ratio_imp( z, delta ) {
	var result;
	if ( floor( delta ) === delta ) {
		if( floor( z ) === z ) {
			// Both z and delta are integers, see if we can just use table lookup
			// of the factorials to get the result:
			if( ( z <= MAX_FACTORIAL ) && ( z + delta <= MAX_FACTORIAL ) ) {
				return factorial( floor(z) - 1 ) / factorial( floor(z+delta) - 1);
			}
		}
		if ( abs(delta) < 20 ) {
			// Delta is a small integer, we can use a finite product:
			if ( delta === 0 ) {
				return 1;
			}
			if ( delta < 0 ) {
				z -= 1;
			}
			result = z;
			while ( 0 !== (delta += 1) ) {
				z -= 1;
				result *= z;
			}
			return result;
		} else {
			result = 1 / z;
			while ( 0 !== (delta -= 1 ) ) {
				z += 1;
				result /= z;
			}
			return result;
		}
	}
	return tgamma_delta_ratio_imp_lanczos( z, delta );
} // end FUNCTION tgamma_delta_ratio_imp()


// EXPORTS //

module.exports = tgamma_delta_ratio_imp;
