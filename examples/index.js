'use strict';

var betainc = require( './../lib' );

for ( var a = 0; a < 5; a++ ) {
	for ( var b = 5; b > 0; b-- ) {
		console.log( 'x: 0.5, \t a: %d, \t b: %d, \t f(x,a,b): %d', a, b, betainc( 0.5, a, b ) );
	}
}
