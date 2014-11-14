//-------|---------|---------|---------|---------|---------|---------|---------|
/*
b_spline.h - operator related routine(s)
*/

#ifndef	B_SPLINE_H
#define	B_SPLINE_H

#include "all.h"
#include "memory.h"
#include "interpolant.h"	//	remove?
#include "polynomial.h"

/*
void compute_B(short degree)
where
	degree is the max degree of the operator in s^2 and
	B is the blurring operator.

To calculate blurring operator, B_p:
c_n = 1/a_0 (b_n - [a_1*c_{n-1} + a_2*c_{n-2} + ... + a_n*c_0])

NOTE:
	* we do not need to store b at all
	* factorials are not calculated anew in every iteration of the loops
		* f keeps a running product that is equivalent to factorial()
	* inner loop only touches columns which will change for given n,k
		* inner loop could be parallelized, although probably not worth it
	* for more details, go to
		http://pmsm.cs.purdue.edu/tiki-index.php?page=omegaprime
		or
		http://en.wikipedia.org/wiki/Formal_power_series#Dividing_series
*/
void compute_blurring_operator(short degree, double* B);

/*
Use O^{-1} = 1 / op = 1 / 1 - D -> D = 1 - O
Then use geometric series
O^{-1} = 1 + D + D^2 + D^3 + ...

NOTE:	O ("oh", not zero) is the input operator
		I ("eye") is the inverse, I = O^{-1}
*/
void compute_operator_inverse(short O_degree, const double* O,
							  short I_degree, double* I);

/*
Convert an operator from delta^2 to shifts: E-2+E

cd -> "central difference" operator (input)
s  -> "shifts" operator (output)

NOTE:
Both arrays should be length cd_degree+1
*/
void convert_delta2_to_shifts(short cd_degree, const double* cd, double* s);

/*
compute omega values
*/
void compute_omega(short p, short n, double* omega);

/*
compute omega' values
*/
void compute_omega_prime(short p, short mu, double* omegap);

#endif

// End of file