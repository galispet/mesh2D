#pragma once

#include <float.h>



/*****************************************************************************/
/*                                                                           */
/*    - Splitter : (DBL_MANT_DIG + 1) >> 1 right bit shift by 1bit			 */
/*		           i.e. the number (54) is divided by 2	(27)				 */
/*				 : 1 << result left bit shift by result (27) bit		     */
/*				 : it results in the number which is represented with 1	in	 */
/*				   binary code in the middle of the mantisa of the double    */
/*                 (the 1 is in the 28-th place in binary: ...00100000...0)  */
/*                                                                           */
/*****************************************************************************/



#define EPS (DBL_EPSILON / 2)

const double epsilonArithmetic	= EPS;
const double splitter			= (1 << ((DBL_MANT_DIG + 1) >> 1)) + 1.0;
const double resultErrorBound	= (3.0 + 8.0 * EPS) * EPS;

const double orientationErrorBoundA = (3.0 + 16.0 * EPS) * EPS;
const double orientationErrorBoundB = (2.0 + 12.0 * EPS) * EPS;
const double orientationErrorBoundC = (9.0 + 64.0 * EPS) * EPS * EPS;

const double inCircumCircleErrorBoundA = (10.0 + 96.0 * EPS) * EPS;
const double inCircumCircleErrorBoundB = (4.0 + 48.0 * EPS) * EPS;
const double inCircumCircleErrorBoundC = (26.0 + 288.0 * EPS) * EPS * EPS;

/*
const double isperrboundA = (16.0 + 224.0 * EPS) * EPS;
const double isperrboundB = (5.0 + 72.0 * EPS) * EPS;
const double isperrboundC = (71.0 + 1408.0 * EPS) * EPS * EPS;
*/



inline void Split(double const a, double & ahi, double & alo) {

	const double c		= splitter * a;
	const double abig	= c - a;

	ahi = c - abig;
	alo = c - ahi;

};
inline double Estimate(int const m, double * e) {

	double Q = e[0];

	for (int i = 1; i < m; i++)
		Q += e[i];

	return Q;

};


inline void Fast_Two_Sum_Tail(double const a, double const b, double & x, double & y) {

	double const bvirtual = x - a;
	y = b - bvirtual;

};
inline void Fast_Two_Sum(double const a, double const b, double & x, double & y) {

	x = a + b;
	Fast_Two_Sum_Tail(a, b, x, y);

};


inline void Two_Sum_Tail(double const a, double const b, double & x, double & y) {

	const double bvirtual = x - a;
	const double avirtual = x - bvirtual;

	const double broundoff = b - bvirtual;
	const double aroundoff = a - avirtual;

	y = aroundoff + broundoff;

};
inline void Two_Sum(double const a, double const b, double & x, double & y) {

	x = a + b;
	Two_Sum_Tail(a, b, x, y);

};


inline void Two_One_Sum(double const a1, double const a0, double const b, double & x2, double & x1, double & x0) {

	double _i;

	Two_Sum(a0, b, _i, x0);
	Two_Sum(a1, _i, x2, x1);

};
inline void Two_Two_Sum(double a1, double a0, double b1, double b0, double & x3, double & x2, double & x1, double & x0) {

	double _j;
	double _0;

	Two_One_Sum(a1, a0, b0, _j, _0, x0);
	Two_One_Sum(_j, _0, b1, x3, x2, x1);

};


inline void Two_Diff_Tail(double const a, double const b, double & x, double & y) {

	const double bvirtual = a - x;
	const double avirtual = x + bvirtual;

	const double broundoff = bvirtual - b;
	const double aroundoff = a - avirtual;

	y = aroundoff + broundoff;

};
inline void Two_Diff(double const a, double const b, double & x, double & y) {

	x = a - b;
	Two_Diff_Tail(a, b, x, y);

};


inline void Two_One_Diff(double const a1, double const a0, double const b, double & x2, double & x1, double & x0) {

	double _i;

	Two_Diff(a0, b, _i, x0);
	Two_Sum(a1, _i, x2, x1);
	// There is really Two_Sum. Check this in Shewchuck paper, if there is really this one;

};
inline void Two_Two_Diff(double a1, double a0, double b1, double b0, double & x3, double & x2, double & x1, double & x0) {

	double _j;
	double _0;

	Two_One_Diff(a1, a0, b0, _j, _0, x0);
	Two_One_Diff(_j, _0, b1, x3, x2, x1);

};


inline void Two_Product_Tail(double const a, double const b, double & x, double & y) {

	double ahi;
	double alo;

	double bhi;
	double blo;

	Split(a, ahi, alo);
	Split(b, bhi, blo);

	const double Error1 = x - (ahi * bhi);
	const double Error2 = Error1 - (alo * bhi);
	const double Error3 = Error2 - (ahi * blo);

	y = alo * blo - Error3;

};
inline void Two_Product(double const a, double const b, double & x, double & y) {

	x = a * b;
	Two_Product_Tail(a, b, x, y);

};
inline void Two_Product_Presplit(double const a, double const b, double const bhi, double const blo, double & x, double & y) {

	x = a * b;

	double ahi;
	double alo;

	Split(a, ahi, alo);

	const double Error1 = x - (ahi * bhi);
	const double Error2 = Error1 - (alo * bhi);
	const double Error3 = Error2 - (ahi * blo);

	y = (alo * blo) - Error3;

};


inline void Square_Tail(double const a, double & x, double & y) {

	double ahi;
	double alo;

	Split(a, ahi, alo);

	const double Error1 = x - (ahi * ahi);
	const double Error3 = Error1 - (ahi + ahi) * alo;

	y = (alo * alo) - Error3;

};
inline void Square(double const a, double & x, double & y) {

	x = a * a;
	Square_Tail(a, x, y);

};


/*****************************************************************************/
/*                                                                           */
/*  fast_expansion_sum_zeroelim()   Sum two expansions, eliminating zero     */
/*                                  components from the output expansion.    */
/*                                                                           */
/*  Sets h = e + f.  See the long version of my paper for details.           */
/*                                                                           */
/*  If round-to-even is used (as with IEEE 754), maintains the strongly      */
/*  nonoverlapping property.  (That is, if e is strongly nonoverlapping, h   */
/*  will be also.)  Does NOT maintain the nonoverlapping or nonadjacent      */
/*  properties.                                                              */
/*                                                                           */
/*****************************************************************************/
int Fast_Expansion_Sum_Zero_Eliminiation(int const eLength, double * e, const int fLength, double * f, double * h) {

	/*****************************************************************************/
	/*                                                                           */
	/*    - INEXACT == volatile													 */
	/*                                                                           */
	/*****************************************************************************/
	double Qnew;
	double hh;
	//double bvirt;

	/*****************************************************************************/
	/*                                                                           */
	/*    - EXACT == nothing													 */
	/*                                                                           */
	/*****************************************************************************/
	double Q;
	//double avirt, bround, around;
	double enow, fnow;

	int eindex, findex, hindex;


	eindex = 0;
	findex = 0;
	hindex = 0;

	enow = e[0];
	fnow = f[0];

	if ((fnow > enow) == (fnow > -enow)) {
		//if (abs(fnow) > abs(enow)) {

		Q = enow;
		enow = e[++eindex];

	}
	else {

		Q = fnow;
		fnow = f[++findex];

	}

	if ((eindex < eLength) && (findex < fLength)) {

		if ((fnow > enow) == (fnow > -enow)) {
			//if (fnow > abs(enow)) {

			Fast_Two_Sum(enow, Q, Qnew, hh);
			enow = e[++eindex];

		}
		else {

			Fast_Two_Sum(fnow, Q, Qnew, hh);
			fnow = f[++findex];

		}

		Q = Qnew;

		if (hh != 0.0)
			h[hindex++] = hh;


		while (eindex < eLength && findex < fLength) {

			if ((fnow > enow) == (fnow > -enow)) {
				//if (fnow > abs(enow)) {

				Two_Sum(Q, enow, Qnew, hh);
				enow = e[++eindex];

			}
			else {

				Two_Sum(Q, fnow, Qnew, hh);
				fnow = f[++findex];

			}

			Q = Qnew;

			if (hh != 0.0)
				h[hindex++] = hh;

		}
	}

	while (eindex < eLength) {

		Two_Sum(Q, enow, Qnew, hh);

		enow = e[++eindex];
		Q = Qnew;

		if (hh != 0.0)
			h[hindex++] = hh;

	}
	while (findex < fLength) {

		Two_Sum(Q, fnow, Qnew, hh);

		fnow = f[++findex];
		Q = Qnew;

		if (hh != 0.0)
			h[hindex++] = hh;

	}

	if ((Q != 0.0) || (hindex == 0))
		h[hindex++] = Q;

	return hindex;

};


/***************************************************************************** /
/*                                                                           */
/*  scale_expansion_zeroelim()   Multiply an expansion by a scalar,          */
/*                               eliminating zero components from the        */
/*                               output expansion.                           */
/*                                                                           */
/*  Sets h = b*e.											                 */
/*                                                                           */
/*****************************************************************************/
int Scale_Expansion_Zero_Elimination(int const elength, double * e, double b, double * h) {


	/*****************************************************************************/
	/*                                                                           */
	/*    - INEXACT == volatile													 */
	/*                                                                           */
	/*****************************************************************************/
	double Q, sum;
	double product1;
	//double bvirt;
	//double c;
	//double abig;

	/*****************************************************************************/
	/*                                                                           */
	/*    - EXACT == nothing													 */
	/*                                                                           */
	/*****************************************************************************/
	double hh;
	double product0;
	double enow;
	//double avirt, bround, around;
	//double ahi, alo;
	double bhi, blo;
	//double error1, error2, error3;

	int hindex;


	Split(b, bhi, blo);
	Two_Product_Presplit(e[0], b, bhi, blo, Q, hh);

	hindex = 0;

	if (hh != 0)
		h[hindex++] = hh;

	for (int eindex = 1; eindex < elength; eindex++) {

		enow = e[eindex];

		Two_Product_Presplit(enow, b, bhi, blo, product1, product0);
		Two_Sum(Q, product0, sum, hh);

		if (hh != 0)
			h[hindex++] = hh;

		Fast_Two_Sum(product1, sum, Q, hh);

		if (hh != 0)
			h[hindex++] = hh;

	}

	if ((Q != 0.0) || (hindex == 0))
		h[hindex++] = Q;

	return hindex;

}