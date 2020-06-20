#pragma once



#include "ExactArithmetic.h"


namespace GeometricKernel {


	template<GeometricPredicatesArithmetic Arithmetic>
	double GeometricPredicates<Arithmetic>::orientation(v_pointer const & a, v_pointer const & b, v_pointer const & c) {


		StatisticsOrientation++;

		double const determinantLeft = (a->x - c->x) * (b->y - c->y);
		double const determinantRight = (a->y - c->y) * (b->x - c->x);

		double const determinant = determinantLeft - determinantRight;


		/*****************************************************************************/
		/*                                                                           */
		/*    - Fast evaluation														 */
		/*                                                                           */
		/*****************************************************************************/
		if (Arithmetic == GeometricPredicatesArithmetic::Fast)
			return determinant;


		/*****************************************************************************/
		/*                                                                           */
		/*    - Robust evaluation													 */
		/*                                                                           */
		/*****************************************************************************/
		StatisticsAdaptOrientation++;

		double determinantSum = 0.0;

		if (determinantLeft > 0.0) {

			if (determinantRight <= 0.0)
				return determinant;
			else
				determinantSum = determinantLeft + determinantRight;

		}
		else if (determinantLeft < 0.0) {

			if (determinantRight >= 0.0)
				return determinant;
			else
				determinantSum = -determinantLeft - determinantRight;

		}
		else {

			return determinant;

		}

		double const errorBound = orientationErrorBoundA * determinantSum;

		if (abs(determinant) >= errorBound)
			return determinant;

		return orientation_adapt(a, b, c, determinantSum);

	};
	template<GeometricPredicatesArithmetic Arithmetic>
	double GeometricPredicates<Arithmetic>::in_circumcircle(t_pointer const & t, v_pointer const & v)  {


		StatisticsInCircle++;

		const double dx = v->x;
		const double dy = v->y;

		double const adx = t->vertices[0]->x - dx;
		double const bdx = t->vertices[1]->x - dx;
		double const cdx = t->vertices[2]->x - dx;

		double const ady = t->vertices[0]->y - dy;
		double const bdy = t->vertices[1]->y - dy;
		double const cdy = t->vertices[2]->y - dy;


		double const bdxcdy = bdx * cdy;
		double const cdxbdy = cdx * bdy;
		double const alift  = adx * adx + ady * ady;

		double const cdxady = cdx * ady;
		double const adxcdy = adx * cdy;
		double const blift  = bdx * bdx + bdy * bdy;

		double const adxbdy = adx * bdy;
		double const bdxady = bdx * ady;
		double const clift  = cdx * cdx + cdy * cdy;

		const double determinant = alift * (bdxcdy - cdxbdy) +
								   blift * (cdxady - adxcdy) +
								   clift * (adxbdy - bdxady);


		/*****************************************************************************/
		/*                                                                           */
		/*    - Fast evaluation														 */
		/*                                                                           */
		/*****************************************************************************/
		if (Arithmetic == GeometricPredicatesArithmetic::Fast)
			return determinant;


		/*****************************************************************************/
		/*                                                                           */
		/*    - Robust evaluation													 */
		/*                                                                           */
		/*****************************************************************************/
		StatisticsAdaptInCircle++;

		const double permanent = alift * (abs(bdxcdy) + abs(cdxbdy)) +
								 blift * (abs(cdxady) + abs(adxcdy)) +
								 clift * (abs(adxbdy) + abs(bdxady));

		const double errorBound = inCircumCircleErrorBoundA * permanent;

		if (abs(determinant) > errorBound)
			return determinant;


		return in_circumcircle_adapt(t, v, permanent);

	};


	template<GeometricPredicatesArithmetic Arithmetic>
	void GeometricPredicates<Arithmetic>::get_edge_intersection(v_pointer const & p1, v_pointer const & p2, v_pointer const & q1, v_pointer const & q2, double & x, double & y) {


		double const p1p2x = p1->x - p2->x;
		double const p1q1x = p1->x - q1->x;
		double const q1q2x = q1->x - q2->x;

		double const p1p2y = p1->y - p2->y;
		double const p1q1y = p1->y - q1->y;
		double const q1q2y = q1->y - q2->y;

		double const t = (q1q2y * p1q1x - q1q2x * p1q1y) / (p1p2x * q1q2y - q1q2x * p1p2y);


		x = p1->x - (t * p1p2x);
		y = p1->y - (t * p1p2y);

	}



	/*****************************************************************************/
	/*                                                                           */
	/*    - Return a positive value if the vertices 'a', 'b', 'c' are oriented   */
	/*      CounterClockwise										             */
	/*                                                                           */
	/*****************************************************************************/
	template<GeometricPredicatesArithmetic Arithmetic>
	double GeometricPredicates<Arithmetic>::orientation_fast(v_pointer const a, v_pointer const b, v_pointer const c) const {

		double const determinant = (a->x - c->x) * (b->y - c->y) - (a->y - c->y) * (b->x - c->x);

		if (abs(determinant) < EpsilonOrientation)
			return 0.0;

		return determinant;

	};
	template<GeometricPredicatesArithmetic Arithmetic>
	double GeometricPredicates<Arithmetic>::orientation_robust(v_pointer const a, v_pointer const b, v_pointer const c) const {


		/*****************************************************************************/
		/*                                                                           */
		/*    - Fast evaluation														 */
		/*                                                                           */
		/*****************************************************************************/
		double const determinantLeft = (a->x - c->x) * (b->y - c->y);
		double const determinantRight = (a->y - c->y) * (b->x - c->x);

		double const determinant = determinantLeft - determinantRight;

		//if (Arithmetic == GeometricPredicatesArithmetic::Fast)
		//	return determinant;

		/*****************************************************************************/
		/*                                                                           */
		/*    - Exact evaluation													 */
		/*                                                                           */
		/*****************************************************************************/
		double determinantSum = 0.0;

		if (determinantLeft > 0.0) {

			if (determinantRight <= 0.0)
				return determinant;
			else
				determinantSum = determinantLeft + determinantRight;

		}
		else if (determinantLeft < 0.0) {

			if (determinantRight >= 0.0)
				return determinant;
			else
				determinantSum = -determinantLeft - determinantRight;

		}
		else {

			return determinant;

		}

		double const errorBound = orientationErrorBoundA * determinantSum;

		if (abs(determinant) >= errorBound)
			return determinant;

		return orientation_adapt(a, b, c, determinantSum);

	};
	template<GeometricPredicatesArithmetic Arithmetic>
	double GeometricPredicates<Arithmetic>::orientation_adapt(v_pointer const a, v_pointer const b, v_pointer const c, double const determinantSum)  {


		int lengthC1, lengthC2, lengthD;

		double determinantLeft_tail;
		double determinantRight_tail;

		double determinant;
		double errorBound;

		double acX_tail;
		double acY_tail;
		double bcX_tail;
		double bcY_tail;

		double B[4], C1[8], C2[12], D[16];
		double u[4];

		double s0, t0;

		double determinantLeft;
		double determinantRight;


		double s1, t1;

		double acX = a->x - c->x;
		double acY = a->y - c->y;
		double bcX = b->x - c->x;
		double bcY = b->y - c->y;





		Two_Product(acX, bcY, determinantLeft, determinantLeft_tail);
		Two_Product(acY, bcX, determinantRight, determinantRight_tail);

		Two_Two_Diff(determinantLeft, determinantLeft_tail, determinantRight, determinantRight_tail, B[3], B[2], B[1], B[0]);

		determinant = Estimate(4, B);
		//determinant = B[0] + B[1] + B[2] + B[3];
		errorBound = orientationErrorBoundB * determinantSum;

		/*****************************************************************************/
		/*                                                                           */
		/*    - First test of the reliability of the sign							 */
		/*                                                                           */
		/*****************************************************************************/
		if (abs(determinant) >= errorBound)
			return determinant;



		Two_Diff_Tail(a->x, c->x, acX, acX_tail);
		Two_Diff_Tail(b->x, c->x, bcX, bcX_tail);
		Two_Diff_Tail(a->y, c->y, acY, acY_tail);
		Two_Diff_Tail(b->y, c->y, bcY, bcY_tail);

		if ((acX_tail == 0.0) && (acY_tail == 0.0) && (bcX_tail == 0.0) && (bcY_tail == 0.0))
			return determinant;


		errorBound   = orientationErrorBoundC * determinantSum + resultErrorBound * abs(determinant);
		determinant += (acX * bcY_tail + bcY * acX_tail)
					 - (acY * bcX_tail + bcX * acY_tail);

		/*****************************************************************************/
		/*                                                                           */
		/*    - Second test of the reliability of the sign							 */
		/*                                                                           */
		/*****************************************************************************/
		if (abs(determinant) >= errorBound)
			return determinant;



		/*****************************************************************************/
		/*                                                                           */
		/*    - Exact evaluation													 */
		/*                                                                           */
		/*****************************************************************************/
		StatisticsExactOrientation++;

		Two_Product(acX_tail, bcY, s1, s0);
		Two_Product(acY_tail, bcX, t1, t0);

		Two_Two_Diff(s1, s0, t1, t0, u[3], u[2], u[1], u[0]);

		lengthC1 = Fast_Expansion_Sum_Zero_Eliminiation(4, B, 4, u, C1);



		Two_Product(acX, bcY_tail, s1, s0);
		Two_Product(acY, bcX_tail, t1, t0);

		Two_Two_Diff(s1, s0, t1, t0, u[3], u[2], u[1], u[0]);

		lengthC2 = Fast_Expansion_Sum_Zero_Eliminiation(lengthC1, C1, 4, u, C2);


		Two_Product(acX_tail, bcY_tail, s1, s0);
		Two_Product(acY_tail, bcX_tail, t1, t0);

		Two_Two_Diff(s1, s0, t1, t0, u[3], u[2], u[1], u[0]);

		lengthD = Fast_Expansion_Sum_Zero_Eliminiation(lengthC2, C2, 4, u, D);


		return D[lengthD - 1];

	};


	/*****************************************************************************/
	/*                                                                           */
	/*    - Return a positive value if the vertex 'v' lies inside the            */
	/*      circum-circle of the triangle							             */
	/*                                                                           */
	/*****************************************************************************/
	template<GeometricPredicatesArithmetic Arithmetic>
	double GeometricPredicates<Arithmetic>::in_circumcircle_fast(t_pointer const t, v_pointer const v) const {

		const double vx = v->x;
		const double vy = v->y;

		double const avX = t->vertices[0]->x - vx;
		double const bvX = t->vertices[1]->x - vx;
		double const cvX = t->vertices[2]->x - vx;

		double const avY = t->vertices[0]->y - vy;
		double const bvY = t->vertices[1]->y - vy;
		double const cvY = t->vertices[2]->y - vy;


		double const bvXcvY = bvX * cvY;
		double const cvXbvY = cvX * bvY;
		double const aLift = avX * avX + avY * avY;

		double const cvXavY = cvX * avY;
		double const avXcvY = avX * cvY;
		double const bLift = bvX * bvX + bvY * bvY;

		double const avXbvY = avX * bvY;
		double const bvXavY = bvX * avY;
		double const cLift = cvX * cvX + cvY * cvY;


		return aLift * (bvXcvY - cvXbvY) + bLift * (cvXavY - avXcvY) + cLift * (avXbvY - bvXavY);

	};
	template<GeometricPredicatesArithmetic Arithmetic>
	double GeometricPredicates<Arithmetic>::in_circumcircle_robust(t_pointer const t, v_pointer const v) const {


		/*****************************************************************************/
		/*                                                                           */
		/*    - Fast evaluation														 */
		/*                                                                           */
		/*****************************************************************************/
		const double vx = v->x;
		const double vy = v->y;

		double const avX = t->vertices[0]->x - vx;
		double const bvX = t->vertices[1]->x - vx;
		double const cvX = t->vertices[2]->x - vx;

		double const avY = t->vertices[0]->y - vy;
		double const bvY = t->vertices[1]->y - vy;
		double const cvY = t->vertices[2]->y - vy;


		double const bvXcvY = bvX * cvY;
		double const cvXbvY = cvX * bvY;
		double const aLift = avX * avX + avY * avY;

		double const cvXavY = cvX * avY;
		double const avXcvY = avX * cvY;
		double const bLift = bvX * bvX + bvY * bvY;

		double const avXbvY = avX * bvY;
		double const bvXavY = bvX * avY;
		double const cLift = cvX * cvX + cvY * cvY;

		const double determinant = aLift * (bvXcvY - cvXbvY)
			+ bLift * (cvXavY - avXcvY)
			+ cLift * (avXbvY - bvXavY);

		//if (Arithmetic == GeometricPredicatesArithmetic::Fast)
		//	return determinant;


		/*****************************************************************************/
		/*                                                                           */
		/*    - Exact evaluation													 */
		/*                                                                           */
		/*****************************************************************************/
		const double permanent = (abs(bvXcvY) + abs(cvXbvY)) * aLift
			+ (abs(cvXavY) + abs(avXcvY)) * bLift
			+ (abs(avXbvY) + abs(bvXavY)) * cLift;

		const double errorBound = inCircumCircleErrorBoundA * permanent;

		if (abs(determinant) > errorBound)
			return determinant;


		return in_circumcircle_adapt(t, v, permanent);

	};
	template<GeometricPredicatesArithmetic Arithmetic>
	double GeometricPredicates<Arithmetic>::in_circumcircle_adapt(t_pointer const t, v_pointer const p, double const permanent) {


		/*****************************************************************************/
		/*                                                                           */
		/*    - INEXACT == volatile													 */
		/*                                                                           */
		/*****************************************************************************/
		double dx, dy;
		double adx, bdx, cdx;
		double ady, bdy, cdy;

		double bdxcdy1, cdxbdy1, cdxady1, adxcdy1, adxbdy1, bdxady1;
		double adxadx1, adyady1, bdxbdx1, bdybdy1, cdxcdx1, cdycdy1;

		double bc3, ca3, ab3;
		double aa3, bb3, cc3;
		double ti1, tj1;
		double u3, v3;

		double abtt3, bctt3, catt3;

		double negate;



		/*****************************************************************************/
		/*                                                                           */
		/*    - EXACT == nothing													 */
		/*                                                                           */
		/*****************************************************************************/
		double determinant;
		double errorBound;

		double adxtail, bdxtail, cdxtail, adytail, bdytail, cdytail;

		double bdxcdy0, cdxbdy0, cdxady0, adxcdy0, adxbdy0, bdxady0;
		double adxadx0, adyady0, bdxbdx0, bdybdy0, cdxcdx0, cdycdy0;

		double axbc[8], axxbc[16], aybc[8], ayybc[16], adet[32];
		double bxca[8], bxxca[16], byca[8], byyca[16], bdet[32];
		double cxab[8], cxxab[16], cyab[8], cyyab[16], cdet[32];
		double abdet[64];
		double fin1[1152], fin2[1152];

		double temp8[8], temp16a[16], temp16b[16], temp16c[16];
		double temp32a[32], temp32b[32], temp48[48], temp64[64];
		double axtbb[8], axtcc[8], aytbb[8], aytcc[8];
		double bxtaa[8], bxtcc[8], bytaa[8], bytcc[8];
		double cxtaa[8], cxtbb[8], cytaa[8], cytbb[8];
		double axtbc[8], aytbc[8], bxtca[8], bytca[8], cxtab[8], cytab[8];
		double axtbct[16], aytbct[16], bxtcat[16], bytcat[16], cxtabt[16], cytabt[16];
		double axtbctt[8], aytbctt[8], bxtcatt[8];
		double bytcatt[8], cxtabtt[8], cytabtt[8];

		double abt[8], bct[8], cat[8];
		double abtt[4], bctt[4], catt[4];


		double bc[4], ca[4], ab[4];
		double aa[4], bb[4], cc[4];
		double ti0, tj0;
		double u[4], v[4];

		double * finnow, *finother, *finswap;


		int axbclen, axxbclen, aybclen, ayybclen, alen;
		int bxcalen, bxxcalen, bycalen, byycalen, blen;
		int cxablen, cxxablen, cyablen, cyyablen, clen;
		int ablen;
		int finlength;

		int temp8len, temp16alen, temp16blen, temp16clen;
		int temp32alen, temp32blen, temp48len, temp64len;
		int axtbblen, axtcclen, aytbblen, aytcclen;
		int bxtaalen, bxtcclen, bytaalen, bytcclen;
		int cxtaalen, cxtbblen, cytaalen, cytbblen;
		int axtbclen, aytbclen, bxtcalen, bytcalen, cxtablen, cytablen;
		int axtbctlen, aytbctlen, bxtcatlen, bytcatlen, cxtabtlen, cytabtlen;
		int axtbcttlen, aytbcttlen, bxtcattlen, bytcattlen, cxtabttlen, cytabttlen;

		int abtlen, bctlen, catlen;
		int abttlen, bcttlen, cattlen;





		dx = p->x;
		dy = p->y;

		adx = t->vertices[0]->x - dx;
		bdx = t->vertices[1]->x - dx;
		cdx = t->vertices[2]->x - dx;

		ady = t->vertices[0]->y - dy;
		bdy = t->vertices[1]->y - dy;
		cdy = t->vertices[2]->y - dy;



		Two_Product(bdx, cdy, bdxcdy1, bdxcdy0);
		Two_Product(cdx, bdy, cdxbdy1, cdxbdy0);
		Two_Two_Diff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0, bc[3], bc[2], bc[1], bc[0]);

		axbclen = Scale_Expansion_Zero_Elimination(4, bc, adx, axbc);
		axxbclen = Scale_Expansion_Zero_Elimination(axbclen, axbc, adx, axxbc);
		aybclen = Scale_Expansion_Zero_Elimination(4, bc, ady, aybc);
		ayybclen = Scale_Expansion_Zero_Elimination(aybclen, aybc, ady, ayybc);

		alen = Fast_Expansion_Sum_Zero_Eliminiation(axxbclen, axxbc, ayybclen, ayybc, adet);

		Two_Product(cdx, ady, cdxady1, cdxady0);
		Two_Product(adx, cdy, adxcdy1, adxcdy0);
		Two_Two_Diff(cdxady1, cdxady0, adxcdy1, adxcdy0, ca[3], ca[2], ca[1], ca[0]);

		bxcalen = Scale_Expansion_Zero_Elimination(4, ca, bdx, bxca);
		bxxcalen = Scale_Expansion_Zero_Elimination(bxcalen, bxca, bdx, bxxca);
		bycalen = Scale_Expansion_Zero_Elimination(4, ca, bdy, byca);
		byycalen = Scale_Expansion_Zero_Elimination(bycalen, byca, bdy, byyca);

		blen = Fast_Expansion_Sum_Zero_Eliminiation(bxxcalen, bxxca, byycalen, byyca, bdet);

		Two_Product(adx, bdy, adxbdy1, adxbdy0);
		Two_Product(bdx, ady, bdxady1, bdxady0);
		Two_Two_Diff(adxbdy1, adxbdy0, bdxady1, bdxady0, ab[3], ab[2], ab[1], ab[0]);

		cxablen = Scale_Expansion_Zero_Elimination(4, ab, cdx, cxab);
		cxxablen = Scale_Expansion_Zero_Elimination(cxablen, cxab, cdx, cxxab);
		cyablen = Scale_Expansion_Zero_Elimination(4, ab, cdy, cyab);
		cyyablen = Scale_Expansion_Zero_Elimination(cyablen, cyab, cdy, cyyab);

		clen = Fast_Expansion_Sum_Zero_Eliminiation(cxxablen, cxxab, cyyablen, cyyab, cdet);

		ablen = Fast_Expansion_Sum_Zero_Eliminiation(alen, adet, blen, bdet, abdet);
		finlength = Fast_Expansion_Sum_Zero_Eliminiation(ablen, abdet, clen, cdet, fin1);



		determinant = Estimate(finlength, fin1);

		errorBound = inCircumCircleErrorBoundB * permanent;

		if (abs(determinant) >= errorBound)
			return determinant;


		Two_Diff_Tail(t->vertices[0]->x, dx, adx, adxtail);
		Two_Diff_Tail(t->vertices[0]->y, dy, ady, adytail);
		Two_Diff_Tail(t->vertices[1]->x, dx, bdx, bdxtail);
		Two_Diff_Tail(t->vertices[1]->y, dy, bdy, bdytail);
		Two_Diff_Tail(t->vertices[2]->x, dx, cdx, cdxtail);
		Two_Diff_Tail(t->vertices[2]->y, dy, cdy, cdytail);



		if ((adxtail == 0.0) && (bdxtail == 0.0) && (cdxtail == 0.0) && (adytail == 0.0) && (bdytail == 0.0) && (cdytail == 0.0))
			return determinant;


		errorBound = inCircumCircleErrorBoundC * permanent + resultErrorBound * abs(determinant);

		double const determinantPart1 = (adx * adx + ady * ady) * ((bdx * cdytail + cdy * bdxtail) - (bdy * cdxtail + cdx * bdytail)) + 2.0 * (adx * adxtail + ady * adytail) * (bdx * cdy - bdy * cdx);
		double const determinantPart2 = (bdx * bdx + bdy * bdy) * ((cdx * adytail + ady * cdxtail) - (cdy * adxtail + adx * cdytail)) + 2.0 * (bdx * bdxtail + bdy * bdytail) * (cdx * ady - cdy * adx);
		double const determinantPart3 = (cdx * cdx + cdy * cdy) * ((adx * bdytail + bdy * adxtail) - (ady * bdxtail + bdx * adytail)) + 2.0 * (cdx * cdxtail + cdy * cdytail) * (adx * bdy - ady * bdx);

		determinant += determinantPart1 + determinantPart2 + determinantPart3;

		if (abs(determinant) >= errorBound)
			return determinant;


		/*****************************************************************************/
		/*                                                                           */
		/*    - Exact evaluation													 */
		/*                                                                           */
		/*****************************************************************************/
		StatisticsExactInCircle++;

		finnow = fin1;
		finother = fin2;

		if ((bdxtail != 0.0) || (bdytail != 0.0) || (cdxtail != 0.0) || (cdytail != 0.0)) {

			Square(adx, adxadx1, adxadx0);
			Square(ady, adyady1, adyady0);

			Two_Two_Sum(adxadx1, adxadx0, adyady1, adyady0, aa[3], aa[2], aa[1], aa[0]);
			//aa[3] = a3;

		}
		if ((cdxtail != 0.0) || (cdytail != 0.0) || (adxtail != 0.0) || (adytail != 0.0)) {

			Square(bdx, bdxbdx1, bdxbdx0);
			Square(bdy, bdybdy1, bdybdy0);

			Two_Two_Sum(bdxbdx1, bdxbdx0, bdybdy1, bdybdy0, bb[3], bb[2], bb[1], bb[0]);
			//bb[3] = bb3;

		}
		if ((adxtail != 0.0) || (adytail != 0.0) || (bdxtail != 0.0) || (bdytail != 0.0)) {

			Square(cdx, cdxcdx1, cdxcdx0);
			Square(cdy, cdycdy1, cdycdy0);

			Two_Two_Sum(cdxcdx1, cdxcdx0, cdycdy1, cdycdy0, cc[3], cc[2], cc[1], cc[0]);
			//cc[3] = cc3;

		}

		if (adxtail != 0.0) {

			axtbclen = Scale_Expansion_Zero_Elimination(4, bc, adxtail, axtbc);
			temp16alen = Scale_Expansion_Zero_Elimination(axtbclen, axtbc, 2.0 * adx, temp16a);

			axtcclen = Scale_Expansion_Zero_Elimination(4, cc, adxtail, axtcc);
			temp16blen = Scale_Expansion_Zero_Elimination(axtcclen, axtcc, bdy, temp16b);

			axtbblen = Scale_Expansion_Zero_Elimination(4, bb, adxtail, axtbb);
			temp16clen = Scale_Expansion_Zero_Elimination(axtbblen, axtbb, -cdy, temp16c);

			temp32alen = Fast_Expansion_Sum_Zero_Eliminiation(temp16alen, temp16a, temp16blen, temp16b, temp32a);
			temp48len = Fast_Expansion_Sum_Zero_Eliminiation(temp16clen, temp16c, temp32alen, temp32a, temp48);
			finlength = Fast_Expansion_Sum_Zero_Eliminiation(finlength, finnow, temp48len, temp48, finother);

			finswap = finnow;
			finnow = finother;
			finother = finswap;

		}
		if (adytail != 0.0) {

			aytbclen = Scale_Expansion_Zero_Elimination(4, bc, adytail, aytbc);
			temp16alen = Scale_Expansion_Zero_Elimination(aytbclen, aytbc, 2.0 * ady, temp16a);

			aytbblen = Scale_Expansion_Zero_Elimination(4, bb, adytail, aytbb);
			temp16blen = Scale_Expansion_Zero_Elimination(aytbblen, aytbb, cdx, temp16b);

			aytcclen = Scale_Expansion_Zero_Elimination(4, cc, adytail, aytcc);
			temp16clen = Scale_Expansion_Zero_Elimination(aytcclen, aytcc, -bdx, temp16c);

			temp32alen = Fast_Expansion_Sum_Zero_Eliminiation(temp16alen, temp16a, temp16blen, temp16b, temp32a);
			temp48len = Fast_Expansion_Sum_Zero_Eliminiation(temp16clen, temp16c, temp32alen, temp32a, temp48);
			finlength = Fast_Expansion_Sum_Zero_Eliminiation(finlength, finnow, temp48len, temp48, finother);

			finswap = finnow;
			finnow = finother;
			finother = finswap;

		}
		if (bdxtail != 0.0) {

			bxtcalen = Scale_Expansion_Zero_Elimination(4, ca, bdxtail, bxtca);
			temp16alen = Scale_Expansion_Zero_Elimination(bxtcalen, bxtca, 2.0 * bdx, temp16a);

			bxtaalen = Scale_Expansion_Zero_Elimination(4, aa, bdxtail, bxtaa);
			temp16blen = Scale_Expansion_Zero_Elimination(bxtaalen, bxtaa, cdy, temp16b);

			bxtcclen = Scale_Expansion_Zero_Elimination(4, cc, bdxtail, bxtcc);
			temp16clen = Scale_Expansion_Zero_Elimination(bxtcclen, bxtcc, -ady, temp16c);

			temp32alen = Fast_Expansion_Sum_Zero_Eliminiation(temp16alen, temp16a, temp16blen, temp16b, temp32a);
			temp48len = Fast_Expansion_Sum_Zero_Eliminiation(temp16clen, temp16c, temp32alen, temp32a, temp48);
			finlength = Fast_Expansion_Sum_Zero_Eliminiation(finlength, finnow, temp48len, temp48, finother);

			finswap = finnow;
			finnow = finother;
			finother = finswap;

		}
		if (bdytail != 0.0) {

			bytcalen = Scale_Expansion_Zero_Elimination(4, ca, bdytail, bytca);
			temp16alen = Scale_Expansion_Zero_Elimination(bytcalen, bytca, 2.0 * bdy, temp16a);

			bytcclen = Scale_Expansion_Zero_Elimination(4, cc, bdytail, bytcc);
			temp16blen = Scale_Expansion_Zero_Elimination(bytcclen, bytcc, adx, temp16b);

			bytaalen = Scale_Expansion_Zero_Elimination(4, aa, bdytail, bytaa);
			temp16clen = Scale_Expansion_Zero_Elimination(bytaalen, bytaa, -cdx, temp16c);

			temp32alen = Fast_Expansion_Sum_Zero_Eliminiation(temp16alen, temp16a, temp16blen, temp16b, temp32a);
			temp48len = Fast_Expansion_Sum_Zero_Eliminiation(temp16clen, temp16c, temp32alen, temp32a, temp48);
			finlength = Fast_Expansion_Sum_Zero_Eliminiation(finlength, finnow, temp48len, temp48, finother);

			finswap = finnow;
			finnow = finother;
			finother = finswap;

		}
		if (cdxtail != 0.0) {

			cxtablen = Scale_Expansion_Zero_Elimination(4, ab, cdxtail, cxtab);
			temp16alen = Scale_Expansion_Zero_Elimination(cxtablen, cxtab, 2.0 * cdx, temp16a);

			cxtbblen = Scale_Expansion_Zero_Elimination(4, bb, cdxtail, cxtbb);
			temp16blen = Scale_Expansion_Zero_Elimination(cxtbblen, cxtbb, ady, temp16b);

			cxtaalen = Scale_Expansion_Zero_Elimination(4, aa, cdxtail, cxtaa);
			temp16clen = Scale_Expansion_Zero_Elimination(cxtaalen, cxtaa, -bdy, temp16c);

			temp32alen = Fast_Expansion_Sum_Zero_Eliminiation(temp16alen, temp16a, temp16blen, temp16b, temp32a);
			temp48len = Fast_Expansion_Sum_Zero_Eliminiation(temp16clen, temp16c, temp32alen, temp32a, temp48);
			finlength = Fast_Expansion_Sum_Zero_Eliminiation(finlength, finnow, temp48len, temp48, finother);

			finswap = finnow;
			finnow = finother;
			finother = finswap;

		}
		if (cdytail != 0.0) {

			cytablen = Scale_Expansion_Zero_Elimination(4, ab, cdytail, cytab);
			temp16alen = Scale_Expansion_Zero_Elimination(cytablen, cytab, 2.0 * cdy, temp16a);

			cytaalen = Scale_Expansion_Zero_Elimination(4, aa, cdytail, cytaa);
			temp16blen = Scale_Expansion_Zero_Elimination(cytaalen, cytaa, bdx, temp16b);

			cytbblen = Scale_Expansion_Zero_Elimination(4, bb, cdytail, cytbb);
			temp16clen = Scale_Expansion_Zero_Elimination(cytbblen, cytbb, -adx, temp16c);

			temp32alen = Fast_Expansion_Sum_Zero_Eliminiation(temp16alen, temp16a, temp16blen, temp16b, temp32a);
			temp48len = Fast_Expansion_Sum_Zero_Eliminiation(temp16clen, temp16c, temp32alen, temp32a, temp48);
			finlength = Fast_Expansion_Sum_Zero_Eliminiation(finlength, finnow, temp48len, temp48, finother);

			finswap = finnow;
			finnow = finother;
			finother = finswap;

		}
		if ((adxtail != 0.0) || (adytail != 0.0)) {

			if ((bdxtail != 0.0) || (bdytail != 0.0) || (cdxtail != 0.0) || (cdytail != 0.0)) {

				Two_Product(bdxtail, cdy, ti1, ti0);
				Two_Product(bdx, cdytail, tj1, tj0);
				Two_Two_Sum(ti1, ti0, tj1, tj0, u[3], u[2], u[1], u[0]);
				//u[3] = u3;
				negate = -bdy;

				Two_Product(cdxtail, negate, ti1, ti0);
				negate = -bdytail;

				Two_Product(cdx, negate, tj1, tj0);
				Two_Two_Sum(ti1, ti0, tj1, tj0, v[3], v[2], v[1], v[0]);
				//v[3] = v3;

				bctlen = Fast_Expansion_Sum_Zero_Eliminiation(4, u, 4, v, bct);

				Two_Product(bdxtail, cdytail, ti1, ti0);
				Two_Product(cdxtail, bdytail, tj1, tj0);
				Two_Two_Diff(ti1, ti0, tj1, tj0, bctt[3], bctt[2], bctt[1], bctt[0]);
				//bctt[3] = bctt3;

				bcttlen = 4;

			}
			else {

				bct[0] = 0.0;
				bctlen = 1;

				bctt[0] = 0.0;
				bcttlen = 1;

			}

			if (adxtail != 0.0) {

				temp16alen = Scale_Expansion_Zero_Elimination(axtbclen, axtbc, adxtail, temp16a);
				axtbctlen = Scale_Expansion_Zero_Elimination(bctlen, bct, adxtail, axtbct);
				temp32alen = Scale_Expansion_Zero_Elimination(axtbctlen, axtbct, 2.0 * adx, temp32a);
				temp48len = Fast_Expansion_Sum_Zero_Eliminiation(temp16alen, temp16a, temp32alen, temp32a, temp48);
				finlength = Fast_Expansion_Sum_Zero_Eliminiation(finlength, finnow, temp48len, temp48, finother);

				finswap = finnow;
				finnow = finother;
				finother = finswap;

				if (bdytail != 0.0) {

					temp8len = Scale_Expansion_Zero_Elimination(4, cc, adxtail, temp8);
					temp16alen = Scale_Expansion_Zero_Elimination(temp8len, temp8, bdytail, temp16a);
					finlength = Fast_Expansion_Sum_Zero_Eliminiation(finlength, finnow, temp16alen, temp16a, finother);

					finswap = finnow;
					finnow = finother;
					finother = finswap;

				}
				if (cdytail != 0.0) {

					temp8len = Scale_Expansion_Zero_Elimination(4, bb, -adxtail, temp8);
					temp16alen = Scale_Expansion_Zero_Elimination(temp8len, temp8, cdytail, temp16a);
					finlength = Fast_Expansion_Sum_Zero_Eliminiation(finlength, finnow, temp16alen, temp16a, finother);

					finswap = finnow;
					finnow = finother;
					finother = finswap;

				}

				temp32alen = Scale_Expansion_Zero_Elimination(axtbctlen, axtbct, adxtail, temp32a);
				axtbcttlen = Scale_Expansion_Zero_Elimination(bcttlen, bctt, adxtail, axtbctt);
				temp16alen = Scale_Expansion_Zero_Elimination(axtbcttlen, axtbctt, 2.0 * adx, temp16a);
				temp16blen = Scale_Expansion_Zero_Elimination(axtbcttlen, axtbctt, adxtail, temp16b);
				temp32blen = Fast_Expansion_Sum_Zero_Eliminiation(temp16alen, temp16a, temp16blen, temp16b, temp32b);
				temp64len = Fast_Expansion_Sum_Zero_Eliminiation(temp32alen, temp32a, temp32blen, temp32b, temp64);
				finlength = Fast_Expansion_Sum_Zero_Eliminiation(finlength, finnow, temp64len, temp64, finother);

				finswap = finnow;
				finnow = finother;
				finother = finswap;

			}
			if (adytail != 0.0) {

				temp16alen = Scale_Expansion_Zero_Elimination(aytbclen, aytbc, adytail, temp16a);
				aytbctlen = Scale_Expansion_Zero_Elimination(bctlen, bct, adytail, aytbct);
				temp32alen = Scale_Expansion_Zero_Elimination(aytbctlen, aytbct, 2.0 * ady, temp32a);
				temp48len = Fast_Expansion_Sum_Zero_Eliminiation(temp16alen, temp16a, temp32alen, temp32a, temp48);
				finlength = Fast_Expansion_Sum_Zero_Eliminiation(finlength, finnow, temp48len, temp48, finother);

				finswap = finnow;
				finnow = finother;
				finother = finswap;

				temp32alen = Scale_Expansion_Zero_Elimination(aytbctlen, aytbct, adytail, temp32a);
				aytbcttlen = Scale_Expansion_Zero_Elimination(bcttlen, bctt, adytail, aytbctt);
				temp16alen = Scale_Expansion_Zero_Elimination(aytbcttlen, aytbctt, 2.0 * ady, temp16a);
				temp16blen = Scale_Expansion_Zero_Elimination(aytbcttlen, aytbctt, adytail, temp16b);
				temp32blen = Fast_Expansion_Sum_Zero_Eliminiation(temp16alen, temp16a, temp16blen, temp16b, temp32b);
				temp64len = Fast_Expansion_Sum_Zero_Eliminiation(temp32alen, temp32a, temp32blen, temp32b, temp64);
				finlength = Fast_Expansion_Sum_Zero_Eliminiation(finlength, finnow, temp64len, temp64, finother);

				finswap = finnow;
				finnow = finother;
				finother = finswap;

			}
		}
		if ((bdxtail != 0.0) || (bdytail != 0.0)) {

			if ((cdxtail != 0.0) || (cdytail != 0.0) || (adxtail != 0.0) || (adytail != 0.0)) {

				Two_Product(cdxtail, ady, ti1, ti0);
				Two_Product(cdx, adytail, tj1, tj0);
				Two_Two_Sum(ti1, ti0, tj1, tj0, u[3], u[2], u[1], u[0]);
				//u[3] = u3;
				negate = -cdy;

				Two_Product(adxtail, negate, ti1, ti0);
				negate = -cdytail;

				Two_Product(adx, negate, tj1, tj0);
				Two_Two_Sum(ti1, ti0, tj1, tj0, v[3], v[2], v[1], v[0]);
				//v[3] = v3;

				catlen = Fast_Expansion_Sum_Zero_Eliminiation(4, u, 4, v, cat);

				Two_Product(cdxtail, adytail, ti1, ti0);
				Two_Product(adxtail, cdytail, tj1, tj0);
				Two_Two_Diff(ti1, ti0, tj1, tj0, catt[3], catt[2], catt[1], catt[0]);
				//catt[3] = catt3;

				cattlen = 4;

			}
			else {

				cat[0] = 0.0;
				catlen = 1;

				catt[0] = 0.0;
				cattlen = 1;

			}

			if (bdxtail != 0.0) {

				temp16alen = Scale_Expansion_Zero_Elimination(bxtcalen, bxtca, bdxtail, temp16a);
				bxtcatlen = Scale_Expansion_Zero_Elimination(catlen, cat, bdxtail, bxtcat);
				temp32alen = Scale_Expansion_Zero_Elimination(bxtcatlen, bxtcat, 2.0 * bdx, temp32a);
				temp48len = Fast_Expansion_Sum_Zero_Eliminiation(temp16alen, temp16a, temp32alen, temp32a, temp48);
				finlength = Fast_Expansion_Sum_Zero_Eliminiation(finlength, finnow, temp48len, temp48, finother);

				finswap = finnow;
				finnow = finother;
				finother = finswap;

				if (cdytail != 0.0) {

					temp8len = Scale_Expansion_Zero_Elimination(4, aa, bdxtail, temp8);
					temp16alen = Scale_Expansion_Zero_Elimination(temp8len, temp8, cdytail, temp16a);
					finlength = Fast_Expansion_Sum_Zero_Eliminiation(finlength, finnow, temp16alen, temp16a, finother);

					finswap = finnow;
					finnow = finother;
					finother = finswap;

				}
				if (adytail != 0.0) {

					temp8len = Scale_Expansion_Zero_Elimination(4, cc, -bdxtail, temp8);
					temp16alen = Scale_Expansion_Zero_Elimination(temp8len, temp8, adytail, temp16a);
					finlength = Fast_Expansion_Sum_Zero_Eliminiation(finlength, finnow, temp16alen, temp16a, finother);

					finswap = finnow;
					finnow = finother;
					finother = finswap;

				}

				temp32alen = Scale_Expansion_Zero_Elimination(bxtcatlen, bxtcat, bdxtail, temp32a);
				bxtcattlen = Scale_Expansion_Zero_Elimination(cattlen, catt, bdxtail, bxtcatt);
				temp16alen = Scale_Expansion_Zero_Elimination(bxtcattlen, bxtcatt, 2.0 * bdx, temp16a);
				temp16blen = Scale_Expansion_Zero_Elimination(bxtcattlen, bxtcatt, bdxtail, temp16b);
				temp32blen = Fast_Expansion_Sum_Zero_Eliminiation(temp16alen, temp16a, temp16blen, temp16b, temp32b);
				temp64len = Fast_Expansion_Sum_Zero_Eliminiation(temp32alen, temp32a, temp32blen, temp32b, temp64);
				finlength = Fast_Expansion_Sum_Zero_Eliminiation(finlength, finnow, temp64len, temp64, finother);

				finswap = finnow;
				finnow = finother;
				finother = finswap;

			}
			if (bdytail != 0.0) {

				temp16alen = Scale_Expansion_Zero_Elimination(bytcalen, bytca, bdytail, temp16a);
				bytcatlen = Scale_Expansion_Zero_Elimination(catlen, cat, bdytail, bytcat);
				temp32alen = Scale_Expansion_Zero_Elimination(bytcatlen, bytcat, 2.0 * bdy, temp32a);
				temp48len = Fast_Expansion_Sum_Zero_Eliminiation(temp16alen, temp16a, temp32alen, temp32a, temp48);
				finlength = Fast_Expansion_Sum_Zero_Eliminiation(finlength, finnow, temp48len, temp48, finother);

				finswap = finnow;
				finnow = finother;
				finother = finswap;

				temp32alen = Scale_Expansion_Zero_Elimination(bytcatlen, bytcat, bdytail, temp32a);
				bytcattlen = Scale_Expansion_Zero_Elimination(cattlen, catt, bdytail, bytcatt);
				temp16alen = Scale_Expansion_Zero_Elimination(bytcattlen, bytcatt, 2.0 * bdy, temp16a);
				temp16blen = Scale_Expansion_Zero_Elimination(bytcattlen, bytcatt, bdytail, temp16b);
				temp32blen = Fast_Expansion_Sum_Zero_Eliminiation(temp16alen, temp16a, temp16blen, temp16b, temp32b);
				temp64len = Fast_Expansion_Sum_Zero_Eliminiation(temp32alen, temp32a, temp32blen, temp32b, temp64);
				finlength = Fast_Expansion_Sum_Zero_Eliminiation(finlength, finnow, temp64len, temp64, finother);

				finswap = finnow;
				finnow = finother;
				finother = finswap;

			}
		}
		if ((cdxtail != 0.0) || (cdytail != 0.0)) {

			if ((adxtail != 0.0) || (adytail != 0.0) || (bdxtail != 0.0) || (bdytail != 0.0)) {

				Two_Product(adxtail, bdy, ti1, ti0);
				Two_Product(adx, bdytail, tj1, tj0);
				Two_Two_Sum(ti1, ti0, tj1, tj0, u[3], u[2], u[1], u[0]);
				//u[3] = u3;
				negate = -ady;

				Two_Product(bdxtail, negate, ti1, ti0);
				negate = -adytail;

				Two_Product(bdx, negate, tj1, tj0);
				Two_Two_Sum(ti1, ti0, tj1, tj0, v[3], v[2], v[1], v[0]);
				//v[3] = v3;

				abtlen = Fast_Expansion_Sum_Zero_Eliminiation(4, u, 4, v, abt);

				Two_Product(adxtail, bdytail, ti1, ti0);
				Two_Product(bdxtail, adytail, tj1, tj0);
				Two_Two_Diff(ti1, ti0, tj1, tj0, abtt[3], abtt[2], abtt[1], abtt[0]);
				//abtt[3] = abtt3;

				abttlen = 4;

			}
			else {

				abt[0] = 0.0;
				abtlen = 1;

				abtt[0] = 0.0;
				abttlen = 1;

			}

			if (cdxtail != 0.0) {

				temp16alen = Scale_Expansion_Zero_Elimination(cxtablen, cxtab, cdxtail, temp16a);
				cxtabtlen = Scale_Expansion_Zero_Elimination(abtlen, abt, cdxtail, cxtabt);
				temp32alen = Scale_Expansion_Zero_Elimination(cxtabtlen, cxtabt, 2.0 * cdx, temp32a);
				temp48len = Fast_Expansion_Sum_Zero_Eliminiation(temp16alen, temp16a, temp32alen, temp32a, temp48);
				finlength = Fast_Expansion_Sum_Zero_Eliminiation(finlength, finnow, temp48len, temp48, finother);

				finswap = finnow;
				finnow = finother;
				finother = finswap;

				if (adytail != 0.0) {

					temp8len = Scale_Expansion_Zero_Elimination(4, bb, cdxtail, temp8);
					temp16alen = Scale_Expansion_Zero_Elimination(temp8len, temp8, adytail, temp16a);
					finlength = Fast_Expansion_Sum_Zero_Eliminiation(finlength, finnow, temp16alen, temp16a, finother);

					finswap = finnow;
					finnow = finother;
					finother = finswap;

				}
				if (bdytail != 0.0) {

					temp8len = Scale_Expansion_Zero_Elimination(4, aa, -cdxtail, temp8);
					temp16alen = Scale_Expansion_Zero_Elimination(temp8len, temp8, bdytail, temp16a);
					finlength = Fast_Expansion_Sum_Zero_Eliminiation(finlength, finnow, temp16alen, temp16a, finother);

					finswap = finnow;
					finnow = finother;
					finother = finswap;

				}

				temp32alen = Scale_Expansion_Zero_Elimination(cxtabtlen, cxtabt, cdxtail, temp32a);
				cxtabttlen = Scale_Expansion_Zero_Elimination(abttlen, abtt, cdxtail, cxtabtt);
				temp16alen = Scale_Expansion_Zero_Elimination(cxtabttlen, cxtabtt, 2.0 * cdx, temp16a);
				temp16blen = Scale_Expansion_Zero_Elimination(cxtabttlen, cxtabtt, cdxtail, temp16b);
				temp32blen = Fast_Expansion_Sum_Zero_Eliminiation(temp16alen, temp16a, temp16blen, temp16b, temp32b);
				temp64len = Fast_Expansion_Sum_Zero_Eliminiation(temp32alen, temp32a, temp32blen, temp32b, temp64);
				finlength = Fast_Expansion_Sum_Zero_Eliminiation(finlength, finnow, temp64len, temp64, finother);

				finswap = finnow;
				finnow = finother;
				finother = finswap;

			}
			if (cdytail != 0.0) {

				temp16alen = Scale_Expansion_Zero_Elimination(cytablen, cytab, cdytail, temp16a);
				cytabtlen = Scale_Expansion_Zero_Elimination(abtlen, abt, cdytail, cytabt);
				temp32alen = Scale_Expansion_Zero_Elimination(cytabtlen, cytabt, 2.0 * cdy, temp32a);
				temp48len = Fast_Expansion_Sum_Zero_Eliminiation(temp16alen, temp16a, temp32alen, temp32a, temp48);
				finlength = Fast_Expansion_Sum_Zero_Eliminiation(finlength, finnow, temp48len, temp48, finother);

				finswap = finnow;
				finnow = finother;
				finother = finswap;

				temp32alen = Scale_Expansion_Zero_Elimination(cytabtlen, cytabt, cdytail, temp32a);
				cytabttlen = Scale_Expansion_Zero_Elimination(abttlen, abtt, cdytail, cytabtt);
				temp16alen = Scale_Expansion_Zero_Elimination(cytabttlen, cytabtt, 2.0 * cdy, temp16a);
				temp16blen = Scale_Expansion_Zero_Elimination(cytabttlen, cytabtt, cdytail, temp16b);
				temp32blen = Fast_Expansion_Sum_Zero_Eliminiation(temp16alen, temp16a, temp16blen, temp16b, temp32b);
				temp64len = Fast_Expansion_Sum_Zero_Eliminiation(temp32alen, temp32a, temp32blen, temp32b, temp64);
				finlength = Fast_Expansion_Sum_Zero_Eliminiation(finlength, finnow, temp64len, temp64, finother);

				finswap = finnow;
				finnow = finother;
				finother = finswap;

			}
		}

		return finnow[finlength - 1];

	};

};