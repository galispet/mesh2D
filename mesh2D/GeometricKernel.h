#pragma once


#include <cmath>


namespace GeometricKernel {

	/*
	double const EPSILON_EQUALITY = 1.e-5;
	double const EPSILON_IN_CIRCUMCIRCLE = 1.e-5;
	double const EPSILON_ON_EDGE = 1.e-3;		// Should be bigger than ORIENTATION -> Otherwise 'COLLINEAR' triangle could be formed => In circum-circle test can be non-symmetric => program failure
	double const EPSILON_SEGMENT_ENCROACH = 1.e-4;
	*/

	#define EpsilonInBall			1e-8
	#define EpsilonOrientation		1e-8
	#define EpsilonOnEdge			0.0


	#define Pi				3.1415926535897932384626433832795028841971693993751058209
	#define SquareRootTwo	1.4142135623730950488016887242096980785696718753769480732

	enum class MarkerVertex						{ Free, Constrained };
	enum class MarkerEdge						{ None, Dirichlet, Neumann };
	enum class MarkerTriangle					{ None, Outside, Inside };
	enum class GeometricPredicatesArithmetic	{ Fast, Robust };

	enum class Location { OnVertex, OnEdge, InTriangle, NotVisible };


	class Vertex;
	class Edge;
	class Triangle;

	typedef Vertex *	v_pointer;
	typedef Edge *		e_pointer;
	typedef Triangle *	t_pointer;


	class Vertex {


		/*****************************************************************************/
		/*                                                                           */
		/*		- Friend classes and functions								         */
		/*                                                                           */
		/*****************************************************************************/

		/*
		friend class Triangle;

		friend class Edge;

		friend class Mesh;

		friend class PlanarStraightLineGraph;

		template<GEOMETRIC_KERNEL GK>
		friend class Triangulation;

		template <GEOMETRIC_KERNEL GK>
		friend class Predicates;
		*/

		/*
		friend bool intersecting_edges(Edge * e, Edge * f);
		friend bool intersecting_edges(Vertex * p1, Vertex * p2, Vertex * q1, Vertex * q2);
		friend bool get_edge_intersection(Vertex * p1, Vertex * p2, Vertex * p3, Vertex * p4, double & x, double & y);
		friend bool in_circle_fast(Triangle * t, Vertex * v);
		friend bool in_circle_robust(Triangle * t, Vertex * v);
		friend bool on_edge_fast(Vertex * a, Vertex * b, Vertex * v);
		friend bool on_edge_fast(Edge * e, Vertex * v);
		friend bool on_edge_robust(Vertex * a, Vertex * b, Vertex * v);
		friend bool on_edge_robust(Edge * e, Vertex * v);
		friend double orientation_fast(Vertex * a, Vertex * b, Vertex * v);
		friend double orientation_robust(Vertex * a, Vertex * b, Vertex * v);
		friend double orientationAdapt(Vertex * a, Vertex * b, Vertex * v, const double detsum);
		friend bool in_circleAdapt(Triangle * t, Vertex * v, const double permanent);
		friend void get_mid_point(Edge * e, double &x_mid, double &y_mid);
		friend bool is_encroached(Vertex * const a, Vertex * const b, Vertex * const v);


		friend bool t_compare_x(Triangle * const t1, Triangle * const t2);
		friend bool t_compare_y(Triangle * const t1, Triangle * const t2);
		friend bool v_compare_x(Vertex * const v1, Vertex * const v2);
		friend bool v_compare_y(Vertex * const v1, Vertex * const v2);
		friend bool e_compare_x(Edge * const e1, Edge * const e2);
		friend bool e_compare_y(Edge * const e1, Edge * const e2);
		*/


	public:



		/*****************************************************************************/
		/*                                                                           */
		/*    - Constructor, copy construcotr, destructor					         */
		/*                                                                           */
		/*****************************************************************************/
		Vertex();
		Vertex(double const X, double const Y);
		Vertex(Vertex const & v);
		Vertex(v_pointer const & v);
		~Vertex();


		/*****************************************************************************/
		/*                                                                           */
		/*		- Class members												         */
		/*                                                                           */
		/*	: x, y																	 */
		/*                                                                           */
		/*		Coordinates of the vertex											 */
		/*                                                                           */
		/*	: index																	 */
		/*                                                                           */
		/*		Global index of the vertex											 */
		/*                                                                           */
		/*	: adjacentTriangle														 */
		/*                                                                           */
		/*		Pointer to unspecified adjacent triangle							 */
		/*                                                                           */
		/*	: marker																 */
		/*                                                                           */
		/*      denotes if the vertex is 'Free' (can be moved or deleted)			 */
		/*		or 'Constrained' (can't be moved or deleted)				         */
		/*                                                                           */
		/*****************************************************************************/
		double x = 0.0;
		double y = 0.0;

		int index = -1;

		t_pointer adjacentTriangle	= NULL;
		MarkerVertex marker			= MarkerVertex::Free;



		/*****************************************************************************/
		/*                                                                           */
		/*    - Class methods												         */
		/*                                                                           */
		/*****************************************************************************/
		void set_coordinates(double const X, double const Y);
		void set_adjacent_triangle(t_pointer const & t);

		bool is_almost_equal(v_pointer const & v) const;


	};

	class Edge {


		/*****************************************************************************/
		/*                                                                           */
		/*		- Friend classes and functions								         */
		/*                                                                           */
		/*****************************************************************************/

		/*
		friend class Vertex;

		friend class Triangle;

		friend class Mesh;

		template<GEOMETRIC_KERNEL GK>
		friend class Triangulation;

		template <GEOMETRIC_KERNEL GK>
		friend class Predicates;
		*/

		/*
		friend bool intersecting_edges(Edge * e, Edge * f);
		friend bool on_edge_fast(Vertex * a, Vertex * b, Vertex * v);
		friend bool on_edge_fast(Edge * e, Vertex * v);
		friend bool on_edge_robust(Edge * e, Vertex * v);
		friend bool is_encroached(Vertex * const a, Vertex * const b, Vertex * const v);
		friend void get_mid_point(Edge * e, double &x_mid, double &y_mid);


		friend bool e_compare_x(Edge * const e1, Edge * const e2);
		friend bool e_compare_y(Edge * const e1, Edge * const e2);
		friend bool marker_compare_neumann(Edge * const e1, Edge * const e2);
		friend bool marker_compare_dirichlet(Edge * const e1, Edge * const e2);
		*/


	public:


		/*****************************************************************************/
		/*                                                                           */
		/*    - Constructor, destructor										         */
		/*                                                                           */
		/*****************************************************************************/
		Edge();
		Edge(v_pointer const & a, v_pointer const & b);
		~Edge();


		/*****************************************************************************/
		/*                                                                           */
		/*		- Class members												         */
		/*                                                                           */
		/*	: vertices																 */
		/*                                                                           */
		/*		Pointers to the vertices of the edge								 */
		/*                                                                           */
		/*	: neighbors																 */
		/*                                                                           */
		/*      Pointers to adjacent triangles on each side of the edge				 */
		/*                                                                           */
		/*	: index																	 */
		/*                                                                           */
		/*		Global index of the edge											 */
		/*                                                                           */
		/*	: isConstrained															 */
		/*                                                                           */
		/*      If the edge is constrained, then it can't be flipped				 */
		/*                                                                           */
		/*	: marker																 */
		/*                                                                           */
		/*      Denotes type of the boundary condition on the edge					 */
		/*                                                                           */
		/*****************************************************************************/
		v_pointer vertices[2]  = { NULL, NULL };
		t_pointer neighbors[2] = { NULL, NULL };

		int	index = -1;
		bool isConstrained = 0;

		MarkerEdge marker = MarkerEdge::None;


		/*****************************************************************************/
		/*                                                                           */
		/*		- Class methods												         */
		/*                                                                           */
		/*****************************************************************************/
		void set_vertices(v_pointer const & a, v_pointer const & b);
		void set_neighbor(t_pointer const & t, unsigned const i);
		t_pointer & set_neighbor(unsigned const i);

		t_pointer   get_neighbor(unsigned const i) const;
		t_pointer & get_neighbor(unsigned const i);

		bool contains(v_pointer const & v);
		bool contains(v_pointer const & a, v_pointer const & b);

		double length();

		v_pointer get_shared_vertex(e_pointer const & e);

	};

	class Triangle {


		/*****************************************************************************/
		/*                                                                           */
		/*		- Friend classes and functions								         */
		/*                                                                           */
		/*****************************************************************************/
		//template<GeometricPredicatesArithmetic Arithmetic>
		//friend class IncrementalTriangulation;


		/*
		friend class Mesh;

		template <GEOMETRIC_KERNEL GK>
		friend class Triangulation;

		template <GEOMETRIC_KERNEL GK>
		friend class Predicates;
		*/

		/*
		friend bool in_triangle_fast(t_pointer t, v_pointer v);
		friend bool in_triangle_robust(t_pointer t, v_pointer v);
		friend bool in_circle_fast(t_pointer t, v_pointer v);
		friend bool in_circle_robust(t_pointer t, v_pointer v);
		friend bool in_circleAdapt(t_pointer t, v_pointer v, const double permanent);


		friend bool t_compare_x(t_pointer const t1, t_pointer const t2);
		friend bool t_compare_y(t_pointer const t1, t_pointer const t2);
		*/


	public:


		/*****************************************************************************/
		/*                                                                           */
		/*    - Constructor, destructor										         */
		/*                                                                           */
		/*****************************************************************************/
		Triangle();
		Triangle(v_pointer const & a, v_pointer const & b, v_pointer const & c);
		~Triangle();


		/*****************************************************************************/
		/*                                                                           */
		/*		- Class members												         */
		/*                                                                           */
		/*	: vertices																 */
		/*                                                                           */
		/*		Pointers to vertices of the triangle								 */
		/*                                                                           */
		/*	: edges																	 */
		/*                                                                           */
		/*		Pointers to edges of the triangle									 */
		/*                                                                           */
		/*	: neighbors																 */
		/*                                                                           */
		/*      Adjacent triangles on each side of the edge							 */
		/*                                                                           */
		/*	: index																	 */
		/*                                                                           */
		/*		Global index of the edge											 */
		/*                                                                           */
		/*	: isConstrained															 */
		/*                                                                           */
		/*      If the edge is constrained, then it can't be flipped				 */
		/*                                                                           */
		/*	: marker																 */
		/*                                                                           */
		/*      Mark the triangle as 'Inside' or 'Outside' of the triangulation		 */
		/*      so the 'virus' algorithm will destroy triangles to be deleted        */
		/*                                                                           */
		/*****************************************************************************/
		v_pointer vertices[3]	= { NULL,NULL,NULL };
		e_pointer edges[3]		= { NULL,NULL,NULL };
		t_pointer neighbors[3]	= { NULL,NULL,NULL };

		int index = -1;

		MarkerTriangle marker = MarkerTriangle::None;



		/*****************************************************************************/
		/*                                                                           */
		/*    - Class methods												         */
		/*                                                                           */
		/*****************************************************************************/
		bool contains(v_pointer const & v) const;
		bool contains(v_pointer const & a, v_pointer const & b) const;
		bool contains(e_pointer const & e) const;


		v_pointer get_vertex(unsigned const i) const;
		v_pointer get_vertex(e_pointer const & e) const;
		v_pointer get_vertex_cw(v_pointer const & v)  const;
		v_pointer get_vertex_ccw(v_pointer const & v) const;
		v_pointer get_vertex_except(v_pointer const & a, v_pointer const & b) const;


		void set_edge(e_pointer const & e);
		e_pointer & set_edge(int const i);

		e_pointer	get_edge(int const i) const;
		e_pointer & get_edge(int const i);
		e_pointer	get_edge(v_pointer const & v) const;
		e_pointer & get_edge(e_pointer const & e);
		e_pointer	get_edge(v_pointer const & a, v_pointer const & b) const;
		e_pointer	get_edge_cw(v_pointer const & v) const;
		e_pointer	get_edge_ccw(v_pointer const & v) const;


		void set_neighbor(t_pointer const & t, unsigned i);
		void set_neighbor(t_pointer const & t, v_pointer const & a, v_pointer const & b);
		void set_neighbor(t_pointer const & t, e_pointer const & e);
		void set_neighbor(t_pointer const & t);
		t_pointer & set_neighbor(unsigned const i);
		t_pointer & set_neighbor(e_pointer const & e);


		t_pointer	get_neighbor(unsigned const i) const;
		t_pointer & get_neighbor(unsigned const i);
		t_pointer	get_neighbor(v_pointer const & v) const;
		t_pointer & get_neighbor(v_pointer const & v);
		t_pointer	get_neighbor(v_pointer const & a, v_pointer const & b) const;
		t_pointer	get_neighbor_cw(v_pointer const & v) const;
		t_pointer	get_neighbor_ccw(v_pointer const & v) const;


		int get_vertex_index(v_pointer const & v) const;
		int get_edge_index(e_pointer const & e) const;
		int get_edge_index(v_pointer const & a, v_pointer const & b) const;
		int get_neighbor_index(t_pointer const & t) const;


		void null_vertices();
		void null_edges();
		void null_neighbors();


		v_pointer get_opposite_vertex(t_pointer const & t, v_pointer const & v) const;

		void rotate_triangle_cw(v_pointer const & v, v_pointer const & ov);


		double area() const;
		double circum_radius_squared() const;
		double shortest_edge_squared() const;
		double ratio_squared() const;

		void circum_center(double & xc, double & yc) const;

		bool is_poor_quality(double const angleBoundInDegrees, double const areaBound);


	};



	template<GeometricPredicatesArithmetic Arithmetic>
	class GeometricPredicates {

	public:

		/*****************************************************************************/
		/*                                                                           */
		/*    - Class members												         */
		/*                                                                           */
		/*      Statistics how many times the computations were executed	         */
		/*                                                                           */
		/*****************************************************************************/
		unsigned StatisticsOrientation		= 0;
		unsigned StatisticsInCircle			= 0;

		unsigned StatisticsAdaptOrientation = 0;
		unsigned StatisticsAdaptInCircle	= 0;

		unsigned StatisticsExactOrientation = 0;
		unsigned StatisticsExactInCircle	= 0;


		/*****************************************************************************/
		/*                                                                           */
		/*    - Class methods												         */
		/*                                                                           */
		/*****************************************************************************/
		double orientation(v_pointer const & a, v_pointer const & b, v_pointer const & c);
		double in_circumcircle(t_pointer const & t, v_pointer const & v);

		void get_edge_intersection(v_pointer const & p1, v_pointer const & p2, v_pointer const & q1, v_pointer const & q2, double & x, double & y);


		double orientation_fast(v_pointer const a, v_pointer const b, v_pointer const c) const;
		double orientation_robust(v_pointer const a, v_pointer const b, v_pointer const c) const;
		double orientation_adapt(v_pointer const a, v_pointer const b, v_pointer const c, double const determinantSum) ;

		double in_circumcircle_fast(t_pointer const t, v_pointer const v) const;
		double in_circumcircle_robust(t_pointer const t, v_pointer const v) const;
		double in_circumcircle_adapt(t_pointer const t, v_pointer const p, double const permanent) ;


		/*****************************************************************************/
		/*                                                                           */
		/*		- Class members												         */
		/*                                                                           */
		/*****************************************************************************/


		/*****************************************************************************/
		/*                                                                           */
		/* In circumcircle test														 */
		/*                                                                           */
		/*****************************************************************************/


	};


};


#include "GeometricPrimitives_impl.h"
#include "GeometricPredicates_impl.h"
