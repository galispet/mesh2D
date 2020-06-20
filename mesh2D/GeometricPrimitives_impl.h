#pragma once

#include <iostream>

namespace GeometricKernel {


	/*****************************************************************************/
	/*                                                                           */
	/*    - Constructor, destructor										         */
	/*                                                                           */
	/*****************************************************************************/
	Vertex::Vertex() {

	};
	Vertex::Vertex(Vertex const & v) {

		x = v.x;
		y = v.y;

		index				= v.index;
		marker				= v.marker;
		adjacentTriangle	= v.adjacentTriangle;

	};
	Vertex::Vertex(v_pointer const & v) {

		x = v->x;
		y = v->y;

		index				= v->index;
		marker				= v->marker;
		adjacentTriangle	= v->adjacentTriangle;

	};
	Vertex::Vertex(double const X, double const Y) : x(X), y(Y) {

	};
	Vertex::~Vertex() {

	};


	/*****************************************************************************/
	/*                                                                           */
	/*	: Set the coordinates of the vertices									 */
	/*                                                                           */
	/*****************************************************************************/
	void Vertex::set_coordinates(double const X, double const Y) {

		x = X;
		y = Y;

	};


	/*****************************************************************************/
	/*                                                                           */
	/*	: Set adjacent triangle													 */
	/*                                                                           */
	/*****************************************************************************/
	void Vertex::set_adjacent_triangle(t_pointer const & t) {

		adjacentTriangle = t;

	};


	/*****************************************************************************/
	/*                                                                           */
	/*	: Check, if the vertex lies in the epsilon-Ball of the vertex 'v'	     */
	/*                                                                           */
	/*****************************************************************************/
	bool Vertex::is_almost_equal(v_pointer const & v) const {

		return sqrt((x - v->x)*(x - v->x) + (y - v->y)*(y - v->y)) < EpsilonInBall;

	};


};


namespace GeometricKernel {


	/*****************************************************************************/
	/*                                                                           */
	/*    - Constructor, destructor										         */
	/*                                                                           */
	/*****************************************************************************/
	Edge::Edge() {

	};
	Edge::Edge(v_pointer const & a, v_pointer const & b) {

		vertices[0] = a;
		vertices[1] = b;

	};
	Edge::~Edge() {


	};


	/*****************************************************************************/
	/*                                                                           */
	/*	: Set the edge's vertices												 */
	/*                                                                           */
	/*****************************************************************************/
	void Edge::set_vertices(v_pointer const & a, v_pointer const & b) {

		vertices[0] = a;
		vertices[1] = b;

	};


	/*****************************************************************************/
	/*                                                                           */
	/*	: Return the pointer to the triangle on the side 'i'					 */
	/*                                                                           */
	/*****************************************************************************/
	void Edge::set_neighbor(t_pointer const & t, unsigned const i) {

		neighbors[i] = t;

	};
	t_pointer & Edge::set_neighbor(unsigned const i) {

		return neighbors[i];

	};

	/*****************************************************************************/
	/*                                                                           */
	/*	: Return the pointer to the triangle on the side 'i'					 */
	/*                                                                           */
	/*****************************************************************************/
	t_pointer Edge::get_neighbor(unsigned const i) const {

		return neighbors[i];

	};


	/*****************************************************************************/
	/*                                                                           */
	/*	: Return the reference of the pointer to the triangle on the side 'i'	 */
	/*                                                                           */
	/*****************************************************************************/
	t_pointer & Edge::get_neighbor(unsigned const i) {

		return neighbors[i];

	};


	/*****************************************************************************/
	/*                                                                           */
	/*	: Check if the edge contains vertex 'v'									 */
	/*                                                                           */
	/*****************************************************************************/
	bool Edge::contains(v_pointer const & v) {

		return vertices[0] == v || vertices[1] == v;

	};


	/*****************************************************************************/
	/*                                                                           */
	/*	: Check if the edge contains vertices 'A' and 'B'						 */
	/*                                                                           */
	/*****************************************************************************/
	bool Edge::contains(v_pointer const & a, v_pointer const & b) {

		return contains(a) && contains(b);

	};


	/*****************************************************************************/
	/*                                                                           */
	/*	: Return the length of the edge											 */
	/*                                                                           */
	/*****************************************************************************/
	double Edge::length() {

		double const x0 = vertices[0]->x;
		double const y0 = vertices[0]->y;

		double const x1 = vertices[1]->x;
		double const y1 = vertices[1]->y;

		return sqrt((x0 - x1)*(x0 - x1) + (y0 - y1)*(y1 - y0));

	};


	/*****************************************************************************/
	/*                                                                           */
	/*	: Return pointer to the vertex which is common for this and edge 'e'	 */
	/*                                                                           */
	/*****************************************************************************/
	v_pointer Edge::get_shared_vertex(e_pointer const & e) {

		if (vertices[0] == e->vertices[0] || vertices[0] == e->vertices[1]) return vertices[0];
		if (vertices[1] == e->vertices[1] || vertices[1] == e->vertices[0]) return vertices[1];

		return NULL;

	};


};


namespace GeometricKernel {


	/*****************************************************************************/
	/*                                                                           */
	/*    - Constructor, destructor										         */
	/*                                                                           */
	/*****************************************************************************/
	Triangle::Triangle() {

	};
	Triangle::Triangle(v_pointer const & a, v_pointer const & b, v_pointer const & c) {


		/*****************************************************************************/
		/*                                                                           */
		/*  It is assumed, that a,b,c are ordered CCW						         */
		/*                                                                           */
		/*****************************************************************************/
		vertices[0] = a;
		vertices[1] = b;
		vertices[2] = c;


		/*****************************************************************************/
		/*                                                                           */
		/*    - Make sure the triangle is oriented CCW (CounterClockWise)	         */
		/*                                                                           */
		/*****************************************************************************/
		//if (orientation_robust(a, b, c) < 0.0) {
		//
		//	vertices[0] = a;
		//	vertices[1] = c;
		//	vertices[2] = b;
		//
		//}
		//else {
		//
		//	vertices[0] = a;
		//	vertices[1] = b;
		//	vertices[2] = c;
		//
		//}


		/*****************************************************************************/
		/*                                                                           */
		/*    - For debugger purpose. Check if the triangle is degenerated (det = 0) */
		/*                                                                           */
		/*****************************************************************************/
		//if (orientation_robust(vertices[0], vertices[1], vertices[2]) <= 0.0) {
		//
		//	std::cout << "Big Phucking Problem" << std::endl;
		//
		//	//double o1 = orientation_robust(vertices[0], vertices[1], vertices[2]);
		//	//double o2 = orientation_fast(vertices[0], vertices[1], vertices[2]);
		//	
		//}


	};
	Triangle::~Triangle() {

	};


	/*****************************************************************************/
	/*                                                                           */
	/*	: Check if the triangle contains vertex 'v' or the couple 'a' and 'b'	 */
	/*	: Check if the triangle contains edge 'e'								 */
	/*                                                                           */
	/*****************************************************************************/
	bool Triangle::contains(v_pointer const & v) const {

		return vertices[0] == v || vertices[1] == v || vertices[2] == v;

	};
	bool Triangle::contains(v_pointer const & a, v_pointer const & b) const {

		return contains(a) && contains(b);

	};
	bool Triangle::contains(e_pointer const & e) const {

		return contains(e->vertices[0], e->vertices[1]);

	};


	/*****************************************************************************/
	/*                                                                           */
	/*	: Return the pointer to the 'i-th' vertex								 */
	/*	: Return the pointer to the vertex opposite to the edge 'e'				 */
	/*	: Return the pointer to the vertex located CW/CCW to the vertex 'v'		 */
	/*	: Return the pointer to the remaining vertex out of 'a' and 'b'			 */
	/*                                                                           */
	/*****************************************************************************/
	v_pointer Triangle::get_vertex(unsigned const i) const {

		return vertices[i];

	};
	v_pointer Triangle::get_vertex(e_pointer const & e) const {

		if		(e == edges[0]) return vertices[0];
		else if (e == edges[1]) return vertices[1];
		else if (e == edges[2]) return vertices[2];

		return NULL;

	};
	v_pointer Triangle::get_vertex_cw(v_pointer const & v)  const {

		if		(vertices[0] == v) return vertices[2];
		else if (vertices[1] == v) return vertices[0];
		else if (vertices[2] == v) return vertices[1];

		return NULL;

	};
	v_pointer Triangle::get_vertex_ccw(v_pointer const & v) const {

		if		(vertices[0] == v) return vertices[1];
		else if (vertices[1] == v) return vertices[2];
		else if (vertices[2] == v) return vertices[0];

		return NULL;

	};
	v_pointer Triangle::get_vertex_except(v_pointer const & a, v_pointer const & b) const {

		if ((vertices[0] == a && vertices[1] == b) || (vertices[1] == a && vertices[0] == b)) return vertices[2];
		else if ((vertices[1] == a && vertices[2] == b) || (vertices[2] == a && vertices[1] == b)) return vertices[0];
		else if ((vertices[2] == a && vertices[0] == b) || (vertices[0] == a && vertices[2] == b)) return vertices[1];

		return NULL;

	};


	/*****************************************************************************/
	/*                                                                           */
	/*	: Set edge into 'edges' to the position where respective pair of		 */
	/*    vertices coincide                                                      */
	/*                                                                           */
	/*****************************************************************************/
	void Triangle::set_edge(e_pointer const & e) {

		if (e->contains(vertices[1], vertices[2])) edges[0] = e;
		else if (e->contains(vertices[0], vertices[2])) edges[1] = e;
		else if (e->contains(vertices[0], vertices[1])) edges[2] = e;

		return;

	};
	e_pointer & Triangle::set_edge(int const i) {

		return edges[i];

	};


	/*****************************************************************************/
	/*                                                                           */
	/*	: Return the pointer to the 'i-th' edge									 */
	/*	: Return the reference to the pointer to the 'i-th' edge				 */
	/*	: Return the pointer to the edge which contains vertices 'a' and 'b'	 */
	/*	: Return the pointer to the edge which located CW/CCW to the vertex 'v'	 */
	/*                                                                           */
	/*****************************************************************************/
	e_pointer	Triangle::get_edge(int const i) const {

		return edges[i];

	};
	e_pointer & Triangle::get_edge(int const i) {

		return edges[i];

	};
	e_pointer	Triangle::get_edge(v_pointer const & v) const {

		if (v == vertices[0]) return edges[0];
		else if (v == vertices[1]) return edges[1];
		else if (v == vertices[2]) return edges[2];

		return NULL;

	};
	e_pointer & Triangle::get_edge(e_pointer const & e) {

		if		(e == edges[0])	return edges[0];
		else if (e == edges[1]) return edges[1];
		else if (e == edges[2]) return edges[2];

		return edges[0];

	};
	e_pointer	Triangle::get_edge(v_pointer const & a, v_pointer const & b) const {

		if ((a == vertices[0] && b == vertices[1]) || (a == vertices[1] && b == vertices[0])) return this->edges[2];
		else if ((a == vertices[1] && b == vertices[2]) || (a == vertices[2] && b == vertices[1])) return this->edges[0];
		else if ((a == vertices[2] && b == vertices[0]) || (a == vertices[0] && b == vertices[2])) return this->edges[1];

		return NULL;

	};
	e_pointer	Triangle::get_edge_cw(v_pointer const & v) const {

		if (vertices[0] == v) return edges[1];
		else if (vertices[1] == v) return edges[2];
		else if (vertices[2] == v) return edges[0];

		return NULL;

	};
	e_pointer	Triangle::get_edge_ccw(v_pointer const & v) const {

		if (vertices[0] == v) return edges[2];
		else if (vertices[1] == v) return edges[0];
		else if (vertices[2] == v) return edges[1];

		return NULL;

	};


	/*****************************************************************************/
	/*                                                                           */
	/*	: Set neighboring triangle 't' into 'neighbors' to the 'i-th' position	 */
	/*	: Set neighboring triangle 't' to the position where there is			 */
	/*	  edge consisting of the vertices 'a' and 'b'							 */
	/*	: Set neighboring triangle 't' to the position where the edge 'e' is	 */
	/*	: Set neighboring triangle 't' to the position where the vertices		 */
	/*	  coincide																 */
	/*                                                                           */
	/*****************************************************************************/
	void Triangle::set_neighbor(t_pointer const & t, unsigned i) {

		neighbors[i] = t;

	};
	void Triangle::set_neighbor(t_pointer const & t, v_pointer const & a, v_pointer const & b) {

		if ((a == vertices[2] && b == vertices[1]) || (a == vertices[1] && b == vertices[2]))	neighbors[0] = t;
		else if ((a == vertices[0] && b == vertices[2]) || (a == vertices[2] && b == vertices[0]))	neighbors[1] = t;
		else if ((a == vertices[0] && b == vertices[1]) || (a == vertices[1] && b == vertices[0]))	neighbors[2] = t;

	};
	void Triangle::set_neighbor(t_pointer const & t, e_pointer const & e) {

		if (e == edges[0]) neighbors[0] = t;
		else if (e == edges[1]) neighbors[1] = t;
		else if (e == edges[2]) neighbors[2] = t;

	};
	void Triangle::set_neighbor(t_pointer const & t) {

		if (!t) return;

		if (t->contains(vertices[1], vertices[2])) {

			neighbors[0] = t;
			t->set_neighbor(this, vertices[1], vertices[2]);

		}
		else if (t->contains(vertices[0], vertices[2])) {

			neighbors[1] = t;
			t->set_neighbor(this, vertices[0], vertices[2]);

		}
		else if (t->contains(vertices[0], vertices[1])) {

			neighbors[2] = t;
			t->set_neighbor(this, vertices[0], vertices[1]);

		}

		return;

	};
	t_pointer & Triangle::set_neighbor(unsigned const i) {

		return neighbors[i];

	};
	t_pointer & Triangle::set_neighbor(e_pointer const & e) {

		if		(e == edges[0]) return neighbors[0];
		else if (e == edges[1]) return neighbors[1];
		else if (e == edges[2]) return neighbors[2];

		return neighbors[0];

	};


	/*****************************************************************************/
	/*                                                                           */
	/*	: Return the pointer to the 'i-th' neighboring triangle					 */
	/*	: Return the reference to the 'i-th' neighboring triangle				 */
	/*	: Return the pointer to the neighboring triangle where the respective	 */
	/*	  vertices coincie with the vertices 'a' and 'b'						 */
	/*	: Return the pointer to the neighboring triangle which is located		 */
	/*	  CW/CCW to the vertex 'v'												 */
	/*                                                                           */
	/*****************************************************************************/
	t_pointer	Triangle::get_neighbor(unsigned const i) const {

		return neighbors[i];

	};
	t_pointer & Triangle::get_neighbor(unsigned const i) {

		return neighbors[i];

	};
	t_pointer	Triangle::get_neighbor(v_pointer const & v) const {

		return get_neighbor(get_vertex_index(v));

	};
	t_pointer & Triangle::get_neighbor(v_pointer const & v) {

		return get_neighbor(get_vertex_index(v));

	};
	t_pointer	Triangle::get_neighbor(v_pointer const & a, v_pointer const & b) const {

		if ((a == vertices[2] && b == vertices[1]) || (a == vertices[1] && b == vertices[2]))	return neighbors[0];
		else if ((a == vertices[0] && b == vertices[2]) || (a == vertices[2] && b == vertices[0]))	return neighbors[1];
		else if ((a == vertices[0] && b == vertices[1]) || (a == vertices[1] && b == vertices[0]))	return neighbors[2];

		return NULL;

	};
	t_pointer	Triangle::get_neighbor_cw(v_pointer const & v) const {

		if (vertices[0] == v) return neighbors[1];
		else if (vertices[1] == v) return neighbors[2];
		else if (vertices[2] == v) return neighbors[0];

		return NULL;

	};
	t_pointer	Triangle::get_neighbor_ccw(v_pointer const & v) const {

		if (vertices[0] == v) return neighbors[2];
		else if (vertices[1] == v) return neighbors[0];
		else if (vertices[2] == v) return neighbors[1];

		return NULL;

	};


	/*****************************************************************************/
	/*                                                                           */
	/*	: Return the local index of the vertex 'v'								 */
	/*	: Return the local index of the edge 'e'								 */
	/*	: Return the local index of the edge consisting of the vertices 'a','b'  */
	/*	: Return the local index of the neighbor 't'							 */
	/*                                                                           */
	/*****************************************************************************/
	int Triangle::get_vertex_index(v_pointer const & v) const {

		if (v == vertices[0]) return 0;
		else if (v == vertices[1]) return 1;
		else if (v == vertices[2]) return 2;

		return -1;

	};
	int Triangle::get_edge_index(e_pointer const & e) const {

		if (e == edges[0]) return 0;
		else if (e == edges[1]) return 1;
		else if (e == edges[2]) return 2;

		return -1;

	};
	int Triangle::get_edge_index(v_pointer const & a, v_pointer const & b) const {

		if ((a == vertices[0] && b == vertices[1]) || (a == vertices[1] && b == vertices[0])) return 2;
		else if ((a == vertices[1] && b == vertices[2]) || (a == vertices[2] && b == vertices[1])) return 0;
		else if ((a == vertices[2] && b == vertices[0]) || (a == vertices[0] && b == vertices[2])) return 1;

		return -1;

	};
	int Triangle::get_neighbor_index(t_pointer const & t) const {

		if		(neighbors[0] == t) return 0;
		else if (neighbors[1] == t) return 1;
		else if (neighbors[2] == t) return 2;

		return -1;

	};


	/*****************************************************************************/
	/*                                                                           */
	/*	: Set the elements in their respective container to NULL				 */
	/*                                                                           */
	/*****************************************************************************/
	void Triangle::null_vertices() {

		vertices[0] = NULL;
		vertices[1] = NULL;
		vertices[2] = NULL;

	};
	void Triangle::null_edges() {

		edges[0] = NULL;
		edges[1] = NULL;
		edges[2] = NULL;

	};
	void Triangle::null_neighbors() {

		neighbors[0] = NULL;
		neighbors[1] = NULL;
		neighbors[2] = NULL;

	};


	/*****************************************************************************/
	/*                                                                           */
	/*	: Return the pointer to the opposite vertex of 'v' in the				 */
	/*	  neighboring triangle 't'												 */
	/*	: Return the pointer to the neighboring triangle which is opposite		 */
	/*	  to the vertex 'v'														 */
	/*                                                                           */
	/*****************************************************************************/
	v_pointer Triangle::get_opposite_vertex(t_pointer const & t, v_pointer const & v) const {

		return this->get_vertex_cw(t->get_vertex_cw(v));

	};


	/*****************************************************************************/
	/*                                                                           */
	/*	: rotate pair of neighboring triangles which have opposite  			 */
	/*	  vertices 'v' and 'ov' (opposite vertex)								 */
	/*                                                                           */
	/*****************************************************************************/
	void Triangle::rotate_triangle_cw(v_pointer const & v, v_pointer const & ov) {

		if (v == vertices[0]) {

			vertices[1] = v;
			vertices[0] = vertices[2];
			vertices[2] = ov;

		}
		else if (v == vertices[1]) {

			vertices[2] = v;
			vertices[1] = vertices[0];
			vertices[0] = ov;

		}
		else if (v == vertices[2]) {

			vertices[0] = v;
			vertices[2] = vertices[1];
			vertices[1] = ov;
		}

	};


	/*****************************************************************************/
	/*                                                                           */
	/*	: Return the area of the triangle										 */
	/*	: Return the squared radius of the circumcircle (RCC) of the triangle	 */
	/*	: Return the squared lenght of the shortest edge (SE) in the triangle	 */
	/*	: Return the ratio (RCC) / (SE)											 */
	/*                                                                           */
	/*****************************************************************************/
	double Triangle::area() const {

		v_pointer const a = vertices[0];
		v_pointer const b = vertices[1];
		v_pointer const c = vertices[2];

		return 0.5 * abs((a->x - c->x) * (b->y - c->y) - (a->y - c->y) * (b->x - c->x));

	};
	double Triangle::circum_radius_squared() const {

		double const ax = vertices[0]->x;
		double const ay = vertices[0]->y;

		double const bx = vertices[1]->x;
		double const by = vertices[1]->y;

		double const cx = vertices[2]->x;
		double const cy = vertices[2]->y;

		double const A = bx - ax;
		double const B = by - ay;
		double const C = cx - ax;
		double const D = cy - ay;

		double const H = cx - bx;
		double const I = cy - by;


		/*****************************************************************************/
		/*                                                                           */
		/*	- Squared lengths of the edges |ab| , |bc| , |ca|						 */
		/*                                                                           */
		/*****************************************************************************/
		double const AB = A * A + B * B;
		double const BC = H * H + I * I;
		double const CA = C * C + D * D;

		double const area_squared = 0.25 * (C*I - D * H)*(C*I - D * H);


		/*****************************************************************************/
		/*                                                                           */
		/*	- Formula for circum radius : sqrt( |ab| * |bc| * |ca| ) / ( 4 * Area )	 */
		/*                                                                           */
		/*****************************************************************************/
		return AB * BC * CA / (16.0 * area_squared);

	};
	double Triangle::shortest_edge_squared() const {

		double const ax = vertices[0]->x;
		double const ay = vertices[0]->y;

		double const bx = vertices[1]->x;
		double const by = vertices[1]->y;

		double const cx = vertices[2]->x;
		double const cy = vertices[2]->y;


		double const A = bx - ax;
		double const B = by - ay;
		double const C = cx - ax;
		double const D = cy - ay;

		double const H = cx - bx;
		double const I = cy - by;


		/*****************************************************************************/
		/*                                                                           */
		/*	- Squared lengths of the edges |ab| , |bc| , |ca|						 */
		/*                                                                           */
		/*****************************************************************************/
		double const AB = A * A + B * B;
		double const BC = H * H + I * I;
		double const CA = C * C + D * D;


		/*****************************************************************************/
		/*                                                                           */
		/*	- Squared lenght of the shortest edge									 */
		/*                                                                           */
		/*****************************************************************************/
		return std::fmin(AB, std::fmin(BC, CA));

	};
	double Triangle::ratio_squared() const {


		double const ax = vertices[0]->x;
		double const ay = vertices[0]->y;

		double const bx = vertices[1]->x;
		double const by = vertices[1]->y;

		double const cx = vertices[2]->x;
		double const cy = vertices[2]->y;

		double const A = bx - ax;
		double const B = by - ay;
		double const C = cx - ax;
		double const D = cy - ay;

		double const H = cx - bx;
		double const I = cy - by;

		/*****************************************************************************/
		/*                                                                           */
		/*	- Squared lengths of the edges |ab| , |bc| , |ca|						 */
		/*                                                                           */
		/*****************************************************************************/
		double const AB = A * A + B * B;
		double const BC = H * H + I * I;
		double const CA = C * C + D * D;

		double const area_squared = 0.25 * (C*I - D * H)*(C*I - D * H);


		double const circum_radius_squared = AB * BC * CA / (16.0 * area_squared);
		double const shortest_edge_squared = std::fmin(AB, std::fmin(BC, CA));


		return circum_radius_squared / shortest_edge_squared;

	};


	/*****************************************************************************/
	/*                                                                           */
	/*	: Return the coordinates of the center of the circumcircle				 */
	/*                                                                           */
	/*****************************************************************************/
	void Triangle::circum_center(double & xc, double & yc) const {

		double const ax = vertices[0]->x;
		double const ay = vertices[0]->y;

		double const bx = vertices[1]->x;
		double const by = vertices[1]->y;

		double const cx = vertices[2]->x;
		double const cy = vertices[2]->y;

		double const A = bx - ax;
		double const B = by - ay;
		double const C = cx - ax;
		double const D = cy - ay;

		double const H = cx - bx;
		double const I = cy - by;

		double const E = A * (ax + bx) + B * (ay + by);
		double const F = C * (ax + cx) + D * (ay + cy);

		double const G = 2.0 * (A*I - B * H);


		/*****************************************************************************/
		/*                                                                           */
		/*	- Circumcenter coordinates												 */
		/*                                                                           */
		/*****************************************************************************/
		xc = (D * E - B * F) / G;
		yc = (A * F - C * E) / G;

	};


	/*****************************************************************************/
	/*                                                                           */
	/*	: Check if the triangle meets quality criteria (min angle, max area)	 */
	/*                                                                           */
	/*****************************************************************************/
	bool Triangle::is_poor_quality(double const angleBoundInDegrees, double const areaBound) {


		double const angleBoundInRadians = Pi * angleBoundInDegrees / 180.0;
		double const angleBound = 0.5 / (sin(angleBoundInRadians));

		bool const poorAngle = ratio_squared() > (angleBound*angleBound);
		bool const poorArea = circum_radius_squared() > (areaBound*areaBound);

		if (poorAngle || poorArea)
			return true;

		return false;

	}


};
