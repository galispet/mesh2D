#pragma once


#include "triangulation.h"
#include "enumerators.h"
#include "shapes.h"

#include "gnuplot.h"

#include <vector>
#include <iostream>
#include <fstream>

#include <cmath>


bool t_compare_x(Triangle * const t1, Triangle * const t2) {


	double const t1_min = std::fmin(std::fmin(t1->vertices[0]->x, t1->vertices[1]->x), t1->vertices[2]->x);
	double const t2_min = std::fmin(std::fmin(t2->vertices[0]->x, t2->vertices[1]->x), t2->vertices[2]->x);

	return t1_min < t2_min;

};
bool t_compare_y(Triangle * const t1, Triangle * const t2) {


	double const t1_min = std::fmin(std::fmin(t1->vertices[0]->y, t1->vertices[1]->y), t1->vertices[2]->y);
	double const t2_min = std::fmin(std::fmin(t2->vertices[0]->y, t2->vertices[1]->y), t2->vertices[2]->y);

	return t1_min < t2_min;

};
bool v_compare_x(Vertex * const v1, Vertex * const v2) {


	double const v1_min = v1->x;
	double const v2_min = v2->x;

	return v1_min < v2_min;

};
bool v_compare_y(Vertex * const v1, Vertex * const v2) {


	double const v1_min = v1->y;
	double const v2_min = v2->y;

	return v1_min < v2_min;

};
bool e_compare_x(Edge * const e1, Edge * const e2) {


	double const e1_min = std::fmin(e1->a->x, e1->b->x);
	double const e2_min = std::fmin(e2->a->x, e2->b->x);

	return e1_min < e2_min;

};
bool e_compare_y(Edge * const e1, Edge * const e2) {


	double const e1_min = std::fmin(e1->a->y, e1->b->y);
	double const e2_min = std::fmin(e2->a->y, e2->b->y);

	return e1_min < e2_min;

};


bool marker_compare_neumann(Edge * const e1, Edge * const e2) {


	E_MARKER const e1_m = e1->marker;
	E_MARKER const e2_m = e2->marker;

	if (e1_m == E_MARKER::NEUMANN && e2_m != E_MARKER::NEUMANN)
		return true;
	
	return false;

};
bool marker_compare_dirichlet(Edge * const e1, Edge * const e2) {


	E_MARKER const e1_m = e1->marker;
	E_MARKER const e2_m = e2->marker;

	if (e1_m == E_MARKER::DIRICHLET && e2_m != E_MARKER::DIRICHLET)
		return true;

	return false;

};



class Mesh {


	std::vector<Vertex *> vertices;
	std::vector<Edge *> edges;
	std::vector<Triangle *> triangles;


	unsigned num_vertices = 0;
	unsigned num_edges = 0;
	unsigned num_triangles = 0;

	unsigned num_dirichlet = 0;
	unsigned num_neumann = 0;


	void clear();


	template < GEOMETRIC_KERNEL GK >
	void allocate_new_primitives(Triangulation<GK> & triangulation);
	template < GEOMETRIC_KERNEL GK >
	void set_adjacencies(Triangulation<GK> & triangulation);

public:



	template < GEOMETRIC_KERNEL GK>
	Mesh(Triangulation<GK> & triangulation);
	~Mesh();


	Vertex * get_vertex(unsigned i) const;
	Edge * get_edge(unsigned i) const;
	Triangle * get_triangle(unsigned i) const;


	unsigned get_number_of_vertices() const;
	unsigned get_number_of_edges() const;
	unsigned get_number_of_triangles() const;

	void export_vertices(std::ofstream & stream) const;
	void export_edges(std::ofstream & stream) const;
	void export_triangles(std::ofstream & stream) const;
	void export_triangle(std::ofstream & stream, Triangle * const t) const;


	template<NUMBERING algorithm>
	void apply_numbering();


	void cuthill_mckee();

};






template < GEOMETRIC_KERNEL GK>
Mesh::Mesh(Triangulation<GK> & triangulation) {



	num_vertices = triangulation.get_number_of_vertices();
	num_edges = triangulation.get_number_of_edges();
	num_triangles = triangulation.get_number_of_triangles();

	num_dirichlet = triangulation.get_num_dirichlet_edges();
	num_neumann = triangulation.get_num_neumann_edges();


	// Alocate new vertices, edges, triangles
	allocate_new_primitives(triangulation);

	std::stable_sort(vertices.begin(), vertices.end(), v_compare_x);

	// Set neighbors of triangles, neighbors of edges, adjacent triangle to vertices
	set_adjacencies(triangulation);



	//// Sort with respective coordinates
	std::stable_sort(vertices.begin(), vertices.end(), v_compare_x);
	std::stable_sort(vertices.begin(), vertices.end(), v_compare_y);
	//
	//// Sort with respective coordinates
	std::stable_sort(triangles.begin(), triangles.end(), t_compare_x);
	std::stable_sort(triangles.begin(), triangles.end(), t_compare_y);
	//
	//// Sort with respective Marker. 1. will be Dirichlet, 2. Neumann, 3. NONE. Edges with respective marker are then sorted with respect to x,y coordinates
	std::stable_sort(edges.begin(), edges.end(), marker_compare_dirichlet);
	std::stable_sort(edges.begin(), edges.begin() + num_dirichlet, e_compare_x);
	std::stable_sort(edges.begin(), edges.begin() + num_dirichlet, e_compare_y);
	//
	//std::stable_sort(edges.begin(), edges.end(), marker_compare_neumann);
	std::stable_sort(edges.begin() + num_dirichlet, edges.end(), marker_compare_neumann);
	std::stable_sort(edges.begin() + num_dirichlet, edges.begin() + num_neumann + num_dirichlet, e_compare_x);
	std::stable_sort(edges.begin() + num_dirichlet, edges.begin() + num_neumann + num_dirichlet, e_compare_y);
	//
	std::stable_sort(edges.begin() + num_neumann + num_dirichlet, edges.end(), e_compare_x);
	std::stable_sort(edges.begin() + num_neumann + num_dirichlet, edges.end(), e_compare_y);

	//for (size_t i = 0; i < edges.size(); i++)
	//	if (edges[i]->marker == E_MARKER::NEUMANN)
	//		cout << 1 << endl;
	//	else
	//		cout << 0 << endl;

	std::cout << std::endl;
	cout << "******************** MESH *******************" << endl;

	std::cout << "\nNumber of vertices : ";
	std::cout << vertices.size() << endl;
	std::cout << "\nNumber of edges : ";
	std::cout << edges.size() << endl;
	std::cout << "\nNumber of triangles : ";
	std::cout << triangles.size() << endl;

	std::cout << "Number of constrained edges : ";
	std::cout << std::count_if(edges.begin(), edges.end(), [](auto e) {return e->is_constrained; }) << std::endl;

	std::cout << "Number of constrained vertices : ";
	std::cout << std::count_if(vertices.begin(), vertices.end(), [](auto v) {return v->marker == V_MARKER::CONSTRAINED; }) << std::endl;

	//std::cout << "\nNumber of vertices on the boundary : ";
	//std::cout << std::count_if(is_on_boundary_vertices.begin(), is_on_boundary_vertices.end(), [](auto it) {return (it); }) << std::endl;


	cout << "*********************************************" << endl;


};
Mesh::~Mesh() {

	clear();

};

template < GEOMETRIC_KERNEL GK >
void Mesh::allocate_new_primitives(Triangulation<GK> & triangulation) {


	// Deep Copy of vertices
	for (size_t i = 0; i < num_vertices; i++) {

		Vertex * const v = triangulation.vertices_tri[i];

		double const x = v->x;
		double const y = v->y;

		Vertex * const new_vertex = new Vertex(x, y);

		new_vertex->index = v->index;
		new_vertex->marker = v->marker;

		vertices.push_back(new_vertex);

	}

	int k = 0;
	// Deep Copy of edges
	for (size_t i = 0; i < num_edges; i++) {


		Edge * const e = triangulation.edges_tri[i];

		Vertex * const a = vertices[e->a->index];
		Vertex * const b = vertices[e->b->index];

		Edge * const new_edge = new Edge(a, b);

		new_edge->index = e->index;
		new_edge->is_constrained = e->is_constrained;
		new_edge->marker = e->marker;

		if (e->index == -1)
			k++;

		edges.push_back(new_edge);

	}

	// Deep Copy of triangles
	for (size_t i = 0; i < num_triangles; i++) {

		Triangle * const t = triangulation.triangles_tri[i];

		Vertex * const A = t->vertices[0];
		Vertex * const B = t->vertices[1];
		Vertex * const C = t->vertices[2];

		Vertex * const a = vertices[A->index];
		Vertex * const b = vertices[B->index];
		Vertex * const c = vertices[C->index];

		Triangle * const new_triangle = new Triangle(a, b, c);

		new_triangle->edges[0] = edges[t->edges[0]->index];
		new_triangle->edges[1] = edges[t->edges[1]->index];
		new_triangle->edges[2] = edges[t->edges[2]->index];

		new_triangle->index = t->index;
		new_triangle->marker = t->marker;

		triangles.push_back(new_triangle);

	}



	//std::stable_sort(vertices.begin(), vertices.end(), [](v_pointer v, v_pointer p) {return v->index < p->index; });
	//for (size_t i = 0; i < vertices.size() - 1; i++) {

	//	unsigned const index = vertices[i + 1]->index - vertices[i]->index;
	//	//unsigned const index = vertices[i]->index;

	//	cout << index << endl;

	//	if (index != 1)
	//		cout << "shit ---------------------------------------- \n----------------------------------------" << endl;

	//}


	//std::stable_sort(edges.begin(), edges.end(), [](e_pointer e, e_pointer f) {return e->index < f->index; });
	//for (size_t i = 0; i < edges.size() - 1; i++) {

	//	unsigned const index = edges[i + 1]->index - edges[i]->index;
	//	//unsigned const index = edges[i]->index;

	//	cout << index << endl;

	//	if (index != 1)
	//		cout << "shit ---------------------------------------- \n----------------------------------------" << endl;

	//}

	//std::stable_sort(triangles.begin(), triangles.end(), [](t_pointer t, t_pointer k) {return t->index < k->index; });
	//for (size_t i = 0; i < triangles.size() - 1; i++) {

	//	unsigned const index = triangles[i + 1]->index - triangles[i]->index;
	//	//unsigned const index = triangles[i]->index;

	//	cout << index << endl;

	//	if (index != 1)
	//		cout << "shit ---------------------------------------- \n----------------------------------------" << endl;

	//}



};
template < GEOMETRIC_KERNEL GK >
void Mesh::set_adjacencies(Triangulation<GK> & triangulation) {


	// Set adjacent triangles of vertices
	for (size_t i = 0; i < num_vertices; i++) {

		Triangle * const adjacent_triangle = triangulation.vertices_tri[i]->adjacent_triangle;

		if (!adjacent_triangle) {

			cout << "Invalid vertex. Continue" << endl;
			continue;

		}

		vertices[i]->adjacent_triangle = triangles[adjacent_triangle->index];

	}

	// Set adjacent triangles of edges
	for (size_t i = 0; i < num_edges; i++) {

		Edge * const e = triangulation.edges_tri[i];

		Triangle * const n0 = e->neighbors[0];
		Triangle * const n1 = e->neighbors[1];

		Edge * const f = edges[i];

		if (n0)	f->neighbors[0] = triangles[n0->index];
		else	f->neighbors[0] = NULL;

		if (n1)	f->neighbors[1] = triangles[n1->index];
		else	f->neighbors[1] = NULL;

	}

	// Set adjacent triangles of triangles
	for (size_t i = 0; i < num_triangles; i++) {

		Triangle * const t = triangulation.triangles_tri[i];

		Triangle * const n0 = t->neighbors[0];
		Triangle * const n1 = t->neighbors[1];
		Triangle * const n2 = t->neighbors[2];

		Triangle * const k = triangles[i];

		if (n0)	k->neighbors[0] = triangles[n0->index];
		else	k->neighbors[0] = NULL;

		if (n1)	k->neighbors[1] = triangles[n1->index];
		else	k->neighbors[1] = NULL;

		if (n2)	k->neighbors[2] = triangles[n2->index];
		else	k->neighbors[2] = NULL;

	}

};


Vertex * Mesh::get_vertex(unsigned i) const {

	return vertices[i];

}
Edge * Mesh::get_edge(unsigned i) const {

	return edges[i];

}
Triangle * Mesh::get_triangle(unsigned i) const {

	return triangles[i];

}



unsigned Mesh::get_number_of_vertices() const {
	
	return num_vertices; 

};
unsigned Mesh::get_number_of_edges() const { 
	
	return num_edges;

};
unsigned Mesh::get_number_of_triangles() const {
	
	return num_triangles;

};


void Mesh::export_vertices(std::ofstream & stream) const {

	size_t const n = vertices.size();

	for (size_t i = 0; i < n; i++) {

		double const v_x = vertices[i]->x;
		double const v_y = vertices[i]->y;

		stream << v_x << "  " << v_y << std::endl;

	}

};
void Mesh::export_edges(std::ofstream & stream) const {

	size_t const n = edges.size();


	//// Gnuplot
	for (size_t i = 0; i < n; i++) {

		//if (edges[i]->is_constrained) {

			double const v0_x = edges[i]->a->x;
			double const v0_y = edges[i]->a->y;

			double const v1_x = edges[i]->b->x;
			double const v1_y = edges[i]->b->y;

			stream << v0_x << "  " << v0_y << std::endl;
			stream << v1_x << "  " << v1_y << std::endl << std::endl;

		//}

	}

	//// Matlab
	//for (size_t i = 0; i < n; i++) {

	//	//if (edges[i]->is_constrained) {

	//		double const v0_x = edges[i]->a->x;
	//		double const v0_y = edges[i]->a->y;

	//		double const v1_x = edges[i]->b->x;
	//		double const v1_y = edges[i]->b->y;

	//		stream << v0_x << "  " << v0_y << "  " << v1_x << "  " << v1_y << std::endl;

	//	//}

	//}

};
void Mesh::export_triangles(std::ofstream & stream) const {

	size_t const n = triangles.size();

	for (size_t i = 0; i < n; i++) {

		Triangle * t = triangles[i];

		double const v0_x = t->vertices[0]->x;
		double const v0_y = t->vertices[0]->y;

		double const v1_x = t->vertices[1]->x;
		double const v1_y = t->vertices[1]->y;

		double const v2_x = t->vertices[2]->x;
		double const v2_y = t->vertices[2]->y;

		stream << v0_x << "  " << v0_y << std::endl;
		stream << v1_x << "  " << v1_y << std::endl;
		stream << v2_x << "  " << v2_y << std::endl;
		stream << v0_x << "  " << v0_y << std::endl << std::endl;

	}

};


void Mesh::clear() {


	Vertex	 * v = NULL;
	Edge	 * e = NULL;
	Triangle * t = NULL;


	// Free allocated memory
	for (size_t i = 0; i < vertices.size(); i++) {

		v = vertices[i];
		delete v;

	}
	for (size_t i = 0; i < edges.size(); i++) {

		e = edges[i];
		delete e;

	}
	for (size_t i = 0; i < triangles.size(); i++) {

		t = triangles[i];
		delete t;

	}

	// Clear the vectors of the rubbish
	vertices.clear();
	edges.clear();
	triangles.clear();

};




template<NUMBERING algorithm>
void Mesh::apply_numbering() {


	switch (algorithm) {


	case NUMBERING::CM:

		cuthill_mckee();
		break;

	}


};

void Mesh::cuthill_mckee() {






};

