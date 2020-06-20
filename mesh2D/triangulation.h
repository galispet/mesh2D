#pragma once



#include "GeometricKernel.h"


#include <random>
#include <iostream>
#include <vector>
#include <list>
#include <fstream>


using namespace GeometricKernel;


namespace DelaunayTriangulation {


	class PlanarStraightLineGraph {


		template<GeometricPredicatesArithmetic Arithmetic>
		friend class IncrementalTriangulation;


	private:

		std::vector<v_pointer>	vertices;
		std::vector<unsigned*>	segments;		// 0. index of first vertex  1. index of second vertex  2. edge_marker
		std::vector<MarkerEdge>	segmentsMarker;


		double xMinimum = +DBL_MAX;
		double xMaximum = -DBL_MAX;
		double yMinimum = +DBL_MAX;
		double yMaximum = -DBL_MAX;


	public:

		PlanarStraightLineGraph();
		~PlanarStraightLineGraph();

		template<MarkerVertex vMarker>
		v_pointer insert_vertex(Vertex const & v);		
		template<MarkerEdge eMarker>
		void insert_constraint(v_pointer const & a, v_pointer const & b);
		template<MarkerEdge eMarker>
		void insert_constraint(unsigned const i, unsigned const j);


		v_pointer get_vertex_pointer(unsigned const i);


		int get_number_of_vertices() const;
		int get_number_of_segments() const;

		void clear();

	};


	template <GeometricPredicatesArithmetic Arithmetic>
	class IncrementalTriangulation {


	public:


		/*****************************************************************************/
		/*                                                                           */
		/*	Constructor and destructor		The input is pslg and seed vertices      */
		/*									defining holes							 */
		/*																			 */
		/*****************************************************************************/
		IncrementalTriangulation(PlanarStraightLineGraph const & pslg, std::vector<Vertex> const & seeds);
		~IncrementalTriangulation();


		/*****************************************************************************/
		/*                                                                           */
		/*  Returns the number of the triangulation primitives						 */
		/*																			 */
		/*****************************************************************************/
		unsigned get_number_of_vertices() const;
		unsigned get_number_of_edges() const;
		unsigned get_number_of_triangles() const;



		/*****************************************************************************/
		/*                                                                           */
		/*    - Returns number of edges with 'NEUMANN' marker						 */
		/*    - Returns number of edges with 'DIRICHLET' marker						 */
		/*																			 */
		/*****************************************************************************/
		unsigned get_num_neumann_edges() const;
		unsigned get_num_dirichlet_edges() const;



		/*****************************************************************************/
		/*                                                                           */
		/*    - Write primitive's coordinates to a text-file						 */
		/*																			 */
		/*****************************************************************************/
		void export_vertices(std::string & fileName) const;
		void export_edges(std::string & fileName) const;
		void export_triangles(std::string & fileName) const;



	private:


		/*****************************************************************************/
		/*                                                                           */
		/*  Containers for geometric primitives										 */
		/*                                                                           */
		/*****************************************************************************/
		std::vector<v_pointer> vertices;
		std::vector<e_pointer> edges;
		std::vector<t_pointer> triangles;


		/*****************************************************************************/
		/*                                                                           */
		/*  Vertices and edges of the super triangle box					         */
		/*                                                                           */
		/*****************************************************************************/
		v_pointer superTriangleV0 = NULL;
		v_pointer superTriangleV1 = NULL;
		v_pointer superTriangleV2 = NULL;

		e_pointer superTriangleE0 = NULL;
		e_pointer superTriangleE1 = NULL;
		e_pointer superTriangleE2 = NULL;


		/*****************************************************************************/
		/*                                                                           */
		/*  Geometric predicates	Geometric functions								 */
		/*                                                                           */
		/*  Template parameter Arithmetic provides Fast or Robust floating-point     */
		/*  arithmetic                                                               */
		/*                                                                           */
		/*****************************************************************************/
		GeometricPredicates<Arithmetic> Predicates;


		/*****************************************************************************/
		/*                                                                           */
		/*  Number of vertices, edges and triangles							         */
		/*                                                                           */
		/*****************************************************************************/
		unsigned numberOfVertices	= 0;
		unsigned numberOfEdges		= 0;
		unsigned numberOfTriangles	= 0;

		unsigned numberOfNeumannEdges	= 0;
		unsigned numberOfDirichletEdges = 0;


		/*****************************************************************************/
		/*                                                                           */
		/*  Triangles labled as outside. They are to be deleted (holes, outside)     */
		/*                                                                           */
		/*****************************************************************************/
		std::vector<t_pointer> trianglesOutside;


		/*****************************************************************************/
		/*                                                                           */
		/*  Additional vertices resulting from inserting constraints (segments)	by   */
		/*  intersection method														 */
		/*                                                                           */
		/*****************************************************************************/
		std::vector<v_pointer> segmentIntersectionVertices;











		/*****************************************************************************/
		/*                                                                           */
		/*  Initialize and remove the super triangle box							 */
		/*                                                                           */
		/*****************************************************************************/
		void super_triangle_create(PlanarStraightLineGraph const & pslg);
		void super_triangle_remove();

	
		/*****************************************************************************/
		/*																			 */
		/*		Returns a triangle which should be a good starting triangle          */
		/*		for vertex location.												 */
		/*																			 */
		/*		Randomly picks few vertices from triangulation and compute			 */
		/*		their distance from the sought point. The closest vertex's			 */
		/*		adjacent triangle is returned										 */
		/*																			 */
		/*****************************************************************************/
		t_pointer & startingTriangleGuess(v_pointer const p);


		/*****************************************************************************/
		/*                                                                           */
		/*  locate_vertex		Locate new inserted vertex p in the triangulation	 */
		/*																			 */
		/*	Vertex can be found either in a triangle, on a edge or on an			 */
		/*	existing vertex															 */
		/*																			 */
		/*																			 */
		/*****************************************************************************/
		Location locate_vertex(t_pointer & t, e_pointer & e, v_pointer const p);

		
		/*****************************************************************************/
		/*																			 */
		/*	insert_vertex		Inserts incrementally single vertex into the		 */
		/*						triangulation and returns its pointer				 */
		/*																			 */
		/*	t is the starting triangle from which the vertex is being sought. If	 */
		/*	the t is not provided, then the starting triangle is the last created	 */
		/*																			 */
		/*****************************************************************************/
		v_pointer insert_vertex(v_pointer const & v, t_pointer t = NULL);


		/*****************************************************************************/
		/*                                                                           */
		/*  insert_vertex_into_triangle		Inserts new vertex p into the triangle t */
		/*																			 */
		/*	Vertex falls strictly into the interior of the triangle breaking it		 */
		/*	into three new triangles												 */
		/*																			 */
		/*																			 */
		/*	insert_vertex_into_edge			Inserts new veretx p into the edge e	 */
		/*																			 */
		/*	Vertex falls strictly on one of triangle's three edges creating	four	 */
		/*	new triangles															 */
		/*																			 */
		/*****************************************************************************/
		t_pointer insert_vertex_into_triangle(t_pointer const & t, v_pointer const & p);
		t_pointer insert_vertex_into_edge(e_pointer const & e, v_pointer const & p);

		
		/*****************************************************************************/
		/*                                                                           */
		/*  rotate_triangle_pair	Legalize each edge of the triangle t			 */
		/*                                                                           */
		/*  Legalization is done via rotating edge's adjacent triangles clockwise    */
		/*  by 90° degrees therefore the non-delaunay edge is flipped                */
		/*																			 */
		/*****************************************************************************/
		bool legalize(t_pointer const & t, v_pointer const & v);
		void rotate_triangle_pair(t_pointer const & t, v_pointer const & v, t_pointer const & ot, v_pointer const & ov);


		/*****************************************************************************/
		/*                                                                           */
		/*	locate_edge		Finds starting triangle contaning vertex a from which    */
		/*					the second vertex b of the constrained edge ab is sought.*/
		/*					Marker is the boundary marker of the constrained edge	 */
		/*																			 */
		/*****************************************************************************/
		t_pointer locate_edge(v_pointer const & a, v_pointer const & b, MarkerEdge marker);


		/*****************************************************************************/
		/*                                                                           */
		/*	insert_constraint		Insert constrained edge							 */
		/*                                                                           */
		/*****************************************************************************/
		void insert_constraint(v_pointer const & a, v_pointer const & b, MarkerEdge const & marker);



		/*****************************************************************************/
		/*                                                                           */
		/*	Mark triangles which are condemned to destruction     !m! 666 !m!        */
		/*                                                                           */
		/*****************************************************************************/
		void mark_triangle(t_pointer const & t, int const previousIndex);
		void remove_marked_triangles();


		/*****************************************************************************/
		/*                                                                           */
		/*	insert_constraint_intersection		Insert constrained edge using the	 */
		/*										intersection algorithm				 */
		/*                                                                           */
		/*****************************************************************************/
		void insert_constraint_intersection(v_pointer const a, v_pointer const b, MarkerEdge const marker);
		

		/*****************************************************************************/
		/*                                                                           */
		/*	remove		Removes triangle/edge/veretx from the respective conatiner   */
		/*                                                                           */
		/*  replace		Set the entity t into appropriate position of the container  */
		/*                                                                           */
		/*	clear		Clear all data and free allocated memory					 */
		/*                                                                           */
		/*****************************************************************************/
		template<typename T>
		void remove(std::vector<T> & vec, T t);
		template<typename T>
		void replace(std::vector<T> & vec, T t);
		void clear();


	};


	template <GeometricPredicatesArithmetic Arithmetic>
	class RuppertTriangulation {


	public:


		/*****************************************************************************/
		/*                                                                           */
		/*	Constructor and destructor		The input is pslg and seed vertices      */
		/*									defining holes							 */
		/*																			 */
		/*****************************************************************************/
		RuppertTriangulation(PlanarStraightLineGraph const & pslg, std::vector<Vertex> const & seeds);
		~RuppertTriangulation();


		/*****************************************************************************/
		/*                                                                           */
		/*  Returns the number of the triangulation primitives						 */
		/*																			 */
		/*****************************************************************************/
		unsigned get_number_of_vertices() const;
		unsigned get_number_of_edges() const;
		unsigned get_number_of_triangles() const;



		/*****************************************************************************/
		/*                                                                           */
		/*    - Returns number of edges with 'NEUMANN' marker						 */
		/*    - Returns number of edges with 'DIRICHLET' marker						 */
		/*																			 */
		/*****************************************************************************/
		unsigned get_num_neumann_edges() const;
		unsigned get_num_dirichlet_edges() const;



		/*****************************************************************************/
		/*                                                                           */
		/*    - Write primitive's coordinates to a text-file						 */
		/*																			 */
		/*****************************************************************************/
		void export_vertices(std::string & fileName) const;
		void export_edges(std::string & fileName) const;
		void export_triangles(std::string & fileName) const;



	private:


		/*****************************************************************************/
		/*                                                                           */
		/*  Containers for geometric primitives										 */
		/*                                                                           */
		/*****************************************************************************/
		std::vector<v_pointer> vertices;
		std::vector<e_pointer> edges;
		std::vector<t_pointer> triangles;


		/*****************************************************************************/
		/*                                                                           */
		/*  Vertices and edges of the super triangle box					         */
		/*                                                                           */
		/*****************************************************************************/
		v_pointer superTriangleV0 = NULL;
		v_pointer superTriangleV1 = NULL;
		v_pointer superTriangleV2 = NULL;

		e_pointer superTriangleE0 = NULL;
		e_pointer superTriangleE1 = NULL;
		e_pointer superTriangleE2 = NULL;


		/*****************************************************************************/
		/*                                                                           */
		/*  Geometric predicates	Geometric functions								 */
		/*                                                                           */
		/*  Template parameter Arithmetic provides Fast or Robust floating-point     */
		/*  arithmetic                                                               */
		/*                                                                           */
		/*****************************************************************************/
		GeometricPredicates<Arithmetic> Predicates;


		/*****************************************************************************/
		/*                                                                           */
		/*  Number of vertices, edges and triangles							         */
		/*                                                                           */
		/*****************************************************************************/
		unsigned numberOfVertices = 0;
		unsigned numberOfEdges = 0;
		unsigned numberOfTriangles = 0;

		unsigned numberOfNeumannEdges = 0;
		unsigned numberOfDirichletEdges = 0;


		/*****************************************************************************/
		/*                                                                           */
		/*  Triangles labled as outside. They are to be deleted (holes, outside)     */
		/*                                                                           */
		/*****************************************************************************/
		std::vector<t_pointer> trianglesOutside;


		/*****************************************************************************/
		/*                                                                           */
		/*  Additional vertices resulting from inserting constraints (segments)	by   */
		/*  intersection method														 */
		/*                                                                           */
		/*****************************************************************************/
		std::vector<v_pointer> segmentIntersectionVertices;











		/*****************************************************************************/
		/*                                                                           */
		/*  Initialize and remove the super triangle box							 */
		/*                                                                           */
		/*****************************************************************************/
		void super_triangle_create(PlanarStraightLineGraph const & pslg);
		void super_triangle_remove();


		/*****************************************************************************/
		/*																			 */
		/*		Returns a triangle which should be a good starting triangle          */
		/*		for vertex location.												 */
		/*																			 */
		/*		Randomly picks few vertices from triangulation and compute			 */
		/*		their distance from the sought point. The closest vertex's			 */
		/*		adjacent triangle is returned										 */
		/*																			 */
		/*****************************************************************************/
		t_pointer & startingTriangleGuess(v_pointer const p);


		/*****************************************************************************/
		/*                                                                           */
		/*  locate_vertex		Locate new inserted vertex p in the triangulation	 */
		/*																			 */
		/*	Vertex can be found either in a triangle, on a edge or on an			 */
		/*	existing vertex															 */
		/*																			 */
		/*																			 */
		/*****************************************************************************/
		Location locate_vertex(t_pointer & t, e_pointer & e, v_pointer const p);


		/*****************************************************************************/
		/*																			 */
		/*	insert_vertex		Inserts incrementally single vertex into the		 */
		/*						triangulation and returns its pointer				 */
		/*																			 */
		/*	t is the starting triangle from which the vertex is being sought. If	 */
		/*	the t is not provided, then the starting triangle is the last created	 */
		/*																			 */
		/*****************************************************************************/
		v_pointer insert_vertex(v_pointer const & v, t_pointer t = NULL);


		/*****************************************************************************/
		/*                                                                           */
		/*  insert_vertex_into_triangle		Inserts new vertex p into the triangle t */
		/*																			 */
		/*	Vertex falls strictly into the interior of the triangle breaking it		 */
		/*	into three new triangles												 */
		/*																			 */
		/*																			 */
		/*	insert_vertex_into_edge			Inserts new veretx p into the edge e	 */
		/*																			 */
		/*	Vertex falls strictly on one of triangle's three edges creating	four	 */
		/*	new triangles															 */
		/*																			 */
		/*****************************************************************************/
		t_pointer insert_vertex_into_triangle(t_pointer const & t, v_pointer const & p);
		t_pointer insert_vertex_into_edge(e_pointer const & e, v_pointer const & p);


		/*****************************************************************************/
		/*                                                                           */
		/*  rotate_triangle_pair	Legalize each edge of the triangle t			 */
		/*                                                                           */
		/*  Legalization is done via rotating edge's adjacent triangles clockwise    */
		/*  by 90° degrees therefore the non-delaunay edge is flipped                */
		/*																			 */
		/*****************************************************************************/
		bool legalize(t_pointer const & t, v_pointer const & v);
		void rotate_triangle_pair(t_pointer const & t, v_pointer const & v, t_pointer const & ot, v_pointer const & ov);


		/*****************************************************************************/
		/*                                                                           */
		/*	locate_edge		Finds starting triangle contaning vertex a from which    */
		/*					the second vertex b of the constrained edge ab is sought.*/
		/*					Marker is the boundary marker of the constrained edge	 */
		/*																			 */
		/*****************************************************************************/
		t_pointer locate_edge(v_pointer const & a, v_pointer const & b, MarkerEdge marker);


		/*****************************************************************************/
		/*                                                                           */
		/*	insert_constraint		Insert constrained edge							 */
		/*                                                                           */
		/*****************************************************************************/
		void insert_constraint(v_pointer const & a, v_pointer const & b, MarkerEdge const & marker);



		/*****************************************************************************/
		/*                                                                           */
		/*	Mark triangles which are condemned to destruction     !m! 666 !m!        */
		/*                                                                           */
		/*****************************************************************************/
		void mark_triangle(t_pointer const & t, int const previousIndex);
		void remove_marked_triangles();


		/*****************************************************************************/
		/*                                                                           */
		/*	insert_constraint_intersection		Insert constrained edge using the	 */
		/*										intersection algorithm				 */
		/*                                                                           */
		/*****************************************************************************/
		void insert_constraint_intersection(v_pointer const a, v_pointer const b, MarkerEdge const marker);


		/*****************************************************************************/
		/*                                                                           */
		/*	remove		Removes triangle/edge/veretx from the respective conatiner   */
		/*                                                                           */
		/*  replace		Set the entity t into appropriate position of the container  */
		/*                                                                           */
		/*	clear		Clear all data and free allocated memory					 */
		/*                                                                           */
		/*****************************************************************************/
		template<typename T>
		void remove(std::vector<T> & vec, T t);
		template<typename T>
		void replace(std::vector<T> & vec, T t);
		void clear();


	};

}


#include "IncrementalTriangulation_impl.h"
#include "RuppertTriangulation_impl.h"
#include "PlanarStraightLineGraph_impl.h"