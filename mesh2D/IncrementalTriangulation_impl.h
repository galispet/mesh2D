#pragma once



namespace DelaunayTriangulation {



	/*****************************************************************************/
	/*                                                                           */
	/*		Constructor and destructor											 */
	/*																			 */
	/*****************************************************************************/
	template<GeometricPredicatesArithmetic Arithmetic>
	IncrementalTriangulation<Arithmetic>::IncrementalTriangulation(PlanarStraightLineGraph const & pslg, std::vector<Vertex> const & seeds) {



		size_t const sizePslgVertices = pslg.vertices.size();
		size_t const sizePslgSegments = pslg.segments.size();
		size_t const sizeSeeds		  = seeds.size();

		if (sizePslgVertices > INT_MAX || sizePslgSegments > INT_MAX || sizeSeeds > INT_MAX)
			throw std::overflow_error("Data size of inputed PSLG is larger than INT_MAX");


		unsigned const nv = static_cast<unsigned>(sizePslgVertices);
		unsigned const ns = static_cast<unsigned>(sizePslgSegments);
		unsigned const nh = static_cast<unsigned>(sizeSeeds);


		/*****************************************************************************/
		/*                                                                           */
		/*  Upper bounds for the number of triangles and edges (theory of graphs)	 */
		/*																			 */
		/*****************************************************************************/
		unsigned const expectedNumberOfVertices = (unsigned) floor((3.0 * nv / 2.0));

		vertices	.reserve(expectedNumberOfVertices);
		edges		.reserve(3 * expectedNumberOfVertices - 6);
		triangles	.reserve(2 * expectedNumberOfVertices - 5);


		/*****************************************************************************/
		/*                                                                           */
		/*  Incremental algorithm needs first initial triangle. Super triangle is    */
		/*	the intial triangle and it ensloses all the vertices and segmenets		 */
		/*	which are given in the input pslg										 */
		/*																			 */
		/*****************************************************************************/
		super_triangle_create(pslg);


		/*****************************************************************************/
		/*                                                                           */
		/*  Deep copy of vertices (memory allocation). IncrementalTriangulation will */
		/*	have its own memory pool.												 */
		/*																			 */
		/*	Method insert_vertex() internally push_back the new vertex nv into the	 */
		/*	container, if succeded insertion (non-overlapping vertex)				 */
		/*																			 */
		/*****************************************************************************/
		for (unsigned i = 0; i < nv; i++) {

			v_pointer const	v = pslg.vertices[i];
			v_pointer const w = new Vertex(v);

			insert_vertex(w);

		}


		/*****************************************************************************/
		/*                                                                           */
		/*  Insert constrained segments defined by the user and provided in th pslg. */
		/*																			 */
		/*  Segments of pslg are defined as 3-element unsigned array.				 */
		/*																			 */
		/*	0. Index of the first vertex 											 */
		/*	1. Index of the second vertex 											 */
		/*	2. Edge marker 															 */
		/*																			 */
		/*	Segments end point vertices are marked constrained (they can't be moved) */
		/*																			 */
		/*****************************************************************************/
		for (unsigned i = 0; i < ns; i++) {


			unsigned * const constraint = pslg.segments[i];


			unsigned const index0 = constraint[0];
			unsigned const index1 = constraint[1];

			MarkerEdge const segmentMarker = pslg.segmentsMarker[i];


			v_pointer const a = vertices[index0];
			v_pointer const b = vertices[index1];

			a->marker = MarkerVertex::Constrained;
			b->marker = MarkerVertex::Constrained;

			insert_constraint(a, b, segmentMarker);

		}


		/*****************************************************************************/
		/*                                                                           */
		/*  Mark triangles in holes to be deleted using seed points. These points 	 */
		/*	are located and then a virus is spread until it reaches closed boundary	 */
		/*	defined by segments														 */
		/*																			 */
		/*****************************************************************************/
		for (unsigned i = 0; i < nh; i++) {

			Vertex s			 = seeds[i];
			v_pointer const seed = &s;

			e_pointer e = NULL;
			t_pointer t = triangles.back();
			//t_pointer t = startingTriangleGuess(&seeds[i]);

			locate_vertex(t, e, seed);

			mark_triangle(t, -1);

		}



		/*****************************************************************************/
		/*                                                                           */
		/*  Mark triangles outside of the triangulation boundary. These triangles    */
		/*	are typically attached to super triangle.								 */
		/*																			 */
		/*****************************************************************************/
		mark_triangle(superTriangleV0->adjacentTriangle, -1);


		/*****************************************************************************/
		/*                                                                           */
		/*  Remove triangles and their edges which are marked for removing			 */
		/*																			 */
		/*****************************************************************************/
		remove_marked_triangles();


		/*****************************************************************************/
		/*                                                                           */
		/*  Remove vertices and edges of the super triangle							 */
		/*																			 */
		/*****************************************************************************/
		super_triangle_remove();


		/*****************************************************************************/
		/*                                                                           */
		/*  After deletetion of triangles, edges and vertices, we have to reindex it */
		/*	Indexing is important for efficient	load/remove in data containers		 */
		/*																			 */
		/*****************************************************************************/
		size_t const sizeVertices  = vertices.size();
		size_t const sizeEdges	   = edges.size();
		size_t const sizeTriangles = triangles.size();


		if (sizeVertices > INT_MAX || sizeEdges > INT_MAX || sizeTriangles > INT_MAX) {

			clear();
			throw std::overflow_error("Data size of triangulation primitives is larger than INT_MAX");

		}


		numberOfVertices	= static_cast<unsigned>(sizeVertices);
		numberOfEdges		= static_cast<unsigned>(sizeEdges);
		numberOfTriangles	= static_cast<unsigned>(sizeTriangles);

		numberOfNeumannEdges   = (unsigned)std::count_if(edges.begin(), edges.end(), [](e_pointer const & e) { return e->marker == MarkerEdge::Neumann; });
		numberOfDirichletEdges = (unsigned)std::count_if(edges.begin(), edges.end(), [](e_pointer const & e) { return e->marker == MarkerEdge::Dirichlet; });


		for (unsigned i = 0; i < numberOfVertices; i++)
			vertices[i]->index = i;

		for (unsigned i = 0; i < numberOfEdges; i++)
			edges[i]->index = i;

		for (unsigned i = 0; i < numberOfTriangles; i++)
			triangles[i]->index = i;


	};
	template<GeometricPredicatesArithmetic Arithmetic>
	IncrementalTriangulation<Arithmetic>::~IncrementalTriangulation() {

		clear();

	};


	/*****************************************************************************/
	/*                                                                           */
	/*		Returns number of triangulation primitives							 */
	/*																			 */
	/*****************************************************************************/
	template<GeometricPredicatesArithmetic Arithmetic>
	unsigned IncrementalTriangulation<Arithmetic>::get_number_of_vertices() const {

		return numberOfVertices;

	};
	template<GeometricPredicatesArithmetic Arithmetic>
	unsigned IncrementalTriangulation<Arithmetic>::get_number_of_edges() const {

		return numberOfEdges;

	};
	template<GeometricPredicatesArithmetic Arithmetic>
	unsigned IncrementalTriangulation<Arithmetic>::get_number_of_triangles() const {

		return numberOfTriangles;

	};


	/*****************************************************************************/
	/*                                                                           */
	/*		Returns number of neumann and dirichlet edges						 */
	/*																			 */
	/*****************************************************************************/
	template<GeometricPredicatesArithmetic Arithmetic>
	unsigned IncrementalTriangulation<Arithmetic>::get_num_neumann_edges() const {

		return numberOfNeumannEdges;

	};
	template<GeometricPredicatesArithmetic Arithmetic>
	unsigned IncrementalTriangulation<Arithmetic>::get_num_dirichlet_edges() const {

		return numberOfDirichletEdges;

	};


	/*****************************************************************************/
	/*                                                                           */
	/*		Write primitives coordinates to a .txt file							 */
	/*																			 */
	/*****************************************************************************/
	template<GeometricPredicatesArithmetic Arithmetic>
	void IncrementalTriangulation<Arithmetic>::export_vertices(std::string & fileName) const {


		std::ofstream stream;

		stream.open(fileName);

		for (unsigned i = 0; i < numberOfVertices; i++) {

			double const x = vertices[i]->x;
			double const y = vertices[i]->y;

			stream << x << "  " << y << std::endl;

		}

		stream.close();

	};
	template<GeometricPredicatesArithmetic Arithmetic>
	void IncrementalTriangulation<Arithmetic>::export_edges(std::string & fileName) const {


		std::ofstream stream;

		stream.open(fileName);

		/*****************************************************************************/
		/*                                                                           */
		/*		Gnuplot																 */
		/*																			 */
		/*****************************************************************************/
		for (unsigned i = 0; i < numberOfEdges; i++) {

			e_pointer const e = edges[i];

			v_pointer const a = e->vertices[0];
			v_pointer const b = e->vertices[1];

			double const ax = a->x;
			double const ay = a->y;

			double const bx = b->x;
			double const by = b->y;

			stream << ax << "  " << ay << std::endl;
			stream << bx << "  " << by << std::endl << std::endl;

		}

		/*****************************************************************************/
		/*                                                                           */
		/*		Matlab																 */
		/*																			 */
		/*****************************************************************************/
		/*for (unsigned i = 0; i < numberOfEdges; i++) {

			e_pointer const e = edges[i];

			v_pointer const a = e->vertices[0];
			v_pointer const b = e->vertices[1];

			double const ax = a->x;
			double const ay = a->y;

			double const bx = b->x;
			double const by = b->y;

			stream << ax << "  " << ay << "  " << bx << "  " << by << std::endl;

		}*/


		stream.close();

	};
	template<GeometricPredicatesArithmetic Arithmetic>
	void IncrementalTriangulation<Arithmetic>::export_triangles(std::string & fileName) const {


		std::ofstream stream;

		stream.open(fileName);

		/*****************************************************************************/
		/*                                                                           */
		/*		Gnuplot																 */
		/*																			 */
		/*****************************************************************************/
		for (unsigned i = 0; i < numberOfTriangles; i++) {

			t_pointer const t = triangles[i];

			double const ax = t->vertices[0]->x;
			double const ay = t->vertices[0]->y;

			double const bx = t->vertices[1]->x;
			double const by = t->vertices[1]->y;

			double const cx = t->vertices[2]->x;
			double const cy = t->vertices[2]->y;

			stream << ax << "  " << ay << std::endl;
			stream << bx << "  " << by << std::endl;
			stream << cx << "  " << cy << std::endl;
			stream << ax << "  " << ay << std::endl << std::endl;

		}

		/*****************************************************************************/
		/*                                                                           */
		/*		Matlab																 */
		/*																			 */
		/*****************************************************************************/
		/*for (unsigned i = 0; i < numberOfTriangles; i++) {

			t_pointer const t = triangles[i];

			double const ax = t->vertices[0]->x;
			double const ay = t->vertices[0]->y;

			double const bx = t->vertices[1]->x;
			double const by = t->vertices[1]->y;

			double const cx = t->vertices[2]->x;
			double const cy = t->vertices[2]->y;

			stream << ax << "  " << ay << "  " << bx << "  " << by << "  " << cx << "  " << cy<< "  " << ax << "  " << ay <<  std::endl;

		}*/


		stream.close();

	};







	/*****************************************************************************/
	/*                                                                           */
	/*		Construct and remove super triangle						             */
	/*                                                                           */
	/*****************************************************************************/
	template<GeometricPredicatesArithmetic Arithmetic>
	void IncrementalTriangulation<Arithmetic>::super_triangle_create(PlanarStraightLineGraph const & pslg) {


		/*****************************************************************************/
		/*                                                                           */
		/*		Get maximum and minimum coordinate on x, y axis			             */
		/*                                                                           */
		/*****************************************************************************/
		double const xMin = pslg.xMinimum;
		double const xMax = pslg.xMaximum;
		double const yMin = pslg.yMinimum;
		double const yMax = pslg.yMaximum;

		/*double minX = pslg.vertices[0]->x;
		double minY = pslg.vertices[0]->y;
		double maxX = minX;
		double maxY = minY;

		size_t const nv = pslg.vertices.size();

		for (size_t i = 0; i < nv; i++) {

			v_pointer const v = pslg.vertices[i];

			if (v->x < minX) minX = v->x;
			if (v->y < minY) minY = v->y;
			if (v->x > maxX) maxX = v->x;
			if (v->y > maxY) maxY = v->y;

		}*/

		double const dx = xMax - xMin;
		double const dy = yMax - yMin;
		double const ox = 0.5 * (xMax + xMin);
		double const oy = 0.5 * (yMax + yMin);

		double const d = std::max(dx, dy);
		double const c = 15.0;

		/*****************************************************************************/
		/*                                                                           */
		/*		Super triangle vertices									             */
		/*                                                                           */
		/*****************************************************************************/
		superTriangleV0 = new Vertex(ox - c * d, oy - d);
		superTriangleV1 = new Vertex(ox + c * d, oy - d);
		superTriangleV2 = new Vertex(ox, oy + c * d);

		//st_v0 = new Vertex(-3 * M, -3 * M);
		//st_v1 = new Vertex(3 * M, 0.0);
		//st_v2 = new Vertex(0.0, 3 * M);

		superTriangleV0->marker = MarkerVertex::Constrained;
		superTriangleV1->marker = MarkerVertex::Constrained;
		superTriangleV2->marker = MarkerVertex::Constrained;


		/*****************************************************************************/
		/*                                                                           */
		/*		Super triangle edges									             */
		/*                                                                           */
		/*****************************************************************************/
		superTriangleE0 = new Edge(superTriangleV1, superTriangleV2);
		superTriangleE1 = new Edge(superTriangleV2, superTriangleV0);
		superTriangleE2 = new Edge(superTriangleV0, superTriangleV1);

		superTriangleE0->isConstrained = true;
		superTriangleE1->isConstrained = true;
		superTriangleE2->isConstrained = true;


		/*****************************************************************************/
		/*                                                                           */
		/*		Super triangle											             */
		/*                                                                           */
		/*****************************************************************************/
		t_pointer const superTriangle = new Triangle(superTriangleV0, superTriangleV1, superTriangleV2);


		/*****************************************************************************/
		/*                                                                           */
		/*		Assign edges its edges to appropriate position			             */
		/*                                                                           */
		/*****************************************************************************/
		superTriangle->set_edge(0) = superTriangleE0;
		superTriangle->set_edge(1) = superTriangleE1;
		superTriangle->set_edge(2) = superTriangleE2;


		/*****************************************************************************/
		/*                                                                           */
		/*		Assign neighboring triangles to edges					             */
		/*                                                                           */
		/*****************************************************************************/
		superTriangleE0->set_neighbor(0) = superTriangle;
		superTriangleE0->set_neighbor(1) = NULL;

		superTriangleE1->set_neighbor(0) = superTriangle;
		superTriangleE1->set_neighbor(1) = NULL;

		superTriangleE2->set_neighbor(0) = superTriangle;
		superTriangleE2->set_neighbor(1) = NULL;


		/*****************************************************************************/
		/*                                                                           */
		/*		Initilize triangulation with super triangle				             */
		/*                                                                           */
		/*****************************************************************************/
		triangles.push_back(superTriangle);


		/*****************************************************************************/
		/*                                                                           */
		/*		Indexing is importent. Indeces are copied to new created triangles   */
		/*      from destroyed ones.                                                 */
		/*                                                                           */
		/*****************************************************************************/
		superTriangle->index = 0;

	};
	template<GeometricPredicatesArithmetic Arithmetic>
	void IncrementalTriangulation<Arithmetic>::super_triangle_remove() {


		delete superTriangleV0;
		delete superTriangleV1;
		delete superTriangleV2;

		delete superTriangleE0;
		delete superTriangleE1;
		delete superTriangleE2;

		superTriangleV0 = NULL;
		superTriangleV1 = NULL;
		superTriangleV2 = NULL;

		superTriangleE0 = NULL;
		superTriangleE1 = NULL;
		superTriangleE2 = NULL;

	};


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
	template<GeometricPredicatesArithmetic Arithmetic>
	t_pointer & IncrementalTriangulation<Arithmetic>::startingTriangleGuess(v_pointer const p) {


		int const min = 1;
		int const max = (int)vertices.size() - 2 < min ? min : (int)vertices.size() - 2;
		int const numberOfSampleVertices = floor(cbrt(vertices.size()));


		std::random_device					randomDevice;
		std::mt19937 						randomNumberGenerator(randomDevice());
		std::uniform_int_distribution<int>  uniformDistribution(min, max);


		double const px = p->x;
		double const py = p->y;



		v_pointer closestVertex = vertices.back();

		double vx = closestVertex->x;
		double vy = closestVertex->y;

		double distance = (px - vx)*(px - vx) + (py - vy)*(py - vy);
		double distanceMinimum = distance;


		vx = vertices.front()->x;
		vy = vertices.front()->y;

		distance = (px - vx)*(px - vx) + (py - vy)*(py - vy);

		if (distance < distanceMinimum) {

			closestVertex = vertices.front();
			distanceMinimum = distance;

		}



		for (unsigned i = 0; i < numberOfSampleVertices; i++) {

			int randomInteger = uniformDistribution(randomNumberGenerator) - 1;

			v_pointer const v = vertices[randomInteger];

			vx = v->x;
			vy - v->y;

			distance = (px - vx)*(px - vx) + (py - vy)*(py - vy);

			if (distance < distanceMinimum)
				closestVertex = v;

		}

		return closestVertex->adjacentTriangle;

	};


	/*****************************************************************************/
	/*                                                                           */
	/*    Locate inserted vertex. Return on edge if the vertex was found on      */
	/*	  the edge or return in triangle if the vertex was found in the			 */
	/*    triangle                                                               */
	/*                                                                           */
	/*****************************************************************************/
	template<GeometricPredicatesArithmetic Arithmetic>
	Location IncrementalTriangulation<Arithmetic>::locate_vertex(t_pointer & t, e_pointer & e, v_pointer const p) {


		/*****************************************************************************/
		/*                                                                           */
		/*		Guibas and Stolfi edge walking algorithm                             */
		/*																			 */
		/*****************************************************************************/
		v_pointer org  = t->vertices[0];
		v_pointer dest = t->vertices[1];
		v_pointer apex = t->vertices[2];

		double const px = p->x;
		double const py = p->y;


		/*****************************************************************************/
		/*                                                                           */
		/*		Initialization step			                                         */
		/*																			 */
		/*		Rotate the triangle so that sougth point is on the left of the 		 */
		/*		line defined by the edge e = (org,dest)								 */
		/*																			 */
		/*****************************************************************************/
		double orientAhead = Predicates.orientation(org, dest, p);

		if (orientAhead < 0.0) {

			//do {
			org  = apex;
			dest = t->get_vertex_ccw(org);
			apex = t->get_vertex_cw(org);
			//} while (Predicates.orientation(org, dest, p) < 0.0);

		}
		else if (abs(orientAhead) <= EpsilonOnEdge) {

			double const ox = org->x;
			double const oy = org->y;

			double const dx = dest->x;
			double const dy = dest->y;

			/*****************************************************************************/
			/*																			 */
			/*		Check if the sought veretx lies on the edge e = (org, dest)			 */
			/*		that is in between vertices org and dest. We know that orient = 0	 */
			/*		so this condition is sufficient										 */
			/*																			 */
			/*****************************************************************************/
			if ((ox < px == px < dx) && (oy < py == py < dy)) {

				e = t->get_edge(apex);

				return Location::OnEdge;

			}
		}


		/*****************************************************************************/
		/*																			 */
		/*		Walking step														 */
		/*																			 */
		/*		Walk toward the sought point through triangles holding the vertex p  */
		/*		left of the line defined by the edge e = (org, dest)				 */
		/*																			 */
		/*****************************************************************************/
		while (true) {

			double orientDest = Predicates.orientation(org, apex, p);
			double orientOrg  = Predicates.orientation(apex, dest, p);

			if (orientDest > 0.0) {

				if (orientOrg > 0.0) {

					double const ox = org->x;
					double const oy = org->y;

					double const dx = dest->x;
					double const dy = dest->y;

					double const ax = apex->x;
					double const ay = apex->y;

					if ((ax - px) * (dx - ox) + (ay - py) * (dy - oy) > 0.0) {

						t	 = t->get_neighbor(dest);
						dest = apex;
						apex = t->get_vertex_cw(org);

					}
					else {

						t	 = t->get_neighbor(org);
						org  = apex;
						apex = t->get_vertex_cw(org);

					}
				}
				else {

					t	 = t->get_neighbor(dest);
					dest = apex;
					apex = t->get_vertex_cw(org);

				}
			}
			else {

				if (orientOrg > 0.0) {

					t	 = t->get_neighbor(org);
					org  = apex;
					apex = t->get_vertex_cw(org);

				}
				else {

					if (abs(orientDest) <= EpsilonOnEdge) {

						e = t->get_edge(dest);
						return Location::OnEdge;

					}
					if (abs(orientOrg) <= EpsilonOnEdge) {

						e = t->get_edge(org);
						return Location::OnEdge;

					}

					return Location::InTriangle;

				}
			}
		}

	};


	/*****************************************************************************/
	/*																			 */
	/*		Inserts incrementally single vertex into the triangulation           */
	/*																			 */
	/*****************************************************************************/
	template<GeometricPredicatesArithmetic Arithmetic>
	v_pointer IncrementalTriangulation<Arithmetic>::insert_vertex(v_pointer const & v, t_pointer t) {


		e_pointer e = NULL;
		t_pointer k = t ? t : triangles.back();
		//t_pointer k = vertices.size()<1 ? triangles.back() : startingTriangleGuess(v);


		switch (locate_vertex(k, e, v)) {

			case Location::InTriangle:

				insert_vertex_into_triangle(k, v);
				break;

			case Location::OnEdge:

				insert_vertex_into_edge(e, v);
				break;

		}

		/*****************************************************************************/
		/*                                                                           */
		/*		Add vertex v into triangulation	and assign it an index	             */
		/*																			 */
		/*****************************************************************************/
		v->index = (int) vertices.size();
		vertices.push_back(v);


		return v;

	};


	/*****************************************************************************/
	/*                                                                           */
	/*		Routines which create new triangles based on the location of the     */
	/*		inserted point														 */
	/*																			 */
	/*****************************************************************************/
	template<GeometricPredicatesArithmetic Arithmetic>
	t_pointer IncrementalTriangulation<Arithmetic>::insert_vertex_into_triangle(t_pointer const & t, v_pointer const & p) {


		v_pointer const v0 = t->vertices[0];
		v_pointer const v1 = t->vertices[1];
		v_pointer const v2 = t->vertices[2];

		t_pointer const n0 = t->neighbors[0];
		t_pointer const n1 = t->neighbors[1];
		t_pointer const n2 = t->neighbors[2];

		t_pointer const t0 = new Triangle(p, v1, v2);
		t_pointer const t1 = new Triangle(p, v2, v0);
		t_pointer const t2 = new Triangle(p, v0, v1);

		e_pointer const e0 = new Edge(p, v0);
		e_pointer const e1 = new Edge(p, v1);
		e_pointer const e2 = new Edge(p, v2);


		/*****************************************************************************/
		/*                                                                           */
		/*		Assign new edges to new triangles						             */
		/*																			 */
		/*****************************************************************************/
		t0->set_edge(0) = t->edges[0];
		t0->set_edge(1) = e2;
		t0->set_edge(2) = e1;

		t1->set_edge(0) = t->edges[1];
		t1->set_edge(1) = e0;
		t1->set_edge(2) = e2;

		t2->set_edge(0) = t->edges[2];
		t2->set_edge(1) = e1;
		t2->set_edge(2) = e0;


		/*****************************************************************************/
		/*                                                                           */
		/*		Assign neighboring triangles to edges					             */
		/*																			 */
		/*****************************************************************************/
		t->edges[0]->set_neighbor(0) = t0;
		t->edges[0]->set_neighbor(1) = n0;

		t->edges[1]->set_neighbor(0) = t1;
		t->edges[1]->set_neighbor(1) = n1;

		t->edges[2]->set_neighbor(0) = t2;
		t->edges[2]->set_neighbor(1) = n2;


		/*****************************************************************************/
		/*                                                                           */
		/*		Set indeces of new triangles							             */
		/*																			 */
		/*****************************************************************************/
		int const nt = (int)triangles.size();

		t0->index = t->index;
		t1->index = nt;
		t2->index = nt + 1;


		/*****************************************************************************/
		/*                                                                           */
		/*		Set indeces of new edges								             */
		/*																			 */
		/*****************************************************************************/
		int const ne = (int)edges.size();

		e0->index = ne;
		e1->index = ne + 1;
		e2->index = ne + 2;


		/*****************************************************************************/
		/*                                                                           */
		/*		Set pointer to any vertex's p adjacent triangle. Also there is a     */
		/*		possibility that some vertices have pointers to the deleted			 */
		/*		triangles t0 and t1. We must renew the adjacency.					 */
		/*                                                                           */
		/*****************************************************************************/
		p->adjacentTriangle = t0;

		v0->adjacentTriangle = t1;
		v1->adjacentTriangle = t2;
		v2->adjacentTriangle = t0;


		/*****************************************************************************/
		/*                                                                           */
		/*		Remove and delete splitted triangle consisting of t0, t1, t2 from    */
		/*		the triangulation													 */
		/*																			 */
		/*****************************************************************************/
		//triangles[t->index] = t0;
		remove(triangles, t);
		delete t;



		/*****************************************************************************/
		/*                                                                           */
		/*		Assign neighboring triangles to edges					             */
		/*																			 */
		/*****************************************************************************/
		e0->set_neighbor(0) = t2;
		e0->set_neighbor(1) = t1;

		e1->set_neighbor(0) = t2;
		e1->set_neighbor(1) = t0;

		e2->set_neighbor(0) = t0;
		e2->set_neighbor(1) = t1;


		/*****************************************************************************/
		/*                                                                           */
		/*		Assign neighborhoods of triangles						             */
		/*																			 */
		/*****************************************************************************/
		t0->set_neighbor(t1);
		t0->set_neighbor(n0);

		t1->set_neighbor(t2);
		t1->set_neighbor(n1);

		t2->set_neighbor(t0);
		t2->set_neighbor(n2);


		/*****************************************************************************/
		/*                                                                           */
		/*		Add new triangles to triangulation						             */
		/*																			 */
		/*****************************************************************************/
		triangles.push_back(t0);
		triangles.push_back(t1);
		triangles.push_back(t2);


		/*****************************************************************************/
		/*                                                                           */
		/*		Add new edges to triangulation							             */
		/*																			 */
		/*****************************************************************************/
		edges.push_back(e0);
		edges.push_back(e1);
		edges.push_back(e2);


		/*****************************************************************************/
		/*                                                                           */
		/*		Perform edge-flip algorithm to legalize non-Delaunay edges           */
		/*																			 */
		/*****************************************************************************/
		legalize(t0, p);
		legalize(t1, p);
		legalize(t2, p);

		return t0;

	};
	template<GeometricPredicatesArithmetic Arithmetic>
	t_pointer IncrementalTriangulation<Arithmetic>::insert_vertex_into_edge(e_pointer const & e, v_pointer const & p) {


		t_pointer const t0 = e->get_neighbor(0);
		t_pointer const t1 = e->get_neighbor(1);

		int const v0_index = t0->get_edge_index(e);
		int const v1_index = t1->get_edge_index(e);

		v_pointer const v0 = t0->get_vertex(v0_index);
		v_pointer const v1 = t1->get_vertex(v1_index);

		v_pointer const w0 = t0->get_vertex_ccw(v0);
		v_pointer const w1 = t0->get_vertex_cw(v0);


		t_pointer const n0 = t0->get_neighbor_ccw(v0);
		t_pointer const n1 = t0->get_neighbor_cw(v0);
		t_pointer const n2 = t1->get_neighbor_ccw(v1);
		t_pointer const n3 = t1->get_neighbor_cw(v1);


		bool const edgeIsConstrainedFlag = e->isConstrained;
		MarkerEdge const edgeMarker = e->marker;


		t_pointer const k0 = new Triangle(p, v0, w0);
		t_pointer const k1 = new Triangle(p, w1, v0);
		t_pointer const k2 = new Triangle(p, v1, w1);
		t_pointer const k3 = new Triangle(p, w0, v1);

		e_pointer const e0 = new Edge(p, w0);
		e_pointer const e1 = new Edge(p, v0);
		e_pointer const e2 = new Edge(p, w1);
		e_pointer const e3 = new Edge(p, v1);


		/*****************************************************************************/
		/*                                                                           */
		/*		Assign new edges to new triangles						             */
		/*																			 */
		/*****************************************************************************/
		k0->set_edge(0) = t0->get_edge_ccw(v0);
		k0->set_edge(1) = e0;
		k0->set_edge(2) = e1;

		k1->set_edge(0) = t0->get_edge_cw(v0);
		k1->set_edge(1) = e1;
		k1->set_edge(2) = e2;

		k2->set_edge(0) = t1->get_edge_ccw(v1);
		k2->set_edge(1) = e2;
		k2->set_edge(2) = e3;

		k3->set_edge(0) = t1->get_edge_cw(v1);
		k3->set_edge(1) = e3;
		k3->set_edge(2) = e0;


		/*****************************************************************************/
		/*                                                                           */
		/*		Copy constraint flag to new edges which were formed by splitting     */
		/*		the edge e with the vertex p										 */
		/*																			 */
		/*****************************************************************************/
		e0->isConstrained = edgeIsConstrainedFlag;
		e2->isConstrained = edgeIsConstrainedFlag;

		e0->marker = edgeMarker;
		e2->marker = edgeMarker;


		/*****************************************************************************/
		/*                                                                           */
		/*		Assign neighboring triangles to edges					             */
		/*																			 */
		/*****************************************************************************/
		e0->set_neighbor(0) = k3;
		e0->set_neighbor(1) = k0;

		e1->set_neighbor(0) = k0;
		e1->set_neighbor(1) = k1;

		e2->set_neighbor(0) = k1;
		e2->set_neighbor(1) = k2;

		e3->set_neighbor(0) = k2;
		e3->set_neighbor(1) = k3;

		k0->get_edge(0)->set_neighbor(0) = k0;
		k0->get_edge(0)->set_neighbor(1) = n0;

		k1->get_edge(0)->set_neighbor(0) = k1;
		k1->get_edge(0)->set_neighbor(1) = n1;

		k2->get_edge(0)->set_neighbor(0) = k2;
		k2->get_edge(0)->set_neighbor(1) = n2;

		k3->get_edge(0)->set_neighbor(0) = k3;
		k3->get_edge(0)->set_neighbor(1) = n3;


		/*****************************************************************************/
		/*                                                                           */
		/*		Assign neighborhoods of triangles						             */
		/*																			 */
		/*****************************************************************************/
		k0->set_neighbor(k1);
		k0->set_neighbor(n0);

		k1->set_neighbor(k2);
		k1->set_neighbor(n1);

		k2->set_neighbor(k3);
		k2->set_neighbor(n2);

		k3->set_neighbor(k0);
		k3->set_neighbor(n3);


		/*****************************************************************************/
		/*                                                                           */
		/*		Set pointer to any vertex's p adjacent triangle. Also there is a     */
		/*		possibility that some vertices have pointers to the deleted			 */
		/*		triangles t0 and t1. We must renew the adjacency.					 */
		/*                                                                           */
		/*****************************************************************************/
		p->adjacentTriangle = k0;

		v0->adjacentTriangle = k1;
		w1->adjacentTriangle = k2;
		v1->adjacentTriangle = k3;
		w0->adjacentTriangle = k0;


		/*****************************************************************************/
		/*                                                                           */
		/*		Remove and delete splitted triangles and edge from the				 */
		/*		triangulation														 */
		/*																			 */
		/*****************************************************************************/
		remove(edges, e);
		delete e;

		remove(triangles, t0);
		delete t0;

		remove(triangles, t1);
		delete t1;


		/*****************************************************************************/
		/*                                                                           */
		/*		Add new triangles to triangulation						             */
		/*																			 */
		/*****************************************************************************/
		triangles.push_back(k0);
		triangles.push_back(k1);
		triangles.push_back(k2);
		triangles.push_back(k3);


		/*****************************************************************************/
		/*                                                                           */
		/*		Add new edges to triangulation							             */
		/*																			 */
		/*****************************************************************************/
		edges.push_back(e0);
		edges.push_back(e1);
		edges.push_back(e2);
		edges.push_back(e3);


		/*****************************************************************************/
		/*                                                                           */
		/*		Perform edge-flip algorithm to legalize non-Delaunay edges           */
		/*																			 */
		/*****************************************************************************/
		legalize(k0, p);
		legalize(k1, p);
		legalize(k2, p);
		legalize(k3, p);

		return k0;

	};


	/*****************************************************************************/
	/*                                                                           */
	/*		Legalization of edges (edge flip)									 */
	/*																			 */
	/*****************************************************************************/
	template<GeometricPredicatesArithmetic Arithmetic>
	bool IncrementalTriangulation<Arithmetic>::legalize(t_pointer const & t, v_pointer const & v) {


		unsigned const edgeToFlip = t->get_vertex_index(v);


		/*****************************************************************************/
		/*                                                                           */
		/*		If the edge to flip is constrained, do not flip it					 */
		/*																			 */
		/*****************************************************************************/
		if (t->edges[edgeToFlip]->isConstrained)
			return false;


		/*****************************************************************************/
		/*                                                                           */
		/*		Get the opposite triangle											 */
		/*																			 */
		/*****************************************************************************/
		t_pointer const ot = t->neighbors[edgeToFlip];


		/*****************************************************************************/
		/*                                                                           */
		/*		When hole or super triangle, ot may be NULL							 */
		/*																			 */
		/*		This should never happen, because super triangle has constrained	 */
		/*		edges, so the function returns above. Holes are defined through		 */
		/*		constrained edges, so the function should return above also			 */
		/*																			 */
		/*****************************************************************************/
		if (!ot) return false;


		/*****************************************************************************/
		/*                                                                           */
		/*		Get the opposite vertex in the opposite triangle					 */
		/*																			 */
		/*****************************************************************************/
		v_pointer const ov = ot->get_opposite_vertex(t, v);


		/*****************************************************************************/
		/*                                                                           */
		/*		Check if the vertex is in the circum circle of the opposite 		 */
		/*      triangle                                                             */
		/*																			 */
		/*****************************************************************************/
		bool const inCircumCircle = Predicates.in_circumcircle(ot, v) > 0.0;


		if (inCircumCircle) {


			/*****************************************************************************/
			/*                                                                           */
			/*		Rotate the triangle pair (quadrilateral) in clockwise direction		 */
			/*																			 */
			/*****************************************************************************/
			rotate_triangle_pair(t, v, ot, ov);


			/*****************************************************************************/
			/*                                                                           */
			/*		Rotating triangles can make vertices's pointers (except v) of		 */
			/*		the formed quadrilateral point to different triangle, which			 */
			/*		does not contain respective vertices								 */
			/*																			 */
			/*****************************************************************************/
			ov->adjacentTriangle = ot;
			t->get_vertex_cw(v)->adjacentTriangle = t;
			ot->get_vertex_cw(ov)->adjacentTriangle = ot;


			/*****************************************************************************/
			/*                                                                           */
			/*		Check Delaunay property of the edges which are newely visible from   */
			/*		the inserted vertex v (as a consequence of the flipping)			 */
			/*																			 */
			/*****************************************************************************/
			legalize(t, v);
			legalize(ot, v);


			/*****************************************************************************/
			/*                                                                           */
			/*		Return the status that the legalization of the edge was accomplished */
			/*																			 */
			/*****************************************************************************/
			return true;

		}


		/*****************************************************************************/
		/*                                                                           */
		/*		Return the status that the legalization didn't take place			 */
		/*																			 */
		/*****************************************************************************/
		return false;

	};
	template<GeometricPredicatesArithmetic Arithmetic>
	void IncrementalTriangulation<Arithmetic>::rotate_triangle_pair(t_pointer const & t, v_pointer const & v, t_pointer const & ot, v_pointer const & ov) {



		/*****************************************************************************/
		/*                                                                           */
		/*		Save neighbors of triangles											 */
		/*																			 */
		/*****************************************************************************/
		t_pointer const n0 = t->get_neighbor_ccw(v);
		t_pointer const n1 = t->get_neighbor_cw(v);
		t_pointer const n2 = ot->get_neighbor_ccw(ov);
		t_pointer const n3 = ot->get_neighbor_cw(ov);


		/*****************************************************************************/
		/*                                                                           */
		/*		Save edges of triangles												 */
		/*																			 */
		/*****************************************************************************/
		e_pointer const edgeToFlip = t->edges[t->get_vertex_index(v)];

		e_pointer const e0 = t->get_edge_ccw(v);
		e_pointer const e1 = t->get_edge_cw(v);
		e_pointer const e2 = ot->get_edge_ccw(ov);
		e_pointer const e3 = ot->get_edge_cw(ov);


		/*****************************************************************************/
		/*                                                                           */
		/*		Flip edge by rotating neighboring triangles	clock-wise				 */
		/*																			 */
		/*****************************************************************************/
		t->rotate_triangle_cw(v, ov);
		ot->rotate_triangle_cw(ov, v);

		edgeToFlip->set_vertices(v, ov);


		/*****************************************************************************/
		/*                                                                           */
		/*		Assign edges to triangles											 */
		/*																			 */
		/*****************************************************************************/
		t->set_edge(edgeToFlip);
		t->set_edge(e2);
		t->set_edge(e1);

		ot->set_edge(edgeToFlip);
		ot->set_edge(e0);
		ot->set_edge(e3);


		/*****************************************************************************/
		/*                                                                           */
		/*		Assign neighboring triangles to edges								 */
		/*																			 */
		/*****************************************************************************/
		e0->set_neighbor(0) = ot;
		e0->set_neighbor(1) = n0;

		e1->set_neighbor(0) = t;
		e1->set_neighbor(1) = n1;

		e2->set_neighbor(0) = t;
		e2->set_neighbor(1) = n2;

		e3->set_neighbor(0) = ot;
		e3->set_neighbor(1) = n3;

		edgeToFlip->set_neighbor(0) = t;
		edgeToFlip->set_neighbor(1) = ot;


		/*****************************************************************************/
		/*                                                                           */
		/*		Assign neighborhood													 */
		/*																			 */
		/*****************************************************************************/
		t->null_neighbors();
		ot->null_neighbors();

		if (n0) ot->set_neighbor(n0);
		if (n1) t->set_neighbor(n1);
		if (n2) t->set_neighbor(n2);
		if (n3) ot->set_neighbor(n3);

		t->set_neighbor(ot);

	};


	/*****************************************************************************/
	/*                                                                           */
	/*		inserts constrained edge using intersections						 */
	/*																			 */
	/*****************************************************************************/
	template<GeometricPredicatesArithmetic Arithmetic>
	void IncrementalTriangulation<Arithmetic>::insert_constraint(v_pointer const & a, v_pointer const & b, MarkerEdge const & marker) {



		/*****************************************************************************/
		/*                                                                           */
		/*		Check if the algorithm did its job by comparing starting and		 */
		/*		ending vertices to be joined by edge segment						 */
		/*																			 */
		/*****************************************************************************/
		if (a == b) return;


		/*****************************************************************************/
		/*                                                                           */
		/*      Initialization step                                                  */
		/*                                                                           */
		/*		Circulate through triangles around the point a so that the vertex b	 */
		/*		lies strictly on the right of the line defined by vertices r and l	 */
		/*																			 */
		/*		We can use the fact, that vertex a is already in the triangulation	 */
		/*		so it is assigned a triangle, which contains vertex a				 */
		/*																			 */
		/*****************************************************************************/
		t_pointer t = a->adjacentTriangle;

		v_pointer right = t->get_vertex_ccw(a);
		v_pointer left  = t->get_vertex_cw(a);

		/*while (true) {

			double const orientationRight = Predicates.orientation(right, a, b);
			double const orientationLeft  = Predicates.orientation(left, a, b);

			if (orientationRight <= 0.0 && orientationLeft >= 0)
				break;

			right = left;
			t = t->get_neighbor_cw(a);
			left = t->get_vertex_cw(a);

		}*/


		if (Predicates.orientation(right, a, b) <= 0.0) {

			while (Predicates.orientation(left, a, b) <= 0.0) {

				right = left;
				t	  = t->get_neighbor_cw(a);
				left  = t->get_vertex_cw(a);

			}
		}
		else {

			do {

				left  = right;
				t	  = t->get_neighbor_ccw(a);
				right = t->get_vertex_ccw(a);

			} while (Predicates.orientation(right, a, b) > 0.0);
		}


		/*****************************************************************************/
		/*                                                                           */
		/*		Check if the edge is already in triangulation. Specifically			 */
		/*		and edge of the triangle t (vertex b is either r or l)				 */
		/*																			 */
		/*****************************************************************************/
		if (left == b || right == b) {


			e_pointer const e = t->get_edge(a, b);

			e->isConstrained = true;
			e->marker		 = marker;

			return;

		}




		/*****************************************************************************/
		/*                                                                           */
		/*		Constrained edge is situated between these two points				 */
		/*																			 */
		/*****************************************************************************/
		v_pointer const r = right;
		v_pointer const l = left;


		if (abs(Predicates.orientation(a, b, r)) <= EpsilonOnEdge) {


			/*****************************************************************************/
			/*                                                                           */
			/*		Denote r as constrained. This vertex must not be removed.			 */
			/*																			 */
			/*****************************************************************************/
			r->marker = MarkerVertex::Constrained;


			/*****************************************************************************/
			/*                                                                           */
			/*		Constrained edge is going through the vertex r therefore part		 */
			/*		of the constrained edge is already in the triangulation. We must	 */
			/*		set it constrained and set its boundary marker.						 */
			/*																			 */
			/*****************************************************************************/
			e_pointer const f = t->get_edge(l);

			f->marker		 = marker;
			f->isConstrained = true;


			/*****************************************************************************/
			/*                                                                           */
			/*		The remaining part of the constrained edge must start at vertex r.	 */
			/*		Start over the algorithm with new starting vertex  a = r			 */
			/*																			 */
			/*****************************************************************************/
			insert_constraint(r, b, marker);

			return;

		}
		else if (abs(Predicates.orientation(a, b, l)) <= EpsilonOnEdge) {


			/*****************************************************************************/
			/*                                                                           */
			/*		Denote l as constrained. This vertex must not be removed.			 */
			/*																			 */
			/*****************************************************************************/
			l->marker = MarkerVertex::Constrained;


			/*****************************************************************************/
			/*                                                                           */
			/*		Constrained edge is going through the vertex l therefore part		 */
			/*		of the constrained edge is already in the triangulation. We must	 */
			/*		set it constrained and set its boundary marker.						 */
			/*																			 */
			/*****************************************************************************/
			e_pointer const f = t->get_edge(r);

			f->marker		 = marker;
			f->isConstrained = true;


			/*****************************************************************************/
			/*                                                                           */
			/*		The remaining part of the constrained edge must start at vertex l.	 */
			/*		Start over the algorithm with new starting vertex  a = l			 */
			/*																			 */
			/*****************************************************************************/
			insert_constraint(l, b, marker);

			return;

		}




		/*****************************************************************************/
		/*                                                                           */
		/*		An edge between r and l vertices. That is edge opposite to vertex a  */
		/*																			 */
		/*****************************************************************************/
		e_pointer const e = t->get_edge(a);


		/*****************************************************************************/
		/*                                                                           */
		/*		Neighboring triangle to the triangle t through the edge e			 */
		/*																			 */
		/*****************************************************************************/
		t_pointer ot		= t->get_neighbor(a);
		v_pointer const oa  = ot->get_vertex(e);






		/*****************************************************************************/
		/*                                                                           */
		/*		Get coordinates of the edges intersection: constrained edge x (l,r)  */
		/*																			 */
		/*****************************************************************************/
		double xi;
		double yi;

		Predicates.get_edge_intersection(a, b, l, r, xi, yi);


		/*****************************************************************************/
		/*                                                                           */
		/*		Insert new vertex into the intersection point. ot is unspecified	 */
		/*		triangle containing intersectionVertex								 */
		/*																			 */
		/*****************************************************************************/
		v_pointer const intersectionVertex = new Vertex(xi, yi);

		vertices					.push_back(intersectionVertex);
		segmentIntersectionVertices	.push_back(intersectionVertex);

		ot = insert_vertex_into_edge(e, intersectionVertex);


		/*****************************************************************************/
		/*                                                                           */
		/*		Denote intersectionVertex as constrained. This vertex must not be 	 */
		/*		displaced or removed. Set the vertex global index					 */
		/*																			 */
		/*****************************************************************************/
		intersectionVertex->marker = MarkerVertex::Constrained;
		intersectionVertex->index  = static_cast<int>(vertices.size());


		/*****************************************************************************/
		/*                                                                           */
		/*		Circulate through triangles around intersectionVertex so that the 	 */
		/*		resulting triangle contains starting vertex a. Then we can denote	 */
		/*		new created edges (segments) as constrained							 */
		/*																			 */
		/*		These were created when inserting intersectionVertex into an edge e	 */
		/*																			 */
		/*****************************************************************************/
		while (!ot->contains(a))
			ot = ot->get_neighbor_ccw(intersectionVertex);

		/*
		while (true) {

			if (!ot->contains(a))	ot = ot->get_neighbor_ccw(intersectionVertex);
			else					break;

		}
		*/


		/*****************************************************************************/
		/*                                                                           */
		/*		From the resulting triangle ot we are able to denote subparts of 	 */
		/*		the constrained edge as constrained and set their boundary marker	 */
		/*																			 */
		/*****************************************************************************/
		e_pointer f = ot->get_edge(intersectionVertex, a);

		f->isConstrained = true;
		f->marker = marker;



		/*****************************************************************************/
		/*                                                                           */
		/*		The vertex oa in the triangle ot is on the constrained edge. This	 */
		/*		vertex must be already in pslg and thus it is constrained			 */
		/*																			 */
		/*****************************************************************************/
		if (abs(Predicates.orientation(a, b, oa)) <= EpsilonOnEdge) {



			/*****************************************************************************/
			/*                                                                           */
			/*		Circulate through triangles around intersectionVertex so that the 	 */
			/*		resulting triangle contains vertex oa. Then we can denote new		 */
			/*		created edges (segments) as constrained								 */
			/*																			 */
			/*		These were created when inserting intersectionVertex into an edge e	 */
			/*																			 */
			/*****************************************************************************/
			while (!ot->contains(oa))
				ot = ot->get_neighbor_ccw(intersectionVertex);


			/*****************************************************************************/
			/*                                                                           */
			/*		From the resulting triangle ot we are able to denote subparts of 	 */
			/*		the constrained edge as constrained and set their boundary marker	 */
			/*																			 */
			/*****************************************************************************/
			f = ot->get_edge(intersectionVertex, oa);

			f->isConstrained = true;
			f->marker		 = marker;


			/*****************************************************************************/
			/*                                                                           */
			/*		Constrained edge is passing through the opposite vertex oa. Vertex 	 */
			/*		can be ending point b. If it is not further subdivision is needed	 */
			/*																			 */
			/*****************************************************************************/
			insert_constraint(oa, b, marker);

			return;

		}


		/*****************************************************************************/
		/*                                                                           */
		/*		This should be typical situation when the constrained edge doesn't   */
		/*		doesn't pass through a vertex except the intersectionVertex			 */
		/*                                                                           */
		/*****************************************************************************/
		insert_constraint(intersectionVertex, b, marker);

	};


	/*****************************************************************************/
	/*                                                                           */
	/*    Mark triangles condemned to destruction !n!, 666					     */
	/*                                                                           */
	/*****************************************************************************/
	template<GeometricPredicatesArithmetic Arithmetic>
	void IncrementalTriangulation<Arithmetic>::mark_triangle(t_pointer const & t, int const previousIndex) {


		t->marker = MarkerTriangle::Outside;

		for (unsigned i = 0; i < 3; i++) {


			/*****************************************************************************/
			/*                                                                           */
			/*    Came here from this triangle. It is already marked so skip it		     */
			/*                                                                           */
			/*****************************************************************************/
			if (i == previousIndex) continue;


			/*****************************************************************************/
			/*                                                                           */
			/*    This is new neighbor triangle where I now go spread virus				 */
			/*                                                                           */
			/*****************************************************************************/
			t_pointer const neighbor = t->neighbors[i];


			/*****************************************************************************/
			/*                                                                           */
			/*    If we are on the edge of super triangle, there are no neighbors		 */
			/*                                                                           */
			/*****************************************************************************/
			if (!neighbor) continue;


			/*****************************************************************************/
			/*                                                                           */
			/*    Check if the neighbor is already marked or the edge is constrained     */
			/*    and in that case we must not cross this edge                           */
			/*                                                                           */
			/*****************************************************************************/
			if (neighbor->marker == MarkerTriangle::Outside || t->edges[i]->isConstrained)
				continue;


			/*****************************************************************************/
			/*                                                                           */
			/*    Recurrently enter into another iteration								 */
			/*                                                                           */
			/*****************************************************************************/
			mark_triangle(neighbor, neighbor->get_neighbor_index(t));

		}


		/*****************************************************************************/
		/*                                                                           */
		/*     Insert this triangle into container of marked triangles				 */
		/*                                                                           */
		/*****************************************************************************/
		trianglesOutside.push_back(t);

	};
	template<GeometricPredicatesArithmetic Arithmetic>
	void IncrementalTriangulation<Arithmetic>::remove_marked_triangles() {


		std::vector<v_pointer> verticesToRemove;


		for (size_t k = 0; k < trianglesOutside.size(); k++) {


			t_pointer const t = trianglesOutside[k];

			/*****************************************************************************/
			/*                                                                           */
			/*		Delete edges which are not constrained inside hole				     */
			/*                                                                           */
			/*****************************************************************************/
			for (unsigned i = 0; i < 3; i++) {


				e_pointer const e = t->edges[i];

				/*****************************************************************************/
				/*                                                                           */
				/*		If edge is already deleted from other iteration then next edge	     */
				/*                                                                           */
				/*****************************************************************************/
				if (!e) continue;

				t_pointer const neighbor = t->neighbors[i];

				v_pointer const a = e->vertices[0];
				v_pointer const b = e->vertices[1];


				if (!e->isConstrained) {


					/*****************************************************************************/
					/*                                                                           */
					/*		The edge e is being destroyed. We need to cancel all pointers and	 */
					/*      references from neighboring triangles to this edge. It's neccessary  */
					/*      because in some point it can be tested for NULL (above). Otherwise   */
					/*      there could be some pointer pointing to rubbish (program can crash)  */
					/*                                                                           */
					/*****************************************************************************/
					neighbor->get_edge(e) = NULL;

					remove(edges, e);
					delete e;


					/*****************************************************************************/
					/*                                                                           */
					/*		Vertex inside a hole is pushed into	verticesToRemove container so    */
					/*      that it is uniquely present in the container                         */
					/*                                                                           */
					/*****************************************************************************/
					bool ain = false;
					bool bin = false;

					for (size_t j = 0; j < verticesToRemove.size(); j++) {

						if (a->marker != MarkerVertex::Constrained)
							if (a == verticesToRemove[j]) ain = true;
						if (b->marker != MarkerVertex::Constrained)
							if (b == verticesToRemove[j]) bin = true;

						if (ain && bin) break;

					}

					if (a->marker != MarkerVertex::Constrained)
						if (!ain) verticesToRemove.push_back(a);
					if (b->marker != MarkerVertex::Constrained)
						if (!bin) verticesToRemove.push_back(b);

				}
				else if (neighbor) {

					neighbor->set_neighbor(e) = NULL;

					e->neighbors[0] = neighbor;
					e->neighbors[1] = NULL;

					a->adjacentTriangle = neighbor;
					b->adjacentTriangle = neighbor;

				}

			}

			remove(triangles, t);
			delete t;

		}


		/*****************************************************************************/
		/*                                                                           */
		/*		Test on some examples what is faster							     */
		/*                                                                           */
		/*****************************************************************************/
		//std::sort(verticesToRemove.begin(), verticesToRemove.end());
		//std::unique(verticesToRemove.begin(), verticesToRemove.end());


		/*****************************************************************************/
		/*                                                                           */
		/*		Delete vertices inside hole										     */
		/*                                                                           */
		/*****************************************************************************/
		for (size_t i = 0; i < verticesToRemove.size(); i++) {

			v_pointer const v = verticesToRemove[i];

			remove(vertices, v);
			delete v;

		}


		/*****************************************************************************/
		/*                                                                           */
		/*		Re-index vertices because some of theme were deleted and therefore   */
		/*      there would be gaps between numbering (some indeces wouldn't exist)  */
		/*                                                                           */
		/*****************************************************************************/
		if (!verticesToRemove.empty())
			for (size_t i = 0; i < vertices.size(); i++)
				vertices[i]->index = i;


		verticesToRemove.clear();
		trianglesOutside.clear();

		verticesToRemove.shrink_to_fit();
		trianglesOutside.shrink_to_fit();

	};


	/*****************************************************************************/
	/*                                                                           */
	/*	remove		Removes triangle / edge / vertex from the container			 */
	/*				   															 */
	/*  clear		Clear all data												 */
	/*                                                                           */
	/*****************************************************************************/
	template<GeometricPredicatesArithmetic Arithmetic> template<typename T>
	void IncrementalTriangulation<Arithmetic>::remove(std::vector<T> & vec, T element) {


		/*****************************************************************************/
		/*                                                                           */
		/*		It depends whether it is faster to lookup from the beginning or		 */
		/*		to lookup from the back for the element t							 */
		/*                                                                           */
		/*****************************************************************************/
		for (auto it = vec.rbegin(); it != vec.rend(); it++) {


			if (*it == element) {

				vec.erase(it.base() - 1);
				break;

			}
		}

	};
	template<GeometricPredicatesArithmetic Arithmetic> template<typename T>
	void IncrementalTriangulation<Arithmetic>::replace(std::vector<T> & vec, T element) {

		vec[element->index] = element;

	};
	template<GeometricPredicatesArithmetic Arithmetic>
	void IncrementalTriangulation<Arithmetic>::clear() {


		/*****************************************************************************/
		/*                                                                           */
		/*		Free all allocated memory											 */
		/*																			 */
		/*****************************************************************************/
		for (size_t i = 0; i < vertices.size(); i++)
			delete vertices[i];

		for (size_t i = 0; i < edges.size(); i++)
			delete edges[i];

		for (size_t i = 0; i < triangles.size(); i++)
			delete triangles[i];


		/*****************************************************************************/
		/*                                                                           */
		/*		Clear the vectors of the rubbish and free memory 					 */
		/*																			 */
		/*****************************************************************************/
		vertices	.clear();
		edges		.clear();
		triangles	.clear();

		vertices	.shrink_to_fit();
		edges		.shrink_to_fit();
		triangles	.shrink_to_fit();


		/*****************************************************************************/
		/*                                                                           */
		/*		These containers contain already deleted primitives. No need to 	 */
		/*		delete theme again (it would only lead to a program failure)		 */
		/*																			 */
		/*****************************************************************************/
		trianglesOutside			.clear();
		segmentIntersectionVertices	.clear();

		trianglesOutside			.shrink_to_fit();
		segmentIntersectionVertices	.shrink_to_fit();

	};



};
