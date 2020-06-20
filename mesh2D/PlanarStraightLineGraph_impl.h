#pragma once



namespace DelaunayTriangulation {


	PlanarStraightLineGraph::PlanarStraightLineGraph() {

	};
	PlanarStraightLineGraph::~PlanarStraightLineGraph() {

		clear();

	};

	template<MarkerVertex vMarker>
	v_pointer PlanarStraightLineGraph::insert_vertex(Vertex const & v) {


		size_t const nv = vertices.size();

		/*****************************************************************************/
		/*                                                                           */
		/*		Check for duplicate vertex											 */
		/*																			 */
		/*****************************************************************************/
		for (size_t i = 0; i < nv; i++)
			if (v.is_almost_equal(vertices[i]))
				return vertices[i];

		
		/*****************************************************************************/
		/*                                                                           */
		/*		Allocate memory in the pslg instance for the vertex					 */
		/*																			 */
		/*****************************************************************************/
		v_pointer const newVertex = new Vertex(v.x, v.y);


		/*****************************************************************************/
		/*                                                                           */
		/*		Get the maximal dimension of the pslg								 */
		/*																			 */
		/*****************************************************************************/
		if		(v.x > xMaximum) xMaximum = v.x;
		else if (v.x < xMinimum) xMinimum = v.x;
		if		(v.y > yMaximum) yMaximum = v.y;
		else if (v.y < yMinimum) yMinimum = v.y;


		/*****************************************************************************/
		/*                                                                           */
		/*		Set the index and marker of the vertex								 */
		/*																			 */
		/*****************************************************************************/
		newVertex->index  = (int) vertices.size();
		newVertex->marker = vMarker;


		/*****************************************************************************/
		/*                                                                           */
		/*		Insert vertex into the container									 */
		/*																			 */
		/*****************************************************************************/
		vertices.push_back(newVertex);


		/*****************************************************************************/
		/*                                                                           */
		/*		Return pointer to the vertex so the User can work with it (create	 */
		/*		constraints etc.)													 */
		/*																			 */
		/*****************************************************************************/
		return newVertex;

	};
	template<MarkerEdge eMarker>
	void PlanarStraightLineGraph::insert_constraint(v_pointer const & a, v_pointer const & b) {


		unsigned const i = a->index;
		unsigned const j = b->index;

		if (i == j) 
			throw "input vertices are duplicates";

		unsigned * const constraint = new unsigned[2];

		constraint[0] = i;
		constraint[1] = j;

		segments		.push_back(constraint);
		segmentsMarker	.push_back(eMarker);

	};
	template<MarkerEdge eMarker>
	void PlanarStraightLineGraph::insert_constraint(unsigned const i, unsigned const j) {


		if (i == j || i >= vertices.size() || j >= vertices.size())
			throw "invalid input (constrained edge)";

		unsigned * const constraint = new unsigned[2];

		constraint[0] = i;
		constraint[1] = j;

		segments		.push_back(constraint);
		segmentsMarker	.push_back(eMarker);

	};


	v_pointer PlanarStraightLineGraph::get_vertex_pointer(unsigned const i) {

		return vertices[i];

	};


	int PlanarStraightLineGraph::get_number_of_vertices() const {

		return (int)vertices.size();

	};
	int PlanarStraightLineGraph::get_number_of_segments() const {

		return (int)segments.size();

	};

	void PlanarStraightLineGraph::clear() {


		size_t const nv = vertices.size();
		size_t const ns = segments.size();

		for (size_t i = 0; i < nv; i++)
			delete vertices[i];

		for (size_t i = 0; i < ns; i++)
			delete segments[i];


		vertices		.clear();
		segments		.clear();
		segmentsMarker	.clear();

		vertices		.shrink_to_fit();
		segments		.shrink_to_fit();
		segmentsMarker	.shrink_to_fit();

	};

}