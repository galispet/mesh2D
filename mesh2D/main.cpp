

/*****************************************************************************/
/*                                                                           */
/*	- To Consider															 */
/*                                                                           */
/*****************************************************************************/
// -Precomputing determinants can be made with additional data struture, that will be deleted after triangulation = > no additional memory is consumed
// - Implement 'infinite bounding triangle'
// - My idiotic solution : if the input vertice is outside bounding triangle->make the triangle inflate.We have to check coordinates of each inserted point : (
// - Verify my triangulation with matlab's -> generate same points and compare 
// - Barycentric refinment(subdivision)
// - Try to normalize given PSLG onto[0, 1] x[0, 1] for better accuracy.new_x = (old_x - min_x) / d_max, new_y = (old_y - min_y) / d_max, d_max = max{ x_max - x_min,y_max - y_min }
// - make laplacian smoothing similiar to the book of 'Jiri Blazek - Computational Fluid Dynamics_ Principles and Applications (2015, Butterworth-Heinemann)' p. 394 / 388
// - Barycentric smoothing, cotangents wights smoothing


/*****************************************************************************/
/*                                                                           */
/*	- To Check																 */
/*                                                                           */
/*****************************************************************************/
// - In GeometricKernel::GeometricPredicates in_circle_robust ->at the bottom there is conditioi for zero determinant (in shewchuck, there is not) Check wether both with it and without gives the same result
// - Function Two_One_Diff There is really Two_Sum. Check this in Shewchuck paper, if there is really this one;
// - Be careful in ShewChuk's Fast expasion sum etc. There are static keywords. Better watch out so I don't break something
// - Discarded the check for CCW position of vertices from constructing triangle. Hope it will not break
// - In insert_contraint chcek the validity of initilizing step and compare those two algorithm I have there


/*****************************************************************************/
/*                                                                           */
/*	- To do																	 */
/*                                                                           */
/*****************************************************************************/
// - Make all the methods in trinagle, vertex, edge pass the argument as reference, so no copies are performed !!!!! It should be much faster
// - Get rid of  erase vector method when removing triangles,edges,vertice. It so so slow. Every time, it shift all the elements resp. reallocate all the elements !!! Use my method replace instead
// - Run tests for orientaion and incircle so that the adapt is uses. because there are variables which may/should be declared as volatile, because of some 
//   optimization issues. INEXACT == volatile in Shewchuk triangle. those varibale are u3-u[3], B3-B[3] etc.
// - Implement on_edge_adapt
// - Try to make the conatiner with primitives as lists and do some benchmarks in comparison with vector
// - Some additional data structure in geometric primitives (such as pointer to the neighboring triangle in the case of the vertex) maybe put in separte container in the triangulation
// - function 'refine_mesh(max_number_of_steiner_points = finite/as it needs)' -> parameter max_number_of_steiner_points input in constructor, or later? If later there can be no estimates for number of triangles
// - Walking algorithm starts not from the back, but from random set of vertices. The vertices adjacent trinagle whose veretx is the closest to the new point is selected
// - Do not create edges immedietelly but as the last step. It will simplify everything (should). And should increae speed


#include "GeometricKernel.h"
#include "Triangulation.h"

//#include "mesh.h"
#include "gnuplot.h"


#include <iostream>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <random>
#include <vector>
#include <iomanip>


enum Generation { Random, Structured };

template<Generation gen>
void generate_vertices(std::vector<GeometricKernel::Vertex > & v);




const double a1 = 0.0;
const double b1 = 10.0;

const double a2 = 0.0;
const double b2 = 10.0;


const int N1 = 70;
const int N2 = 80;

const int N_random = 10000;


using namespace DelaunayTriangulation;


int main() {


	PlanarStraightLineGraph pslg;

	std::vector<Vertex> vertices;
	std::vector<Vertex> seeds;


	/*****************************************************************************/
	/*                                                                           */
	/*		Generate and insert vertices										 */
	/*																			 */
	/*****************************************************************************/
	//generate_vertices<Random>(vertices);
	generate_vertices<Structured>(vertices);

	for (size_t i = 0; i < vertices.size(); i++)
		pslg.insert_vertex<MarkerVertex::Constrained>(vertices[i]);


	/*****************************************************************************/
	/*                                                                           */
	/*		Insert constraints													 */
	/*																			 */
	/*****************************************************************************/
	v_pointer const va = pslg.get_vertex_pointer(0);
	v_pointer const vb = pslg.get_vertex_pointer(N2 - 1);
	v_pointer const vc = pslg.get_vertex_pointer(N1*N2 - 1);
	v_pointer const vd = pslg.get_vertex_pointer(N1*N2 - 1 - (N2 - 1));
	v_pointer const ve = pslg.get_vertex_pointer(0);

	pslg.insert_constraint<MarkerEdge::Dirichlet>(va, vb);
	pslg.insert_constraint<MarkerEdge::Dirichlet>(vb, vc);
	pslg.insert_constraint<MarkerEdge::Dirichlet>(vc, vd);
	pslg.insert_constraint<MarkerEdge::Dirichlet>(vd, ve);

	unsigned const n1 = 100;
	double const r1 = 2.5;

	for (unsigned i = 0; i < n1; i++) {

		double const u1 = i * 2.0 * Pi / n1;
		double const u2 = (i + 1) * 2.0 * Pi / n1;

		v_pointer const a = pslg.insert_vertex<MarkerVertex::Free>(Vertex(r1*sin(u1) + 7.0, r1*cos(u1) + 5.0));
		v_pointer const b = pslg.insert_vertex<MarkerVertex::Free>(Vertex(r1*sin(u2) + 7.0, r1*cos(u2) + 5.0));

		pslg.insert_constraint<MarkerEdge::Dirichlet>(a, b);

	}

	seeds.push_back(Vertex(5.0, 5.0));



	unsigned const n2 = 100;
	double const r2 = 1.5;

	for (unsigned i = 0; i < n2; i++) {

		double const u1 = i * 2.0 * Pi / n1;
		double const u2 = (i + 1) * 2.0 * Pi / n1;

		v_pointer const a = pslg.insert_vertex<MarkerVertex::Free>(Vertex(r2*sin(u1) + 2.0, r2*cos(u1) + 5.0));
		v_pointer const b = pslg.insert_vertex<MarkerVertex::Free>(Vertex(r2*sin(u2) + 2.0, r2*cos(u2) + 5.0));

		pslg.insert_constraint<MarkerEdge::Dirichlet>(a, b);

	}

	seeds.push_back(Vertex(2.0, 5.0));


	/*****************************************************************************/
	/*                                                                           */
	/*		Compute Constrained triangulation									 */
	/*																			 */
	/*****************************************************************************/
	clock_t const begin = clock();
	IncrementalTriangulation<GeometricPredicatesArithmetic::Robust> cdt(pslg, seeds);
	clock_t const end = clock();

	std::cout << std::endl << "CPU clocks : " << (end - begin) << std::endl << std::endl;


	/*****************************************************************************/
	/*                                                                           */
	/*		Export triangulation primitives										 */
	/*																			 */
	/*****************************************************************************/
	std::string verticesFile	= "C:\\Users\\pgali\\Desktop\\verticesPoints.txt";
	std::string edgesFile		= "C:\\Users\\pgali\\Desktop\\edgesPoints.txt";
	std::string trianglesFile	= "C:\\Users\\pgali\\Desktop\\trianglesPoints.txt";

	std::cout << std::endl << "Exporting files ...\n" << std::endl;

	cdt.export_vertices(verticesFile);
	cdt.export_edges(edgesFile);
	cdt.export_triangles(trianglesFile);


	//GNUPLOT gnuplot;
	//gnuplot("plot 'C:\\Users\\pgali\\Desktop\\trianglesPoints.txt' with lines");
	////gnuplot("plot 'C:\\Users\\pgali\\Desktop\\triangles_txt_files\\trianglesPoints.txt' with lines");



	system("pause");

	return 0; 


	/*

	//v_pointer va = pslg.insert_vertex<V_MARKER::CONSTRAINED>(Vertex(a1, a2));
	//v_pointer vb = pslg.insert_vertex<V_MARKER::CONSTRAINED>(Vertex(a1, b2));
	//v_pointer vc = pslg.insert_vertex<V_MARKER::CONSTRAINED>(Vertex(b1, b2));
	//v_pointer vd = pslg.insert_vertex<V_MARKER::CONSTRAINED>(Vertex(b1, a2));

	//v_pointer vaa = pslg.insert_vertex<V_MARKER::CONSTRAINED>(Vertex(-5.0, -20.0));
	//v_pointer vbb = pslg.insert_vertex<V_MARKER::CONSTRAINED>(Vertex(-5.0, -15.0));
	//v_pointer vcc = pslg.insert_vertex<V_MARKER::CONSTRAINED>(Vertex(5.0, -15.0));
	//v_pointer vdd = pslg.insert_vertex<V_MARKER::CONSTRAINED>(Vertex(5.0, -20.0));


	//v_pointer hhhhh = pslg.insert_vertex<V_MARKER::CONSTRAINED>(Vertex(-5.0, -11.0));

	//pslg.insert_constraint<E_MARKER::DIRICHLET>(hhhhh, vb);

	//v_pointer va = pslg.get_vertex(0);
	//v_pointer vb = pslg.get_vertex(N2 - 1);
	//v_pointer vc = pslg.get_vertex(N1*N2 - 1);
	//v_pointer vd = pslg.get_vertex(N1*N2 - 1 - (N2 - 1));

	//pslg.insert_constraint<E_MARKER::DIRICHLET>(va, vb);
	//pslg.insert_constraint<E_MARKER::NEUMANN>(vb, vc);
	//pslg.insert_constraint<E_MARKER::NEUMANN>(vc, vd);
	//pslg.insert_constraint<E_MARKER::NEUMANN>(vd, va);

	//pslg.insert_constraint<E_MARKER::NEUMANN>(vaa, vbb);
	//pslg.insert_constraint<E_MARKER::NEUMANN>(vbb, vcc);
	//pslg.insert_constraint<E_MARKER::NEUMANN>(vcc, vdd);
	//pslg.insert_constraint<E_MARKER::NEUMANN>(vaa, vdd);
	

	// These are for creating small input angles and thus it serves to implement concentric shells + rules of Miller, Pav, and Walkington + off-center circumcenter
	//pslg.insert_constraint<E_MARKER::NEUMANN>(va, vaa);
	//pslg.insert_constraint<E_MARKER::DIRICHLET>(vd, vdd);

	*/


	//std::vector<Vertex > hole1;
	//std::vector<Vertex > hole2;
	//std::vector<Vertex > hole3;
	//std::vector<Vertex > hole4;
	//std::vector<Vertex > hole5;
	//
	//
	//unsigned n1 = 5;
	//unsigned n2 = 30;
	//unsigned n3 = 30;
	//unsigned n4 = 40;
	//unsigned n5 = 30;
	//
	//double r1 = 3.0;
	//double r2 = 2.5;
	//double r3 = 4.0;
	//double r4 = 5.0;
	//double r5 = 2.5;
	//
	//double x1 = -6.0;
	//double y1 = -5.0;
	//
	//double x2 = 0.0;
	//double y2 = 0.0;
	//
	//double x3 = +7.0;
	//double y3 = +7.0;
	//
	//double x4 = 9.0;
	//double y4 = -8.0;
	//
	//double x5 = 9.0;
	//double y5 = -8.0;
	//
	//Vertex * const seed1 = new Vertex(x1, y1);
	//Vertex * const seed2 = new Vertex(x2, y2);
	//Vertex * const seed3 = new Vertex(x3, y3);
	//Vertex * const seed4 = new Vertex(10.0, -3.6);
	////Vertex * const seed4 = new Vertex(x4, y4);
	//Vertex * const seed5 = new Vertex(x5, y5);

	/*
	//for (unsigned j = 0; j < n1; j++) {

	//	double const theta = j * 2 * Pi / (n1 - 0);

	//	hole1.push_back(Vertex(x1 + r1 * sin(theta), y1 + r1 * cos(theta)));

	//}

	//for (unsigned j = 0; j < n2; j++) {

	//	double const theta = j * 2 * Pi / (n2 - 0);

	//	hole2.push_back(Vertex(x2 + r2 * sin(theta), y2 + r2 * cos(theta)));;

	//}

	//for (unsigned j = 0; j < n3; j++) {

	//	double const theta = j * 2 * Pi / (n3 + 0);

	//	hole3.push_back(Vertex(x3 + r3 * sin(theta), y3 + (r3 - 3.0) * cos(theta)));

	//}

	//for (unsigned j = 0; j < n4; j++) {

	//	double const theta = j * 2 * Pi / (n4 - 0);

	//	hole4.push_back(Vertex(x4 + r4 * sin(theta), y4 + r4 * cos(theta)));

	//}

	////hole4.push_back(Vertex(0.0, 20.0));

	//for (unsigned j = 0; j < n5; j++) {

	//	double const theta = j * 2 * Pi / (n5 - 0);

	//	hole5.push_back(Vertex(x5 + r5 * sin(theta), y5 + r5 * cos(theta)));

	//}



	//------- AIRFOIL -------/
	//
	//// airfoil : 200, airfoil3 : 160

	//std::ifstream inFile;
	////inFile.open("C:\\Users\\pgali\\Desktop\\verticesPoints.txt");
	//inFile.open("C:\\Users\\pgali\\Desktop\\airfoil3.txt");
	//
	//std::vector<Vertex> airfoil;

	//double x, y;
	//char delim;

	//for (unsigned i = 0; i < 94; i++) {

	//	inFile >> x >> delim >> y;

	//	Vertex v(20*x - 10.0, 60*y);

	//	airfoil.push_back(v);

	//}

	//inFile.close();
	//
*/


	//tri.insert_vertex(vertices.begin(), vertices.end());

	//Vertex_handle va = tri.get_vertex(0);
	//Vertex_handle vb = tri.get_vertex(N2 - 1);
	//Vertex_handle vc = tri.get_vertex(N1*N2 - 1);
	//Vertex_handle vd = tri.get_vertex(N1*N2 - 1 - (N2 - 1));
	//Vertex_handle ve = tri.get_vertex(0);
		
	//tri.insert_individual_constraint<INSERT_CONSTRAINT::INTERSECTION, E_MARKER::NEUMANN>(va, vb);
	//tri.insert_individual_constraint<INSERT_CONSTRAINT::INTERSECTION, E_MARKER::DIRICHLET>(vb, vc);
	//tri.insert_individual_constraint<INSERT_CONSTRAINT::INTERSECTION, E_MARKER::DIRICHLET>(vc, vd);
	//tri.insert_individual_constraint<INSERT_CONSTRAINT::INTERSECTION, E_MARKER::NEUMANN>(vd, va);


	//tri.insert_individual_constraint(vertices[0], vertices[N1*N2 - 1]);
	//tri.insert_individual_constraint(vertices[N2 - 1], vertices[N1*N2 - 1 - (N2 - 1)]);
	

	//tri.insert_sequential_constraint<INSERT_CONSTRAINT::INTERSECTION, POLYGON_TYPE::CLOSED, E_MARKER::NEUMANN>(airfoil.begin(), airfoil.end());

	//tri.insert_sequential_constraint<INSERT_CONSTRAINT::INTERSECTION, POLYGON_TYPE::CLOSED, E_MARKER::NEUMANN>(hole1.begin(), hole1.end());
	//tri.insert_sequential_constraint<INSERT_CONSTRAINT::INTERSECTION, POLYGON_TYPE::CLOSED, E_MARKER::NEUMANN>(hole2.begin(), hole2.end());
	//tri.insert_sequential_constraint<INSERT_CONSTRAINT::INTERSECTION, POLYGON_TYPE::CLOSED, E_MARKER::NEUMANN>(hole3.begin(), hole3.end());
	//tri.insert_sequential_constraint<INSERT_CONSTRAINT::INTERSECTION, POLYGON_TYPE::CLOSED, E_MARKER::NEUMANN>(hole4.begin(), hole4.end());
	//tri.insert_sequential_constraint<INSERT_CONSTRAINT::INTERSECTION, POLYGON_TYPE::CLOSED, E_MARKER::NEUMANN>(hole5.begin(), hole5.end());


	//std::vector<Vertex *> additional_vertices = tri.get_new_vertices_handles();

	//tri.legalize_vertices(additional_vertices.begin(), additional_vertices.end());
	//tri.legalize_vertices();

	//cout << "smoothing" << endl;
	//tri.apply_smoothing<SMOOTHING::LAPLACIAN>(40);

	//double const angle = 32.0; // Degrees

	//tri.refinement_ruppert(angle / 180 * Pi);
	//tri.apply_smoothing<SMOOTHING::LAPLACIAN>(40);

};


template<Generation gen>
void generate_vertices(std::vector<GeometricKernel::Vertex > & v) {


	if (gen == Generation::Random) {


		std::random_device					rd;				// obtain a random number from hardware
		std::mt19937						eng(rd());		// seed the generator - Merssene twister
		std::uniform_real_distribution<>	distrX(a1, b1);	// define the range
		std::uniform_real_distribution<>	distrY(a2, b2);	// define the range

		for (int i = 0; i < N_random; i++) {

			double const x = distrX(eng);
			double const y = distrY(eng);

			v.push_back(Vertex(x, y));

		}
	}
	else if (gen == Generation::Structured) {


		for (int i = 0; i < N1; i++) {

			double const x = a1 + i * (b1 - a1) / (N1 - 1);

			for (int j = 0; j < N2; j++) {

				double const y = a2 + j * (b2 - a2) / (N2 - 1);

				v.push_back(Vertex(x, y));

			}
		}

		
		//for (int i = 0; i < N1; i++) {
		//
		//	
		//	double r = 0.2 + i * 10.0 / (N1 - 1);
		//
		//	for (int j = 0; j < N2; j++) {
		//
		//		double theta = j*2*3.14 / (N2);
		//
		//		v.push_back(new Vertex(r*sin(theta), r*cos(theta)));
		//
		//	}
		//}
		

	}

};
