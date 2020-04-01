

#include "shapes.h"

#include "PSLG.h"
#include "triangulation.h"

#include "mesh.h"
#include "gnuplot.h"


#include <iostream>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <random>

#include <vector>

#include <iomanip>



GEOMETRIC_KERNEL const GK = GEOMETRIC_KERNEL::INEXACT_PREDICATES_INEXACT_CONSTRUCT;
CHECK_KERNEL const CK = CHECK_KERNEL::CHECK_DUPLICATE;




enum GENERATION { RANDOM, STRUCTURED };

template<GENERATION gen>
void generate_vertices(std::vector<Vertex > & v);




//const double a1 = -0.2;
//const double b1 = 1.2;
//
//const double a2 = -0.05;
//const double b2 = 0.2;

const double a1 = -25.0;
const double b1 = 25.0;

const double a2 = -25.0;
const double b2 = 25.0;

const int N1 = 40;// 87;
const int N2 = 40;// 157;

const int N_random = 200;




int main() {


	PlanarStraightLineGraph pslg;



	/*

	std::vector<Vertex> vertices;
	generate_vertices<RANDOM>(vertices);

	for (size_t i = 0; i < vertices.size(); i++)
		pslg.insert_vertex<V_MARKER::FREE>(vertices[i]);


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

	//unsigned const n1 = 100;

	//double r1 = 11.0;

	//for (unsigned i = 0; i < n1; i++) {

	//	double u1 = i * 2.0 * Pi / n1;
	//	double u2 = (i + 1) * 2.0 * Pi / n1;


	//	v_pointer a = pslg.insert_vertex<V_MARKER::FREE>(Vertex(r1*sin(u1), 1.0*r1*cos(u1)));
	//	v_pointer b = pslg.insert_vertex<V_MARKER::FREE>(Vertex(r1*sin(u2), 1.0*r1*cos(u2)));

	//	pslg.insert_constraint<E_MARKER::DIRICHLET>(a, b);

	//}

	//unsigned const n = 50;

	//double r = 6.0;

	//for (unsigned i = 0; i < n; i++) {

	//	double u1 = i * 2.0 * Pi / n;
	//	double u2 = (i + 1) * 2.0 * Pi / n;

	//	double r = 6.0;


	//	v_pointer a = pslg.insert_vertex<V_MARKER::FREE>(Vertex(r*sin(u1), 1.5*r*cos(u1)));
	//	v_pointer b = pslg.insert_vertex<V_MARKER::FREE>(Vertex(r*sin(u2), 1.5*r*cos(u2)));

	//	pslg.insert_constraint<E_MARKER::DIRICHLET>(a, b);

	//}



	//v_pointer ve = pslg.insert_vertex<V_MARKER::FREE>(Vertex(a1 + 5.0, a2 + 5.0));
	//v_pointer vf = pslg.insert_vertex<V_MARKER::FREE>(Vertex(a1 + 5.0, b2 - 5.0));
	//v_pointer vg = pslg.insert_vertex<V_MARKER::FREE>(Vertex(b1 - 5.0, b2 - 5.0));
	//v_pointer vh = pslg.insert_vertex<V_MARKER::FREE>(Vertex(b1 - 5.0, a1 + 5.0));

	//pslg.insert_constraint<E_MARKER::NEUMANN>(ve, vf);
	//pslg.insert_constraint<E_MARKER::NEUMANN>(vf, vg);
	//pslg.insert_constraint<E_MARKER::NEUMANN>(vg, vh);
	//pslg.insert_constraint<E_MARKER::NEUMANN>(vh, ve);

	//pslg.insert_constraint<E_MARKER::NEUMANN>(va, vg);


	//pslg.insert_constraint<E_MARKER::DIRICHLET>(va, vb);
	//pslg.insert_constraint<E_MARKER::DIRICHLET>(vb, vc);
	//pslg.insert_constraint<E_MARKER::DIRICHLET>(vc, vd);
	//pslg.insert_constraint<E_MARKER::DIRICHLET>(vd, va);


	std::vector<Vertex> seeds;

	//seeds.push_back(Vertex(8.5, 0.1));
	//seeds.push_back(Vertex(0.0, -16.0));



	clock_t begin = clock();

	Triangulation<GK> tri(pslg, seeds);

	//tri.refinement_ruppert<REFINEMENT_PRIORITY::WORST>(32.0, 0.0005);


	clock_t end = clock();

	Mesh mesh(tri);


	cout << "*************** Triangulation ***************" << endl;
	cout << tri.get_number_of_vertices() << endl;
	cout << tri.get_number_of_edges() << endl;
	cout << tri.get_number_of_triangles() << endl;
	cout << "*********************************************" << endl;


	std::ofstream vertices_txt;
	std::ofstream edges_txt;
	std::ofstream triangles_txt;


	std::cout << std::endl << "Creating files" << std::endl;

	std::cout << "Vertices..." << std::endl;
	vertices_txt.open("C:\\Users\\pgali\\Desktop\\verticesPoints.txt");
	mesh.export_vertices(vertices_txt);
	vertices_txt.close();

	std::cout << "Triangles..." << std::endl;
	triangles_txt.open("C:\\Users\\pgali\\Desktop\\trianglesPoints.txt");
	mesh.export_triangles(triangles_txt);
	triangles_txt.close();

	std::cout << "Edges..." << std::endl << std::endl;
	edges_txt.open("C:\\Users\\pgali\\Desktop\\edgesPoints.txt");
	mesh.export_edges(edges_txt);
	edges_txt.close();
	*/



	std::vector<Vertex > vertices;
	std::vector<Vertex> seeds;

	//std::vector<Vertex > hole1;
	//std::vector<Vertex > hole2;
	//std::vector<Vertex > hole3;
	//std::vector<Vertex > hole4;
	//std::vector<Vertex > hole5;

	generate_vertices<RANDOM>(vertices);
	//generate_vertices<STRUCTURED>(vertices);


	//unsigned n1 = 5;
	//unsigned n2 = 30;
	//unsigned n3 = 30;
	//unsigned n4 = 40;
	//unsigned n5 = 30;

	//double r1 = 3.0;
	//double r2 = 2.5;
	//double r3 = 4.0;
	//double r4 = 5.0;
	//double r5 = 2.5;

	//double x1 = -6.0;
	//double y1 = -5.0;

	//double x2 = 0.0;
	//double y2 = 0.0;

	//double x3 = +7.0;
	//double y3 = +7.0;

	//double x4 = 9.0;
	//double y4 = -8.0;

	//double x5 = 9.0;
	//double y5 = -8.0;

	//Vertex * const seed1 = new Vertex(x1, y1);
	//Vertex * const seed2 = new Vertex(x2, y2);
	//Vertex * const seed3 = new Vertex(x3, y3);
	//Vertex * const seed4 = new Vertex(10.0, -3.6);
	////Vertex * const seed4 = new Vertex(x4, y4);
	//Vertex * const seed5 = new Vertex(x5, y5);


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



	///******* AIRFOIL *******/
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

	/***********************/

	for (size_t i = 0; i < vertices.size(); i++)
		pslg.insert_vertex<V_MARKER::FREE>(vertices[i]);


	Triangulation<GK> tri(pslg, seeds);

	std::ofstream vertices_txt;
	std::ofstream edges_txt;
	std::ofstream triangles_txt;


	std::cout << std::endl << "Creating files" << std::endl;

	std::cout << "Vertices..." << std::endl;
	vertices_txt.open("C:\\Users\\pgali\\Desktop\\verticesPoints.txt");
	tri.export_vertices(vertices_txt);
	vertices_txt.close();

	std::cout << "Triangles..." << std::endl;
	triangles_txt.open("C:\\Users\\pgali\\Desktop\\trianglesPoints.txt");
	tri.export_triangles(triangles_txt);
	triangles_txt.close();

	std::cout << "Edges..." << std::endl << std::endl;
	edges_txt.open("C:\\Users\\pgali\\Desktop\\edgesPoints.txt");
	tri.export_edges(edges_txt);
	edges_txt.close();



	//for (size_t i = 0; i < vertices.size(); i++)
	//	tri.insert_vertex(vertices[i]);

	//tri.insert_vertex(vertices.begin(), vertices.end());





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

	//cout << "legalizing" << endl;
	//tri.legalize_vertices();


	//tri.make_individual_hole(seed1);
	////tri.make_individual_hole(seed2);
	////tri.make_individual_hole(seed3);
	////tri.make_individual_hole(seed4);
	//////tri.make_individual_hole(Vertex(3.91, 7.13));
	////////tri.make_individual_hole(Vertex(3.55, 10.53));
	//tri.make_individual_hole(Vertex(-5.0, 5.0));
	//////tri.make_individual_hole(Vertex(8.0, 0.82));
	//////tri.make_individual_hole(Vertex(6.0, 3.7));





	//Mesh mesh(tri);



	//GNUPLOT gnuplot;

	//std::ofstream triangle_txt;

	//unsigned const nt = mesh.get_number_of_triangles();

	//for (size_t i = 0; i < nt; i++) {

	//	triangle_txt.open("C:\\Users\\pgali\\Desktop\\triangles_txt_files\\trianglesPoints.txt");

	//	for (size_t j = 0; j < i + 1; j++)
	//		mesh.export_triangle(triangle_txt, mesh.get_triangle(j));

	//	triangle_txt.close();

	//	gnuplot("plot 'C:\\Users\\pgali\\Desktop\\triangles_txt_files\\trianglesPoints.txt' with lines");

	//}







	//std::ofstream vertices_txt;
	//std::ofstream edges_txt;
	//std::ofstream triangles_txt;


	//std::cout << std::endl << "Creating files" << std::endl;
	//
	//std::cout << "Vertices..." << std::endl;
	//vertices_txt.open("C:\\Users\\pgali\\Desktop\\verticesPoints.txt");
	//mesh.export_vertices(vertices_txt);
	//vertices_txt.close();

	//std::cout << "Triangles..." << std::endl;
	//triangles_txt.open("C:\\Users\\pgali\\Desktop\\trianglesPoints.txt");
	//mesh.export_triangles(triangles_txt);
	//triangles_txt.close();
	//
	//std::cout << "Edges..." << std::endl << std::endl;
	//edges_txt.open("C:\\Users\\pgali\\Desktop\\edgesPoints.txt");
	//mesh.export_edges(edges_txt);
	//edges_txt.close();







	//std::vector<Edge> edges = tri.get_edges();
	//std::vector<Vertex_handle> vert = tri.get_vertices_handles();

	//std::cout << std::endl;
	//std::cout << "/************* Debugging info *************/" << std::endl << std::endl;


	//std::cout << "Number of constrained edges + 3 of super triangle: ";
	//std::cout << std::count_if(edges.begin(), edges.end(), [](Edge e) {return e.is_constrained; }) << std::endl;

	//std::cout << "Number of Neumann edges : ";
	//std::cout << tri.get_num_neumann_edges() << std::endl;

	//std::cout << "Number of Dirichlet edges : ";
	//std::cout << tri.get_num_dirichlet_edges() << std::endl;

	//std::cout << "Number of additonal vertices : ";
	//std::cout << tri.get_new_vertices().size() << std::endl;

	//std::cout << "Number of deleted vertices inside holes : ";
	//std::cout << tri.num_deleted_vertices << std::endl;

	//std::cout << "\nNumber of constrained vertices : ";
	//std::cout << std::count_if(vert.begin(), vert.end(), [](auto it) {return (it)->is_constrained == true; }) << std::endl;

	////std::cout << "num_neumann edges : " << tri.num_neumann << std::endl;

	//std::cout << "\n" <<  "/******************************************/" << std::endl << std::endl;

	//tri.clear();

	//vertices.clear();

	//hole1.clear();
	//hole2.clear();
	//hole3.clear();
	//hole4.clear();
	//hole5.clear();





	std::cout << std::endl;
	//std::cout << "CPU clocks : " << (end - begin) << std::endl << std::endl;
	system("pause");






	return 0; 

};


template<GENERATION gen> 
void generate_vertices(std::vector<Vertex > & v) {


	if (gen == GENERATION::RANDOM) {


		std::random_device					rd;				// obtain a random number from hardware
		std::mt19937						eng(rd());		// seed the generator - Merssene twister
		std::uniform_real_distribution<>	distrX(a1, b1);	// define the range
		std::uniform_real_distribution<>	distrY(a2, b2);	// define the range

		double x;
		double y;

		for (int i = 0; i < N_random; i++) {

			x = distrX(eng);
			y = distrY(eng);

			v.push_back(Vertex(x, y));

		}

	}
	else if (gen == GENERATION::STRUCTURED) {

		double x;
		double y;

		for (int i = 0; i < N1; i++) {

			x = a1 + i * (b1 - a1) / (N1 - 1);

			for (int j = 0; j < N2; j++) {

				y = a2 + j * (b2 - a2) / (N2 - 1);

				v.push_back(Vertex(x, y));

			}
		}

		/*
		for (int i = 0; i < N1; i++) {

			
			double r = 0.2 + i * 10.0 / (N1 - 1);

			for (int j = 0; j < N2; j++) {

				double theta = j*2*3.14 / (N2);

				v.push_back(new Vertex(r*sin(theta), r*cos(theta)));

			}
		}
		*/

	}

};
