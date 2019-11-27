// Generator class definition.
// Member functions defined in generator.cpp
#ifndef GENERATOR_H
#define GENERATOR_H
#include "gensolver.h" // definition of Gensolver class
//#include "gurobi_c++.h" // definition of GUROBI header file
//#include "mosek.h" // definition of MOSEK header file
#include <vector>
#include <array>

using namespace std;

class Node; // forward declaration

class Generator {

public:
	Generator( int, Node *, Gensolver &, int ); // constructor
	~Generator(); // destructor
	int getGenID(); // returns the ID of the Generator
	int getGenNodeID(); // returns the ID of the Node to which the Generator is connected
	void gpowerangleMessage( double, double[], double[], double[], double[], double[], double[] ); // Real power and angle iterate with CVXGEN
	/*static bool genOptGurobi(); // Real power and angle iterate with GUROBI (dense)
	void genOptGurobiSparse( double, double[], double[], double[], double[], double[], double[] ); // Real power and angle iterate with GUROBI (sparse)
	void genOptMosek( double, double[], double[], double[], double[], double[], double[] ); // Real power and angle iterate with MOSEK*/
	void setGenData(); // Function to set values of generator data
	double genPower(); // Function to return the value of the power iterate
	double calcPtilde( int ); // Calculates the difference between Power iterate and average power
	double calcPavInit() const; // gets the Ptilde before iterations start from node
	double calcPavInitc( int ) const; // gets the Ptilde before iterations start for the contingency case
	double calcvtilde( int ) const; // Calculates the difference between v and average v
	double objectiveGen(); // Calculates the objective function value after each iteration of ADMM
	double getu( int ) const; // Gets the value of price of real power from node
	double calcThetatilde( int ); // calculates the difference between voltage angle and average voltage angle
	double getv( int ); // Gets the value of price of voltage angle 
	double getPmax(); // Gets the maximum power generating limit
	double getPmin(); // Gets the minimum power generating limit
	double getquadCoeff(); // Gets the Quadratic cost coefficient
	double getlinCoeff(); // Gets the Linear cost coefficient
	double getconstCoeff(); // Gets the Constant coefficient

private:
	int genID; // Generator object id number
	int deviceNature; // 1 for Generating device, 0 otherwise
	double Pg; // Power iterates
	double *Thetag; // Angle iterates
	int contCountGen; // Total number of contingency scenarios
	Node *connNodegPtr; // connection node object
	double *v; // Voltage angle constraint Lagrange Multiplier
	Gensolver genSolver; // Generator Solver object

}; //end class Generator

#endif // GENERATOR_H
	
