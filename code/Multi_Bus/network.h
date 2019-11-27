// Network class definition.
#ifndef NETWORK_H
#define NETWORK_H

// include the definitions of classes, generator, load, loadcont, transmission line, tlinecont, and node
#include "generator.h"
#include "load.h"
#include "loadcont.h"
#include "transl.h"
#include "tlinecont.h"
#include "node.h"
#include <vector>

using namespace std;

class Network {

public:
	Network( int, int, int ); // constructor
	~Network(); // destructor
	void setNetworkVariables( int ); // sets variables of the network
	void runSimulation(int, double[], double[], int); // runs the SCOPF simulation
	int retContCount(): // returns the number of contingency scenarios
	int getGenNumber(); // returns the number of Generators in the network
	double *getPowPrev(); // returns what I think about previous dispatch interval generators
	double *getPowNext(); // returns what I think about next door fellow
	double *getPowSelf(); // returns the values of what I think about myself

private:
	// Define the Network
	int networkID; // ID number of the network instance
	int scenarioIndex; // ID number of the contingency scenario, corresponding to which the network instance is generated, 0 for base case
	int nextUpcoming; // Takes values 1 or 0, to indicate if the network instance belongs to the next to upcoming interval or upcoming interval, respectively
	int genNumber, genFields; // Number of Generators & Fields
	int loadNumber, loadFields; // Number of Loads & Fields
	int translNumber, translFields; // Number of Transmission Lines & Fields
	int deviceTermCount; // Number of device terminals
	int nodeNumber; // Number of Nodes	
	double Rho; // ADMM Tuning Parameter
	int contingencyCount; // Total number of contingency scenarios
	int solverChoice; // Choose a solver from amongst CVXGEN, GUROBI, or MOSEK
	bool Verbose; // Parameter to decide whether to display intermediate results or not
	vector< double > pSelfBeleif; // Belief about the generator MW output of the generators in this dispatch interval from the previous APP iteration
	vector< double > pPrevBeleif; // Belief about the generator MW output of the generators in previous dispatch interval from the previous APP iteration
	vector< double > pNextBeleif; // Belief about the generator MW output of the generators in next dispatch interval from the previous APP iteration
	double pSelfBuffer[100];
	double pPrevBuffer[100];
	double pNextBuffer[100];
	double divConvMWPU; // Divisor, which is set to 100 for all other systems, except two bus system, for which it is set to 1

	// Create vectors of Generators, Loads, Loads under contingency, Transmission lines, Transmission lines under contingency, and Nodes
	vector< Generator > genObject;
	vector< Load > loadObject;
	vector< loadContingency > loadContObject;
	vector< transmissionLine > translObject;
	vector< tlineContingency > tlineContObject;
	vector< Node > nodeObject;

}; // end class Network

#endif // NETWORK_H
