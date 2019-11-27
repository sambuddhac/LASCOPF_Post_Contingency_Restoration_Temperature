// Main function for implementing APMP Algorithm for the LASCOPF for Post-Contingency Restoration in One Interval case in serial mode
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <cstring>
using namespace std;

#include "network.h" // Network class definition

int main() // function main begins program execution
{
	int netID; // Network ID number to indicate the type of the system with specifying the number of buses
	vector< Network* > futureNetVector; // Vector of second future dispatch interval network objects

	cout << "\nEnter the number of nodes to initialize the network. (Allowed choices are 2, 3, 5, 14, 30, 48, 57, 118, and 300 Bus IEEE Test Bus Systems as of now. So, please restrict yourself to one of these)\n";
	cin >> netID;

	cout << endl << "\n*** NETWORK INITIALIZATION STAGE BEGINS ***\n" << endl << endl;
	Network* network = new Network( netID, 0, 0 ); // create network object corresponding to the upcoming dispatch interval
	int numberOfCont = network->retContCount(); // gets the number of contingency scenarios in the variable numberOfCont
	futureNetVector.push_back( network ); // push to the vector of future network instances 
	for ( int i = 0; i <= numberOfCont; ++i ) {
		Network* network = new Network( netID, i, 1 ); // create the network instances for the future dispatch intervals, which includes one base case and as many networks as the number of contingency scenarios
		futureNetVector.push_back( network ); // push to the vector of future network instances
	}

	cout << "\n*** NETWORK INITIALIZATION STAGE ENDS ***\n" << endl;

	int numberOfGenerators = futureNetVector[0]->getGenNumber(); // get the number of generators in the system
	int iterCountAPP = 1; // Iteration counter for APP coarse grain decomposition algorithm
	double alphaAPP = 100.0; // APP Parameter/Path-length
	int consLagDim = 2*(numberOfCont+1)*numberOfGnenerators; // Dimension of the vectors of APP Lagrange Multipliers and Power Generation Consensus	
	double lambdaAPP[consLagDim]; // Array of APP Lagrange Multipliers for achieving consensus among the values of power generated, as guessed by different intervals
	double powDiff[consLagDim]; // Array of lack of consensus between generation values, as guessed by different intervals
	for ( int i = 0; i < consLagDim; ++i ) {
		lambdaAPP[i] = 0.0; // Initialize lambdaAPP for the first iteration of APP and ADMM-PMP
		powDiff[i] = 0.0; // Initialize powDiff for the first iteration of APP and ADMM-PMP
	}
	double finTol = 1000.0; // Initial Guess of the Final tolerance of the APP iteration/Stopping criterion
	string outputAPPFileName = "resultAPP.txt";
	ofstream matrixResultAPPOut( outputAPPFileName, ios::out ); // create a new file result.txt to output the results

	// exit program if unable to create file
	if ( !matrixResultAPPOut ) {
		cerr << "File could not be opened" << endl;
		exit( 1 );
	}
	
	cout << endl << "\n*** APMP ALGORITHM BASED LASCOPF FOR POST-CONTINGENCY RESTORATION IN ONE INTERVAL SIMULATION (SERIAL IMPLEMENTATION) BEGINS ***\n" << endl << endl;
	cout << endl << "\n*** SIMULATION IN PROGRESS; PLEASE DON'T CLOSE ANY WINDOW OR OPEN ANY OUTPUT FILE YET ... ***\n" << endl << endl;
	matrixResultAPPOut << endl << "\nInitial Value of the Tolerance to kick-start the APP outer iterations= " << finTol << "\n" << endl << endl;
	matrixResultAPPOut << "APP Itearion Count" << "\t" << "APP Tolerance" << "\n";	
	clock_t start_s = clock(); // begin keeping track of the time
	cout << endl << "\n*** APMP ALGORITHM BASED LASCOPF FOR POST-CONTINGENCY RESTORATION IN ONE INTERVAL (SERIAL IMPLEMENTATION) BEGINS ***\n" << endl << endl;
	cout << endl << "\n*** SIMULATION IN PROGRESS; PLEASE DON'T CLOSE ANY WINDOW OR OPEN ANY OUTPUT FILE YET ... ***\n" << endl << endl;

//*********************************************AUXILIARY PROBLEM PRINCIPLE (APP) COARSE GRAINED DECOMPOSITION COMPONENT******************************************************//	
	do { // APP Coarse grain iterations start
	//for ( iterCountAPP = 1; iterCountAPP <= 100; ++iterCountAPP ) {
		for ( int netSimCount = 0; netSimCount <= (numberOfCont+1); ++netSimCount )
			futureNetVector[netSimCount]->runSimulation(lambdaAPP, powDiff, netSimCount); // start simulation
		for ( int i = 0; i <= numberOfCont; ++i ) {
			for ( int j = 0; j < numberOfGenerators; ++j ) {
				powDiff[2*i*numberOfGenerators+j]=*(futureNetVector[i]->getPowSelf()+j)-*(futureNetVector[i+1]->getPowPrev()+j); // what I think about myself Vs. what next door fellow thinks about me
				powDiff[(2*i+1)*numberOfGenerators+j]=*(futureNetVector[i]->getPowNext()+j)-*(futureNetVector[i+1]->getPowSelf()+j); // what I think about next door fellow Vs. what next door fellow thinks about himself
			}
		}
		for ( int i = 0; i <= numberOfCont; ++i ) {
			for ( int j = 0; j < numberOfGenerators; ++j ) {
				lambdaAPP[2*i*numberOfGenerators+j] = lambdaAPP[2*i*numberOfGenerators+j] + alphaAPP * (powDiff[2*i*numberOfGenerators+j]); // what I think about myself Vs. what next door fellow thinks about me
				lambdaAPP[(2*i+1)*numberOfGenerators+j] = lambdaAPP[(2*i+1)*numberOfGenerators+j] + alphaAPP * (powDiff[(2*i+1)*numberOfGenerators+j]); // what I think about next door fellow Vs. what next door fellow thinks about himself
			}
		}
		++iterCountAPP; // increment the APP iteration counter
		double tolAPP = 0.0;
		for ( int i = 0; i < consLagDim; ++i ) {
			tolAPP = tolAPP + pow(powDiff[i], 2);
		}
		finTol = sqrt(tolAPP);
		matrixResultAPPOut << (iterCountAPP-1) << "\t" << finTol << "\n";
	} while (finTol>=0.005); //Check the termination criterion of the APP iterations
	//}
//****************************************END OF AUXILIARY PROBLEM PRINCIPLE (APP) COARSE GRAINED DECOMPOSITION COMPONENT******************************************************//

	cout << "\n*** LASCOPF FOR POST-CONTINGENCY RESTORATION IN ONE INTERVAL SIMULATION ENDS ***\n" << endl;
	clock_t stop_s = clock();  // end
	matrixResultAPPOut << "\nExecution time (s): " << static_cast<double>( stop_s - start_s ) / CLOCKS_PER_SEC << endl;
	cout << "\nExecution time (s): " << static_cast<double>( stop_s - start_s ) / CLOCKS_PER_SEC << endl;

	for ( int i = 0; i <= (numberOfCont+1); ++i ) {
		delete futureNetVector[i]; // push to the vector of future network instances
	}

	return 0; // indicates successful program termination

} // end main
