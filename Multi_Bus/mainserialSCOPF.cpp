// Main function for implementing APMP Algorithm for the LASCOPF for Post-Contingency Restoration in Multiple Intervals, Limiting Line Temperature case in serial mode
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

#include "network.h" // Network class definition

int main() // function main begins program execution
{
	int netID; // Network ID number to indicate the type of the system with specifying the number of buses
	int dispatchIntervals; // Number of dispatch intervals for the look-ahead simulation/Maximum Restoration Duration (\Gamma_{MRD})
	int numRND; // Number of dispatch intervals for Restoration to Normal Duration (\Gamma_{RND})
	vector< Network > futureNetVector; // Vector of future dispatch interval network objects

	cout << "\nEnter the number of nodes to initialize the network. (Allowed choices are 3, 5, 14, 30, 48, 57, 118, and 300 Bus IEEE Test Bus Systems as of now. So, please restrict yourself to one of these)\n";
	cin >> netID;

	cout << "\nEnter the number of Maximum Restoration Duration (\Gamma_{MRD}).\n";
	cin >> dispatchIntervals;
	cout << "\nEnter the number of dispatch intervals for Restoration to Normal Duration (\Gamma_{RND}).\n";
	cin >> numRND;
	int numRSD = dispatchIntervals - numRND; // Number of dispatch intervals for Restoration to Secure Duration (\Gamma_{RSD})
	cout << endl << "\n*** NETWORK INITIALIZATION STAGE BEGINS ***\n" << endl << endl;
	Network network( netID, 0, 0 ); // create network object corresponding to the upcoming dispatch interval
	int numberOfCont = network.retContCount(); // gets the number of contingency scenarios in the variable numberOfCont
	futureNetVector.push_back( network ); // push to the vector of future network instances 
	for ( int i = 0; i <= numberOfCont; ++i ) {
		for ( int j = 1; j <= dispatchIntervals; ++j ) {
			Network futureNet( netID, i, j ); // create the network instances for the future dispatch intervals, which includes one base case and as many networks as the number of contingency scenarios
			futureNetVector.push_back( futureNet ); // push to the vector of future network instances
		}
	}

	cout << "\n*** NETWORK INITIALIZATION STAGE ENDS ***\n" << endl;

	int numberOfGenerators = futureNetVector[0].getGenNumber(); // get the number of generators in the system
	int iterCountAPP = 1; // Iteration counter for APP coarse grain decomposition algorithm	
	int lambdaAPP[2*(numberOfCont+1)*numberOfGnenerators]; // Array of APP Lagrange Multipliers for achieving consensus among the values of power generated, as guessed by different intervals
	double powDiff[2*(numberOfCont+1)*numberOfGnenerators]; // Array of lack of consensus between generation values, as guessed by different intervals
	double finTol = 1000.0; //Initial Guess of the Final tolerance of the APP iteration/Stopping criterion
	cout << endl << "\n*** APMP ALGORITHM LASCOPF FOR POST-CONTINGENCY RESTORATION IN MULTIPLE INTERVALS WITH LINE TEMPERATURE LIMITING SIMULATION (SERIAL IMPLEMENTATION) BEGINS ***\n" << endl << endl;
	cout << endl << "\n*** SIMULATION IN PROGRESS; PLEASE DON'T CLOSE ANY WINDOW OR OPEN ANY OUTPUT FILE YET ... ***\n" << endl << endl;

//*********************************************AUXILIARY PROBLEM PRINCIPLE (APP) COARSE GRAINED DECOMPOSITION COMPONENT******************************************************//	
	do { // APP Coarse grain iterations start
		network.runSimulation(lambdaAPP, powDiff, 0); // start simulation
		for ( int netSimCount = 0; netSimCount <= numberOfCont; ++netSimCount )
			futureNetVector[netSimCount].runSimulation(lambdaAPP, powDiff, netSimCount); // start simulation
		for ( int i = 0; i <= numberOfCont; ++i ) {
			for ( int j = 0; j < numberOfGenerators; ++j ) {
				powDiff[2*i*numberOfGenerators+j]=futureNetVector[i].getPowDiffSelf(); // what I think about myself Vs. what next door fellow thinks about me
				powDiff[(2*i+1)*numberOfGenerators+j]=futureNetVector[i].getPowDiffNext(); // what I think about next door fellow Vs. what next door fellow thinks about himself
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
		for ( int i = 0; i < 2*(numberOfCont+1)*numberOfGnenerators; ++i ) {
			tolAPP = tolAPP + pow(powDiff[i], 2);
		}
		finTol = sqrt(tolAPP);
	} while (finTol>=0.005); //Check the termination criterion of the APP iterations
//****************************************END OF AUXILIARY PROBLEM PRINCIPLE (APP) COARSE GRAINED DECOMPOSITION COMPONENT******************************************************//

	cout << "\n*** LASCOPF FOR POST-CONTINGENCY RESTORATION IN MULTIPLE INTERVALS WITH LINE TEMPERATURE LIMITING SIMULATION ENDS ***\n" << endl;

	return 0; // indicates successful program termination

} // end main
