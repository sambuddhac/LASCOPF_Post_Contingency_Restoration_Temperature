// Main function for implementing APMP Algorithm for the LASCOPF for Post-Contingency Restoration Controlling Line Temperature case in serial mode
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <ctime>
#include <cmath>
#include <fstream>
#include <cstring>
#include "gurobi_c++.h"
#include "supernetwork.h" // Super Network class definition
using namespace std;

int main() // function main begins program execution
{
	int netID; // Network ID number to indicate the type of the system with specifying the number of buses/nodes
	int solverChoice; // Solver choice among CVXGEN-ADMM-PMP+APP fully distributed, GUROBI-ADMM-PMP+APP fully distributed, or GUROBI APP half distributed
	int dummyIntervalChoice; // The choice pertaining to whether to include a dummy interval at the start or not (Inclusion of a dummy interval may speed up convergence and/or improve accuracy of solution). Enter 1 to include and 0 to not include
	int RNDIntervals; // Number of dispatch intervals for Restoration to Normal Rating
	int RSDIntervals; // Number of dispatch intervals for Restoration to Secure Rating
	int nextChoice; // The index which decides whether to include the ramping constraint for the last supernetwork for considering the next imaginary interval or not
	int contSolverAccuracy; // Switch to select between whether an extensive/exhaustive (and presumably more accurate) solver for contingency scenarios is desired, or just a simpler one is desired; 1 for former, 0 for latter 
	int last = 0; // flag to indicate the last interval; last = 0, for dispatch interval that is not the last one; last = 1, for the last interval
	vector< superNetwork* > futureNetVector; // Vector of future look-ahead dispatch interval supernetwork objects
	cout << "\nEnter the number of nodes to initialize the network. (Allowed choices are 2, 3, 5, 14, 30, 48, 57, 118, and 300 Bus IEEE Test Bus Systems as of now. So, please restrict yourself to one of these)\n";
	cin >> netID;
	cout << "\nEnter the switch value to select between whether an extensive/exhaustive (and presumably more accurate) solver for contingency scenarios is desired, or just a simpler one is desired; 1 for former, 0 for latter\n";
	cin >> contSolverAccuracy;
	cout << "\nEnter the choice of the solver for SCOPF of each dispatch interval, 1 for GUROBI-APMP(ADMM/PMP+APP), 2 for CVXGEN-APMP(ADMM/PMP+APP), 3 for GUROBI APP Coarse Grained, 4 for centralized GUROBI SCOPF" << endl;
	cin >> solverChoice;
	cout << "\nEnter the choice pertaining to whether you want to consider the ramping constraint to the next interval, for the last interval: 0 for not considering and 1 for considering" << endl;
	cin >> nextChoice;
	int setRhoTuning; // parameter to select adaptive rho, fixed rho, and type of adaptive rho
	if ((solverChoice==1) || (solverChoice==2)) { // APMP Fully distributed, Bi-layer (N-1) SCOPF Simulation 
		cout << "Enter the tuning mode; Enter 1 for maintaining Rho * primTol = dualTol; 2 for primTol = dualTol; anything else for Adaptive Rho (with mode-1 being implemented for the first 3000 iterations and then Rho is held constant).\n" << endl;
		cin >> setRhoTuning;
	}
	else
		setRhoTuning = 0; // Otherwise, if we aren't using ADMM-PMP, Rho tuning is unnecessary, 0 is a dummy value
	cout << "\nEnter the choice pertaining to whether to include a dummy interval at the start or not (Inclusion of a dummy interval may speed up convergence and/or improve accuracy of solution). Enter 1 to include and 0 to not include\n";
	cin >> dummyIntervalChoice;
	cout << "\nEnter the number of look-ahead dispatch intervals for restoring line flows to within normal long-term ratings.\n";
	cin >> RNDIntervals;
	cout << "\nEnter the number of furthermore look-ahead dispatch intervals for making the system secure w.r.t. next set of contingencies.\n";
	cin >> RSDIntervals;

	cout << endl << "\n*** SUPERNETWORK INITIALIZATION STAGE BEGINS ***\n" << endl << endl;
	GRBEnv* environmentGUROBI = new GRBEnv("GUROBILogFile.log"); // GUROBI Environment object for storing the different optimization models
	superNetwork* supernet = new superNetwork( netID, solverChoice, setRhoTuning, 0, 0, 0, 0, nextChoice, dummyIntervalChoice, contSolverAccuracy, 0 ); // create the network instances for the future dummy zero dispatch intervals
	int numberOfCont = supernet->retContCount(); // gets the number of contingency scenarios in the variable numberOfCont
	futureNetVector.push_back( supernet ); // push to the vector of future network instances 
	superNetwork* supernet1 = new superNetwork( netID, solverChoice, setRhoTuning, 0, 0, 1, 0, nextChoice, dummyIntervalChoice, contSolverAccuracy, 0 ); // create the network instances for the future upcoming dispatch intervals
	futureNetVector.push_back( supernet1 ); // push to the vector of future network instances 
	for ( int i = 0; i <= numberOfCont; ++i ) {
		for ( int j = 1; j < RNDIntervals; ++j ) {
			int lineOutaged = 0; // the serial number of transmission line outaged in any scenario: default value is zero
			if (i > 0) // for the post-contingency scenarios
				lineOutaged = futureNetVector[0]->indexOfLineOut(i); // gets the serial number of transmission line outaged in this scenario 
			superNetwork* supernet = new superNetwork( netID, solverChoice, setRhoTuning, i, j, 2, last, nextChoice, dummyIntervalChoice, contSolverAccuracy, lineOutaged ); // create the network instances for the future next-to-upcoming-dispatch intervals for pos-contingency cases
			futureNetVector.push_back( supernet ); // push to the vector of future network instances
		}
		for ( int j = 0; j <= RSDIntervals; ++j ) {
			int lineOutaged = 0; // the serial number of transmission line outaged in any scenario: default value is zero
			if (i > 0) // for the post-contingency scenarios
				lineOutaged = futureNetVector[0]->indexOfLineOut(i); // gets the serial number of transmission line outaged in this scenario 
			if (j==RSDIntervals) // set the flag to 1 to indicate the last interval
				last = 1; // set the flag to 1 to indicate the last interval
			superNetwork* supernet = new superNetwork( netID, solverChoice, setRhoTuning, i, (j+RNDIntervals), 2, last, nextChoice, dummyIntervalChoice, contSolverAccuracy, lineOutaged ); // create the network instances for the future next-to-upcoming-dispatch intervals for pos-contingency cases
			futureNetVector.push_back( supernet ); // push to the vector of future network instances
		}
	}
	cout << "\n*** SUPERNETWORK INITIALIZATION STAGE ENDS ***\n" << endl;

	int numberOfGenerators = futureNetVector[0]->getGenNumber(); // get the number of generators in the system
	int numberOfLines = futureNetVector[0]->getTransNumber(); // get the number of remaining transmission lines in the system
	int iterCountAPP = 1; // Iteration counter for APP coarse grain decomposition algorithm
	double alphaAPP = 100.0; // APP Parameter/Path-length
	int consLagDim; // Dimension of the vectors of APP Lagrange Multipliers and Power Generation Consensus
	int consLineLagDim; // Dimension of the vectors of APP Lagrange Multipliers and Line Flow consensus for (RND-1) intervals for temperature limiting
	if (dummyIntervalChoice==1)
		consLagDim = 2*((numberOfCont+1)*(RNDIntervals+RSDIntervals)+1)*numberOfGenerators; // Dimension of the vectors of APP Lagrange Multipliers and Power Generation Consensus
	else
		consLagDim = 2*((numberOfCont+1)*(RNDIntervals+RSDIntervals))*numberOfGenerators; // Dimension of the vectors of APP Lagrange Multipliers and Power Generation Consensus	
	consLineLagDim = (RNDIntervals-1)*numberOfLine*(numberOfCont+1); // Dimension of the vectors of APP Lagrange Multipliers and remaining transmission lines in the system-flow Consensus
	double lambdaAPP[consLagDim]; // Array of APP Lagrange Multipliers for achieving consensus among the values of power generated, as guessed by different intervals
	double powDiff[consLagDim]; // Array of lack of consensus between generation values, as guessed by different intervals
	double lambdaAPPLine[consLineLagDim]; // Array of APP Lagrange Multipliers for achieving consensus among the values of line flows, as guessed by different intervals
	double powDiffLine[consLineLagDim]; // Array of lack of consensus between line flows, as guessed by different intervals
	int supernetNum; // Number of supernetworks considered for computation (whether dummy interval is included or not)
	int supernetNumNext; // Number of future supernetworks about which generation belief are held by the existing supernetworks
	int supernetLineNumNext; // Number of future supernetworks about which line flow beliefs are held by the existing supernetworks
	if (dummyIntervalChoice==1) {
		supernetNum=(numberOfCont+1)*(RNDIntervals+RSDIntervals)+2;
		supernetNumNext=(numberOfCont+1)*(RNDIntervals+RSDIntervals+1)+1;
	}
	else {
		supernetNum=(numberOfCont+1)*(RNDIntervals+RSDIntervals)+1;
		supernetNumNext=(numberOfCont+1)*(RNDIntervals+RSDIntervals+1);
	}
	supernetLineNumNext=(numberOfCont+1)*numberOfLines*(RNDIntervals-1);
	double powerSelfGen[supernetNum*numberOfGenerators]; // what I think about myself
	double powerNextBel[supernetNumNext*numberOfGenerators]; // what I think about next door fellow
	double powerPrevBel[supernetNum*numberOfGenerators]; // what I think about previous door fellow
	double powerNextFlowBel[supernetLineNumNext]; // what I think about flows for next door fellow
	double powerSelfFlowBel[supernetLineNumNext]; // what I think about flows for myself (only look-ahead intervals 1 to (RNDIntervals-1))
	for ( int i = 0; i < consLagDim; ++i ) {
		lambdaAPP[i] = 0.0; // Initialize lambdaAPP for the first iteration of APP and ADMM-PMP
		powDiff[i] = 0.0; // Initialize powDiff for the first iteration of APP and ADMM-PMP
	}
	for ( int i = 0; i < consLineLagDim; ++i ) {
		lambdaAPPLine[i] = 0.0; // Initialize lambdaAPPLine for the first iteration of APP and ADMM-PMP
		powDiffLine[i] = 0.0; // Initialize powDiffLine for the first iteration of APP and ADMM-PMP
	}
	// Initializing the self belief, next belief, and previous beliefs about MW generated by a warm start with the respective generation values of last realized dispatch
	for ( int i = 0; i < supernetNum; ++i ) {
		for ( int j = 0; j < numberOfGenerators; ++j ) {
			powerSelfGen[i*numberOfGenerators+j] = *((futureNetVector[0])->getPowPrev()+j); // Use 0.0 if warm start is not desired
			if (i==0) {
				powerPrevBel[i*numberOfGenerators+j] = *((futureNetVector[0])->getPowPrev()+j); // Actual value of previous interval dispatch for the first interval
			}
			else {
				powerPrevBel[i*numberOfGenerators+j] = *((futureNetVector[0])->getPowPrev()+j); // Use 0.0 if warm start is not desired
			}
		}
	}
	for ( int i = 0; i < supernetNumNext; ++i ) {
		for ( int j = 0; j < numberOfGenerators; ++j ) {
			powerNextBel[i*numberOfGenerators+j] = *((futureNetVector[0])->getPowPrev()+j); // Use 0.0 if warm start is not desired
		}
	}
	for ( int i = 0; i <= numberOfCont; ++i ) {
		for ( int k = 1; k < RNDIntervals; ++k ) {
			for ( int j = 0; j < numberOfLines; ++j ) {
				powerNextFlowBel[i*(RNDIntervals-1)*numberOfLines+(k-1)*numberOfLines+j] = 0.0; // Difficult to warm start so just assume 0
				powerSelfFlowBel[i*(RNDIntervals-1)*numberOfLines+(k-1)*numberOfLines+j] = 0.0; // Difficult to warm start so just assume 0
			}
		}
	}
	double finTol = 1000.0; // Initial Guess of the Final tolerance of the APP iteration/Stopping criterion
	double finTolDelayed = 1000.0; // Initial Guess of the Final tolerance delayed of the APP iteration/Stopping criterion
	string outputAPPFileName;
	if (solverChoice==1)
		outputAPPFileName = "/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/LASCOPF_Post_Contingency_Restoration_Temperature/output/ADMM_PMP_GUROBI/resultOuterAPP-SCOPF.txt";
	if (solverChoice==2)
		outputAPPFileName = "/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/LASCOPF_Post_Contingency_Restoration_Temperature/output/ADMM_PMP_CVXGEN/resultOuterAPP-SCOPF.txt";
	if (solverChoice==3)
		outputAPPFileName = "/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/LASCOPF_Post_Contingency_Restoration_Temperature/output/APP_Quasi_Decent_GUROBI/resultOuterAPP-SCOPF.txt";
	if (solverChoice==4)
		outputAPPFileName = "/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/LASCOPF_Post_Contingency_Restoration_Temperature/output/APP_GUROBI_Centralized_SCOPF/resultOuterAPP-SCOPF.txt";
	ofstream matrixResultAPPOut( outputAPPFileName, ios::out ); // create a new file to output the results
	// exit program if unable to create file
	if ( !matrixResultAPPOut ) {
		cerr << "File could not be opened" << endl;
		exit( 1 );
	}
	
	matrixResultAPPOut << endl << "\n*** APMP ALGORITHM BASED LASCOPF FOR POST CONTINGENCY RESTORATION CONTROLLING LINE TEMPERATURE SIMULATION (SERIAL IMPLEMENTATION) SUPERNETWORK LAYER BEGINS ***\n" << endl << endl;
	matrixResultAPPOut << endl << "\n*** SIMULATION IN PROGRESS; PLEASE DON'T CLOSE ANY WINDOW OR OPEN ANY OUTPUT FILE YET ... ***\n" << endl << endl;	
	matrixResultAPPOut << endl << "\nInitial Value of the Tolerance to kick-start the APP outer iterations= " << finTol << "\n" << endl << endl;
	matrixResultAPPOut << "APP Iteration Count" << "\t" << "APP Tolerance" << "\n";	
	cout << endl << "\n*** APMP ALGORITHM BASED LASCOPF FOR POST CONTINGENCY RESTORATION CONTROLLING LINE TEMPERATURE SIMULATION (SERIAL IMPLEMENTATION) SUPERNETWORK LAYER BEGINS ***\n" << endl << endl;
	cout << endl << "\n*** SIMULATION IN PROGRESS; PLEASE DON'T CLOSE ANY WINDOW OR OPEN ANY OUTPUT FILE YET ... ***\n" << endl << endl;

//*********************************************AUXILIARY PROBLEM PRINCIPLE (APP) COARSE GRAINED DECOMPOSITION COMPONENT******************************************************//
	vector<double> largestSuperNetTimeVec; // vector largest value of the computational time in a particular outer APP iteration for any supernetwork
	vector<double> singleSuperNetTimeVec; // vector of the computational times in a particular outer APP iteration for all supernetworks
	largestSuperNetTimeVec.clear(); // clear for upcoming iteration
	double actualSuperNetTime = 0; // Initialize the supernetwork computational time
	clock_t start_s = clock(); // begin keeping track of the time
	do { // APP Coarse grain iterations start
	//for ( iterCountAPP = 1; iterCountAPP <= 100; ++iterCountAPP ) {
		singleSuperNetTimeVec.clear(); // clear for upcoming iteration
		if (dummyIntervalChoice==1) { // Outermost APP layer with a dummy zero interval at the beginning
			for ( int netSimCount = 0; netSimCount <= (numberOfCont+1)*(RNDIntervals+RSDIntervals)+1; ++netSimCount ) {
				if (netSimCount == 0)
					cout << "\nStart of " << iterCountAPP << " -th Outermost APP iteration for dummy zero dispatch interval" << endl;
				if (netSimCount == 1)
					cout << "\nStart of " << iterCountAPP << " -th Outermost APP iteration for " << netSimCount << " -th dispatch interval" << endl;
				else
					cout << "\nStart of " << iterCountAPP << " -th Outermost APP iteration for second dispatch interval for " << netSimCount-2 << " -th post-contingency scenario" << endl;
				futureNetVector[netSimCount]->runSimulation(iterCountAPP, lambdaAPP, powDiff, powerSelfGen, powerNextBel, powerPrevBel, lambdaAPPLine, powDiffLine, powerSelfFlowBel, powerNextFlowBel, environmentGUROBI); // start simulation
				double singleSuperNetTime = futureNetVector[netSimCount]->getvirtualNetExecTime(); // get the computational time for each supernetwork under the assumption of nested and complete parallelism of each generator optimization, within each coarse grain optimization in the supernetworks
				actualSuperNetTime += singleSuperNetTime; // Actual time
				singleSuperNetTimeVec.push_back(singleSuperNetTime); // Vector of all independent supernet solve times
			}
			double largestSuperNetTime = *max_element(singleSuperNetTimeVec.begin(), singleSuperNetTimeVec.end()); // get the laziest solve-time for this iteration
			largestSuperNetTimeVec.push_back(largestSuperNetTime); // vector of all te laziest supernet calculations over all iterations
			// Calculate the power generation opinions and disagreements between the different dispatch interval coarse grains
			for ( int i = 0; i <= (numberOfCont+1)*(RNDIntervals+RSDIntervals)+1; ++i ) {
				if (i==0) {
					for ( int j = 0; j < numberOfGenerators; ++j ) {
						powDiff[2*i*numberOfGenerators+j]=*(futureNetVector[i]->getPowSelf()+j)-*(futureNetVector[i+1]->getPowPrev()+j); // what I think about myself Vs. what next door fellow thinks about me
						powerSelfGen[i*numberOfGenerators+j]=*(futureNetVector[i]->getPowSelf()+j); // what I think about myself
						powerNextBel[i*numberOfGenerators+j]=*(futureNetVector[i]->getPowNext(0,i)+j); // what I think about next door fellow
						powerPrevBel[i*numberOfGenerators+j]=*(futureNetVector[i]->getPowPrev()+j); // what I think about previous interval
						powDiff[(2*i+1)*numberOfGenerators+j]=*(futureNetVector[i]->getPowNext(0,i)+j)-*(futureNetVector[i+1]->getPowSelf()+j); // what I think about next door fellow Vs. what next door fellow thinks about himself
					}
				}
				else {
					for ( int j = 0; j < numberOfGenerators; ++j ) {
						powerSelfGen[i*numberOfGenerators+j]=*(futureNetVector[i]->getPowSelf()+j); // what I think about myself
						if (i == 1) {
							for (int continCounter=0; continCounter <= numberOfCont; ++continCounter) {
								powerNextBel[(i+continCounter)*numberOfGenerators+j]=*(futureNetVector[i]->getPowNext(continCounter, i)+j); // what I think about next door fellow
								powDiff[2*(i+continCounter)*numberOfGenerators+j]=*(futureNetVector[i]->getPowSelf()+j)-*(futureNetVector[i+continCounter+1]->getPowPrev()+j); // what I think about myself Vs. what next door fellow thinks about me
								powDiff[(2*(i+continCounter)+1)*numberOfGenerators+j]=*(futureNetVector[i]->getPowNext(continCounter, i)+j)-*(futureNetVector[i+continCounter+1]->getPowSelf()+j); // what I think about next door fellow Vs. what next door fellow thinks about himself
							}
						}
						else {
							powerNextBel[(i+numberOfCont)*numberOfGenerators+j]=*(futureNetVector[i]->getPowNext(0, i)+j); // what I think about next door fellow
							for (int continCounter=0; continCounter <= numberOfCont; ++continCounter) { // Inefficient: Better, should run only for the particular value of continCounter for the particular i
								if (i!=((continCounter+1)*(RNDIntervals+RSDIntervals)+1)) { // Make sure the last supernetworks for any post-cont scenario are left out
									powDiff[2*(i+numberOfCont)*numberOfGenerators+j]=*(futureNetVector[i]->getPowSelf()+j)-*(futureNetVector[i+1]->getPowPrev()+j); // what I think about myself Vs. what next door fellow thinks about me
									powDiff[(2*(i+numberOfCont)+1)*numberOfGenerators+j]=*(futureNetVector[i]->getPowNext(0, i)+j)-*(futureNetVector[i+1]->getPowSelf()+j); // what I think about next door fellow Vs. what next door fellow thinks about himself
								}
							}
						}
						powerPrevBel[i*numberOfGenerators+j]=*(futureNetVector[i]->getPowPrev()+j); // what I think about previous interval
					}
					for ( int j = 0; j < numberOfLines; ++j ) {
						if (i == 1) {
							for (int continCounter=0; continCounter <= numberOfCont; ++continCounter) {
								for ( int k = 1; k < RNDIntervals; ++k ) {
									powerNextFlowBel[continCounter*(RNDIntervals-1)*numberOfLines+(k-1)*numberOfLines+j]=*(futureNetVector[i]->getPowFlowNext(continCounter, i, k)+j); // what I think about next door fellow
									powDiffLine[continCounter*(RNDIntervals-1)*numberOfLines+(k-1)*numberOfLines+j]=*(futureNetVector[i]->getPowFlowNext(continCounter, i, k)+j)-*(futureNetVector[2+continCounter*(RNDIntervals+RSDIntervals)+(k-1)]->getPowFlowSelf()+j);
								}
							}
						}
						else {
							for (int continCounter=0; continCounter <= numberOfCont; ++continCounter) {
								for ( int k = 1; k < RNDIntervals; ++k ) {
									if (i==2+continCounter*(RNDIntervals+RSDIntervals)+(k-1))
										powerSelfFlowBel[continCounter*(RNDIntervals-1)*numberOfLines+(k-1)*numberOfLines+j]=*(futureNetVector[i]->getPowFlowSelf()+j);
								}
							}
						}
					}
				}
			}
			// Tuning the alphaAPP
			if ( ( iterCountAPP > 5 ) && ( iterCountAPP <= 10 ) )
				alphaAPP = 75.0;
			if ( ( iterCountAPP > 10 ) && ( iterCountAPP <= 15 ) )
				alphaAPP = 50;
			if ( ( iterCountAPP > 15 ) && ( iterCountAPP <= 20 ) )
				alphaAPP = 25;
			if ( ( iterCountAPP > 20 ) )
				alphaAPP = 10;
			// Update power disagreement Lagrange Multipliers
			for ( int i = 0; i < consLagDim; ++i ) {
				lambdaAPP[i] = lambdaAPP[i] + alphaAPP * (powDiff[i]);
			}
			for ( int i = 0; i < consLineLagDim; ++i ) {
				lambdaAPPLine[i] = lambdaAPPLine[i] + alphaAPP * (powDiffLine[i]);
			}
			//++iterCountAPP; // increment the APP iteration counter
			double tolAPP = 0.0;
			double tolAPPDelayed = 0.0; // APP tolerance, excluding the first (dummy) interval
			matrixResultAPPOut << (iterCountAPP-1) << "\t";
			for ( int i = 0; i < consLagDim; ++i ) {
				tolAPP = tolAPP + pow(powDiff[i], 2);
				if (i>=2*numberOfGenerators)
					tolAPPDelayed = tolAPPDelayed + pow(powDiff[i], 2);
				matrixResultAPPOut << powDiff[i] << "\t";
			}
			for ( int i = 0; i < consLineLagDim; ++i ) {
				tolAPP = tolAPP + pow(powDiffLine[i], 2);
				tolAPPDelayed = tolAPPDelayed + pow(powDiffLine[i], 2);
				matrixResultAPPOut << powDiffLine[i] << "\t";
			}
			matrixResultAPPOut << "\n";
			finTol = sqrt(tolAPP);
			finTolDelayed = sqrt(tolAPPDelayed);
			matrixResultAPPOut << (iterCountAPP-1) << "\t" << finTol << "\t" << finTolDelayed <<"\n";
			++iterCountAPP;
			cout << "\nFinal Value of Outer APP Tolerance " << finTol << "\nAnd Final Value of Outer APP Delayed Tolerance " << finTolDelayed << endl;
			finTol = finTolDelayed; // Assign finTolDelayed to finTol in order for checking on the condition at the end of the loop
		}
		else { // Outermost APP layer without a dummy zero interval at the beginning
			for ( int netSimCount = 0; netSimCount < (numberOfCont+1)*(RNDIntervals+RSDIntervals)+1; ++netSimCount ) {
				if (netSimCount == 0)
					cout << "\nStart of " << iterCountAPP << " -th Outermost APP iteration for " << netSimCount+1 << " -th dispatch interval" << endl;
				else
					cout << "\nStart of " << iterCountAPP << " -th Outermost APP iteration for second dispatch interval for " << netSimCount-1 << " -th post-contingency scenario" << endl;
				futureNetVector[netSimCount+1]->runSimulation(iterCountAPP, lambdaAPP, powDiff, powerSelfGen, powerNextBel, powerPrevBel, lambdaAPPLine, powDiffLine, powerSelfFlowBel, powerNextFlowBel, environmentGUROBI); // start simulation
				double singleSuperNetTime = futureNetVector[netSimCount+1]->getvirtualNetExecTime(); // get the computational time for each supernetwork under the assumption of nested and complete parallelism of each generator optimization, within each coarse grain optimization in the supernetworks
				actualSuperNetTime += singleSuperNetTime; // Actual time
				singleSuperNetTimeVec.push_back(singleSuperNetTime); // Vector of all independent supernet solve times
			}
			double largestSuperNetTime = *max_element(singleSuperNetTimeVec.begin(), singleSuperNetTimeVec.end()); // get the laziest solve-time for this iteration
			largestSuperNetTimeVec.push_back(largestSuperNetTime); // vector of all te laziest supernet calculations over all iterations
			// Calculate the power generation opinions and disagreements between the different dispatch interval coarse grains
			for ( int i = 0; i < (numberOfCont+1)*(RNDIntervals+RSDIntervals)+1; ++i ) {
				if (i==0) {
					for ( int j = 0; j < numberOfGenerators; ++j ) {
						powerSelfGen[i*numberOfGenerators+j]=*(futureNetVector[i+1]->getPowSelf()+j); // what I think about myself
						for (int continCounter=0; continCounter <= numberOfCont; ++continCounter) {
							powerNextBel[(i+continCounter)*numberOfGenerators+j]=*(futureNetVector[i+1]->getPowNext(continCounter, (i+1))+j); // what I think about next door fellow
							powDiff[2*(i+continCounter)*numberOfGenerators+j]=*(futureNetVector[i+1]->getPowSelf()+j)-*(futureNetVector[i+continCounter+2]->getPowPrev()+j); // what I think about myself Vs. what next door fellow thinks about me
							powDiff[(2*(i+continCounter)+1)*numberOfGenerators+j]=*(futureNetVector[i+1]->getPowNext(continCounter, (i+1))+j)-*(futureNetVector[i+continCounter+2]->getPowSelf()+j); // what I think about next door fellow Vs. what next door fellow thinks about himself
						}
						powerPrevBel[i*numberOfGenerators+j]=*(futureNetVector[i+1]->getPowPrev()+j); // what I think about previous interval
					}
				}
				else {
					for ( int j = 0; j < numberOfGenerators; ++j ) {
						powerSelfGen[i*numberOfGenerators+j]=*(futureNetVector[i+1]->getPowSelf()+j); // what I think about myself
						powerNextBel[(i+numberOfCont)*numberOfGenerators+j]=*(futureNetVector[i+1]->getPowNext(0, (i+1))+j); // what I think about next door fellow
						powerPrevBel[i*numberOfGenerators+j]=*(futureNetVector[i+1]->getPowPrev()+j); // what I think about previous interval
						powerNextBel[(i+numberOfCont)*numberOfGenerators+j]=*(futureNetVector[i]->getPowNext(0, i)+j); // what I think about next door fellow
						for (int continCounter=0; continCounter <= numberOfCont; ++continCounter) { // Inefficient: Better, should run only for the particular value of continCounter for the particular i
							if (i!=((continCounter+1)*(RNDIntervals+RSDIntervals))) { // Make sure the last supernetworks for any post-cont scenario are left out
								powDiff[2*(i+numberOfCont)*numberOfGenerators+j]=*(futureNetVector[i+1]->getPowSelf()+j)-*(futureNetVector[i+2]->getPowPrev()+j); // what I think about myself Vs. what next door fellow thinks about me
								powDiff[(2*(i+numberOfCont)+1)*numberOfGenerators+j]=*(futureNetVector[i+1]->getPowNext(0, (i+1))+j)-*(futureNetVector[i+2]->getPowSelf()+j); // what I think about next door fellow Vs. what next door fellow thinks about himself
							}
						}
					}
				}
				for ( int j = 0; j < numberOfLines; ++j ) {
					if (i == 0) {
						for (int continCounter=0; continCounter <= numberOfCont; ++continCounter) {
							for ( int k = 1; k < RNDIntervals; ++k ) {
								powerNextFlowBel[continCounter*(RNDIntervals-1)*numberOfLines+(k-1)*numberOfLines+j]=*(futureNetVector[i+1]->getPowFlowNext(continCounter, (i+1), k)+j); // what I think about next door fellow
								powDiffLine[continCounter*(RNDIntervals-1)*numberOfLines+(k-1)*numberOfLines+j]=*(futureNetVector[i+1]->getPowFlowNext(continCounter, (i+1), k)+j)-*(futureNetVector[2+continCounter*(RNDIntervals+RSDIntervals)+(k-1)]->getPowFlowSelf()+j);
							}
						}
					}
					else {
						for (int continCounter=0; continCounter <= numberOfCont; ++continCounter) {
							for ( int k = 1; k < RNDIntervals; ++k ) {
								if (i==1+continCounter*(RNDIntervals+RSDIntervals)+(k-1))
									powerSelfFlowBel[continCounter*(RNDIntervals-1)*numberOfLines+(k-1)*numberOfLines+j]=*(futureNetVector[i+1]->getPowFlowSelf()+j);
							}
						}
					}
				}
			}
			// Tuning the alphaAPP
			if ( ( iterCountAPP > 5 ) && ( iterCountAPP <= 10 ) )
				alphaAPP = 75.0;
			if ( ( iterCountAPP > 10 ) && ( iterCountAPP <= 15 ) )
				alphaAPP = 50;
			if ( ( iterCountAPP > 15 ) && ( iterCountAPP <= 20 ) )
				alphaAPP = 25;
			if ( ( iterCountAPP > 20 ) )
				alphaAPP = 10;
			// Update power disagreement Lagrange Multipliers
			for ( int i = 0; i < consLagDim; ++i ) {
				lambdaAPP[i] = lambdaAPP[i] + alphaAPP * (powDiff[i]);
			}
			for ( int i = 0; i < consLineLagDim; ++i ) {
				lambdaAPPLine[i] = lambdaAPPLine[i] + alphaAPP * (powDiffLine[i]);
			}
			//++iterCountAPP; // increment the APP iteration counter
			double tolAPP = 0.0;
			matrixResultAPPOut << (iterCountAPP-1) << "\t";
			for ( int i = 0; i < consLagDim; ++i ) {
				tolAPP = tolAPP + pow(powDiff[i], 2);
				matrixResultAPPOut << powDiff[i] << "\t";
			}
			for ( int i = 0; i < consLineLagDim; ++i ) {
				tolAPP = tolAPP + pow(powDiffLine[i], 2);
				matrixResultAPPOut << powDiffLine[i] << "\t";
			}
			matrixResultAPPOut << "\n";
			finTol = sqrt(tolAPP);
			matrixResultAPPOut << (iterCountAPP-1) << "\t" << finTol << "\n";
			++iterCountAPP;
			cout << "\nFinal Value of Outer APP Tolerance " << finTol << endl;
		}
	} while (finTol>=0.005); //Check the termination criterion of the APP iterations
	//}
//****************************************END OF AUXILIARY PROBLEM PRINCIPLE (APP) COARSE GRAINED DECOMPOSITION COMPONENT******************************************************//

	clock_t stop_s = clock();  // end
	cout << "\n*** LASCOPF FOR POST-CONTINGENCY RESTORATION CONTROLLING LINE TEMPERATURE SIMULATION SUPERNETWORK LAYER ENDS ***\n" << endl;
	matrixResultAPPOut << "\nExecution Outermost layer time (s): " << static_cast<double>( stop_s - start_s ) / CLOCKS_PER_SEC << endl;
	matrixResultAPPOut << "\nVirtual Outermost layer Execution time (s): " << static_cast<double>( stop_s - start_s ) / CLOCKS_PER_SEC  - actualSuperNetTime + accumulate(largestSuperNetTimeVec.begin(), largestSuperNetTimeVec.end(), 0.0)<< endl;
	cout << "\nVirtual Outermost layer Execution time (s): " << static_cast<double>( stop_s - start_s ) / CLOCKS_PER_SEC  - actualSuperNetTime + accumulate(largestSuperNetTimeVec.begin(), largestSuperNetTimeVec.end(), 0.0);
	cout << "\nExecution time (s): " << static_cast<double>( stop_s - start_s ) / CLOCKS_PER_SEC << endl;
	delete environmentGUROBI; // Free the memory of the GUROBI environment object
	for ( int i = 0; i <= (numberOfCont+1)*(RNDIntervals+RSDIntervals)+1; ++i ) {
		delete futureNetVector[i]; // Free the memory of future network instances
	}
	return 0; // indicates successful program termination

} // end main
