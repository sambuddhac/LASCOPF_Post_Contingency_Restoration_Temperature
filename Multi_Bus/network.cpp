// Member functions for class Network
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cstring>
// include definitions for classes generator, load, transmission line, network and node
#include "generator.h"
#include "transl.h"
#include "tlinecont.h"
#include "load.h"
#include "loadcont.h"
#include "node.h"
#include "network.h"
//#include "gurobi_c++.h"
//#include "mosek.h"
using namespace std;

Network::Network( int val, int countOfScenario, int countOfInterval )
	: networkID( val ), // constructor begins; initialize networkID, Rho, scenarioIndex, dispatch interval count, and contingencyCount through constructor initializer list
	  Rho( 100.0 ),
	  intervalCount( countOfInterval ),
	  scenarioIndex( countOfScenario ),
	  contingencyCount( 0 )
{
	// Initializes the number and fields of Transmission lines, Generators, Loads, Nodes, and Device Terminals.
	translNumber = 0;
	translFields = 0;
	genNumber = 0;
	genFields = 0;
	loadNumber = 0; 
	loadFields = 0;
	deviceTermCount = 0;
	nodeNumber = 0;
	setNetworkVariables( networkID ); // sets the variables of the network

} // end constructor

// destructor
Network::~Network()
{	
	cout << "\nNetwork instance: " << networkID << " for this simulation destroyed. You can now open the output files to view the results of the simulation.\n" << endl;

} // end destructor

void Network::setNetworkVariables( int networkID ) // Function setNetworkVariables starts to initialize the parameters and variables
{
	Verbose = false; // "false" disables intermediate result display. If you want, make it "true"
	do {

		nodeNumber = networkID; // set the number of nodes of the network
		cout << "\nEnter the choice of the solver. Enter 1 for MOSEK, 2 for GUROBI dense, and 3 for GUROBI sparse, 4 for CVXGEN Custom Solver." << endl;
		cin >> solverChoice; // Enter the value of the choice of the solver
		char genFile[ 15 ];
		char tranFile[ 15 ];
		char loadFile[ 15 ];
		
		switch ( nodeNumber ) { // switch structure determines which system to simulate for

			case 14: // 14 Bus case
				
				strcpy( genFile, "Gen14.txt" );			
				strcpy( tranFile, "Tran14.txt" );
				strcpy( loadFile, "Load14.txt" );				
				break; // exit switch

			case 30: // 30 Bus case

				strcpy( genFile, "Gen30.txt" );			
				strcpy( tranFile, "Tran30.txt" );
				strcpy( loadFile, "Load30.txt" );	
				break; // exit switch

			case 57: // 57 Bus case

				strcpy( genFile, "Gen57.txt" );			
				strcpy( tranFile, "Tran57.txt" );
				strcpy( loadFile, "Load57.txt" );	
				break; // exit switch

			case 118: // 118 Bus case

				strcpy( genFile, "Gen118.txt" );			
				strcpy( tranFile, "Tran118.txt" );
				strcpy( loadFile, "Load118.txt" );	
				break; // exit switch

			case 300: // 300 Bus case

				strcpy( genFile, "Gen300.txt" );			
				strcpy( tranFile, "Tran300.txt" );
				strcpy( loadFile, "Load300.txt" );	
				break; // exit switch

			case 3: // 3 Bus case

				strcpy( genFile, "Gen3A.txt" );			
				strcpy( tranFile, "Tran3A.txt" );
				strcpy( loadFile, "Load3A.txt" );	
				break; // exit switch

			case 5: // 5 Bus case

				strcpy( genFile, "Gen5.txt" );			
				strcpy( tranFile, "Tran5.txt" );
				strcpy( loadFile, "Load5.txt" );	
				break; // exit switch

			
			default: // catch all other entries

				cout << "Sorry, invalid case. Can't do simulation at this moment.\n" << endl;
				break; // exit switch

		} // end switch

		// Transmission Lines
		ifstream matrixFirstFile( tranFile, ios::in ); // ifstream constructor opens the file of Transmission lines

		// exit program if ifstream could not open file
		if ( !matrixFirstFile ) {
			exit( 1 );
		} // end if
	
		matrixFirstFile >> translNumber >> translFields; // get the dimensions of the Transmission line matrix
		double matrixTran[ translNumber ][ translFields ]; // Transmission line matrix
		for ( int i = 0; i < translNumber; ++i ) {
			for ( int j = 0; j < translFields; ++j ) {
				matrixFirstFile >> matrixTran[ i ][ j ]; // read the Transmission line matrix
			}
		}

		// Count the total number of contingency scenarios
		for ( int k = 0; k < translNumber; ++k ) {
			contingencyCount += matrixTran[ k ][ 4 ]; // count the number of contingency scenarios
		}
		//*cout << "\nThe total number of contingency scenarios considerred is : " << contingencyCount << endl;
		//*contingencyCount = 0; // Uncomment this statement for purposes of Base-Case/OPF Simulation (Comment out for SCOPF Simulation)
		// Nodes
		for ( int l = 0; l < nodeNumber; ++l ) {

			Node nodeInstance( l + 1, contingencyCount ); // creates nodeInstance object with ID l + 1

			nodeObject.push_back( nodeInstance ); // pushes the nodeInstance object into the vector

		} // end initialization for Nodes

	
		// Resume Creation of Transmission Lines
		int lineContingencyIndex = 0; // variable to keep track of the contingency scenario corresponding to which a line should be outaged
		for ( int k = 0; k < translNumber; ++k ) {
			int tNodeID1, tNodeID2; // node object IDs to which the particular transmission line object is connected
			do {
				//*cout << "Stuck while creating nodes of transmission line: " << ( k + 1 ) << endl;
				//node IDs of the node objects to which this transmission line is connected.
				tNodeID1 = matrixTran[ k ][ 0 ]; //From end
				tNodeID2 = matrixTran[ k ][ 1 ]; //To end
			} while ( ( tNodeID1 <= 0 ) || ( tNodeID1 > nodeNumber ) || ( tNodeID2 <= 0 ) || 
				( tNodeID2 > nodeNumber ) || ( tNodeID1 == tNodeID2) ); // validity check
			double resT, reacT, ptMax, ptMin; // Parameters for Transmission Line
			double contingencyIndex; // Parameter indicating consideration for contingency
			do {
				//*cout << "Stuck while creating transmission line: " << ( k + 1 ) << endl;
				//Resistance:
				resT = matrixTran[ k ][ 2 ];
				//Reactance:
				reacT = matrixTran[ k ][ 3 ];
				//values of maximum allowable power flow on line in the forward and reverse direction:
				//Forward direction:
				ptMax = matrixTran[ k ][ 5 ] / 100;
				ptMin = -ptMax; //Reverse direction
				contingencyIndex = matrixTran[ k ][ 4 ]; // indicates whether the line needs to be considered for contingency analysis or not
			} while ( ( resT < 0 ) || ( reacT <= 0 ) || ( ptMax <= ptMin ) || ( ( contingencyIndex != 1.0 ) && ( contingencyIndex != 0.0 ) ) ); // check the bounds and validity of the parameter values

			// creates transLineInstance object with ID k + 1
			transmissionLine transLineInstance( k + 1, &nodeObject[ tNodeID1 - 1 ], &nodeObject[ tNodeID2 - 1 ], ptMax, reacT, resT ); 
			for ( int counter = 0; counter < contingencyCount; ++counter ) {
				// creates the transmission line contingency object
				if ( ( contingencyIndex == 1.0 ) && ( lineContingencyIndex == counter ) ) {
					tlineContingency tlineContInstance( k + 1, &nodeObject[ tNodeID1 - 1 ], &nodeObject[ tNodeID2 - 1 ], pow( 10.0, -100.0 ), pow( 10.0, 4.0 ), pow( 10.0, 4.0 ), contingencyIndex, ( counter + 1 ) ); // Instantiate the outaged line for contingency
					tlineContInstance.reduceNodeCount(); // Reduce the node connectivity corresponding to the outaged line
					tlineContObject.push_back( tlineContInstance ); // creates the vector of transmission lines corr. to contingency
					//*cout << "\nTransmission line # " << k + 1 << " is being considerred for contingency analysis in scenario # " << counter + 1 << endl;
				}
				else {
					tlineContingency tlineContInstance( k + 1, &nodeObject[ tNodeID1 - 1 ], &nodeObject[ tNodeID2 - 1 ], ptMax, reacT, resT, contingencyIndex, ( counter + 1 ) ); // Instantiate the non-outaged line for contingency
					tlineContObject.push_back( tlineContInstance ); // creates the vector of transmission lines corr. to contingency
				}

			}

			translObject.push_back( transLineInstance ); // pushes the transLineInstance object into the vector

			if ( contingencyIndex == 1.0 )
				++lineContingencyIndex; // increment the line contingency index by one to match to the next contingency scenario

		} // end initialization for Transmission Lines

		// Generators
		ifstream matrixSecondFile( genFile, ios::in ); // ifstream constructor opens the file of Generators

		// exit program if ifstream could not open file
		if ( !matrixSecondFile ) {
			cerr << "\nFile for Generators could not be opened\n" << endl;
			exit( 1 );
		} // end if
	
		matrixSecondFile >> genNumber >> genFields; // get the dimensions of the Generator matrix
		double matrixGen[ genNumber ][ genFields ]; // Generator matrix
		for ( int i = 0; i < genNumber; ++i ) {
			for ( int j = 0; j < genFields; ++j ) {
				matrixSecondFile >> matrixGen[ i ][ j ]; // read the Generator matrix
			}
		}
		
		// Create Generators
		for ( int i = 0; i < genNumber; ++i ) {
			int gNodeID; // node object ID to which the particular generator object is connected
			do {
				gNodeID = matrixGen[ i ][ 0 ];
			} while ( ( gNodeID <= 0 ) || ( gNodeID > nodeNumber ) ); // validity check

			double c2, c1, c0, PgMax, PgMin; // Parameters for Generator
			do {
				//Quadratic Coefficient: 
				c2 = matrixGen[ i ][ 1 ] * 10000;
				//Linear coefficient: 
				c1 = matrixGen[ i ][ 2 ] * 100;
				//Constant term: 
				c0 = matrixGen[ i ][ 3 ];
				//Maximum Limit: 
				PgMax = matrixGen[ i ][ 4 ] / 100;
				//Minimum Limit: 
				PgMin = matrixGen[ i ][ 5 ] / 100;
			} while ( (c2 < 0 ) || ( c1 < 0 ) || ( PgMax <= 0 ) || ( PgMin < 0 ) || ( PgMax <= PgMin ) ); 
			// check the bounds and validity of the parameter values

			Gensolver genParam( contingencyCount, c2, c1, c0, PgMax, PgMin ); // Instantiate the copy constructor for the generator solver object
	
			Generator generatorInstance( i + 1, &nodeObject[ gNodeID - 1 ], genParam, contingencyCount ); // creates generatorInstance object with ID number i + 1

			genObject.push_back( generatorInstance ); // pushes the generatorInstance object into the vector

		} // end initialization for Generators
	
		// Loads
		ifstream matrixThirdFile( loadFile, ios::in ); // ifstream constructor opens the file of Loads

		// exit program if ifstream could not open file
		if ( !matrixThirdFile ) {
			cerr << "\nFile for Loads could not be opened\n" << endl;
			exit( 1 );
		} // end if
	
		matrixThirdFile >> loadNumber >> loadFields; // get the dimensions of the Load matrix
		double matrixLoad[ loadNumber ][ loadFields ]; // Load matrix
		for ( int i = 0; i < loadNumber; ++i ) {
			for ( int j = 0; j < loadFields; ++j ) {
				matrixThirdFile >> matrixLoad[ i ][ j ]; // read the Load matrix
			}
		}
		// Create Loads
		for ( int j = 0; j < loadNumber; ++j ) {
			//cout << "\nEnter the parameters of the " << j + 1 << " -th Load:\n";
			int lNodeID; // node object ID to which the particular load object is connected
			do {
				//node ID of the node object to which this load object is connected.
				lNodeID = matrixLoad[ j ][ 0 ]; 
			} while ( ( lNodeID <= 0 ) || ( lNodeID > nodeNumber ) ); // validity check

			double P_Load; // Parameters for Load
			do {
				//value of allowable power consumption capability of load with a negative sign to indicate consumption:
				//Power Consumption:
				P_Load = matrixLoad[ j ][ 1 ] / 100;
			} while ( -P_Load <= 0 ); // check the bounds and validity of the parameter values

			Load loadInstance( j + 1, &nodeObject[ lNodeID - 1 ], P_Load ); // creates loadInstance object object with ID number j + 1
			for ( int counter = 0; counter < contingencyCount; ++counter ) {
				loadContingency loadContInstance( j + 1, &nodeObject[ lNodeID - 1 ], P_Load, ( counter + 1 ) ); // creates load contingency object
				loadContObject.push_back( loadContInstance ); // creates the vector of transmission lines corr. to contingency
			}

			loadObject.push_back( loadInstance ); // pushes the loadInstance object into the vector

		} // end initialization for Loads

	} while ( (genNumber <= 0 ) || ( nodeNumber <= 0 ) || ( loadNumber <= 0 ) || ( translNumber <= 0 ) || (genFields <= 0 ) || ( loadFields <= 0 ) || ( translFields <= 0 ) );
	// check the bounds and validity of the parameter values
	
	deviceTermCount = genNumber + loadNumber + 2 * translNumber; // total number of device-terminals

} // end setNetworkVariables function

// runSimulation function definition
void Network::runSimulation() //Function runSimulation begins
{
	// Declaration of intermerdiate variables and parameters for running the simulation
	int iteration_count = 1; // iteration counter

	double dualTol = 1.0; // initialize the dual tolerance
	double primalTol; // primal tolerance
        double PrimalTol[ contingencyCount + 1 ]; // Array of primal tolerances for different base-case and contingency scenarios
        double DualTol[ contingencyCount + 1 ]; // Array of dual tolerances for the different base-case and contingency scenarios
        double PrimalTolSplit[ nodeNumber ][ contingencyCount + 1 ]; // 2-D Array of primal tolerances for different base-case and contingency scenarios and nodes
        double DualTolSplit[ nodeNumber ][ contingencyCount + 1 ]; // 2-D Array of dual tolerances for the different base-case and contingency scenarios and nodes
	double ptolsq = 0.0; // initialize the primal tolerance square
	
	vector< int > iterationGraph; // vector of iteration counts
	vector< double > primTolGraph; // vector of primal tolerance
	vector< double > dualTolGraph; // vector of dual tolerance
	vector< double > objectiveValue; // vector of objective function values

	int bufferIndex; // index of the buffer to store past values of voltage iterate, power and angle iterate

	double V_avg[ nodeNumber ][ contingencyCount + 1 ]; // array of average node angle imbalance price from last to last iterate in contingency case
	double vBuffer1[ nodeNumber ][ contingencyCount + 1 ]; // intermediate buffer for average node angle price from last to last iterate in contingency case
	double vBuffer2[ nodeNumber ][ contingencyCount + 1 ]; // intermediate buffer for average node angle price from last iterate in contingency case

	double angleBuffer[ nodeNumber ][ contingencyCount + 1 ]; // buffer for average node voltage angles from present iterate in contingency case
	double angleBuffer1[ nodeNumber ][ contingencyCount + 1 ]; // buffer for average node voltage angles from last iterate in contingency case
	double angtildeBuffer[ deviceTermCount ][ contingencyCount + 1 ]; // Thetatilde from present iterate in contingency case

	double powerBuffer[ deviceTermCount ][ contingencyCount + 1 ]; // Ptilde from present iterate in contingency case
	double powerBuffer1[ deviceTermCount ][ contingencyCount + 1 ]; // Ptilde from last iterate in contingency case
	double pavBuffer[ nodeNumber ][ contingencyCount + 1 ]; // Pav from present iterate in contingency case
	double ptildeinitBuffer[ deviceTermCount ][ contingencyCount + 1 ]; // Ptilde before iterations begin in contingency case

	double uPrice[ deviceTermCount ][ contingencyCount + 1 ]; // u parameter from previous iteration in contingency case
	double vPrice[ deviceTermCount ][ contingencyCount + 1 ]; // v parameter from previous iteration in contingency case
	double LMP[ nodeNumber ][ contingencyCount + 1 ]; // vector of LMPs in contingency case


	double Rho1 = 100.0; // Previous value of Rho from previous iteration
	double W, Wprev; // Present and previous values of W for the PID controller for modifying Rho
	double lambdaAdap = 0.01; // Parameter of the Proportional (P) controller for adjusting the ADMM tuning parameter
	double muAdap = 0.0; // Parameter of the Derivative (D) controller for adjusting the ADMM tuning parameter
        double xiAdap = 0.002; // Parameter of the Integral (I) controller for adjusting the ADMM tuning parameter
        double controllerSum = 0.0; // Integral term of the PID controller
	int setTuning; // parameter to select adaptive rho, fixed rho, and type of adaptive rho

	// Set the type of tuning
	cout << "Enter the tuning mode; Enter 1 for maintaining Rho * primTol = dualTol; 2 for primTol = dualTol; 3 for partly variable and partly fixed Rho; Anything else for Adaptive Rho (with mode-1 being implemented for the first 3000 iterations and then Rho is held constant).\n" << endl;
	cin >> setTuning;

	/*GRBEnv* env = 0; // GUROBI environment for using the GUROBI Solver
	MSKenv_t menv = NULL; 
	MSKtask_t mtask = NULL; 
	MSKrescodee res = MSK_RES_OK; */
	
	// Calculation of initial value of Primal Tolerance before the start of the iterations
	vector< Load >::iterator loadIterator;	
	for ( loadIterator = loadObject.begin(); loadIterator != loadObject.end(); loadIterator++ ) {
		ptolsq = ptolsq + pow( loadIterator->pinitMessage(), 2.0 ); // calls the node to divide by the number of devices connected in base case
	}
	
	vector< loadContingency >::iterator loadContIterator;	
	for ( loadContIterator = loadContObject.begin(); loadContIterator != loadContObject.end(); loadContIterator++ ) {
		ptolsq = ptolsq + pow( loadContIterator->pinitMessageCont(), 2.0 ); // calls the node to divide by the number of devices connected in contingency cases
	}
	primalTol = sqrt( ptolsq ); // initial value of primal tolerance to kick-start the iterations
        dualTol = Rho1 * primalTol; 
	
	// Calculation of initial value of Ptilde before the iterations start
	vector< Generator >::iterator generatorIterator;	
	for ( generatorIterator = genObject.begin(); generatorIterator != genObject.end(); generatorIterator++ ) {
		bufferIndex = generatorIterator->getGenID() - 1;
		ptildeinitBuffer[ bufferIndex ][ 0 ] = -generatorIterator->calcPavInit();
		for ( int contin = 1; contin <= contingencyCount; ++contin ) {
			ptildeinitBuffer[ bufferIndex ][ contin ] = -generatorIterator->calcPavInitc( contin );
		}
	}
	
	loadContIterator = loadContObject.begin();
	for ( loadIterator = loadObject.begin(); loadIterator != loadObject.end(); loadIterator++ ) {
		bufferIndex = genNumber + ( loadIterator->getLoadID() - 1 );
		ptildeinitBuffer[ bufferIndex ][ 0 ] = loadIterator->calcPavInit();
		for ( int contin = 1; contin <= contingencyCount; ++contin ) {
			//cout << "Parameters Initialized" << endl;
			ptildeinitBuffer[ bufferIndex ][ contin ] = loadContIterator->calcPavInitc( contin );
			++loadContIterator;
		}
	}
	
	int temptrans1 = 0; // counter to make sure that two values of Ptilde are accounted for each line
	vector< transmissionLine >::iterator translIterator;	
	vector< tlineContingency >::iterator translContIterator;
	translContIterator = tlineContObject.begin();
	for ( translIterator = translObject.begin(); translIterator != translObject.end(); translIterator++ ) {
		bufferIndex = genNumber + loadNumber + ( translIterator->getTranslID() - 1 ) + temptrans1;
		ptildeinitBuffer[ bufferIndex ][ 0 ] = -translIterator->calcPavInit1(); // Ptilde corresponding to 'from' end
		ptildeinitBuffer[ ( bufferIndex + 1 ) ][ 0 ] = -translIterator->calcPavInit2(); // Ptilde corresponding to 'to' end
		for ( int contin = 1; contin <= contingencyCount; ++contin ) {
			ptildeinitBuffer[ bufferIndex ][ contin ] = -translContIterator->calcPavInitc1( contin );
			ptildeinitBuffer[ ( bufferIndex + 1 ) ][ contin ] = -translContIterator->calcPavInitc2( contin );
			++translContIterator;
		}
		temptrans1++;
	}

	ofstream matrixResultOut( "result.txt", ios::out ); // create a new file result.txt to output the results
	
	// exit program if unable to create file
	if ( !matrixResultOut ) {
		cerr << "File could not be opened" << endl;
		exit( 1 );
	}
	
	if ( Verbose ) {
		matrixResultOut << "\nThe initial value of primal tolerance to kick-start iterations is: " << primalTol << "\nThe initial value of dual tolerance to kick-start iterations is: " << dualTol << endl;	
	}

	clock_t start_s = clock(); // begin keeping track of the time

	// Starting of the ADMM Based Proximal Message Passing Algorithm Iterations
	while( ( primalTol >= 0.01 ) || ( dualTol >= 0.01 ) ) {
	//for ( iteration_count = 1; ( iteration_count < 5001 ); iteration_count++ ) {
		// ( primalTol >= 0.001 ) && ( dualTol >= 0.001 ) // ( iteration_count <= 122 )
	
		if ( Verbose ) {
			matrixResultOut << "\nThe value of primal tolerance before this iteration is: " << primalTol << "\nThe value of dual tolerance before this iteration is: " << dualTol << endl;
			matrixResultOut << "\n**********Start of " << iteration_count << " -th iteration***********\n";
		}
				
		// Recording data for plotting graphs
		
		iterationGraph.push_back( iteration_count ); // stores the iteration count to be graphed
		primTolGraph.push_back( primalTol ); // stores the primal tolerance value to be graphed
		dualTolGraph.push_back( dualTol ); // stores the dual tolerance value to be graphed
		//Initialize the average node angle imbalance price (v) vector from last to last interation, V_avg
		//**if ( iteration_count <= 2 ) {
			for ( int i = 0; i < nodeNumber; i++ ) {
				V_avg[ i ][ 0 ] = 0.0; // initialize to zero for the first and second iterations in base case
				for ( int contin = 1; contin <= contingencyCount; ++contin ) 
					V_avg[ i ][ contin ] = 0.0; // initialize to zero for the 1st and 2nd iteration in contingency case
			}
		//**}
		//**else {
			//**for ( int j = 0; j < nodeNumber; j++ )
				//**V_avg[ j ] = vBuffer1[ j ]; // initialize to the average node v from last to last iteration for 3rd iteration on
		
		//**}
		// Initialize average v, average theta, ptilde, average P before the start of a particular iteration
		if ( iteration_count >= 2 ) {
			for ( int contin = 0; contin <= contingencyCount; ++contin )
				angleBuffer1[ 0 ][ contin ] = 0.0; // set the first node as the slack node, the average voltage angle is always zero
			for ( int i = 1; i < nodeNumber; i++ ) {
				//**vBuffer1[ i ] = vBuffer2[ i ]; // Save to vBuffer1, the average v from last iteration for use in next iteration
				angleBuffer1[ i ][ 0 ] = angleBuffer[ i ][ 0 ]; // Save to angleBuffer1, the average node voltage angle from last iteration in base case
				for ( int contin = 1; contin <= contingencyCount; ++contin ) 
					angleBuffer1[ i ][ contin ] = angleBuffer[ i ][ contin ]; // Save to angleBuffer1, the average node voltage angle from last iteration in contingency case
			}

			for ( int j = 0; j < deviceTermCount; j++ ) {
				powerBuffer1[ j ][ 0 ] = powerBuffer[ j ][ 0 ]; // Save to powerBuffer1, the Ptilde for each device term. from last itern in base case
				for ( int contin = 1; contin <= contingencyCount; ++contin )
					powerBuffer1[ j ][ contin ] = powerBuffer[ j ][ contin ]; // Save to powerBuffer1, the Ptilde for each device term. from last itern in base case
			}

		}
		
		else {
			Wprev = 0.0; // for the first iteration
			for ( int i = 0; i < nodeNumber; i++ ) {
			
				angleBuffer1[ i ][ 0 ] = 0.0; // Set average node voltage angle to zero for 1st iteration in base case
				for ( int contin = 1; contin <= contingencyCount; ++contin )
					angleBuffer1[ i ][ contin ] = 0.0; // Set average node voltage angle to zero for 1st iteration in contingency case
			}

			vector< Node >::iterator nodeIterator;
			for ( nodeIterator = nodeObject.begin(); nodeIterator != nodeObject.end(); nodeIterator++ ) {
				bufferIndex = nodeIterator->getNodeID() - 1;
				pavBuffer[ bufferIndex ][ 0 ] = nodeIterator->devpinitMessage(); // Average node power injection before 1st iteration in base case
				for ( int contin = 1; contin <= contingencyCount; ++contin )
					pavBuffer[ bufferIndex ][ contin ] = nodeIterator->devpinitMessageCont( contin ); // Average node power injection before 1st iteration in contingency case
			}
			for ( int j = 0; j < deviceTermCount; j++ ) {
				powerBuffer1[ j ][ 0 ] = ptildeinitBuffer[ j ][ 0 ]; // Save to powerBuffer1, the Ptilde before the 1st iteration in base case
				for ( int contin = 1; contin <= contingencyCount; ++contin )
					powerBuffer1[ j ][ contin ] = ptildeinitBuffer[ j ][ contin ]; // Save to powerBuffer1, the Ptilde before the 1st iteration in contingency case
			}
		}

		// Distributed Optimizations; Generators' Opt. Problems
		double calcObjective = 0.0;	// initialize the total generator cost for this iteration
		for ( generatorIterator = genObject.begin(); generatorIterator != genObject.end(); generatorIterator++ ) {
			double Pgit[ contingencyCount + 1 ], PowerPrice[ contingencyCount + 1 ], APrice[ contingencyCount + 1 ]; // Generator Power, Power Price, & Angle Price iterates from last iterations
			bufferIndex = generatorIterator->getGenID() - 1;
			int gnid = generatorIterator->getGenNodeID() - 1; // gets the ID number of connection node
			if ( iteration_count > 1 ) { // If 2nd or higher iterations, initialize to previous iterate values
				Pgit[ 0 ] = generatorIterator->genPower();
				PowerPrice[ 0 ] = uPrice[ bufferIndex ][ 0 ]; 
				if ( gnid == 0 ) 
					APrice[ 0 ] = 0.0; // Consider node-1 as the slack node, the angle price is zero always
				else
					APrice[ 0 ] = vPrice[ bufferIndex ][ 0 ];
				for ( int contin = 1; contin <= contingencyCount; ++contin ) {
					Pgit[ contin ] = Pgit[ 0 ]; // To expedite, not doing the function call every time
					PowerPrice[ contin ] = uPrice[ bufferIndex ][ contin ];
					if ( gnid == 0 )
						APrice[ contin ] = 0.0;
					else
						APrice[ contin ] = vPrice[ bufferIndex ][ contin ];
				}
			}
			else { // If 1st iteration, initialize to zero
				for ( int contin = 0; contin <= contingencyCount; ++contin ) {
					Pgit[ contin ] = 0.0;
					PowerPrice[ contin ] = 0.0;
					APrice[ contin ] = 0.0;
				} 
			}
			if ( Verbose ) {
				matrixResultOut << "\nStarting of Generator Optimization Iterations for Generator " << bufferIndex + 1 << "\n";
				for ( int contin = 0; contin <= contingencyCount; ++contin ) {
					matrixResultOut << "Previous power iterate (MW)\n" << Pgit[ contin ] << "\nPrevious average power (MW) for this node\n" << pavBuffer[ gnid ][ contin ] << "\nPrevious power price ($/MWh, LMP)\n" << PowerPrice[ contin ] << "\nAngle price from last to last iterate\n" << V_avg[ gnid ][ contin ] << "\nAngle value from last iterate\n" << angleBuffer1[ gnid ][ contin ] << "\nPrevious angle price\n" << APrice[ contin ] << endl;
				}
			}
			/*if ( solverChoice == 1 ) { // For MOSEK

				generatorIterator->genOptMosek( Rho, Pgit, pavBuffer[ gnid ], PowerPrice, V_avg[ gnid ], angleBuffer1[ gnid ], APrice ); // Solve the Optimization Problem
			}
			if ( solverChoice == 2 ) { // For GUROBI dense
				  try {
    					env = new GRBEnv();
    					double c[] = {1, 1, 0};
    					double  Q[3][3] = {{1, 1, 0}, {0, 1, 1}, {0, 0, 1}};
    					double  A[2][3] = {{1, 2, 3}, {1, 1, 0}};
    					char    sense[] = {'>', '<'};
    					double  rhs[]   = {generatorIterator->getPmin(), generatorIterator->getPmax()};
    					double  lb[]    = {0, 0, 0};
    					bool    success;
    					double  objval, sol[3];

    					success = generatorIterator->genOptGurobi(env, 2, 3, c, &Q[0][0], &A[0][0], sense, rhs,
                             		lb, NULL, NULL, sol, &objval);

    					cout << "x: " << sol[0] << " y: " << sol[1] << " z: " << sol[2] << endl;

  				} catch(GRBException e) {
    				cout << "Error code = " << e.getErrorCode() << endl;
    				cout << e.getMessage() << endl;
  				} catch(...) {
    				cout << "Exception during optimization" << endl;
  				}

  				delete env;
			}
			if ( solverChoice == 3 ) { // for GUROBI Sparse
				generatorIterator->genOptGurobiSparse( Rho, Pgit, pavBuffer[ gnid ], PowerPrice, V_avg[ gnid ], angleBuffer1[ gnid ], APrice ); // Solve the Optimization Problem
			}*/
			if ( solverChoice == 4 ) { // For CVXGEN Custom Solver
				generatorIterator->gpowerangleMessage( Rho, Pgit, pavBuffer[ gnid ], PowerPrice, V_avg[ gnid ], angleBuffer1[ gnid ], APrice ); // Solve the Optimization Problem
			}
			calcObjective = calcObjective + generatorIterator->objectiveGen(); // calculate the total objective after this iteration
		}

		//vector< Load >::const_iterator loadIterator;	// Distributed Optimizations; Loads' Optimization Problems
		for (loadIterator = loadObject.begin(); loadIterator != loadObject.end(); ++loadIterator) {
			double APrice, PPrice; // Load Power Price and Angle Price from last iterations
			bufferIndex = genNumber + ( loadIterator->getLoadID() - 1 );
			int lnid = loadIterator->getLoadNodeID() - 1; // gets ID number of connection node
			if ( iteration_count > 1 ) { // If 2nd or higher iterations, initialize to previous iterate values
				if ( lnid == 0 )
					APrice = 0.0;
				else
					APrice = vPrice[ bufferIndex ][ 0 ];
				PPrice = uPrice[ bufferIndex ][ 0 ];
			}
			else 
				APrice = 0.0; // If 1st iteration, initialize to zero
			if ( Verbose ) {
				matrixResultOut << "\nStarting of Load Optimization Iterations for Load " << loadIterator->getLoadID() << "\n for base-case scenario: " << "\n";
				matrixResultOut << "\nAngle price from last to last iterate\n" << V_avg[ lnid ][ 0 ] << "\nAngle value from last iterate\n" << angleBuffer1[ lnid ][ 0 ] << "\nPrevious angle price\n" << APrice << endl;
			}
			loadIterator->lpowerangleMessage( Rho, V_avg[ lnid ][ 0 ], angleBuffer1[ lnid ][ 0 ], APrice, 0 ); // Solve the Optimization Problem
		}

		//vector< loadContingency >::const_iterator loadContIterator;	// Distributed Optimizations; Loads' Optimization Problems
		for ( loadContIterator = loadContObject.begin(); loadContIterator != loadContObject.end(); loadContIterator++ ) {
			double APrice, PPrice; // Load Power Price and Angle Price from last iterations
			bufferIndex = genNumber + ( loadContIterator->getLoadID() - 1 );
			int lnid = loadContIterator->getLoadNodeID() - 1; // gets ID number of connection node
			int continCounter = loadContIterator->getLoadContCounter(); // gets the contingency scenario count of the load
			if ( iteration_count > 1 ) { // If 2nd or higher iterations, initialize to previous iterate values
				if ( lnid == 0 ) 
					APrice = 0.0;
				else
					APrice = vPrice[ bufferIndex ][ continCounter ];
				PPrice = uPrice[ bufferIndex ][ continCounter ];
			}
			else 
				APrice = 0.0; // If 1st iteration, initialize to zero
			if ( Verbose ) {
				matrixResultOut << "\nStarting of Load Optimization Iterations for Load " << loadContIterator->getLoadID() << "\n for contingency scenario: " << continCounter << "\n";
				matrixResultOut << "\nAngle price from last to last iterate\n" << V_avg[ lnid ][ continCounter ] << "\nAngle value from last iterate\n" << angleBuffer1[ lnid ][ continCounter ] << "\nPrevious angle price\n" << APrice << endl;
			}
			loadContIterator->lpowerangleMessage( Rho, V_avg[ lnid ][ continCounter ], angleBuffer1[ lnid ][ continCounter ], APrice, continCounter ); // Solve the Optimization Problem
		}

		//vector< transmissionLine >::const_iterator translIterator;// Distributed Optimizations; TLine' Optimization Problems
		int temptrans2 = 0;	
		for ( translIterator = translObject.begin(); translIterator != translObject.end(); translIterator++ ) {
			double Ptit1, Ptit2, PowerPrice1, PowerPrice2, APrice1, APrice2; // Tline Power, Power price, Angle price at both ends
			bufferIndex = genNumber + loadNumber + ( translIterator->getTranslID() - 1 ) + temptrans2;
			int tnid1 = translIterator->getTranslNodeID1() - 1; // gets ID number of first conection node
			int tnid2 = translIterator->getTranslNodeID2() - 1; // gets ID number of second connection node
			if (iteration_count > 1 ) { // If 2nd or higher iterations, initialize to previous iterate values
				Ptit1 = translIterator->translPower1();
				Ptit2 = translIterator->translPower2();
				PowerPrice1 = uPrice[ bufferIndex ][ 0 ];
				PowerPrice2 = uPrice[ ( bufferIndex + 1 ) ][ 0 ];
				if ( tnid1 == 0 )
					APrice1 = 0.0;
				else
					APrice1 = vPrice[ bufferIndex ][ 0 ];
				if ( tnid2 == 0 )
					APrice2 = 0.0;
				else
					APrice2 = vPrice[ ( bufferIndex + 1 ) ][ 0 ];
			}
			else { // If 1st iteration, initialize to zero
				Ptit1 = 0.0;
				Ptit2 = 0.0;
				PowerPrice1 = 0.0;
				PowerPrice2 = 0.0;
				APrice1 = 0.0;
				APrice2 = 0.0;
			}
			if ( Verbose ) {
				matrixResultOut << "\nStarting of Transmission Line Optimization Iterations for Transmission line " << translIterator->getTranslID() << "\n for base-case scenario: " << "\n";
				matrixResultOut << "Previous power iterate (MW) for end-1\n" << Ptit1 << "\nPrevious average power (MW) for end-1\n" << pavBuffer[ tnid1 ][ 0 ] << "\nPrevious power price ($/MWh, LMP) for end-1\n" << PowerPrice1 << "\nAngle price from last to last iterate for end-1\n" << V_avg[ tnid1 ][ 0 ] << "\nAngle value from last iterate for end-1\n" << angleBuffer1[ tnid1 ][ 0 ] << "\nPrevious angle price for end-1\n" << APrice1 << "\nPrevious power iterate (MW) for end-2\n" << Ptit2 << "\nPrevious average power (MW) for end-2\n" << pavBuffer[ tnid2 ][ 0 ] << "\nPrevious power price ($/MWh) for end-2\n" << PowerPrice2 << "\nAngle price from last to last iterate for end-2\n" << V_avg[ tnid2 ][ 0 ] << "\nAngle value from last iterate for end-2\n" << angleBuffer1[ tnid2 ][ 0 ] << "\nPrevious angle price for end-2\n" << APrice2 << endl;
			}				
			translIterator->tpowerangleMessage( Rho, Ptit1, pavBuffer[ tnid1 ][ 0 ], PowerPrice1, V_avg[ tnid1 ][ 0 ], angleBuffer1[ tnid1 ][ 0 ], APrice1, Ptit2, pavBuffer[ tnid2 ][ 0 ], PowerPrice2, V_avg[ tnid2 ][ 0 ], angleBuffer1[ tnid2 ][ 0 ], APrice2, 0 ); // Solve the Opt. Problem
			temptrans2++; 
		}

		//vector< tlineContingency >::const_iterator translContIterator;// Distributed Optimizations; TLine' Optimization Problems
		int temptrans3 = 0;
		int temptranscount3 = 0;	
		for (translContIterator = tlineContObject.begin(); translContIterator != tlineContObject.end(); ++translContIterator) {
			double Ptit1, Ptit2, PowerPrice1, PowerPrice2, APrice1, APrice2; // Tline Power, Power price, Angle price at both ends
			bufferIndex = genNumber + loadNumber + ( translContIterator->getTranslID() - 1 ) + temptrans3;
			int tnid1 = translContIterator->getTranslNodeID1() - 1; // gets ID number of first conection node
			int tnid2 = translContIterator->getTranslNodeID2() - 1; // gets ID number of second connection node
			int continCounter = translContIterator->getTranslContCounter(); // gets the contingency scenario count of the line
			if (iteration_count > 1 ) { // If 2nd or higher iterations, initialize to previous iterate values
				Ptit1 = translContIterator->translPower1();
				Ptit2 = translContIterator->translPower2();
				PowerPrice1 = uPrice[ bufferIndex ][ continCounter ];
				PowerPrice2 = uPrice[ ( bufferIndex + 1 ) ][ continCounter ];
				if ( tnid1 == 0 )
					APrice1 = 0.0;
				else
					APrice1 = vPrice[ bufferIndex ][ continCounter ];
				if ( tnid2 == 0 )
					APrice2 = 0.0;
				else
					APrice2 = vPrice[ ( bufferIndex + 1 ) ][ continCounter ];
			}
			else { // If 1st iteration, initialize to zero
				Ptit1 = 0.0;
				Ptit2 = 0.0;
				PowerPrice1 = 0.0;
				PowerPrice2 = 0.0;
				APrice1 = 0.0;
				APrice2 = 0.0;
			}
			if ( Verbose ) {
				matrixResultOut << "\nStarting of Transmission Line Optimization Iterations for Transmission line " << translContIterator->getTranslID() << "\n for contingency scenario: " << continCounter << "\n";
				matrixResultOut << "Previous power iterate (MW) for end-1\n" << Ptit1 << "\nPrevious average power (MW) for end-1\n" << pavBuffer[ tnid1 ][ continCounter ] << "\nPrevious power price ($/MWh, LMP) for end-1\n" << PowerPrice1 << "\nAngle price from last to last iterate for end-1\n" << V_avg[ tnid1 ][ continCounter ] << "\nAngle value from last iterate for end-1\n" << angleBuffer1[ tnid1 ][ continCounter ] << "\nPrevious angle price for end-1\n" << APrice1 << "\nPrevious power iterate (MW) for end-2\n" << Ptit2 << "\nPrevious average power (MW) for end-2\n" << pavBuffer[ tnid2 ][ continCounter ] << "\nPrevious power price ($/MWh) for end-2\n" << PowerPrice2 << "\nAngle price from last to last iterate for end-2\n" << V_avg[ tnid2 ][ continCounter ] << "\nAngle value from last iterate for end-2\n" << angleBuffer1[ tnid2 ][ continCounter ] << "\nPrevious angle price for end-2\n" << APrice2 << endl;
			}				
			translContIterator->tpowerangleMessage( Rho, Ptit1, pavBuffer[ tnid1 ][ continCounter ], PowerPrice1, V_avg[ tnid1 ][ continCounter ], angleBuffer1[ tnid1 ][ continCounter ], APrice1, Ptit2, pavBuffer[ tnid2 ][ continCounter ], PowerPrice2, V_avg[ tnid2 ][ continCounter ], angleBuffer1[ tnid2 ][ continCounter ], APrice2, continCounter ); // Solve the Opt. Problem
			temptranscount3++;
			if ( temptranscount3 == contingencyCount ) 
				temptrans3++; 
		}
		
		if ( Rho <= 10000 ) {
		 if ( setTuning == 1 ) {
			W = ( Rho1 ) * ( primalTol / dualTol ) - 1; // Definition of W for adaptive Rho with Rho1 * primalTol = dualTol
		 }
		 else {
			if ( setTuning == 2 ) {
				W = ( primalTol / dualTol ) - 1; // Definition of W for adaptive Rho with primalTol = dualTol
			}
			else {
	 			if ( setTuning == 3 ) {
					if ( iteration_count <= 1000 ) {
						W = ( Rho1 ) * ( primalTol / dualTol ) - 1; // Definition of W for adaptive Rho with Rho1 * primalTol = dualTol
					}
					if ( ( iteration_count > 1000 ) && ( iteration_count <= 3000 ) ) {
						W = ( Rho1 ) * ( primalTol / dualTol ) - 1; // Definition of W for adaptive Rho with Rho1 * primalTol = dualTol
						xiAdap = 0.0;
					}
					else 	{
						W = 0.0; // Definition of W for fixed Rho
						xiAdap = 0.0;
					}
				}
				else 
					W = 0.0;
			}
		 }
                } else {
                 W = 0.0;
		 xiAdap = 0.0;
                }
		// Calculation of Adaptive Rho
		controllerSum = controllerSum + W;
		Rho1 = Rho; // Store previous Rho
		Rho = ( Rho1 ) * ( exp( ( lambdaAdap * W ) + ( muAdap * ( W - Wprev ) ) + ( xiAdap * controllerSum  ) ) ); // Next iterate value of Rho
		Wprev = W; // Buffering
		//*cout << "\nThe value of Rho and W after iteration " << iteration_count << " are " << Rho << " and " << W << endl;
		if ( ( iteration_count >= 2900 ) && ( iteration_count <= 2910 ) ) {
			cout << "\nThe values of Primal and Dual Tolerances are " << primalTol << " and " << dualTol << endl;
		}
	
		if ( Verbose ) {	
			matrixResultOut << "\n*********Starting of Gather Operation************\n";
		}
		vector< Node >::iterator nodeIterator; // Distributed Optimizations; Nodes' Optimization Problem; Gather Operation
		for ( nodeIterator = nodeObject.begin(); nodeIterator != nodeObject.end(); nodeIterator++ ) {
			bufferIndex = nodeIterator->getNodeID() - 1;
			//**vBuffer2[ bufferIndex ] = ( Rho1 / Rho ) * ( nodeIterator->vavMessage() ); // Gather & Calculate average v after present iteration/node
			if ( bufferIndex == 0 )
				angleBuffer[ bufferIndex ][ 0 ] = 0.0; // consider node 1 as slack node; average voltage angle always zero
			else
				angleBuffer[ bufferIndex ][ 0 ] = nodeIterator->ThetaavMessageCont(); // Calculate average angle after present iteration/node at base case
			//*cout << "\nFor base-case scenario: " << "\nNode Number: " << bufferIndex + 1 /*<< "\nV_avg = " << vBuffer2[ bufferIndex ] */<< "\nTheta_avg = " << angleBuffer[ bufferIndex ][ 0 ] << endl;
			pavBuffer[ bufferIndex ][ 0 ] = nodeIterator->PavMessageCont(); // Calculate average power after present iteration/node at base case
			if ( Verbose ) {
				matrixResultOut << "\nFor base-case scenario: " << "\nNode Number: " << bufferIndex + 1 /*<< "\nV_avg = " << vBuffer2[ bufferIndex ] */<< "\nTheta_avg = " << angleBuffer[ bufferIndex ][ 0 ] << "\nP_avg = " << pavBuffer[ bufferIndex ][ 0 ] << endl;
			}
			for ( int contin = 1; contin <= contingencyCount; ++contin ) {
				if ( bufferIndex == 0 )
					angleBuffer[ bufferIndex ][ contin ] = 0.0; // consider node 1 as slack node; average voltage angle always zero
				else
					angleBuffer[ bufferIndex ][ contin ] = nodeIterator->ThetaavMessageCont( contin ); // Calculate average angle after present iteration/node at contingency
				pavBuffer[ bufferIndex ][ contin ] = nodeIterator->PavMessageCont( contin ); // Calculate average power after present iteration/node at contingency
				//*cout << "\nFor contingency scenario: " << contin << "\nNode Number: " << bufferIndex + 1 /*<< "\nV_avg = " << vBuffer2[ bufferIndex ] */<< "\nTheta_avg = " << angleBuffer[ bufferIndex ][ contin ] << endl;
				if ( Verbose ) {
					matrixResultOut << "\nFor contingency scenario: " << contin << "\nNode Number: " << bufferIndex + 1 /*<< "\nV_avg = " << vBuffer2[ bufferIndex ] */<< "\nTheta_avg = " << angleBuffer[ bufferIndex ][ contin ] << "\nP_avg = " << pavBuffer[ bufferIndex ][ contin ] << endl;
				}
			}
		}

		if ( Verbose ) {
			matrixResultOut << "\n*******Starting of Broadcast Operation*******\n";
		}
		// vector< Generator >::const_iterator generatorIterator;	// Broadcast to Generators
		for ( generatorIterator = genObject.begin(); generatorIterator != genObject.end(); generatorIterator++ ) {
			bufferIndex = generatorIterator->getGenID() - 1;
			if ( Verbose ) {
				matrixResultOut << "\n***Generator: " << bufferIndex + 1 << " results***\n" << endl;
			}
			for ( int contin = 0; contin <= contingencyCount; ++contin ) {
				powerBuffer[ bufferIndex ][ contin ] = generatorIterator->calcPtilde( contin );
				uPrice[ bufferIndex ][ contin ] = ( Rho1 / Rho ) * ( generatorIterator->getu( contin ) );
				angtildeBuffer[ bufferIndex ][ contin ] = generatorIterator->calcThetatilde( contin );
				//*cout << "\nPrimal residual component for angle after iteration number: " << iteration_count << ", terminal n.: " << bufferIndex << ", and contingency count: " << contin << " is: " << angtildeBuffer[ bufferIndex ][ contin ] << endl;
				//*cout << "\nFor contingency scenario: " << contin << " ***Generator: " << bufferIndex + 1 << "\nThetatilde at the end of this iteration is: " << angtildeBuffer[ bufferIndex ][ contin ] << endl;
				//generatorIterator->calcvtilde();
				vPrice[ bufferIndex ][ contin ] = ( Rho1 / Rho ) * ( generatorIterator->getv( contin ) );
				if ( Verbose ) {
					matrixResultOut << "\nFor contingency scenario: " << contin << "\nPower price after this iteration ($/MWh, LMP) is: " << ( Rho / 100 ) * uPrice[ bufferIndex ][ contin ] << "\nAngle price after this iteration is: " << ( Rho ) * vPrice[ bufferIndex ][ contin ] << "\nPtilde after this iteration is: " << powerBuffer[ bufferIndex ][ contin ] << "\nThetatilde at the end of this iteration is: " << angtildeBuffer[ bufferIndex ][ contin ] << endl;
				}
			}
		}

		// Broadcast to Loads
		for ( loadIterator = loadObject.begin(); loadIterator != loadObject.end(); loadIterator++ ) {
			bufferIndex = genNumber + ( loadIterator->getLoadID() - 1 );
			if ( Verbose ) {
				matrixResultOut << "\n***Load: " << loadIterator->getLoadID() << " for base-case scenario results***\n" << endl;
			}
			powerBuffer[ bufferIndex ][ 0 ] = loadIterator->calcPtilde();
			uPrice[ bufferIndex ][ 0 ] = ( Rho1 / Rho ) * ( loadIterator->getu() );
			angtildeBuffer[ bufferIndex ][ 0 ] = loadIterator->calcThetatilde();
			//*cout << "\nPrimal residual component for angle after iteration number: " << iteration_count << ", terminal n.: " << bufferIndex << ", and contingency count: " << 0 << " is: " << angtildeBuffer[ bufferIndex ][ 0 ] << endl;
			//loadIterator->calcvtilde();
			vPrice[ bufferIndex ][ 0 ] = ( Rho1 / Rho ) * ( loadIterator->getv() );
			if ( Verbose ) {
				matrixResultOut << "\nPower price after this iteration ($/MWh, LMP) is: " << ( Rho / 100 ) * uPrice[ bufferIndex ][ 0 ] << "\nAngle price after this iteration is: " << ( Rho ) * vPrice[ bufferIndex ][ 0 ] << "\nPtilde after this iteration is: " << powerBuffer[ bufferIndex ][ 0 ] << "\nThetatilde at the end of this iteration is: " << angtildeBuffer[ bufferIndex ][ 0 ] << endl;
			}
		}

		// Broadcast to Load-Contingency scenarios
		for ( loadContIterator = loadContObject.begin(); loadContIterator != loadContObject.end(); loadContIterator++ ) {
			bufferIndex = genNumber + ( loadContIterator->getLoadID() - 1 );
			int continCounter = loadContIterator->getLoadContCounter(); // gets the contingency scenario count of the load
			if ( Verbose ) {
				matrixResultOut << "\n***Load: " << loadContIterator->getLoadID() << " for contingency scenario: " << continCounter << " results***\n" << endl;
			}
			powerBuffer[ bufferIndex ][ continCounter ] = loadContIterator->calcPtilde( continCounter );
			uPrice[ bufferIndex ][ continCounter ] = ( Rho1 / Rho ) * ( loadContIterator->getu( continCounter ) );
			angtildeBuffer[ bufferIndex ][ continCounter ] = loadContIterator->calcThetatilde( continCounter );
			//*cout << "\nPrimal residual component for angle after iteration number: " << iteration_count << ", terminal n.: " << bufferIndex << ", and contingency count: " << continCounter << " is: " << angtildeBuffer[ bufferIndex ][ continCounter ] << endl;
			//loadIterator->calcvtilde();
			vPrice[ bufferIndex ][ continCounter ] = ( Rho1 / Rho ) * ( loadContIterator->getv( continCounter ) );
			if ( Verbose ) {
				matrixResultOut << "\nPower price after this iteration ($/MWh, LMP) is: " << ( Rho / 100 ) * uPrice[ bufferIndex ][ continCounter ] << "\nAngle price after this iteration is: " << ( Rho ) * vPrice[ bufferIndex ][ continCounter ] << "\nPtilde after this iteration is: " << powerBuffer[ bufferIndex ][ continCounter ] << "\nThetatilde at the end of this iteration is: " << angtildeBuffer[ bufferIndex ][ continCounter ] << endl;
			}
		}

		int temptrans = 0; // temporary count of transmission lines to account for both the ends // Broadcast to Transmission Lines
		for (translIterator = translObject.begin(); translIterator != translObject.end(); ++translIterator) {
			bufferIndex = genNumber + loadNumber + ( translIterator->getTranslID() - 1 ) + temptrans;
			if ( Verbose ) {
				matrixResultOut << "\n***Transmission Line: " << translIterator->getTranslID() << " for base-case scenario results***\n" << endl;
			}
			powerBuffer[ bufferIndex ][ 0 ] = translIterator->calcPtilde1();
			uPrice[ bufferIndex ][ 0 ] = ( Rho1 / Rho ) * ( translIterator->getu1() );
			angtildeBuffer[ bufferIndex ][ 0 ] = translIterator->calcThetatilde1();
			//*cout << "\nPrimal residual component for angle after iteration number: " << iteration_count << ", terminal n.: " << bufferIndex << ", and contingency count: " << 0 << " is: " << angtildeBuffer[ bufferIndex ][ 0 ] << endl;
			//translIterator->calcvtilde1();
			vPrice[ bufferIndex ][ 0 ] = ( Rho1 / Rho ) * ( translIterator->getv1() );
			powerBuffer[ ( bufferIndex + 1 ) ][ 0 ] = translIterator->calcPtilde2();
			uPrice[ ( bufferIndex + 1 ) ][ 0 ] = ( Rho1 / Rho ) * ( translIterator->getu2() );
			angtildeBuffer[ ( bufferIndex + 1 ) ][ 0 ] = translIterator->calcThetatilde2();
			//*cout << "\nPrimal residual component for angle after iteration number: " << iteration_count << ", terminal n.: " << ( bufferIndex + 1 ) << ", and contingency count: " << 0 << " is: " << angtildeBuffer[ ( bufferIndex + 1 ) ][ 0 ] << endl;
			//translIterator->calcvtilde2();
			vPrice[ ( bufferIndex + 1 ) ][ 0 ] = ( Rho1 / Rho ) * ( translIterator->getv2() );
			temptrans++;
			if ( Verbose ) {
				matrixResultOut << "\nPower price ($/MWh, LMP at end-1) after this iteration is: " << ( Rho / 100 ) * uPrice[ bufferIndex ][ 0 ] << "\nAngle price (end-1) after this iteration is: " << ( Rho ) * vPrice[ bufferIndex ][ 0 ] << "\nPtilde (end-1) after this iteration is: " << powerBuffer[ bufferIndex ][ 0 ] << "\nThetatilde (end-1) at the end of this iteration is: " << angtildeBuffer[ bufferIndex ][ 0 ] << "\nPower price ($/MWh, LMP at end-2) after this iteration is: " << ( Rho / 100 ) * uPrice[ ( bufferIndex + 1 ) ][ 0 ] << "\nAngle price (end-2) after this iteration is: " << ( Rho ) * vPrice[ ( bufferIndex + 1 ) ][ 0 ] << "\nPtilde (end-2) after this iteration is: " << powerBuffer[ ( bufferIndex + 1 ) ][ 0 ] << "\nThetatilde (end-2)  at the end of this iteration is: " << angtildeBuffer[ ( bufferIndex + 1 ) ][ 0 ] << endl;
			}
		}

		int temptrans4 = 0; // temporary count of transmission lines to account for both the ends // Broadcast to Transmission Lines-Contingency scenarios
		int temptranscount4 = 0;
		for (translContIterator = tlineContObject.begin(); translContIterator != tlineContObject.end(); ++translContIterator) {
			bufferIndex = genNumber + loadNumber + ( translContIterator->getTranslID() - 1 ) + temptrans4;
			int continCounter = translContIterator->getTranslContCounter(); // gets the contingency scenario count of the line
			if ( Verbose ) {
				matrixResultOut << "\n***Transmission Line: " << translContIterator->getTranslID() << " for contingency scenario: " << continCounter << " results***\n" << endl;
			}
			powerBuffer[ bufferIndex ][ continCounter ] = translContIterator->calcPtilde1( continCounter );
			uPrice[ bufferIndex ][ continCounter ] = ( Rho1 / Rho ) * ( translContIterator->getu1( continCounter ) );
			angtildeBuffer[ bufferIndex ][ continCounter ] = translContIterator->calcThetatilde1( continCounter );
			//*cout << "\nPrimal residual component for angle after iteration number: " << iteration_count << ", terminal n.: " << bufferIndex << ", and contingency count: " << continCounter << " is: " << angtildeBuffer[ bufferIndex ][ continCounter ] << endl;
			//translIterator->calcvtilde1();
			vPrice[ bufferIndex ][ continCounter ] = ( Rho1 / Rho ) * ( translContIterator->getv1( continCounter ) );
			powerBuffer[ ( bufferIndex + 1 ) ][ continCounter ] = translContIterator->calcPtilde2( continCounter );
			uPrice[ ( bufferIndex + 1 ) ][ continCounter ] = ( Rho1 / Rho ) * ( translContIterator->getu2( continCounter ) );
			angtildeBuffer[ ( bufferIndex + 1 ) ][ continCounter ] = translContIterator->calcThetatilde2( continCounter );
			//*cout << "\nPrimal residual component for angle after iteration number: " << iteration_count << ", terminal n.: " << ( bufferIndex + 1 ) << ", and contingency count: " << continCounter << " is: " << angtildeBuffer[ ( bufferIndex + 1 ) ][ continCounter ] << endl;
			//translIterator->calcvtilde2();
			vPrice[ ( bufferIndex + 1 ) ][ continCounter ] = ( Rho1 / Rho ) * ( translContIterator->getv2( continCounter ) );
			temptranscount4++;
			if ( temptranscount4 == contingencyCount )
				temptrans4++;
			if ( Verbose ) {
				matrixResultOut << "\nPower price ($/MWh, LMP at end-1) after this iteration is: " << ( Rho / 100 ) * uPrice[ bufferIndex ][ continCounter ] << "\nAngle price (end-1) after this iteration is: " << ( Rho ) * vPrice[ bufferIndex ][ continCounter ] << "\nPtilde (end-1) after this iteration is: " << powerBuffer[ bufferIndex ][ continCounter ] << "\nThetatilde (end-1) at the end of this iteration is: " << angtildeBuffer[ bufferIndex ][ continCounter ] << "\nPower price ($/MWh, LMP at end-2) after this iteration is: " << ( Rho / 100 ) * uPrice[ ( bufferIndex + 1 ) ][ continCounter ] << "\nAngle price (end-2) after this iteration is: " << ( Rho ) * vPrice[ ( bufferIndex + 1 ) ][ continCounter ] << "\nPtilde (end-2) after this iteration is: " << powerBuffer[ ( bufferIndex + 1 ) ][ continCounter ] << "\nThetatilde (end-2)  at the end of this iteration is: " << angtildeBuffer[ ( bufferIndex + 1 ) ][ continCounter ] << endl;
			}
		}

		//if ( ( iteration_count >= 100 ) && ( ( ( iteration_count % 100 ) == 0 ) || ( iteration_count == MAX_ITER - 1 ) ) ) {
		int i = 0;
		for ( nodeIterator = nodeObject.begin(); nodeIterator != nodeObject.end(); nodeIterator++ ) {
			for ( int j = 0; j <= contingencyCount; ++j ) {
				LMP[ i ][ j ] = ( Rho / 100 ) * nodeIterator->uMessage( j ); // record the LMP values; rescaled and converted to $/MWh
				//nodeIterator->reset(); // reset the node variables that need to start from zero in the next iteration
			}
			++i;	
		}
		//++first;
		//}

		/*for ( int j = 0; j < deviceTermCount; j++ ) {
			for ( int beta = 0; beta <= contingencyCount; ++beta ) {
				cout << "\nBefore reset, Primal residual component for angle after iteration number: " << iteration_count << ", terminal n.: " << j << ", and contingency count: " << beta << " is: " << angtildeBuffer[ j ][ beta ] << endl;
			}
		}*/
		for ( nodeIterator = nodeObject.begin(); nodeIterator != nodeObject.end(); nodeIterator++ ) {
			nodeIterator->reset(); 
		}

		/*for ( int j = 0; j < deviceTermCount; j++ ) {
			for ( int beta = 0; beta <= contingencyCount; ++beta ) {
				cout << "\nAfter reset, Primal residual component for angle after iteration number: " << iteration_count << ", terminal n.: " << j << ", and contingency count: " << beta << " is: " << angtildeBuffer[ j ][ beta ] << endl;
			}
		}*/
		// Calculation of Primal Residual, primalTol at the end of this particular iteration
		double primsum = 0.0;
		for ( int i = 0; i < nodeNumber; i++ ) {
			for ( int alpha = 0; alpha <= contingencyCount; ++alpha ) {
				primsum = primsum + ( pavBuffer[ i ][ alpha ] ) * ( pavBuffer[ i ][ alpha ] );
				//*cout << "\nPrimal residual component for power after iteration number: " << iteration_count << ", node n.: " << i << ", and contingency count: " << alpha << " is: " << pavBuffer[ i ][ alpha ] << endl;
			}
		}
		//*int PrimSum = primsum;
		//*cout << "\nPrimal residual for power after iteration number: " << iteration_count << " is: " << primsum << endl;
		for ( int j = 0; j < deviceTermCount; j++ ) {
			for ( int beta = 0; beta <= contingencyCount; ++beta ) {
				//*cout << "\nPrimal residual component for angle after iteration number: " << iteration_count << ", terminal n.: " << j << ", and contingency count: " << beta << " is: " << angtildeBuffer[ j ][ beta ] << endl;
				primsum = primsum + ( angtildeBuffer[ j ][ beta ] ) * ( angtildeBuffer[ j ][ beta ] );
			}
		}
		//*cout << "\nPrimal residual for angle after iteration number: " << iteration_count << " is: " << primsum - PrimSum << endl;
		primalTol = sqrt( primsum );
		if ( Verbose ) {
			matrixResultOut << "\nPrimal Residual at the end of this iteration is: " << "\t" << primalTol << endl;
		}
		
		// Calculation of Dual Residual, dualTol at the end of this particular iteration
		double sum = 0.0;
		if ( iteration_count > 1 ) {
			for ( int k = 0; k < deviceTermCount; k++ ) {
				for ( int beta = 0; beta <= contingencyCount; ++beta ) {
					sum = sum + ( ( powerBuffer[ k ][ beta ] - powerBuffer1[ k ][ beta ] ) ) * ( ( powerBuffer[ k ][ beta ] - powerBuffer1[ k ][ beta ] ) ); 
					//matrixResultOut << "\npowerBuffer: " << powerBuffer[ k ] << "\npowerBuffer1: " << powerBuffer1[ k ] << endl;
					//*cout << "\nDual residual component for power after iteration number: " << iteration_count << ", terminal n.: " << k << ", and contingency count: " << beta << " is: " << powerBuffer[ k ][ beta ] - powerBuffer1[ k ][ beta ] << endl;
				}
			}
			//*int Sum = sum;
			//*cout << "\nDual residual for power after iteration number: " << iteration_count << " is: " << Sum << endl;
			for ( int i = 0; i < nodeNumber; i++ ) {
				for ( int alpha = 0; alpha <= contingencyCount; ++alpha ) {
					sum = sum + ( ( angleBuffer[ i ][ alpha ] - angleBuffer1[ i ][ alpha ] ) ) * ( ( angleBuffer[ i ][ alpha ] - angleBuffer1[ i ][ alpha ] ) );
					//matrixResultOut << "\nangleBuffer: " << angleBuffer[ i ] << "\nangleBuffer1: " << angleBuffer1[ i ] << endl;
					//*cout << "\nDual residual component for angle after iteration number: " << iteration_count << ", node n.: " << i + 1 << ", and contingency count: " << alpha << " is: " << angleBuffer[ i ][ alpha ] - angleBuffer1[ i ][ alpha ] << endl;
				}
			}
			//*cout << "\nDual residual for angle after iteration number: " << iteration_count << " is: " << sum - Sum << endl;
		}
		else {
			for ( int i = 0; i < nodeNumber; i++ ) {
				for ( int alpha = 0; alpha <= contingencyCount; ++alpha ) 
					sum = sum + ( ( angleBuffer[ i ][ alpha ] ) ) * ( ( angleBuffer[ i ][ alpha ] ) ); 
			}
			//*int Sum = sum;
			//*cout << "\nDual residual for angle after iteration number: " << iteration_count << " is: " << Sum << endl;
			for ( int k = 0; k < deviceTermCount; k++ ) {
				for ( int beta = 0; beta <= contingencyCount; ++beta )
					sum = sum + ( ( powerBuffer[ k ][ beta ] - ptildeinitBuffer[ k ][ beta ] ) ) * ( ( powerBuffer[ k ][ beta ] - ptildeinitBuffer[ k ][ beta ] ) );
			}
			//*cout << "\nDual residual for power after iteration number: " << iteration_count << " is: " << sum - Sum << endl;
		}
		
		dualTol = ( Rho1 ) * sqrt( sum );
		if ( Verbose ) {
			matrixResultOut << sqrt( sum ) << endl;
			matrixResultOut << "\nDual Residual at the end of this iteration is: " << "\t" << dualTol << endl;
			matrixResultOut << "\nObjective value at the end of this iteration is ($): " << "\t" << calcObjective << endl;
			matrixResultOut << "\n****************End of " << iteration_count << " -th iteration***********\n";
		}
		objectiveValue.push_back( calcObjective ); // record the objective values

		iteration_count++;

	} // end of one iteration
	
	clock_t stop_s = clock();  // end
	matrixResultOut << "\nExecution time (s): " << static_cast<double>( stop_s - start_s ) / CLOCKS_PER_SEC << endl;
	matrixResultOut << "\nLast value of dual residual / Rho = " << dualTol / Rho1 << endl;
	matrixResultOut << "\nLast value of primal residual = " << primalTol << endl;
	matrixResultOut << "\nLast value of Rho = " << Rho1 << endl;
	matrixResultOut << "\nLast value of dual residual = " << dualTol << endl;
	matrixResultOut << "\nTotal Number of Iterations = " << iteration_count - 1 << endl;
	cout << "\nExecution time (s): " << static_cast<double>( stop_s - start_s ) / CLOCKS_PER_SEC << endl;

	/**PRINT MW**/
	ofstream devProdOut( "powerResult.txt", ios::out ); // create a new file powerResult.txt to output the results
	
	// exit program if unable to create file
	if ( !devProdOut ) {
		cerr << "File could not be opened" << endl;
		exit( 1 );
	}
	
	devProdOut << "Gen#" << "\t" << "Conn." << "\t" << "MW" << endl;
	for ( generatorIterator = genObject.begin(); generatorIterator != genObject.end(); generatorIterator++ ) {
		devProdOut << generatorIterator->getGenID() << "\t" << generatorIterator->getGenNodeID() << "\t" <<  generatorIterator->genPower() * 100 << endl;
	}
	devProdOut << "T.line#" << "\t" << "From" << "\t" << "To" << "\t" << "From MW" << "\t" << "To MW" << endl;
	for ( translIterator = translObject.begin(); translIterator != translObject.end(); translIterator++ ) {
		devProdOut << translIterator->getTranslID() << "\t" << translIterator->getTranslNodeID1() << "\t" << translIterator->getTranslNodeID2() << "\t" << translIterator->translPower1() * 100 << "\t" << translIterator->translPower2() * 100 << endl;
	}

	devProdOut << "T.line#" << "\t" << "Contingency Scenario" << "\t" << "From" << "\t" << "To" << "\t" << "From MW" << "\t" << "To MW" << endl;
	for (translContIterator = tlineContObject.begin(); translContIterator != tlineContObject.end(); ++translContIterator) {
		devProdOut << translContIterator->getTranslID() << "\t" << translContIterator->getTranslContCounter() << "\t" << translContIterator->getTranslNodeID1() << "\t" << translContIterator->getTranslNodeID2() << "\t" << translContIterator->translPower1() * 100 << "\t" << translContIterator->translPower2() * 100 << endl;
	}

	
	/**PRINT ITERATION COUNTS**/
	// create a new file itresult.txt to output the Iteration Count values
	ofstream iterationResultOut( "itresult.txt", ios::out ); 
	
	// exit program if unable to create file
	if ( !iterationResultOut ) {
		cerr << "File could not be opened" << endl;
		exit( 1 );
	}
	
	iterationResultOut << "\nIteration Count: " << endl;
	vector< int >::iterator iterationCountIterator; 
	for ( iterationCountIterator = iterationGraph.begin(); iterationCountIterator != iterationGraph.end(); iterationCountIterator++ ) {
		iterationResultOut << *iterationCountIterator << endl;
	}

	/**PRINT LMPs**/
	ofstream lmpResultOut( "LMPresult.txt", ios::out ); // create a new file itresult.txt to output the results
	
	// exit program if unable to create file
	if ( !lmpResultOut ) {
		cerr << "File could not be opened" << endl;
		exit( 1 );
	}
	
	lmpResultOut << "\nLocational Marginal Prices for Real Power at nodes ($/MWh): " << endl;
	
	//for ( int j = 0; j < firstIndex; ++j ) {
		//lmpResultOut << "After " << ( j + 1 ) * 100 << " iterations, LMPs are:" << endl;
		for ( int i = 0; i < nodeNumber; ++i ) {
			double Price = 0.0; // Initialize the LMP of the node
			for ( int j = 0; j <= contingencyCount; ++j ) {
				Price += LMP[ i ][ j ];
				lmpResultOut << i + 1 << "\t" << "Contingency Scenario: " << "\t" << j << "\t" << LMP[ i ][ j ] << endl; // print the LMP values
			}
			lmpResultOut << "\nNode : " << "\t" << i + 1 << "\t" << " LMP : " << "\t" << Price << "\t" << " $/MWh" << endl;
		}
	//}

	/**PRINT OBJECTIVE VALUES**/
	// create a new file objective.txt to output the Objective Function value results
	ofstream objectiveResultOut( "objective.txt", ios::out ); 
	
	// exit program if unable to create file
	if ( !objectiveResultOut ) {
		cerr << "File could not be opened" << endl;
		exit( 1 );
	}
	
	objectiveResultOut << "\nObjective value: " << endl;
	vector< double >::iterator objectiveIterator; 
	for ( objectiveIterator = objectiveValue.begin(); objectiveIterator != objectiveValue.end(); objectiveIterator++ )  {
		objectiveResultOut << *objectiveIterator << endl;
	}
	matrixResultOut << "\nLast value of Objective = " << *(objectiveIterator-1) << endl;

	/**PRINT PRIMAL RESIDUAL**/	
	// create a new file primresult.txt to output the Primal Residual results
	ofstream primalResultOut( "primresult.txt", ios::out ); 
	
	// exit program if unable to create file
	if ( !primalResultOut ) {
		cerr << "File could not be opened" << endl;
		exit( 1 );
	}
	
	primalResultOut << "\nPrimal Residual: " << endl;
	vector< double >::iterator primalToleranceIterator;
	for ( primalToleranceIterator = primTolGraph.begin(); primalToleranceIterator != primTolGraph.end(); primalToleranceIterator++ )  {
		primalResultOut << *primalToleranceIterator << endl;
	}

	/**PRINT DUAL RESIDUAL**/	
	// create a new file dualresult.txt to output the Dual Residual results
	ofstream dualResultOut( "dualresult.txt", ios::out ); 
	
	// exit program if unable to create file
	if ( !dualResultOut ) {
		cerr << "File could not be opened" << endl;
		exit( 1 );
	}
	
	dualResultOut << "\nDual Residual: " << endl;
	vector< double >::iterator dualToleranceIterator;
	for ( dualToleranceIterator = dualTolGraph.begin(); dualToleranceIterator != dualTolGraph.end(); dualToleranceIterator++ )  		
	{
		dualResultOut << *dualToleranceIterator << endl;
	}
} // end runSimulation

	




	

	
