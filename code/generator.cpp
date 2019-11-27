// Member functions for class Generator.
#include <iostream>
#include <iomanip>
// include Generator class definition from generator.h
#include "generator.h"
// Include definition of Node class 
#include "node.h"
// Include Generator solver class defintion
#include "gensolverFirst.h" // definition of Gensolver class for base case scenario first interval
#include "gensolverFirstBase.h" // definition of Gensolver class for base case scenario dummy zeroth interval in case of Dummy zero (i.e. gensolver for zeroth interval)
#include "gensolverFirstCont.h" // definition of Gensolver class for contingency case scenario first interval
#include "gensolverFirstDZ.h" // definition of Gensolver class for base case scenario first interval with dummy zero interval
#include "gensolverFirstDZCont.h" // definition of Gensolver class for contingency case scenario first interval with dummy zero interval
#include "gensolverSecondBase.h" // definition of Gensolver class for base case scenario for second interval
#include "gensolverSecondCont.h" // definition of Gensolver class for contingency case scenario for second interval
#include "gensolverCont.h" // definition of Gensolver class for contingency scenario
#include "gurobi_c++.h" // includes definition of the GUROBI solver header file

using namespace std;

Generator::Generator( int idOfGen, int interval, int lastFlag, int contScenarioCount, int PCScenarioCount, int baseCont, int dummyZero, int accuracy, Node *nodeConng, GensolverFirstBase &paramOfGenFirstBase, GensolverFirstDZ &paramOfGenDZBase, GensolverFirst &paramOfGenFirst, GensolverSecondBase &paramOfGenSecondBase, GensolverFirstCont &paramOfGenFirstCont, GensolverFirstDZCont &paramOfGenDZCont, GensolverSecondCont &paramOfGenSecondCont, GensolverCont &paramOfGenCont, int countOfContingency, int genTotal ) // constructor begins
	: genID( idOfGen ),
	  numberOfGenerators(genTotal),
	  dispatchInterval(interval),
	  flagLast(lastFlag),
	  dummyZeroIntFlag(dummyZero),
	  contSolverAccuracy(accuracy),
	  scenarioContCount( contScenarioCount ),
	  postContScenCount( PCScenarioCount ),
	  baseContScenario( baseCont ),
	  connNodegPtr( nodeConng ),
	  genSolverFirstBase( paramOfGenFirstBase ),
	  genSolverDZBase( paramOfGenDZBase ),
	  genSolverFirst( paramOfGenFirst ),
	  genSolverSecondBase( paramOfGenSecondBase ),
	  genSolverDZCont( paramOfGenDZCont ),
	  genSolverFirstCont( paramOfGenFirstCont ),
	  genSolverSecondCont( paramOfGenSecondCont ),
	  genSolverCont( paramOfGenCont ),
	  contCountGen( countOfContingency )
{
	//cout << "\nInitializing the parameters of the generator with ID: " << genID << endl;
	connNodegPtr->setgConn(idOfGen); // increments the generation connection variable to node
	PgenPrev=genSolverFirstBase.getPgPrev();
	setGenData(); // calls setGenData member function to set the parameter values

} // constructor ends

Generator::~Generator() // destructor
{
	//cout << "\nThe generator object having ID " << genID << " have been destroyed.\n";

} // end of destructor

int Generator::getGenID() // function getGenID begins
{
	return genID; // returns the ID of the generator object
} // end of getGenID function

int Generator::getGenNodeID() // function getGenNodeID begins
{
	return connNodegPtr->getNodeID(); // returns the ID number of the node to which the generator object is connected
} // end of getGenNodeID function

void Generator::setGenData() // start setGenData function
{
	Pg = 0.0; // Initialize power iterate
	PgenNextPtr = NULL;
	Thetag = 0.0; // Initialize angle iterate
	v = 0.0; // Initialize the Lagrange multiplier corresponding voltage angle constraint to zero
	
} // end of setGenData function

void Generator::gpowerangleMessage(int outerAPPIt, int  APPItCount, double gsRho, double Pgenprev, double Pgenavg, double Powerprice, double Angpriceavg, double Angavg, double Angprice, double PgenPrevAPP, double PgenAPP, double PgenAPPInner, double PgenNextAPP[], double AAPPExternal, double BAPPExternal[], double DAPPExternal[], double LambAPP1External[], double LambAPP2External[], double LambAPP3External, double LambAPP4External, double BAPP[], double LambAPP1[]) //const // function gpowerangleMessage begins
{
	double BAPPNew[contCountGen];
	double LambdaAPPNew[contCountGen];
	double BAPPExtNew[contCountGen+1];
	double DAPPExtNew[contCountGen+1];
	double LambdaAPP1ExtNew[contCountGen+1];
	double LambdaAPP2ExtNew[contCountGen+1];
	double PgNextAPPNew[contCountGen+1];
	for (int counterCont = 0; counterCont < contCountGen; ++counterCont) {
		BAPPNew[counterCont]=0; 
		LambdaAPPNew[counterCont]=0;
	}
	for (int counterCont = 0; counterCont <= contCountGen; ++counterCont) {
		BAPPExtNew[counterCont]=0; 
		LambdaAPP1ExtNew[counterCont]=0;
		DAPPExtNew[counterCont]=0; 
		LambdaAPP2ExtNew[counterCont]=0;
		PgNextAPPNew[counterCont]=0;
	}
	if ( baseContScenario == 0 ) { // Use the solver for base cases
		if ( dummyZeroIntFlag == 1 ) { // If dummy zero interval is considered
			if ((dispatchInterval==0) && (flagLast==0)) { // For the dummy zeroth interval
				for (int counterCont = 0; counterCont < contCountGen; ++counterCont) {
					BAPPNew[counterCont]=BAPP[counterCont*numberOfGenerators+(genID-1)]; 
					LambdaAPPNew[counterCont]=LambAPP1[counterCont*numberOfGenerators+(genID-1)];
				}
				BAPPExtNew[0]=*BAPPExternal; 
				LambdaAPP1ExtNew[0]=*LambAPP1External;
				DAPPExtNew[0]=*DAPPExternal; 
				LambdaAPP2ExtNew[0]=*LambAPP2External;
				PgNextAPPNew[0]=PgenNextAPP[(genID-1)];
				genSolverFirstBase.mainsolve( outerAPPIt, APPItCount, gsRho, Pgenprev, Pgenavg, Powerprice, Angpriceavg, Angavg, Angprice, PgenAPP, PgenAPPInner, PgNextAPPNew[0], BAPPExtNew[0], DAPPExtNew[0], LambAPP1ExtNew[0], LambAPP2ExtNew[0], BAPPNew, LambdaAPPNew ); // calls the Generator optimization solver
				Pg = genSolverFirstBase.getPSol(); // get the Generator Power iterate
				PgenNext = genSolverFirstBase.getPNextSol();
				PgenPrev = genSolverFirstBase.getPgPrev();
				Thetag = *(genSolverFirstBase.getThetaPtr());
			}
			if ((dispatchInterval!=0) && (flagLast==0)) { // For the first interval
				for (int counterCont = 0; counterCont < contCountGen; ++counterCont) {
					BAPPNew[counterCont]=BAPP[counterCont*numberOfGenerators+(genID-1)]; 
					LambdaAPPNew[counterCont]=LambAPP1[counterCont*numberOfGenerators+(genID-1)];
				}
				for (int counterCont = 0; counterCont <= contCountGen; ++counterCont) {
					BAPPExtNew[counterCont]=BAPPExternal[counterCont]; 
					LambdaAPP1ExtNew[counterCont]=LambAPP1External[counterCont];
					DAPPExtNew[counterCont]=DAPPExternal[counterCont]; 
					LambdaAPP2ExtNew[counterCont]=LambAPP2External[counterCont];
					PgNextAPPNew[counterCont]=PgenNextAPP[counterCont];
				}
				genSolverDZBase.mainsolve( outerAPPIt, APPItCount, gsRho, Pgenprev, Pgenavg, Powerprice, Angpriceavg, Angavg, Angprice, PgenAPP, PgenAPPInner, PgNextAPPNew, PgenPrevAPP, AAPPExternal, BAPPExtNew, DAPPExtNew, LambdaAPP1ExtNew, LambdaAPP2ExtNew, LambAPP3External, LambAPP4External, BAPPNew, LambdaAPPNew ); // calls the Generator optimization solver
				Pg = genSolverDZBase.getPSol(); // get the Generator Power iterate
				PgenNextPtr = genSolverDZBase.getPNextSol();
				PgenPrev = genSolverDZBase.getPPrevSol();
				Thetag = *(genSolverDZBase.getThetaPtr());
			}
			if ((dispatchInterval!=0) && (flagLast==1)) { // For the second (or, in this case, the last) interval
				for (int counterCont = 0; counterCont < contCountGen; ++counterCont) {
					BAPPNew[counterCont]=BAPP[counterCont*numberOfGenerators+(genID-1)]; 
					LambdaAPPNew[counterCont]=LambAPP1[counterCont*numberOfGenerators+(genID-1)];
				}
				BAPPExtNew[0]=-(*BAPPExternal);
				DAPPExtNew[0]=*DAPPExternal;
				genSolverSecondBase.mainsolve( outerAPPIt, APPItCount, gsRho, Pgenprev, Pgenavg, Powerprice, Angpriceavg, Angavg, Angprice, PgenAPP, PgenAPPInner, PgenPrevAPP, AAPPExternal, BAPPExtNew[0], LambAPP3External, LambAPP4External,  BAPPNew, LambdaAPPNew ); // calls the Generator optimization solver
				Pg = genSolverSecondBase.getPSol(); // get the Generator Power iterate
				PgenNext = genSolverSecondBase.getPNextSol();
				PgenPrev = genSolverSecondBase.getPPrevSol();
				Thetag = *(genSolverSecondBase.getThetaPtr());
			}
		}
		if ( dummyZeroIntFlag == 0 ) { // If dummy zero interval is not considered
			if ((dispatchInterval!=0) && (flagLast==0)) { // For the first interval
				for (int counterCont = 0; counterCont < contCountGen; ++counterCont) {
					BAPPNew[counterCont]=BAPP[counterCont*numberOfGenerators+(genID-1)]; 
					LambdaAPPNew[counterCont]=LambAPP1[counterCont*numberOfGenerators+(genID-1)];
				}
				for (int counterCont = 0; counterCont <= contCountGen; ++counterCont) {
					BAPPExtNew[counterCont]=BAPPExternal[counterCont*numberOfGenerators+(genID-1)]; 
					LambdaAPP1ExtNew[counterCont]=LambAPP1External[counterCont*numberOfGenerators+(genID-1)];
					DAPPExtNew[counterCont]=DAPPExternal[counterCont*numberOfGenerators+(genID-1)]; 
					LambdaAPP2ExtNew[counterCont]=LambAPP2External[counterCont*numberOfGenerators+(genID-1)];
				}
				genSolverFirst.mainsolve( outerAPPIt, APPItCount, gsRho, Pgenprev, Pgenavg, Powerprice, Angpriceavg, Angavg, Angprice, PgenAPP, PgenAPPInner, PgNextAPPNew, BAPPExtNew, DAPPExtNew, LambdaAPP1ExtNew, LambdaAPP2ExtNew, BAPPNew, LambdaAPPNew ); // calls the Generator optimization solver
				Pg = genSolverFirst.getPSol(); // get the Generator Power iterate
				PgenNextPtr = genSolverFirst.getPNextSol();
				PgenPrev = genSolverFirst.getPgPrev();
				Thetag = *(genSolverFirst.getThetaPtr());
			}
			if ((dispatchInterval!=0) && (flagLast==1)) { // For the second (or, in this case, the last) interval
				for (int counterCont = 0; counterCont < contCountGen; ++counterCont) {
					BAPPNew[counterCont]=BAPP[counterCont*numberOfGenerators+(genID-1)]; 
					LambdaAPPNew[counterCont]=LambAPP1[counterCont*numberOfGenerators+(genID-1)];
				}
				BAPPExtNew[0]=-(*BAPPExternal);
				DAPPExtNew[0]=*DAPPExternal;
				genSolverSecondBase.mainsolve( outerAPPIt, APPItCount, gsRho, Pgenprev, Pgenavg, Powerprice, Angpriceavg, Angavg, Angprice, PgenAPP, PgenAPPInner, PgenPrevAPP, AAPPExternal, BAPPExtNew[0], LambAPP3External, LambAPP4External,  BAPPNew, LambdaAPPNew ); // calls the Generator optimization solver
				Pg = genSolverSecondBase.getPSol(); // get the Generator Power iterate
				PgenNext = genSolverSecondBase.getPNextSol();
				PgenPrev = genSolverSecondBase.getPPrevSol();
				Thetag = *(genSolverSecondBase.getThetaPtr());
			}
		}
	}
	else { // If a contingency scenario
		if ( dummyZeroIntFlag == 1 ) { // If dummy zero interval is considered
			if ((dispatchInterval==0) && (flagLast==0)) { // For the dummy zeroth interval
				genSolverCont.mainsolve( APPItCount, gsRho, Pgenprev, Pgenavg, Powerprice, Angpriceavg, Angavg, Angprice, PgenAPPInner, -BAPP[(scenarioContCount-1)*numberOfGenerators+(genID-1)], LambAPP1[(scenarioContCount-1)*numberOfGenerators+(genID-1)] ); // calls the Generator optimization solver
				Pg = genSolverCont.getPSol(); // get the Generator Power iterate
				Thetag = *(genSolverCont.getThetaPtr());
				//*cout << "\nThiterate from generator: " << *(Thetag+i) << endl;
			}
			if ((dispatchInterval!=0) && (flagLast==0)) { // For the first interval
				if (contSolverAccuracy == 0) {// If the exhaustive calculation for contingency scenarios is not desired
					genSolverCont.mainsolve( APPItCount, gsRho, Pgenprev, Pgenavg, Powerprice, Angpriceavg, Angavg, Angprice, PgenAPPInner, -BAPP[(scenarioContCount-1)*numberOfGenerators+(genID-1)], LambAPP1[(scenarioContCount-1)*numberOfGenerators+(genID-1)] ); // calls the Generator optimization solver
					Pg = genSolverCont.getPSol(); // get the Generator Power iterate
					Thetag = *(genSolverCont.getThetaPtr());
					//*cout << "\nThiterate from generator: " << *(Thetag+i) << endl;
				}
				if (contSolverAccuracy != 0) {// If the exhaustive calculation for contingency scenarios not desired
					for (int counterCont = 0; counterCont <= contCountGen; ++counterCont) {
						BAPPExtNew[counterCont]=BAPPExternal[counterCont]; 
						LambdaAPP1ExtNew[counterCont]=LambAPP1External[counterCont];
						DAPPExtNew[counterCont]=DAPPExternal[counterCont]; 
						LambdaAPP2ExtNew[counterCont]=LambAPP2External[counterCont];
						PgNextAPPNew[counterCont]=PgenNextAPP[counterCont];
					}
					genSolverDZCont.mainsolve( gsRho, Pgenprev, Pgenavg, Powerprice, Angpriceavg, Angavg, Angprice, PgenAPP, PgenAPPInner, PgNextAPPNew, PgenPrevAPP, AAPPExternal, BAPPExtNew, DAPPExtNew, LambdaAPP1ExtNew, LambdaAPP2ExtNew, LambAPP3External, LambAPP4External, -BAPP[(scenarioContCount-1)*numberOfGenerators+(genID-1)], LambAPP1[(scenarioContCount-1)*numberOfGenerators+(genID-1)] ); // calls the Generator optimization solver
					Pg = genSolverDZCont.getPSol(); // get the Generator Power iterate
					Thetag = *(genSolverDZCont.getThetaPtr());
					//*cout << "\nThiterate from generator: " << *(Thetag+i) << endl;
				}
			}
			if ((dispatchInterval!=0) && (flagLast==1)) { // For the second (or, in this case, the last) interval
				if (contSolverAccuracy == 0) {// If the exhaustive calculation for contingency scenarios is not desired
					genSolverCont.mainsolve( APPItCount, gsRho, Pgenprev, Pgenavg, Powerprice, Angpriceavg, Angavg, Angprice, PgenAPPInner, -BAPP[(scenarioContCount-1)*numberOfGenerators+(genID-1)], LambAPP1[(scenarioContCount-1)*numberOfGenerators+(genID-1)] ); // calls the Generator optimization solver
					Pg = genSolverCont.getPSol(); // get the Generator Power iterate
					Thetag = *(genSolverCont.getThetaPtr());
					//*cout << "\nThiterate from generator: " << *(Thetag+i) << endl;
				}
				if (contSolverAccuracy != 0) {// If the exhaustive calculation for contingency scenarios is desired
					BAPPExtNew[0]=-(*BAPPExternal);
					DAPPExtNew[0]=*DAPPExternal;
					genSolverSecondCont.mainsolve( gsRho, Pgenprev, Pgenavg, Powerprice, Angpriceavg, Angavg, Angprice, PgenAPP, PgenAPPInner, PgenPrevAPP, AAPPExternal, BAPPExtNew[0], LambAPP3External, LambAPP4External, -BAPP[(scenarioContCount-1)*numberOfGenerators+(genID-1)], LambAPP1[(scenarioContCount-1)*numberOfGenerators+(genID-1)] ); // calls the Generator optimization solver
					Pg = genSolverSecondCont.getPSol(); // get the Generator Power iterate
					Thetag = *(genSolverSecondCont.getThetaPtr());
					//*cout << "\nThiterate from generator: " << *(Thetag+i) << endl;
				}
			}
		}
		if ( dummyZeroIntFlag == 0 ) { // If dummy zero interval is not considered
			if ((dispatchInterval==0) && (flagLast==0)) { // For the dummy zeroth interval **/ Will not be used in this case**/
				genSolverCont.mainsolve( APPItCount, gsRho, Pgenprev, Pgenavg, Powerprice, Angpriceavg, Angavg, Angprice, PgenAPPInner, -BAPP[(scenarioContCount-1)*numberOfGenerators+(genID-1)], LambAPP1[(scenarioContCount-1)*numberOfGenerators+(genID-1)] ); // calls the Generator optimization solver
				Pg = genSolverCont.getPSol(); // get the Generator Power iterate
				Thetag = *(genSolverCont.getThetaPtr());
				//*cout << "\nThiterate from generator: " << *(Thetag+i) << endl;
			}
			if ((dispatchInterval!=0) && (flagLast==0)) { // For the first interval
				if (contSolverAccuracy == 0) {// If the exhaustive calculation for contingency scenarios is not desired
					genSolverCont.mainsolve( APPItCount, gsRho, Pgenprev, Pgenavg, Powerprice, Angpriceavg, Angavg, Angprice, PgenAPPInner, -BAPP[(scenarioContCount-1)*numberOfGenerators+(genID-1)], LambAPP1[(scenarioContCount-1)*numberOfGenerators+(genID-1)] ); // calls the Generator optimization solver
					Pg = genSolverCont.getPSol(); // get the Generator Power iterate
					Thetag = *(genSolverCont.getThetaPtr());
					//*cout << "\nThiterate from generator: " << *(Thetag+i) << endl;
				}
				if (contSolverAccuracy != 0) {// If the exhaustive calculation for contingency scenarios not desired
					for (int counterCont = 0; counterCont <= contCountGen; ++counterCont) {
						BAPPExtNew[counterCont]=BAPPExternal[counterCont*numberOfGenerators+(genID-1)]; 
						LambdaAPP1ExtNew[counterCont]=LambAPP1External[counterCont*numberOfGenerators+(genID-1)];
						DAPPExtNew[counterCont]=DAPPExternal[counterCont*numberOfGenerators+(genID-1)]; 
						LambdaAPP2ExtNew[counterCont]=LambAPP2External[counterCont*numberOfGenerators+(genID-1)];
					}
					genSolverFirstCont.mainsolve( gsRho, Pgenprev, Pgenavg, Powerprice, Angpriceavg, Angavg, Angprice, PgenAPP, PgenAPPInner, PgNextAPPNew, BAPPExtNew, DAPPExtNew, LambdaAPP1ExtNew, LambdaAPP2ExtNew, -BAPP[(scenarioContCount-1)*numberOfGenerators+(genID-1)], LambAPP1[(scenarioContCount-1)*numberOfGenerators+(genID-1)] ); // calls the Generator optimization solver
					Pg = genSolverFirstCont.getPSol(); // get the Generator Power iterate
					Thetag = *(genSolverFirstCont.getThetaPtr());
					//*cout << "\nThiterate from generator: " << *(Thetag+i) << endl;
				}
			}
			if ((dispatchInterval!=0) && (flagLast==1)) { // For the second (or, in this case, the last) interval
				if (contSolverAccuracy == 0) {// If the exhaustive calculation for contingency scenarios is not desired
					genSolverCont.mainsolve( APPItCount, gsRho, Pgenprev, Pgenavg, Powerprice, Angpriceavg, Angavg, Angprice, PgenAPPInner, -BAPP[(scenarioContCount-1)*numberOfGenerators+(genID-1)], LambAPP1[(scenarioContCount-1)*numberOfGenerators+(genID-1)] ); // calls the Generator optimization solver
					Pg = genSolverCont.getPSol(); // get the Generator Power iterate
					Thetag = *(genSolverCont.getThetaPtr());
					//*cout << "\nThiterate from generator: " << *(Thetag+i) << endl;
				}
				if (contSolverAccuracy != 0) {// If the exhaustive calculation for contingency scenarios is desired
					BAPPExtNew[0]=-(*BAPPExternal);
					DAPPExtNew[0]=*DAPPExternal;
					genSolverSecondCont.mainsolve( gsRho, Pgenprev, Pgenavg, Powerprice, Angpriceavg, Angavg, Angprice, PgenAPP, PgenAPPInner, PgenPrevAPP, AAPPExternal, BAPPExtNew[0], LambAPP3External, LambAPP4External, -BAPP[(scenarioContCount-1)*numberOfGenerators+(genID-1)], LambAPP1[(scenarioContCount-1)*numberOfGenerators+(genID-1)] ); // calls the Generator optimization solver
					Pg = genSolverSecondCont.getPSol(); // get the Generator Power iterate
					Thetag = *(genSolverSecondCont.getThetaPtr());
					//*cout << "\nThiterate from generator: " << *(Thetag+i) << endl;
				}
			}
		}			
	}
	connNodegPtr->powerangleMessage( Pg, v, Thetag ); // passes to node object the corresponding iterates of power, angle, v, and number of scenarios
} // function gpowerangleMessage ends

void Generator::gpowerangleMessageGUROBI(int outerAPPIt, int  APPItCount, double gsRho, double Pprevit, double Pnetavg, double uprev, double vprevavg, double Aprevavg, double vprev, double PgenPrevAPP, double PgenAPP, double PgenAPPInner, double PgenNextAPP[], double AAPPExternal, double BAPPExternal[], double DAPPExternal[], double LambAPP1External[], double LambAPP2External[], double LambAPP3External, double LambAPP4External, double BAPP[], double LambAPP1[], GRBEnv* environmentGUROBI) //const // function gpowerangleMessage begins
{
	// CREATION OF THE MIP SOLVER INSTANCE //
        int dimRow; // Total number of rows of the A matrix (number of structural constraints of the QP): first term for the upper generation limit, the next term for the lower generation limit
        int dimCol; // Total number of columns of the QP (number of Decision Variables) first term to account for power generation MW outputs, second term for voltage phase angles for generation node
	// Base-Case with Dummy Zero Interval
	if ( ( dummyZeroIntFlag == 1 ) && ( baseContScenario == 0 ) && (dispatchInterval==0) && (flagLast==0) ){ 
        	dimRow = 6; // Total number of rows of the A matrix (number of structural constraints of the QP): first term for the upper generation limit, the next term for the lower generation limit
        	dimCol = 3; // Total number of columns of the QP (number of Decision Variables) first term to account for power generation MW outputs, second term for voltage phase angles for generation node
	}
	if ( ( dummyZeroIntFlag == 1 ) && ( baseContScenario == 0 ) && (dispatchInterval!=0) && (flagLast==0) ){ 
        	dimRow = 4+2*(1+contCountGen); // Total number of rows of the A matrix (number of structural constraints of the QP): first term for the upper generation limit, the next term for the lower generation limit
        	dimCol = 4+contCountGen; // Total number of columns of the QP (number of Decision Variables) first term to account for power generation MW outputs, second term for voltage phase angles for generation node
	}
	if ( ( dummyZeroIntFlag == 1 ) && ( baseContScenario == 0 ) && (dispatchInterval!=0) && (flagLast==1) ){ 
        	dimRow = 6; // Total number of rows of the A matrix (number of structural constraints of the QP): first term for the upper generation limit, the next term for the lower generation limit
        	dimCol = 3; // Total number of columns of the QP (number of Decision Variables) first term to account for power generation MW outputs, second term for voltage phase angles for generation node
	}
	// Base-Case without Dummy Zero Interval
	if ( ( dummyZeroIntFlag == 0 ) && ( baseContScenario == 0 ) && (dispatchInterval!=0) && (flagLast==0) ){ 
        	dimRow = 4+2*(1+contCountGen); // Total number of rows of the A matrix (number of structural constraints of the QP): first term for the upper generation limit, the next term for the lower generation limit
        	dimCol = 3+contCountGen; // Total number of columns of the QP (number of Decision Variables) first term to account for power generation MW outputs, second term for voltage phase angles for generation node
	}
	if ( ( dummyZeroIntFlag == 0 ) && ( baseContScenario == 0 ) && (dispatchInterval!=0) && (flagLast==1) ){ 
        	dimRow = 6; // Total number of rows of the A matrix (number of structural constraints of the QP): first term for the upper generation limit, the next term for the lower generation limit
        	dimCol = 3; // Total number of columns of the QP (number of Decision Variables) first term to account for power generation MW outputs, second term for voltage phase angles for generation node
	}
	// Contingency Scenario-Case
	if ( baseContScenario != 0 ) {
		if (contSolverAccuracy == 0){ 
        		dimRow = 2; // Total number of rows of the A matrix (number of structural constraints of the QP): first term for the upper generation limit, the next term for the lower generation limit
        		dimCol = 2; // Total number of columns of the QP (number of Decision Variables) first term to account for power generation MW outputs, second term for voltage phase angles for generation node
		}
		if ( (contSolverAccuracy == 1)&&( dummyZeroIntFlag == 1 ) && (dispatchInterval==0) && (flagLast==0)){ 
        		dimRow = 2; // Total number of rows of the A matrix (number of structural constraints of the QP): first term for the upper generation limit, the next term for the lower generation limit
        		dimCol = 2; // Total number of columns of the QP (number of Decision Variables) first term to account for power generation MW outputs, second term for voltage phase angles for generation node
		}
		if ( (contSolverAccuracy == 1)&&( dummyZeroIntFlag == 1 ) && (dispatchInterval!=0) && (flagLast==0)){ 
        		dimRow = 4+2*(1+contCountGen); // Total number of rows of the A matrix (number of structural constraints of the QP): first term for the upper generation limit, the next term for the lower generation limit
        		dimCol = 4+contCountGen; // Total number of columns of the QP (number of Decision Variables) first term to account for power generation MW outputs, second term for voltage phase angles for generation node
		}
		if ( (contSolverAccuracy == 1)&&( dummyZeroIntFlag == 1 ) && (dispatchInterval!=0) && (flagLast==1)){ 
        		dimRow = 6; // Total number of rows of the A matrix (number of structural constraints of the QP): first term for the upper generation limit, the next term for the lower generation limit
        		dimCol = 3; // Total number of columns of the QP (number of Decision Variables) first term to account for power generation MW outputs, second term for voltage phase angles for generation node
		}
		if ( (contSolverAccuracy == 1)&&( dummyZeroIntFlag == 0 ) && (dispatchInterval!=0) && (flagLast==0)){ 
        		dimRow = 4+2*(1+contCountGen); // Total number of rows of the A matrix (number of structural constraints of the QP): first term for the upper generation limit, the next term for the lower generation limit
        		dimCol = 3+contCountGen; // Total number of columns of the QP (number of Decision Variables) first term to account for power generation MW outputs, second term for voltage phase angles for generation node
		}
		if ( (contSolverAccuracy == 1)&&( dummyZeroIntFlag == 0 ) && (dispatchInterval!=0) && (flagLast==1)){ 
        		dimRow = 6; // Total number of rows of the A matrix (number of structural constraints of the QP): first term for the upper generation limit, the next term for the lower generation limit
        		dimCol = 3; // Total number of columns of the QP (number of Decision Variables) first term to account for power generation MW outputs, second term for voltage phase angles for generation node
		}
	}
	// Instantiate GUROBI Problem model
	GRBModel *modelGenQP = new GRBModel(*environmentGUROBI);
    	modelGenQP->set(GRB_StringAttr_ModelName, "assignment");
	modelGenQP->set(GRB_IntParam_OutputFlag, 0);
	GRBVar decvar[dimCol+1];
	double z; // variable to store the objective value

	// SPECIFICATION OF PROBLEM PARAMETERS //
	// Dummy Decision Variable //
	decvar[0] = modelGenQP->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
	//Decision Variable Definitions, Bounds, and Objective Function Co-efficients//
	int colCount = 1;
	if ( dummyZeroIntFlag == 1 ) {
		if ( ( baseContScenario == 0 ) && (dispatchInterval==0) ){ 
			//Columns corresponding to Power Generation continuous variables for different generators//
			decvar[colCount] = modelGenQP->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
			++colCount;
			//Columns corresponding to Power Generation continuous variables for different generators for next interval//
			decvar[colCount] = modelGenQP->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
			++colCount;
			//Columns corresponding to Voltage Phase Angles continuous variables for different nodes//
			decvar[colCount] = modelGenQP->addVar((0), (44/7), 0.0, GRB_CONTINUOUS);
		}
		if ( ( baseContScenario == 0 ) && (dispatchInterval!=0) && (flagLast==0) ){ 
			//Columns corresponding to Power Generation continuous variables for different generators//
			decvar[colCount] = modelGenQP->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
			++colCount;
			for (int counterCont = 0; counterCont <= contCountGen; ++counterCont) {
				//Columns corresponding to Power Generation continuous variables for different generators for next interval//
				decvar[colCount] = modelGenQP->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
				++colCount;
			}
			//Columns corresponding to Power Generation continuous variables for different generators for previous interval//
			decvar[colCount] = modelGenQP->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
			++colCount;
			//Columns corresponding to Voltage Phase Angles continuous variables for different nodes//
			decvar[colCount] = modelGenQP->addVar((0), (44/7), 0.0, GRB_CONTINUOUS);
		}
		if ( ( baseContScenario == 0 ) && (dispatchInterval!=0) && (flagLast==1) ){ 
			//Columns corresponding to Power Generation continuous variables for different generators//
			decvar[colCount] = modelGenQP->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
			++colCount;
			//Columns corresponding to Power Generation continuous variables for different generators for previous interval//
			decvar[colCount] = modelGenQP->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
			++colCount;
			//Columns corresponding to Voltage Phase Angles continuous variables for different nodes//
			decvar[colCount] = modelGenQP->addVar((0), (44/7), 0.0, GRB_CONTINUOUS);
		}
	}
	if ( dummyZeroIntFlag == 0 ) {
		if ( ( baseContScenario == 0 ) && (dispatchInterval!=0) && (flagLast==0) ){ 
			//Columns corresponding to Power Generation continuous variables for different generators//
			decvar[colCount] = modelGenQP->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
			++colCount;
			for (int counterCont = 0; counterCont <= contCountGen; ++counterCont) {
				//Columns corresponding to Power Generation continuous variables for different generators for next interval//
				decvar[colCount] = modelGenQP->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
				++colCount;
			}
			//Columns corresponding to Voltage Phase Angles continuous variables for different nodes//
			decvar[colCount] = modelGenQP->addVar((0), (44/7), 0.0, GRB_CONTINUOUS);
		}
		if ( ( baseContScenario == 0 ) && (dispatchInterval!=0) && (flagLast==1) ){ 
			//Columns corresponding to Power Generation continuous variables for different generators//
			decvar[colCount] = modelGenQP->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
			++colCount;
			//Columns corresponding to Power Generation continuous variables for different generators for previous interval//
			decvar[colCount] = modelGenQP->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
			++colCount;
			//Columns corresponding to Voltage Phase Angles continuous variables for different nodes//
			decvar[colCount] = modelGenQP->addVar((0), (44/7), 0.0, GRB_CONTINUOUS);
		}
	}
	if (( baseContScenario != 0 )&&(contSolverAccuracy == 0)){ 
		//Columns corresponding to Power Generation continuous variables for different generators//
		decvar[colCount] = modelGenQP->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
		++colCount;
		//Columns corresponding to Voltage Phase Angles continuous variables for different nodes//
		decvar[colCount] = modelGenQP->addVar((0), (44/7), 0.0, GRB_CONTINUOUS);
	}
	if (( baseContScenario != 0 )&&(contSolverAccuracy == 1)){ 
		if ( dummyZeroIntFlag == 1 ) {
			if (dispatchInterval==0) { 
				//Columns corresponding to Power Generation continuous variables for different generators//
				decvar[colCount] = modelGenQP->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
				++colCount;
				//Columns corresponding to Voltage Phase Angles continuous variables for different nodes//
				decvar[colCount] = modelGenQP->addVar((0), (44/7), 0.0, GRB_CONTINUOUS);
			}
			if ( (dispatchInterval!=0) && (flagLast==0) ){ 
				//Columns corresponding to Power Generation continuous variables for different generators//
				decvar[colCount] = modelGenQP->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
				++colCount;
				for (int counterCont = 0; counterCont <= contCountGen; ++counterCont) {
					//Columns corresponding to Power Generation continuous variables for different generators for next interval//
					decvar[colCount] = modelGenQP->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
					++colCount;
				}
				//Columns corresponding to Power Generation continuous variables for different generators for previous interval//
				decvar[colCount] = modelGenQP->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
				++colCount;
				//Columns corresponding to Voltage Phase Angles continuous variables for different nodes//
				decvar[colCount] = modelGenQP->addVar((0), (44/7), 0.0, GRB_CONTINUOUS);
			}
			if ( (dispatchInterval!=0) && (flagLast==1) ){ 
				//Columns corresponding to Power Generation continuous variables for different generators//
				decvar[colCount] = modelGenQP->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
				++colCount;
				//Columns corresponding to Power Generation continuous variables for different generators for previous interval//
				decvar[colCount] = modelGenQP->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
				++colCount;
				//Columns corresponding to Voltage Phase Angles continuous variables for different nodes//
				decvar[colCount] = modelGenQP->addVar((0), (44/7), 0.0, GRB_CONTINUOUS);
			}
		}
		if ( dummyZeroIntFlag == 0 ) {
			if ( (dispatchInterval!=0) && (flagLast==0) ){ 
				//Columns corresponding to Power Generation continuous variables for different generators//
				decvar[colCount] = modelGenQP->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
				++colCount;
				for (int counterCont = 0; counterCont <= contCountGen; ++counterCont) {
					//Columns corresponding to Power Generation continuous variables for different generators for next interval//
					decvar[colCount] = modelGenQP->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
					++colCount;
				}
				//Columns corresponding to Voltage Phase Angles continuous variables for different nodes//
				decvar[colCount] = modelGenQP->addVar((0), (44/7), 0.0, GRB_CONTINUOUS);
			}
			if ( (dispatchInterval!=0) && (flagLast==1) ){ 
				//Columns corresponding to Power Generation continuous variables for different generators//
				decvar[colCount] = modelGenQP->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
				++colCount;
				//Columns corresponding to Power Generation continuous variables for different generators for previous interval//
				decvar[colCount] = modelGenQP->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
				++colCount;
				//Columns corresponding to Voltage Phase Angles continuous variables for different nodes//
				decvar[colCount] = modelGenQP->addVar((0), (44/7), 0.0, GRB_CONTINUOUS);
			}
		}
	}
	//Setting Objective//
	GRBQuadExpr obj = 0.0;
	// Objective Contribution from Dummy Decision Variable //
	obj += 0*(decvar[0]);
	colCount = 1;
	double BAPPNew[contCountGen];
	double LambdaAPPNew[contCountGen];
	double BAPPExtNew[contCountGen+1];
	double DAPPExtNew[contCountGen+1];
	double LambdaAPP1ExtNew[contCountGen+1];
	double LambdaAPP2ExtNew[contCountGen+1];
	double PgNextAPPNew[contCountGen+1];
	double Lambda1Sum = 0;
	double PPresB = 0; 
	double BAPPSum = 0;
	double LambdaAPPSum = 0;
	for (int counterCont = 0; counterCont < contCountGen; ++counterCont) {
		BAPPNew[counterCont]=0; 
		LambdaAPPNew[counterCont]=0;
	}
	for (int counterCont = 0; counterCont <= contCountGen; ++counterCont) {
		BAPPExtNew[counterCont]=0; 
		LambdaAPP1ExtNew[counterCont]=0;
		DAPPExtNew[counterCont]=0; 
		LambdaAPP2ExtNew[counterCont]=0;
		PgNextAPPNew[counterCont]=0;
	}
	//Columns corresponding to Power Generation continuous variables for different generators//
	if ( baseContScenario == 0 ) { // Base-Case
		if ( dummyZeroIntFlag == 1 ) { // If there is a dummy-zero interval at the beginning
			if (dispatchInterval==0) { // The dummy zero interval
				for (int counterCont = 0; counterCont < contCountGen; ++counterCont) { // Disagreements and Lagrange Multipliers for contingency scenarios
					BAPPNew[counterCont]=BAPP[counterCont*numberOfGenerators+(genID-1)]; 
					LambdaAPPNew[counterCont]=LambAPP1[counterCont*numberOfGenerators+(genID-1)];
					BAPPSum += BAPPNew[counterCont];
					LambdaAPPSum += LambdaAPPNew[counterCont];
				}
				// Disagreements and Lagrange Multipliers for post-contingency base-cases
				BAPPExtNew[0]=*BAPPExternal; 
				LambdaAPP1ExtNew[0]=*LambAPP1External;
				DAPPExtNew[0]=*DAPPExternal; 
				LambdaAPP2ExtNew[0]=*LambAPP2External;
				PgNextAPPNew[0]=PgenNextAPP[(genID-1)];
				// Objective function component for present estimate of power
				obj += (genSolverFirstBase.getQuadCoeff())*(decvar[colCount])*(decvar[colCount])+(genSolverFirstBase.getLinCoeff())*(decvar[colCount])+(genSolverFirstBase.getConstCoeff())+((genSolverFirstBase.getBeta())/2)*((decvar[colCount])-PgenAPP)*((decvar[colCount])-PgenAPP)+((genSolverFirstBase.getIntBeta())/2)*((decvar[colCount])-PgenAPPInner)*((decvar[colCount])-PgenAPPInner)+(genSolverFirstBase.getIntGamma())*((decvar[colCount])*BAPPSum)+(decvar[colCount])*LambdaAPPSum+(genSolverFirstBase.getGamma())*((decvar[colCount])*(BAPPExtNew[0]))+(decvar[colCount])*(LambdaAPP1ExtNew[0])+(gsRho/2)*(decvar[colCount]-Pprevit+Pnetavg+uprev)*(decvar[colCount]-Pprevit+Pnetavg+uprev);
				++colCount;
				// Objective function component for future estimate of power
				obj += ((genSolverFirstBase.getBeta())/2)*((decvar[colCount])-PgNextAPPNew[0])*((decvar[colCount])-PgNextAPPNew[0])+(genSolverFirstBase.getGamma())*((decvar[colCount])*(DAPPExtNew[0]))+(decvar[colCount])*(LambdaAPP2ExtNew[0]);
			}
			if ( (dispatchInterval!=0) && (flagLast==0) ){ // The first dispatch interval in consideration
				for (int counterCont = 0; counterCont < contCountGen; ++counterCont) { // Disagreements and Lagrange Multipliers for contingency scenarios
					BAPPNew[counterCont]=BAPP[counterCont*numberOfGenerators+(genID-1)]; 
					LambdaAPPNew[counterCont]=LambAPP1[counterCont*numberOfGenerators+(genID-1)];
					BAPPSum += BAPPNew[counterCont];
					LambdaAPPSum += LambdaAPPNew[counterCont];
				}
				for (int counterCont = 0; counterCont <= contCountGen; ++counterCont) { // Disagreements and Lagrange Multipliers for post-contingency base-cases
					BAPPExtNew[counterCont]=BAPPExternal[counterCont]; 
					LambdaAPP1ExtNew[counterCont]=LambAPP1External[counterCont];
					DAPPExtNew[counterCont]=DAPPExternal[counterCont]; 
					LambdaAPP2ExtNew[counterCont]=LambAPP2External[counterCont];
					PgNextAPPNew[counterCont]=PgenNextAPP[counterCont];
					Lambda1Sum += LambdaAPP1ExtNew[counterCont];
					PPresB += BAPPExtNew[counterCont]; 
				}
				// Objective function component for present estimate of power
				obj += (genSolverDZBase.getQuadCoeff())*(decvar[colCount])*(decvar[colCount])+(genSolverDZBase.getLinCoeff())*(decvar[colCount])+(genSolverDZBase.getConstCoeff())+((genSolverDZBase.getBeta())/2)*((decvar[colCount])-PgenAPP)*((decvar[colCount])-PgenAPP)+((genSolverDZBase.getIntBeta())/2)*((decvar[colCount])-PgenAPPInner)*((decvar[colCount])-PgenAPPInner)+(genSolverDZBase.getIntGamma())*((decvar[colCount])*BAPPSum)+(decvar[colCount])*LambdaAPPSum+(genSolverDZBase.getGamma())*((decvar[colCount])*PPresB)+(decvar[colCount])*(Lambda1Sum-LambAPP4External)+(gsRho/2)*(decvar[colCount]-Pprevit+Pnetavg+uprev)*(decvar[colCount]-Pprevit+Pnetavg+uprev);
				++colCount;
				// Objective function component for future estimate of power
				for (int counterCont = 0; counterCont <= contCountGen; ++counterCont) {
					obj += ((genSolverDZBase.getBeta())/2)*((decvar[colCount])-PgNextAPPNew[counterCont])*((decvar[colCount])-PgNextAPPNew[counterCont])+(genSolverDZBase.getGamma())*((decvar[colCount])*DAPPExtNew[counterCont])+(decvar[colCount])*(LambdaAPP2ExtNew[counterCont]);
					++colCount;
				}
				// Objective function component for previous interval estimate of power
				obj += ((genSolverDZBase.getBeta())/2)*((decvar[colCount])-PgenPrevAPP)*((decvar[colCount])-PgenPrevAPP)+(genSolverDZBase.getGamma())*((decvar[colCount])*AAPPExternal)-(decvar[colCount])*LambAPP3External;
			}
			if ( (dispatchInterval!=0) && (flagLast==1) ){ // The second or last dispatch interval in consideration
				for (int counterCont = 0; counterCont < contCountGen; ++counterCont) { // Disagreements and Lagrange Multipliers for contingency scenarios
					BAPPNew[counterCont]=BAPP[counterCont*numberOfGenerators+(genID-1)]; 
					LambdaAPPNew[counterCont]=LambAPP1[counterCont*numberOfGenerators+(genID-1)];
					BAPPSum += BAPPNew[counterCont];
					LambdaAPPSum += LambdaAPPNew[counterCont];
				}
				// Disagreements and Lagrange Multipliers for post-contingency base-cases
				BAPPExtNew[0]=-(*BAPPExternal);
				DAPPExtNew[0]=*DAPPExternal;
				// Objective function component for present estimate of power
				obj += (genSolverSecondBase.getQuadCoeff())*(decvar[colCount])*(decvar[colCount])+(genSolverSecondBase.getLinCoeff())*(decvar[colCount])+(genSolverSecondBase.getConstCoeff())+((genSolverSecondBase.getBeta())/2)*((decvar[colCount])-PgenAPP)*((decvar[colCount])-PgenAPP)+((genSolverSecondBase.getIntBeta())/2)*((decvar[colCount])-PgenAPPInner)*((decvar[colCount])-PgenAPPInner)+(genSolverSecondBase.getIntGamma())*((decvar[colCount])*BAPPSum)+(decvar[colCount])*LambdaAPPSum+(genSolverSecondBase.getGamma())*((decvar[colCount])*(BAPPExtNew[0]))+(decvar[colCount])*(-LambAPP4External)+(gsRho/2)*(decvar[colCount]-Pprevit+Pnetavg+uprev)*(decvar[colCount]-Pprevit+Pnetavg+uprev);
				++colCount;
				// Objective function component for previous interval estimate of power
				obj += ((genSolverSecondBase.getBeta())/2)*((decvar[colCount])-PgenPrevAPP)*((decvar[colCount])-PgenPrevAPP)+(genSolverSecondBase.getGamma())*((decvar[colCount])*AAPPExternal)-(decvar[colCount])*LambAPP3External;
			}	
		}
		if ( dummyZeroIntFlag == 0 ) { // If there is no dummy-zero interval at the beginning
			if ( (dispatchInterval!=0) && (flagLast==0) ){ // The first dispatch interval in consideration 
				for (int counterCont = 0; counterCont < contCountGen; ++counterCont) { // Disagreements and Lagrange Multipliers for contingency scenarios
					BAPPNew[counterCont]=BAPP[counterCont*numberOfGenerators+(genID-1)]; 
					LambdaAPPNew[counterCont]=LambAPP1[counterCont*numberOfGenerators+(genID-1)];
					BAPPSum += BAPPNew[counterCont];
					LambdaAPPSum += LambdaAPPNew[counterCont];
				}
				for (int counterCont = 0; counterCont <= contCountGen; ++counterCont) { // Disagreements and Lagrange Multipliers for post-contingency base-cases
					BAPPExtNew[counterCont]=BAPPExternal[counterCont]; 
					LambdaAPP1ExtNew[counterCont]=LambAPP1External[counterCont];
					DAPPExtNew[counterCont]=DAPPExternal[counterCont]; 
					LambdaAPP2ExtNew[counterCont]=LambAPP2External[counterCont];
					PgNextAPPNew[counterCont]=PgenNextAPP[counterCont];
					Lambda1Sum += LambdaAPP1ExtNew[counterCont];
					PPresB += BAPPExtNew[counterCont]; 
				}
				// Objective function component for present estimate of power
				obj += (genSolverFirst.getQuadCoeff())*(decvar[colCount])*(decvar[colCount])+(genSolverFirst.getLinCoeff())*(decvar[colCount])+(genSolverFirst.getConstCoeff())+((genSolverFirst.getBeta())/2)*((decvar[colCount])-PgenAPP)*((decvar[colCount])-PgenAPP)+((genSolverFirst.getIntBeta())/2)*((decvar[colCount])-PgenAPPInner)*((decvar[colCount])-PgenAPPInner)+(genSolverFirst.getIntGamma())*((decvar[colCount])*BAPPSum)+(decvar[colCount])*LambdaAPPSum+(genSolverFirst.getGamma())*((decvar[colCount])*PPresB)+(decvar[colCount])*(Lambda1Sum-LambAPP4External)+(gsRho/2)*(decvar[colCount]-Pprevit+Pnetavg+uprev)*(decvar[colCount]-Pprevit+Pnetavg+uprev);
				// Objective function component for future estimate of power
				for (int counterCont = 0; counterCont <= contCountGen; ++counterCont) {
					++colCount;
					obj += ((genSolverFirst.getBeta())/2)*((decvar[colCount])-PgNextAPPNew[counterCont])*((decvar[colCount])-PgNextAPPNew[counterCont])+(genSolverFirst.getGamma())*((decvar[colCount])*DAPPExtNew[counterCont])+(decvar[colCount])*(LambdaAPP2ExtNew[counterCont]);
				}
			}
			if ( (dispatchInterval!=0) && (flagLast==1) ){ 
				for (int counterCont = 0; counterCont < contCountGen; ++counterCont) { // Disagreements and Lagrange Multipliers for contingency scenarios
					BAPPNew[counterCont]=BAPP[counterCont*numberOfGenerators+(genID-1)]; 
					LambdaAPPNew[counterCont]=LambAPP1[counterCont*numberOfGenerators+(genID-1)];
					BAPPSum += BAPPNew[counterCont];
					LambdaAPPSum += LambdaAPPNew[counterCont];
				}
				// Disagreements and Lagrange Multipliers for post-contingency base-cases
				BAPPExtNew[0]=-(*BAPPExternal);
				DAPPExtNew[0]=*DAPPExternal;
				// Objective function component for present estimate of power
				obj += (genSolverSecondBase.getQuadCoeff())*(decvar[colCount])*(decvar[colCount])+(genSolverSecondBase.getLinCoeff())*(decvar[colCount])+(genSolverSecondBase.getConstCoeff())+((genSolverSecondBase.getBeta())/2)*((decvar[colCount])-PgenAPP)*((decvar[colCount])-PgenAPP)+((genSolverSecondBase.getIntBeta())/2)*((decvar[colCount])-PgenAPPInner)*((decvar[colCount])-PgenAPPInner)+(genSolverSecondBase.getIntGamma())*((decvar[colCount])*BAPPSum)+(decvar[colCount])*LambdaAPPSum+(genSolverSecondBase.getGamma())*((decvar[colCount])*(BAPPExtNew[0]))+(decvar[colCount])*(-LambAPP4External)+(gsRho/2)*(decvar[colCount]-Pprevit+Pnetavg+uprev)*(decvar[colCount]-Pprevit+Pnetavg+uprev);
				++colCount;
				// Objective function component for previous interval estimate of power
				obj += ((genSolverSecondBase.getBeta())/2)*((decvar[colCount])-PgenPrevAPP)*((decvar[colCount])-PgenPrevAPP)+(genSolverSecondBase.getGamma())*((decvar[colCount])*AAPPExternal)-(decvar[colCount])*LambAPP3External;
			}	
		}
	}
	if ( baseContScenario != 0 ) { // Contingency Scenarios
		if (contSolverAccuracy == 0) {// If the exhaustive calculation for contingency scenarios is not desired
			// Use the solver for first dispatch interval
			BAPPSum = -BAPP[(scenarioContCount-1)*numberOfGenerators+(genID-1)]; 
			LambdaAPPSum = LambAPP1[(scenarioContCount-1)*numberOfGenerators+(genID-1)];
			obj += (genSolverCont.getQuadCoeff())*(decvar[colCount])*(decvar[colCount])+(genSolverCont.getLinCoeff())*(decvar[colCount])+(genSolverCont.getConstCoeff())+((genSolverCont.getBeta())/2)*((decvar[colCount])-PgenAPPInner)*((decvar[colCount])-PgenAPPInner)+(genSolverCont.getGamma())*((decvar[colCount])*BAPPSum)-(decvar[colCount])*LambdaAPPSum+(gsRho/2)*(decvar[colCount]-Pprevit+Pnetavg+uprev)*(decvar[colCount]-Pprevit+Pnetavg+uprev);
		}
		else  { // For exhaustive calculation for contingency scenarios 
			if ( dummyZeroIntFlag == 1 ) { // If there is a dummy-zero interval at the beginning
				if (dispatchInterval==0) { // The dummy zero interval
					// Use the solver for first dispatch interval
					BAPPSum = -BAPP[(scenarioContCount-1)*numberOfGenerators+(genID-1)]; 
					LambdaAPPSum = LambAPP1[(scenarioContCount-1)*numberOfGenerators+(genID-1)];
					// Objective function component for present estimate of power
					obj += (genSolverCont.getQuadCoeff())*(decvar[colCount])*(decvar[colCount])+(genSolverCont.getLinCoeff())*(decvar[colCount])+(genSolverCont.getConstCoeff())+((genSolverCont.getBeta())/2)*((decvar[colCount])-PgenAPPInner)*((decvar[colCount])-PgenAPPInner)+(genSolverCont.getGamma())*((decvar[colCount])*BAPPSum)-(decvar[colCount])*LambdaAPPSum+(gsRho/2)*(decvar[colCount]-Pprevit+Pnetavg+uprev)*(decvar[colCount]-Pprevit+Pnetavg+uprev);
				}
				if ( (dispatchInterval!=0) && (flagLast==0) ){ // The first dispatch interval in consideration
					// Use the solver for first dispatch interval
					BAPPSum = -BAPP[(scenarioContCount-1)*numberOfGenerators+(genID-1)]; 
					LambdaAPPSum = LambAPP1[(scenarioContCount-1)*numberOfGenerators+(genID-1)];
					for (int counterCont = 0; counterCont <= contCountGen; ++counterCont) { // Disagreements and Lagrange Multipliers for post-contingency base-cases
						BAPPExtNew[counterCont]=BAPPExternal[counterCont]; 
						LambdaAPP1ExtNew[counterCont]=LambAPP1External[counterCont];
						DAPPExtNew[counterCont]=DAPPExternal[counterCont]; 
						LambdaAPP2ExtNew[counterCont]=LambAPP2External[counterCont];
						PgNextAPPNew[counterCont]=PgenNextAPP[counterCont];
						Lambda1Sum += LambdaAPP1ExtNew[counterCont];
						PPresB += BAPPExtNew[counterCont]; 
					}
					// Objective function component for present estimate of power
					obj += (genSolverDZBase.getQuadCoeff())*(decvar[colCount])*(decvar[colCount])+(genSolverDZBase.getLinCoeff())*(decvar[colCount])+(genSolverDZBase.getConstCoeff())+((genSolverDZBase.getBeta())/2)*((decvar[colCount])-PgenAPP)*((decvar[colCount])-PgenAPP)+((genSolverDZBase.getIntBeta())/2)*((decvar[colCount])-PgenAPPInner)*((decvar[colCount])-PgenAPPInner)+(genSolverDZBase.getIntGamma())*((decvar[colCount])*BAPPSum)+(decvar[colCount])*LambdaAPPSum+(genSolverDZBase.getGamma())*((decvar[colCount])*PPresB)+(decvar[colCount])*(Lambda1Sum-LambAPP4External)+(gsRho/2)*(decvar[colCount]-Pprevit+Pnetavg+uprev)*(decvar[colCount]-Pprevit+Pnetavg+uprev);
					++colCount;
					// Objective function component for future estimate of power
					for (int counterCont = 0; counterCont <= contCountGen; ++counterCont) {
						obj += ((genSolverDZBase.getBeta())/2)*((decvar[colCount])-PgNextAPPNew[counterCont])*((decvar[colCount])-PgNextAPPNew[counterCont])+(genSolverDZBase.getGamma())*((decvar[colCount])*DAPPExtNew[counterCont])+(decvar[colCount])*(LambdaAPP2ExtNew[counterCont]);
						++colCount;
					}
					// Objective function component for previous interval estimate of power
					obj += ((genSolverDZBase.getBeta())/2)*((decvar[colCount])-PgenPrevAPP)*((decvar[colCount])-PgenPrevAPP)+(genSolverDZBase.getGamma())*((decvar[colCount])*AAPPExternal)-(decvar[colCount])*LambAPP3External;
				}
				if ( (dispatchInterval!=0) && (flagLast==1) ){ // The second or last dispatch interval in consideration
					// Use the solver for first dispatch interval
					BAPPSum = -BAPP[(scenarioContCount-1)*numberOfGenerators+(genID-1)]; 
					LambdaAPPSum = LambAPP1[(scenarioContCount-1)*numberOfGenerators+(genID-1)];
					// Disagreements and Lagrange Multipliers for post-contingency base-cases
					BAPPExtNew[0]=-(*BAPPExternal);
					DAPPExtNew[0]=*DAPPExternal;
					// Objective function component for present estimate of power
					obj += (genSolverSecondBase.getQuadCoeff())*(decvar[colCount])*(decvar[colCount])+(genSolverSecondBase.getLinCoeff())*(decvar[colCount])+(genSolverSecondBase.getConstCoeff())+((genSolverSecondBase.getBeta())/2)*((decvar[colCount])-PgenAPP)*((decvar[colCount])-PgenAPP)+((genSolverSecondBase.getIntBeta())/2)*((decvar[colCount])-PgenAPPInner)*((decvar[colCount])-PgenAPPInner)+(genSolverSecondBase.getIntGamma())*((decvar[colCount])*BAPPSum)+(decvar[colCount])*LambdaAPPSum+(genSolverSecondBase.getGamma())*((decvar[colCount])*(BAPPExtNew[0]))+(decvar[colCount])*(-LambAPP4External)+(gsRho/2)*(decvar[colCount]-Pprevit+Pnetavg+uprev)*(decvar[colCount]-Pprevit+Pnetavg+uprev);
					++colCount;
					// Objective function component for previous interval estimate of power
					obj += ((genSolverSecondBase.getBeta())/2)*((decvar[colCount])-PgenPrevAPP)*((decvar[colCount])-PgenPrevAPP)+(genSolverSecondBase.getGamma())*((decvar[colCount])*AAPPExternal)-(decvar[colCount])*LambAPP3External;
				}	
			}
			if ( dummyZeroIntFlag == 0 ) { // If there is no dummy-zero interval at the beginning
				if ( (dispatchInterval!=0) && (flagLast==0) ){ // The first dispatch interval in consideration 
					// Use the solver for first dispatch interval
					BAPPSum = -BAPP[(scenarioContCount-1)*numberOfGenerators+(genID-1)]; 
					LambdaAPPSum = LambAPP1[(scenarioContCount-1)*numberOfGenerators+(genID-1)];
					for (int counterCont = 0; counterCont <= contCountGen; ++counterCont) { // Disagreements and Lagrange Multipliers for post-contingency base-cases
						BAPPExtNew[counterCont]=BAPPExternal[counterCont]; 
						LambdaAPP1ExtNew[counterCont]=LambAPP1External[counterCont];
						DAPPExtNew[counterCont]=DAPPExternal[counterCont]; 
						LambdaAPP2ExtNew[counterCont]=LambAPP2External[counterCont];
						PgNextAPPNew[counterCont]=PgenNextAPP[counterCont];
						Lambda1Sum += LambdaAPP1ExtNew[counterCont];
						PPresB += BAPPExtNew[counterCont]; 
					}
					// Objective function component for present estimate of power
					obj += (genSolverFirst.getQuadCoeff())*(decvar[colCount])*(decvar[colCount])+(genSolverFirst.getLinCoeff())*(decvar[colCount])+(genSolverFirst.getConstCoeff())+((genSolverFirst.getBeta())/2)*((decvar[colCount])-PgenAPP)*((decvar[colCount])-PgenAPP)+((genSolverFirst.getIntBeta())/2)*((decvar[colCount])-PgenAPPInner)*((decvar[colCount])-PgenAPPInner)+(genSolverFirst.getIntGamma())*((decvar[colCount])*BAPPSum)+(decvar[colCount])*LambdaAPPSum+(genSolverFirst.getGamma())*((decvar[colCount])*PPresB)+(decvar[colCount])*(Lambda1Sum-LambAPP4External)+(gsRho/2)*(decvar[colCount]-Pprevit+Pnetavg+uprev)*(decvar[colCount]-Pprevit+Pnetavg+uprev);
					// Objective function component for future estimate of power
					for (int counterCont = 0; counterCont <= contCountGen; ++counterCont) {
						++colCount;
						obj += ((genSolverFirst.getBeta())/2)*((decvar[colCount])-PgNextAPPNew[counterCont])*((decvar[colCount])-PgNextAPPNew[counterCont])+(genSolverFirst.getGamma())*((decvar[colCount])*DAPPExtNew[counterCont])+(decvar[colCount])*(LambdaAPP2ExtNew[counterCont]);
					}
				}
				if ( (dispatchInterval!=0) && (flagLast==1) ){ 
					// Use the solver for first dispatch interval
					BAPPSum = -BAPP[(scenarioContCount-1)*numberOfGenerators+(genID-1)]; 
					LambdaAPPSum = LambAPP1[(scenarioContCount-1)*numberOfGenerators+(genID-1)];
					// Disagreements and Lagrange Multipliers for post-contingency base-cases
					BAPPExtNew[0]=-(*BAPPExternal);
					DAPPExtNew[0]=*DAPPExternal;
					// Objective function component for present estimate of power
					obj += (genSolverSecondBase.getQuadCoeff())*(decvar[colCount])*(decvar[colCount])+(genSolverSecondBase.getLinCoeff())*(decvar[colCount])+(genSolverSecondBase.getConstCoeff())+((genSolverSecondBase.getBeta())/2)*((decvar[colCount])-PgenAPP)*((decvar[colCount])-PgenAPP)+((genSolverSecondBase.getIntBeta())/2)*((decvar[colCount])-PgenAPPInner)*((decvar[colCount])-PgenAPPInner)+(genSolverSecondBase.getIntGamma())*((decvar[colCount])*BAPPSum)+(decvar[colCount])*LambdaAPPSum+(genSolverSecondBase.getGamma())*((decvar[colCount])*(BAPPExtNew[0]))+(decvar[colCount])*(-LambAPP4External)+(gsRho/2)*(decvar[colCount]-Pprevit+Pnetavg+uprev)*(decvar[colCount]-Pprevit+Pnetavg+uprev);
					++colCount;
					// Objective function component for previous interval estimate of power
					obj += ((genSolverSecondBase.getBeta())/2)*((decvar[colCount])-PgenPrevAPP)*((decvar[colCount])-PgenPrevAPP)+(genSolverSecondBase.getGamma())*((decvar[colCount])*AAPPExternal)-(decvar[colCount])*LambAPP3External;
				}	
			}			
		}
	}
	++colCount;
	//Columns corresponding to Voltage Phase Angles continuous variables for different nodes//
	obj += (gsRho/2)*(decvar[colCount]-vprevavg-Aprevavg+vprev )*(decvar[colCount]-vprevavg-Aprevavg+vprev );

	modelGenQP->setObjective(obj, GRB_MINIMIZE);
	//Row Definitions: Specification of b<=Ax<=b//
	GRBLinExpr lhs[dimRow+1];
	//Row Definitions and Bounds Corresponding to Constraints/
	//Non-Zero entries of A matrix (Constraint/Coefficient matrix entries)//
	// Dummy Constraint //
	lhs[0] = 0*(decvar[0]);
	modelGenQP->addConstr(lhs[0], GRB_EQUAL, 0);
	int rCount = 1; // Initialize the row count
	if ( baseContScenario == 0 ) { // Base-Case
		if ( dummyZeroIntFlag == 1 ) { // If there is a dummy-zero interval at the beginning
			if (dispatchInterval==0) { // The dummy zero interval
				// Coefficients corresponding to lower generation limits
				lhs[rCount] = 0;
				lhs[rCount] += decvar[rCount];
				modelGenQP->addConstr(lhs[rCount] >= (genSolverFirstBase.getPMin()));
				++rCount; // Increment the row count to point to the next generator object
				// Coefficients corresponding to upper generation limits
				lhs[rCount] = 0;
				lhs[rCount] += decvar[rCount - 1];
				modelGenQP->addConstr(lhs[rCount] <= (genSolverFirstBase.getPMax()));
				++rCount; // Increment the row count to point to the next generator object
				// Coefficients corresponding to lower ramp limits for next interval
				lhs[rCount] = 0;
				lhs[rCount] += (decvar[rCount-1]-decvar[rCount-2]);
				modelGenQP->addConstr(lhs[rCount] >= (genSolverFirstBase.getRMin()));
				++rCount; // Increment the row count to point to the next generator object
				// Coefficients corresponding to upper ramp limits for next interval
				lhs[rCount] = 0;
				lhs[rCount] += (decvar[rCount-2]-decvar[rCount-3]);
				modelGenQP->addConstr(lhs[rCount] <= (genSolverFirstBase.getRMax()));
				++rCount; // Increment the row count to point to the next generator object
				// Coefficients corresponding to lower ramp limits for previous interval
				lhs[rCount] = 0;
				lhs[rCount] += (decvar[rCount-4]-genSolverFirstBase.getPgPrev());
				modelGenQP->addConstr(lhs[rCount] >= (genSolverFirstBase.getRMin()));
				++rCount; // Increment the row count to point to the next generator object
				// Coefficients corresponding to upper ramp limits for previous interval
				lhs[rCount] = 0;
				lhs[rCount] += (decvar[rCount-5]-genSolverFirstBase.getPgPrev());
				modelGenQP->addConstr(lhs[rCount] <= (genSolverFirstBase.getRMax()));
			}
			if ( (dispatchInterval!=0) && (flagLast==0) ){ 
				// Coefficients corresponding to lower generation limits
				lhs[rCount] = 0;
				lhs[rCount] += decvar[rCount];
				modelGenQP->addConstr(lhs[rCount] >= (genSolverDZBase.getPMin()));
				++rCount; // Increment the row count to point to the next generator object
				// Coefficients corresponding to upper generation limits
				lhs[rCount] = 0;
				lhs[rCount] += decvar[rCount - 1];
				modelGenQP->addConstr(lhs[rCount] <= (genSolverDZBase.getPMax()));
				++rCount; // Increment the row count to point to the next generator object
				for (int counterCont = 0; counterCont <= contCountGen; ++counterCont) {
					// Coefficients corresponding to lower ramp limits for next interval
					lhs[rCount] = 0;
					lhs[rCount] += (decvar[rCount-(1+counterCont)]-decvar[1]);
					modelGenQP->addConstr(lhs[rCount] >= (genSolverDZBase.getRMin()));
					++rCount; // Increment the row count to point to the next generator object
					// Coefficients corresponding to upper ramp limits for next interval
					lhs[rCount] = 0;
					lhs[rCount] += (decvar[rCount-(2+counterCont)]-decvar[1]);
					modelGenQP->addConstr(lhs[rCount] <= (genSolverDZBase.getRMax()));
					++rCount; // Increment the row count to point to the next generator object
				}
				// Coefficients corresponding to lower ramp limits for previous interval
				lhs[rCount] = 0;
				lhs[rCount] += (decvar[1]-decvar[rCount-(2+contCountGen)]);
				modelGenQP->addConstr(lhs[rCount] >= (genSolverDZBase.getRMin()));
				++rCount; // Increment the row count to point to the next generator object
				// Coefficients corresponding to upper ramp limits for previous interval
				lhs[rCount] = 0;
				lhs[rCount] += (decvar[1]-decvar[rCount-(3+contCountGen)]);
				modelGenQP->addConstr(lhs[rCount] <= (genSolverDZBase.getRMax()));
			}
			if ( (dispatchInterval!=0) && (flagLast==1) ){ 
				// Coefficients corresponding to lower generation limits
				lhs[rCount] = 0;
				lhs[rCount] += decvar[rCount];
				modelGenQP->addConstr(lhs[rCount] >= (genSolverSecondBase.getPMin()));
				++rCount; // Increment the row count to point to the next generator object
				// Coefficients corresponding to upper generation limits
				lhs[rCount] = 0;
				lhs[rCount] += decvar[rCount - 1];
				modelGenQP->addConstr(lhs[rCount] <= (genSolverSecondBase.getPMax()));
				++rCount; // Increment the row count to point to the next generator object
				// Coefficients corresponding to lower ramp limits for next interval
				lhs[rCount] = 0;
				lhs[rCount] += (genSolverSecondBase.getLastIntChoice())*(PgenAPP-decvar[rCount-2]);
				modelGenQP->addConstr(lhs[rCount] >= (genSolverSecondBase.getRMin()));
				++rCount; // Increment the row count to point to the next generator object
				// Coefficients corresponding to upper ramp limits for next interval
				lhs[rCount] = 0;
				lhs[rCount] += (genSolverSecondBase.getLastIntChoice())*(PgenAPP-decvar[rCount-3]);
				modelGenQP->addConstr(lhs[rCount] <= (genSolverSecondBase.getRMax()));
				++rCount; // Increment the row count to point to the next generator object
				// Coefficients corresponding to lower ramp limits for previous interval
				lhs[rCount] = 0;
				lhs[rCount] += (decvar[rCount-4]-decvar[rCount-3]);
				modelGenQP->addConstr(lhs[rCount] >= (genSolverSecondBase.getRMin()));
				++rCount; // Increment the row count to point to the next generator object
				// Coefficients corresponding to upper ramp limits for previous interval
				lhs[rCount] = 0;
				lhs[rCount] += (decvar[rCount-5]-decvar[rCount-4]);
				modelGenQP->addConstr(lhs[rCount] <= (genSolverSecondBase.getRMax()));
			}
		}
		if ( dummyZeroIntFlag == 0 ) { // If there is no dummy-zero interval at the beginning
			if ( (dispatchInterval!=0) && (flagLast==0) ){ 
				// Coefficients corresponding to lower generation limits
				lhs[rCount] = 0;
				lhs[rCount] += decvar[rCount];
				modelGenQP->addConstr(lhs[rCount] >= (genSolverFirst.getPMin()));
				++rCount; // Increment the row count to point to the next generator object
				// Coefficients corresponding to upper generation limits
				lhs[rCount] = 0;
				lhs[rCount] += decvar[rCount - 1];
				modelGenQP->addConstr(lhs[rCount] <= (genSolverFirst.getPMax()));
				++rCount; // Increment the row count to point to the next generator object
				for (int counterCont = 0; counterCont <= contCountGen; ++counterCont) {
					// Coefficients corresponding to lower ramp limits for next interval
					lhs[rCount] = 0;
					lhs[rCount] += (decvar[rCount-(1+counterCont)]-decvar[1]);
					modelGenQP->addConstr(lhs[rCount] >= (genSolverFirst.getRMin()));
					++rCount; // Increment the row count to point to the next generator object
					// Coefficients corresponding to upper ramp limits for next interval
					lhs[rCount] = 0;
					lhs[rCount] += (decvar[rCount-(2+counterCont)]-decvar[1]);
					modelGenQP->addConstr(lhs[rCount] <= (genSolverFirst.getRMax()));
					++rCount; // Increment the row count to point to the next generator object
				}
				// Coefficients corresponding to lower ramp limits for previous interval
				lhs[rCount] = 0;
				lhs[rCount] += (decvar[1]-decvar[rCount-genSolverFirst.getPgPrev()]);
				modelGenQP->addConstr(lhs[rCount] >= (genSolverFirst.getRMin()));
				++rCount; // Increment the row count to point to the next generator object
				// Coefficients corresponding to upper ramp limits for previous interval
				lhs[rCount] = 0;
				lhs[rCount] += (decvar[1]-decvar[rCount-genSolverFirst.getPgPrev()]);
				modelGenQP->addConstr(lhs[rCount] <= (genSolverFirst.getRMax()));
			}
			if ( (dispatchInterval!=0) && (flagLast==1) ){ 
				// Coefficients corresponding to lower generation limits
				lhs[rCount] = 0;
				lhs[rCount] += decvar[rCount];
				modelGenQP->addConstr(lhs[rCount] >= (genSolverSecondBase.getPMin()));
				++rCount; // Increment the row count to point to the next generator object
				// Coefficients corresponding to upper generation limits
				lhs[rCount] = 0;
				lhs[rCount] += decvar[rCount - 1];
				modelGenQP->addConstr(lhs[rCount] <= (genSolverSecondBase.getPMax()));
				++rCount; // Increment the row count to point to the next generator object
				// Coefficients corresponding to lower ramp limits for next interval
				lhs[rCount] = 0;
				lhs[rCount] += (genSolverSecondBase.getLastIntChoice())*(PgenAPP-decvar[rCount-2]);
				modelGenQP->addConstr(lhs[rCount] >= (genSolverSecondBase.getRMin()));
				++rCount; // Increment the row count to point to the next generator object
				// Coefficients corresponding to upper ramp limits for next interval
				lhs[rCount] = 0;
				lhs[rCount] += (genSolverSecondBase.getLastIntChoice())*(PgenAPP-decvar[rCount-3]);
				modelGenQP->addConstr(lhs[rCount] <= (genSolverSecondBase.getRMax()));
				++rCount; // Increment the row count to point to the next generator object
				// Coefficients corresponding to lower ramp limits for previous interval
				lhs[rCount] = 0;
				lhs[rCount] += (decvar[rCount-4]-decvar[rCount-3]);
				modelGenQP->addConstr(lhs[rCount] >= (genSolverSecondBase.getRMin()));
				++rCount; // Increment the row count to point to the next generator object
				// Coefficients corresponding to upper ramp limits for previous interval
				lhs[rCount] = 0;
				lhs[rCount] += (decvar[rCount-5]-decvar[rCount-4]);
				modelGenQP->addConstr(lhs[rCount] <= (genSolverSecondBase.getRMax()));
			}
		}
	}
	if ( baseContScenario != 0 ){ 
		if (contSolverAccuracy == 0) {// If the exhaustive calculation for contingency scenarios is not desired
			// Coefficients corresponding to lower generation limits
			lhs[rCount] = 0;
			lhs[rCount] += decvar[rCount];
			modelGenQP->addConstr(lhs[rCount] >= (genSolverCont.getPMin()));
			++rCount; // Increment the row count to point to the next generator object
			// Coefficients corresponding to upper generation limits
			lhs[rCount] = 0;
			lhs[rCount] += decvar[rCount - 1];
			modelGenQP->addConstr(lhs[rCount] <= (genSolverCont.getPMax()));
		}
		else {
			if ( dummyZeroIntFlag == 1 ) { // If there is a dummy-zero interval at the beginning
				if (dispatchInterval==0) { // The dummy zero interval
					// Coefficients corresponding to lower generation limits
					lhs[rCount] = 0;
					lhs[rCount] += decvar[rCount];
					modelGenQP->addConstr(lhs[rCount] >= (genSolverCont.getPMin()));
					++rCount; // Increment the row count to point to the next generator object
					// Coefficients corresponding to upper generation limits
					lhs[rCount] = 0;
					lhs[rCount] += decvar[rCount - 1];
					modelGenQP->addConstr(lhs[rCount] <= (genSolverCont.getPMax()));
				}
				if ( (dispatchInterval!=0) && (flagLast==0) ){ 
					// Coefficients corresponding to lower generation limits
					lhs[rCount] = 0;
					lhs[rCount] += decvar[rCount];
					modelGenQP->addConstr(lhs[rCount] >= (genSolverDZBase.getPMin()));
					++rCount; // Increment the row count to point to the next generator object
					// Coefficients corresponding to upper generation limits
					lhs[rCount] = 0;
					lhs[rCount] += decvar[rCount - 1];
					modelGenQP->addConstr(lhs[rCount] <= (genSolverDZBase.getPMax()));
					++rCount; // Increment the row count to point to the next generator object
					for (int counterCont = 0; counterCont <= contCountGen; ++counterCont) {
						// Coefficients corresponding to lower ramp limits for next interval
						lhs[rCount] = 0;
						lhs[rCount] += (decvar[rCount-(1+counterCont)]-decvar[1]);
						modelGenQP->addConstr(lhs[rCount] >= (genSolverDZBase.getRMin()));
						++rCount; // Increment the row count to point to the next generator object
						// Coefficients corresponding to upper ramp limits for next interval
						lhs[rCount] = 0;
						lhs[rCount] += (decvar[rCount-(2+counterCont)]-decvar[1]);
						modelGenQP->addConstr(lhs[rCount] <= (genSolverDZBase.getRMax()));
						++rCount; // Increment the row count to point to the next generator object
					}
					// Coefficients corresponding to lower ramp limits for previous interval
					lhs[rCount] = 0;
					lhs[rCount] += (decvar[1]-decvar[rCount-(2+contCountGen)]);
					modelGenQP->addConstr(lhs[rCount] >= (genSolverDZBase.getRMin()));
					++rCount; // Increment the row count to point to the next generator object
					// Coefficients corresponding to upper ramp limits for previous interval
					lhs[rCount] = 0;
					lhs[rCount] += (decvar[1]-decvar[rCount-(3+contCountGen)]);
					modelGenQP->addConstr(lhs[rCount] <= (genSolverDZBase.getRMax()));
				}
				if ( (dispatchInterval!=0) && (flagLast==1) ){ 
					// Coefficients corresponding to lower generation limits
					lhs[rCount] = 0;
					lhs[rCount] += decvar[rCount];
					modelGenQP->addConstr(lhs[rCount] >= (genSolverSecondBase.getPMin()));
					++rCount; // Increment the row count to point to the next generator object
					// Coefficients corresponding to upper generation limits
					lhs[rCount] = 0;
					lhs[rCount] += decvar[rCount - 1];
					modelGenQP->addConstr(lhs[rCount] <= (genSolverSecondBase.getPMax()));
					++rCount; // Increment the row count to point to the next generator object
					// Coefficients corresponding to lower ramp limits for next interval
					lhs[rCount] = 0;
					lhs[rCount] += (genSolverSecondBase.getLastIntChoice())*(PgenAPP-decvar[rCount-2]);
					modelGenQP->addConstr(lhs[rCount] >= (genSolverSecondBase.getRMin()));
					++rCount; // Increment the row count to point to the next generator object
					// Coefficients corresponding to upper ramp limits for next interval
					lhs[rCount] = 0;
					lhs[rCount] += (genSolverSecondBase.getLastIntChoice())*(PgenAPP-decvar[rCount-3]);
					modelGenQP->addConstr(lhs[rCount] <= (genSolverSecondBase.getRMax()));
					++rCount; // Increment the row count to point to the next generator object
					// Coefficients corresponding to lower ramp limits for previous interval
					lhs[rCount] = 0;
					lhs[rCount] += (decvar[rCount-4]-decvar[rCount-3]);
					modelGenQP->addConstr(lhs[rCount] >= (genSolverSecondBase.getRMin()));
					++rCount; // Increment the row count to point to the next generator object
					// Coefficients corresponding to upper ramp limits for previous interval
					lhs[rCount] = 0;
					lhs[rCount] += (decvar[rCount-5]-decvar[rCount-4]);
					modelGenQP->addConstr(lhs[rCount] <= (genSolverSecondBase.getRMax()));
				}
			}
			if ( dummyZeroIntFlag == 0 ) { // If there is no dummy-zero interval at the beginning
				if ( (dispatchInterval!=0) && (flagLast==0) ){ 
					// Coefficients corresponding to lower generation limits
					lhs[rCount] = 0;
					lhs[rCount] += decvar[rCount];
					modelGenQP->addConstr(lhs[rCount] >= (genSolverFirst.getPMin()));
					++rCount; // Increment the row count to point to the next generator object
					// Coefficients corresponding to upper generation limits
					lhs[rCount] = 0;
					lhs[rCount] += decvar[rCount - 1];
					modelGenQP->addConstr(lhs[rCount] <= (genSolverFirst.getPMax()));
					++rCount; // Increment the row count to point to the next generator object
					for (int counterCont = 0; counterCont <= contCountGen; ++counterCont) {
						// Coefficients corresponding to lower ramp limits for next interval
						lhs[rCount] = 0;
						lhs[rCount] += (decvar[rCount-(1+counterCont)]-decvar[1]);
						modelGenQP->addConstr(lhs[rCount] >= (genSolverFirst.getRMin()));
						++rCount; // Increment the row count to point to the next generator object
						// Coefficients corresponding to upper ramp limits for next interval
						lhs[rCount] = 0;
						lhs[rCount] += (decvar[rCount-(2+counterCont)]-decvar[1]);
						modelGenQP->addConstr(lhs[rCount] <= (genSolverFirst.getRMax()));
						++rCount; // Increment the row count to point to the next generator object
					}
					// Coefficients corresponding to lower ramp limits for previous interval
					lhs[rCount] = 0;
					lhs[rCount] += (decvar[1]-decvar[rCount-genSolverFirst.getPgPrev()]);
					modelGenQP->addConstr(lhs[rCount] >= (genSolverFirst.getRMin()));
					++rCount; // Increment the row count to point to the next generator object
					// Coefficients corresponding to upper ramp limits for previous interval
					lhs[rCount] = 0;
					lhs[rCount] += (decvar[1]-decvar[rCount-genSolverFirst.getPgPrev()]);
					modelGenQP->addConstr(lhs[rCount] <= (genSolverFirst.getRMax()));
				}
				if ( (dispatchInterval!=0) && (flagLast==1) ){ 
					// Coefficients corresponding to lower generation limits
					lhs[rCount] = 0;
					lhs[rCount] += decvar[rCount];
					modelGenQP->addConstr(lhs[rCount] >= (genSolverSecondBase.getPMin()));
					++rCount; // Increment the row count to point to the next generator object
					// Coefficients corresponding to upper generation limits
					lhs[rCount] = 0;
					lhs[rCount] += decvar[rCount - 1];
					modelGenQP->addConstr(lhs[rCount] <= (genSolverSecondBase.getPMax()));
					++rCount; // Increment the row count to point to the next generator object
					// Coefficients corresponding to lower ramp limits for next interval
					lhs[rCount] = 0;
					lhs[rCount] += (genSolverSecondBase.getLastIntChoice())*(PgenAPP-decvar[rCount-2]);
					modelGenQP->addConstr(lhs[rCount] >= (genSolverSecondBase.getRMin()));
					++rCount; // Increment the row count to point to the next generator object
					// Coefficients corresponding to upper ramp limits for next interval
					lhs[rCount] = 0;
					lhs[rCount] += (genSolverSecondBase.getLastIntChoice())*(PgenAPP-decvar[rCount-3]);
					modelGenQP->addConstr(lhs[rCount] <= (genSolverSecondBase.getRMax()));
					++rCount; // Increment the row count to point to the next generator object
					// Coefficients corresponding to lower ramp limits for previous interval
					lhs[rCount] = 0;
					lhs[rCount] += (decvar[rCount-4]-decvar[rCount-3]);
					modelGenQP->addConstr(lhs[rCount] >= (genSolverSecondBase.getRMin()));
					++rCount; // Increment the row count to point to the next generator object
					// Coefficients corresponding to upper ramp limits for previous interval
					lhs[rCount] = 0;
					lhs[rCount] += (decvar[rCount-5]-decvar[rCount-4]);
					modelGenQP->addConstr(lhs[rCount] <= (genSolverSecondBase.getRMax()));
				}
			}
		}	
	}
	// RUN THE OPTIMIZATION SIMULATION ALGORITHM //
	modelGenQP->optimize(); // Solves the optimization problem
	int stat = modelGenQP->get(GRB_IntAttr_Status); // Outputs the solution status of the problem 

	// DISPLAY THE SOLUTION DETAILS //
	if (stat == GRB_INFEASIBLE){
		cout << "\nThe solution to the problem is INFEASIBLE." << endl;
		delete modelGenQP; // Free the memory of the GUROBI Problem Model
	} else if (stat == GRB_INF_OR_UNBD) {
		cout << "\nNO FEASIBLE or BOUNDED solution to the problem exists." << endl;
		delete modelGenQP; // Free the memory of the GUROBI Problem Model
	} else if (stat == GRB_UNBOUNDED) {
		cout << "\nThe solution to the problem is UNBOUNDED." << endl;
		delete modelGenQP; // Free the memory of the GUROBI Problem Model
	} else if (stat == GRB_OPTIMAL) {
		//Get the Optimal Objective Value results//
		z = modelGenQP->get(GRB_DoubleAttr_ObjVal);
		// writing results of different variables
		vector<double> x; // Vector for storing decision variable output 
		x.push_back(0); // Initialize the decision Variable vector
		objOpt = 0;
		//Power Generation
		int arrayInd = 1;
		if ( baseContScenario == 0 ) { // Use the solver for base cases
			if ( dummyZeroIntFlag == 1 ) { // If dummy zero interval is considered
				if ((dispatchInterval==0) && (flagLast==0)) { // For the dummy zeroth interval
					x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
					Pg = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator Power iterate
					objOpt += (genSolverFirstBase.getQuadCoeff())*(Pg)*(Pg)+(genSolverFirstBase.getLinCoeff())*(Pg)+(genSolverFirstBase.getConstCoeff());
					++arrayInd;
					x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
					PgenNext = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator Power iterate for next interval	
					// Internal node voltage phase angle variables
					++arrayInd;
					x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
					Thetag = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator voltage angle iterate	
					PgenPrev = genSolverFirstBase.getPgPrev();
				}
				if ( (dispatchInterval!=0) && (flagLast==0) ){ 
					x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
					Pg = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator Power iterate
					objOpt += (genSolverFirstBase.getQuadCoeff())*(Pg)*(Pg)+(genSolverFirstBase.getLinCoeff())*(Pg)+(genSolverFirstBase.getConstCoeff());
					++arrayInd;
					for (int counterCont = 0; counterCont <= contCountGen; ++counterCont) {
						x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
						PgenNextPtr[counterCont] = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator Power iterate for next interval	
						// Internal node voltage phase angle variables
						++arrayInd;
					}
					x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
					PgenPrev = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator Power iterate for previous interval	
					// Internal node voltage phase angle variables
					++arrayInd;
					// Internal node voltage phase angle variables
					x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
					Thetag = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator voltage angle iterate
				}
				if ( (dispatchInterval!=0) && (flagLast==1) ){ 
					x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
					Pg = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator Power iterate
					objOpt += (genSolverFirstBase.getQuadCoeff())*(Pg)*(Pg)+(genSolverFirstBase.getLinCoeff())*(Pg)+(genSolverFirstBase.getConstCoeff());
					++arrayInd;
					x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
					PgenPrev = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator Power iterate for previous interval	
					// Internal node voltage phase angle variables
					++arrayInd;
					// Internal node voltage phase angle variables
					x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
					Thetag = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator voltage angle iterate	
				}
			}
			if ( dummyZeroIntFlag == 0 ) { // If dummy zero interval is not considered
				if ( (dispatchInterval!=0) && (flagLast==0) ){ 
					x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
					Pg = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator Power iterate
					objOpt += (genSolverFirstBase.getQuadCoeff())*(Pg)*(Pg)+(genSolverFirstBase.getLinCoeff())*(Pg)+(genSolverFirstBase.getConstCoeff());
					++arrayInd;
					for (int counterCont = 0; counterCont <= contCountGen; ++counterCont) {
						x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
						PgenNextPtr[counterCont] = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator Power iterate for next interval	
						// Internal node voltage phase angle variables
						++arrayInd;
					}
					// Internal node voltage phase angle variables
					x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
					Thetag = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator voltage angle iterate
				}
				if ( (dispatchInterval!=0) && (flagLast==1) ){ 
					x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
					Pg = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator Power iterate
					objOpt += (genSolverFirstBase.getQuadCoeff())*(Pg)*(Pg)+(genSolverFirstBase.getLinCoeff())*(Pg)+(genSolverFirstBase.getConstCoeff());
					++arrayInd;
					x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
					PgenPrev = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator Power iterate for previous interval	
					// Internal node voltage phase angle variables
					++arrayInd;
					// Internal node voltage phase angle variables
					x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
					Thetag = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator voltage angle iterate	
				}
			}
		}
		if ( baseContScenario != 0 ){ 
			if (contSolverAccuracy == 0) {// If the exhaustive calculation for contingency scenarios is not desired
				x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
				Pg = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator Power iterate
				objOpt += (genSolverFirstBase.getQuadCoeff())*(Pg)*(Pg)+(genSolverFirstBase.getLinCoeff())*(Pg)+(genSolverFirstBase.getConstCoeff());
				++arrayInd;
				// Internal node voltage phase angle variables
				x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
				Thetag = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator voltage angle iterate	
			}
			else {
				if ( dummyZeroIntFlag == 1 ) { // If there is a dummy-zero interval at the beginning
					if (dispatchInterval==0) { // The dummy zero interval
						x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
						Pg = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator Power iterate
						objOpt += (genSolverFirstBase.getQuadCoeff())*(Pg)*(Pg)+(genSolverFirstBase.getLinCoeff())*(Pg)+(genSolverFirstBase.getConstCoeff());
						++arrayInd;
						// Internal node voltage phase angle variables
						x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
						Thetag = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator voltage angle iterate	
					}
					if ( (dispatchInterval!=0) && (flagLast==0) ){ 
						x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
						Pg = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator Power iterate
						objOpt += (genSolverFirstBase.getQuadCoeff())*(Pg)*(Pg)+(genSolverFirstBase.getLinCoeff())*(Pg)+(genSolverFirstBase.getConstCoeff());
						++arrayInd;
						for (int counterCont = 0; counterCont <= contCountGen; ++counterCont) {
							x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
							PgenNextPtr[counterCont] = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator Power iterate for next interval	
							// Internal node voltage phase angle variables
							++arrayInd;
						}
						x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
						PgenPrev = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator Power iterate for previous interval	
						// Internal node voltage phase angle variables
						++arrayInd;
						// Internal node voltage phase angle variables
						x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
						Thetag = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator voltage angle iterate
					}
					if ( (dispatchInterval!=0) && (flagLast==1) ){ 
						x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
						Pg = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator Power iterate
						objOpt += (genSolverFirstBase.getQuadCoeff())*(Pg)*(Pg)+(genSolverFirstBase.getLinCoeff())*(Pg)+(genSolverFirstBase.getConstCoeff());
						++arrayInd;
						x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
						PgenPrev = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator Power iterate for previous interval	
						// Internal node voltage phase angle variables
						++arrayInd;
						// Internal node voltage phase angle variables
						x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
						Thetag = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator voltage angle iterate	
					}
				}
				if ( dummyZeroIntFlag == 0 ) { // If dummy zero interval is not considered
					if ( (dispatchInterval!=0) && (flagLast==0) ){ 
						x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
						Pg = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator Power iterate
						objOpt += (genSolverFirstBase.getQuadCoeff())*(Pg)*(Pg)+(genSolverFirstBase.getLinCoeff())*(Pg)+(genSolverFirstBase.getConstCoeff());
						++arrayInd;
						for (int counterCont = 0; counterCont <= contCountGen; ++counterCont) {
							x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
							PgenNextPtr[counterCont] = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator Power iterate for next interval	
							// Internal node voltage phase angle variables
							++arrayInd;
						}
						// Internal node voltage phase angle variables
						x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
						Thetag = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator voltage angle iterate
					}
					if ( (dispatchInterval!=0) && (flagLast==1) ){ 
						x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
						Pg = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator Power iterate
						objOpt += (genSolverFirstBase.getQuadCoeff())*(Pg)*(Pg)+(genSolverFirstBase.getLinCoeff())*(Pg)+(genSolverFirstBase.getConstCoeff());
						++arrayInd;
						x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
						PgenPrev = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator Power iterate for previous interval	
						// Internal node voltage phase angle variables
						++arrayInd;
						// Internal node voltage phase angle variables
						x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
						Thetag = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator voltage angle iterate	
					}
				}
			}				
		}
		connNodegPtr->powerangleMessage( Pg, v, Thetag ); // passes to node object the corresponding iterates of power, angle and v
		delete modelGenQP; 
	}
} // function gpowerangleMessage ends

double Generator::genPower() //const // function genPower begins
{
	return Pg; // returns the Pg iterate
} // function genPower ends

double Generator::genPowerPrev() //const // function genPower begins
{
	if (dispatchInterval==0)
		return genSolverFirstBase.getPgPrev(); // returns the Pg iterate
	else 
		return PgenPrev;
} // function genPower ends

double Generator::genPowerNext(int nextScen) //const // function genPower begins
{
	if (flagLast==1)
		return Pg; // returns the Pg iterate
	if ((dispatchInterval!=0) && (flagLast==0))
		return PgenNextPtr[nextScen];
	else 
		return PgenNext;
} // function genPower ends

double Generator::objectiveGen() // function objectiveGen begins
{
	if ( baseContScenario == 0 ) { // Use the solver for base cases
		if ( dummyZeroIntFlag == 1 ) { // If dummy zero interval is considered
			if ((dispatchInterval==0) && (flagLast==0)) { // For the dummy zeroth interval
				return genSolverFirstBase.getObj(); //returns the evaluated objective
			}
			if ((dispatchInterval!=0) && (flagLast==0)) { // For the first interval
				return genSolverDZBase.getObj(); //returns the evaluated objective
			}
			if ((dispatchInterval!=0) && (flagLast==1)) { // For the second (or, in this case, the last) interval
				return genSolverSecondBase.getObj(); //returns the evaluated objective
			}
		}
		if ( dummyZeroIntFlag == 0 ) { // If dummy zero interval is considered
			if ((dispatchInterval==0) && (flagLast==0)) { // For the dummy zeroth interval
				return genSolverFirstBase.getObj(); //returns the evaluated objective
			}
			if ((dispatchInterval!=0) && (flagLast==0)) { // For the first interval
				return genSolverFirst.getObj(); //returns the evaluated objective
			}
			if ((dispatchInterval!=0) && (flagLast==1)) { // For the second (or, in this case, the last) interval
				return genSolverSecondBase.getObj(); //returns the evaluated objective
			}
		}
	}
	if ( baseContScenario != 0 ) {
		if ( dummyZeroIntFlag == 1 ) { // If dummy zero interval is considered
			if ((dispatchInterval==0) && (flagLast==0)) { // For the dummy zeroth interval
				return genSolverCont.getObj(); //returns the evaluated objective
			}
			if ((dispatchInterval!=0) && (flagLast==0)) { // For the first interval
				if (contSolverAccuracy == 0) {// If the exhaustive calculation for contingency scenarios is not desired
					return genSolverCont.getObj(); //returns the evaluated objective
				}
				if (contSolverAccuracy != 0) {// If the exhaustive calculation for contingency scenarios not desired
					return genSolverDZCont.getObj();
				}
			}
			if ((dispatchInterval!=0) && (flagLast==1)) { // For the second (or, in this case, the last) interval
				if (contSolverAccuracy == 0) {// If the exhaustive calculation for contingency scenarios is not desired
					return genSolverCont.getObj(); //returns the evaluated objective
				}
				if (contSolverAccuracy != 0) {// If the exhaustive calculation for contingency scenarios is desired
					return genSolverSecondCont.getObj(); //returns the evaluated objective
				}
			}
		}
		if ( dummyZeroIntFlag == 0 ) { // If dummy zero interval is not considered
			if ((dispatchInterval==0) && (flagLast==0)) { // For the dummy zeroth interval **/ Will not be used in this case**/
				return genSolverCont.getObj(); //returns the evaluated objective
			}
			if ((dispatchInterval!=0) && (flagLast==0)) { // For the first interval
				if (contSolverAccuracy == 0) {// If the exhaustive calculation for contingency scenarios is not desired
					return genSolverCont.getObj(); //returns the evaluated objective
				}
				if (contSolverAccuracy != 0) {// If the exhaustive calculation for contingency scenarios not desired
					return genSolverFirstCont.getObj();
				}
			}
			if ((dispatchInterval!=0) && (flagLast==1)) { // For the second (or, in this case, the last) interval
				if (contSolverAccuracy == 0) {// If the exhaustive calculation for contingency scenarios is not desired
					return genSolverCont.getObj(); //returns the evaluated objective
				}
				if (contSolverAccuracy != 0) {// If the exhaustive calculation for contingency scenarios is desired
					return genSolverSecondCont.getObj(); //returns the evaluated objective
				}
			}
		}			
	}
} // function objectiveGen ends

double Generator::objectiveGenGUROBI() // Objective from GUROBI ADMM 
{
	return objOpt; //returns the evaluated objective
} // function objectiveGen ends

double Generator::calcPtilde() //const // function calcPtilde begins
{
	double P_avg = connNodegPtr->PavMessage(); // Gets average power from the corresponding node object
	double Ptilde = Pg - P_avg; // calculates the difference between power iterate and average
	return Ptilde; // returns the difference
} // function calcPtilde ends

double Generator::calcPavInit() const // function calcPavInit begins
{
	return connNodegPtr->devpinitMessage(); // seeks the initial Ptilde from the node
} // function calcPavInit ends

double Generator::getu() const // function getu begins
{
	double u = connNodegPtr->uMessage(); // gets the value of the price corresponding to power balance from node
	//cout << "u: " << u << endl;
	return u; // returns the price
} // function getu ends

double Generator::calcThetatilde() //const // function calcThetatilde begins
{
	//cout << "Thetag: " << Thetag << endl;
	double Theta_avg = connNodegPtr->ThetaavMessage(); // get the average voltage angle at the particular node
	//cout << "Theta_avg: " << Theta_avg << endl;
	double Theta_tilde = Thetag - Theta_avg; // claculate the deviation between the voltage angle of the device and the average
	return Theta_tilde; // return the deviation
} // function calcThetatilde ends

double Generator::calcvtilde() const // function calcvtilde begins
{
	double v_avg = connNodegPtr->vavMessage(); // get the average of the Lagrange multiplier corresponding to voltage angle balance
	//cout << "v_avg: " << v_avg << endl;
	double v_tilde = v - v_avg; // calculate the deviation of the node Lagrange multiplier to the average
	return v_tilde; // return the deviation
} // function calcvtilde ends

double Generator::getv() // function getv begins
{
	//cout << "v_initial: " << v << endl;
	v = v + calcThetatilde(); // Calculate the value of the Lagrange multiplier corresponding to angle constraint
	//cout << "v_final: " << v << endl;
	return v; // Calculate the value of the Lagrange multiplier corresponding to angle constraint
} // function getv ends		
double Generator::getPMax(){return genSolverFirstBase.getPMax();}
double Generator::getPMin(){return genSolverFirstBase.getPMin();}
double Generator::getQuadCoeff(){return genSolverFirstBase.getQuadCoeff();}
double Generator::getLinCoeff(){return genSolverFirstBase.getLinCoeff();}
double Generator::getConstCoeff(){return genSolverFirstBase.getConstCoeff();}
double Generator::getPgenPrev(){return genSolverFirstBase.getPgPrev();}
double Generator::getPgenNext(){return genSolverFirstBase.getPNextSol();}
double Generator::getRMax(){return genSolverFirstBase.getRMax();}
double Generator::getRMin(){return genSolverFirstBase.getRMin();}
