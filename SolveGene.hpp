
#ifndef INSTANCETSP_HPP
#define INSTANCETSP_HPP

#include <vector>
#include <string>
#include <climits>

using namespace std;

class SolveGene {

private:

    //parameters related to the input data
    // Number of Genes
    int nbGenes;
    // Vector/matrix of wi<j between two Gene
    vector<float>  wijMatrix;

    //related to the solution of instance

    vector<int> curSol;
    vector<int> revSol;
	float currentPossibility;
    float revSolutionProb;
    float subSolutionProb;

    vector<int>  bestSolution;
	float bestPossibility;

	//parameters of LocalSearch
	int nbMaxIterLS;
	int nbMaxRestart;

	int nbIterLS;
	bool verboseMode;

      //following variables are used to count each mode moves
    int countSwap = 0;//it's used to count how many moves is swap mode
    int countRev = 0;
    int countInsert =0;
    double ilsHCtime;
    
    int nbIterSimulated;

    //private methods, only used by solving methods
	
	//WARNING: Different Gene rankings, i and j, the wi<j + wj<i = 1 
    float ProbaIbeforeJ(int i,int j) { return wijMatrix[nbGenes * i+j]; };
    float sumProba(int i,int j) { return wijMatrix[nbGenes * (curSol[i]-1)+(curSol[j]-1)]; };
    // version pour ne pas se soucier des bords?
    //    int distSol(int i1,int i2){return distanceMatrix[nbCities * curSol[i1%nbCities]+curSol[i2%nbCities]];};

    void printStatus();
    bool checkFeasibility();

    bool hillClimbingIter(bool swapMoves, bool revMoves, bool insertMoves);//TODO

    void computeReversedPossibility();
    void reverseGenes();

public:

    void testGene();
    //public methods, to call in external test functions

    void updateBestSolution();
    void computeProbability();
    float getCurrentProbability() { return currentPossibility; };


    void setLSparameters(int nbIter, int nbRestart){nbMaxIterLS = nbIter; nbMaxRestart= nbRestart;};
    void setVerboseMode(bool b){verboseMode=true;};
    float getSolutionCost(){updateBestSolution(); return bestPossibility;};
    int getNbIterLS(){return nbIterLS;};
    vector<int> getCurSol(){ return curSol;}
    double getilsHCtime() { return ilsHCtime;}

    //simulated climbing
    int getnbIterSimulated(){ return this->nbIterSimulated;}
    //solving methods
    void hillClimbing(bool swapMoves, bool revMoves, bool insertMoves);

    void naiveInit();
    void greedyInit();

    void randomizedConstr();
    // void randomizedGreedy();
    void graspHC(bool swapMoves, bool revMoves, bool insertMoves);
    void ilsHC(bool swapMoves, bool revMoves, bool insertMoves);

    void simulatedAnnealing();//TODO


    /* Reads and displays instance data. */
    void readInstanceFromFile(const string& fileName);
    void displayInstance();

};

#endif
