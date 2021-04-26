#include "SolveGene.hpp"
#include <iostream>
#include<time.h>
#include <fstream>

using std::cout; using std::ofstream;
using std::endl; using std::string;
using std::fstream;

using namespace std;

int main(int argc, char** argv) {
        if (argc < 2) {
        //cerr << "Usage: tsp inputFile " << endl;
        return 1;
    }
    const char* instanceFile = argv[1];
    SolveGene genepb;
    genepb.setVerboseMode(true);
	genepb.setLSparameters(100,1);

    genepb.readInstanceFromFile(instanceFile);
    //display the instances
    genepb.displayInstance();
    //test naiveInit
    // genepb.naiveInit();
    //test greedyInit
    // genepb.greedyInit();
    // test randomizedConstr
    // genepb.randomizedConstr();
    // test randomizedGreedy
    // genepb.randomizedGreedy();
    clock_t start,end;
    vector<float> objective_val; // objective value
    vector<double> simulAnnealtimes ;
    for (auto i = 0; i < 30; i++)
    {
        genepb.randomizedConstr();
        start = clock();
        genepb.simulatedAnnealing();
        end = clock();
        double calcul_temps = (double)(end-start)/CLOCKS_PER_SEC;
        cout<<"Simulated Annealing takes "<<calcul_temps<<" s"<<endl; //afficher
        objective_val.push_back(genepb.getCurrentProbability());
        simulAnnealtimes.push_back(calcul_temps);
    }
    
    string filename("./Result/Simulated_Annealing/SA_G200_80.txt");
    fstream file_out;

    file_out.open(filename, std::ios_base::out);
    if (!file_out.is_open()) {
        cout << "failed to open " << filename << '\n';
    } else {
        //value
        file_out<<"value: ";
        for(int it=0;it<objective_val.size()-1;it++){
            file_out<< objective_val[it] <<" ";
        }
        file_out<<objective_val[objective_val.size()-1]<<endl;
        file_out<<"simulAnnealtimes: ";
        for (auto it3 = 0; it3 < simulAnnealtimes.size()-1; it3++)
        {
            file_out<< simulAnnealtimes[it3] <<" ";
        }
        file_out<<simulAnnealtimes[simulAnnealtimes.size()-1]<<endl;
        
        cout << "Done Writing!" << endl;

    }
    // clock_t start,end;
    // vector<float> objective_val;
    // vector<int> nbiterLS;
    // vector<double> ilsHCtimes ;
    // for (int i = 0; i < 30; i++)
    // {
        
    // }
    
    /*vector<float> objective_val;
    vector<int> nbiterLS;
    vector<double> ilsHCtimes ;

    for(int i=0;i<30;i++){
        genepb.ilsHC(true,true,true);
        objective_val.push_back(genepb.getCurrentProbability());
        nbiterLS.push_back(genepb.getNbIterLS());
        ilsHCtimes.push_back(genepb.getilsHCtime());
    }
    

    string filename("./Result/rand_G20_20.txt");
    fstream file_out;

    file_out.open(filename, std::ios_base::out);
    if (!file_out.is_open()) {
        cout << "failed to open " << filename << '\n';
    } else {
        //value
        file_out<<"value: ";
        for(int it=0;it<objective_val.size()-1;it++){
            file_out<< objective_val[it] <<" ";
        }
        file_out<<objective_val[objective_val.size()-1]<<endl;
        //nbiterLS
        file_out<<"nbiterLS: ";
        for (auto it2 = 0; it2 <nbiterLS.size()-1 ; it2++)
        {
            file_out<< nbiterLS[it2] <<" ";
        }
        file_out<<nbiterLS[nbiterLS.size()-1]<<endl;
        //ilsHCtimes
        file_out<<"ilsHCtimes: ";
        for (auto it3 = 0; it3 < ilsHCtimes.size()-1; it3++)
        {
            file_out<< ilsHCtimes[it3] <<" ";
        }
        file_out<<ilsHCtimes[ilsHCtimes.size()-1]<<endl;
        
        cout << "Done Writing!" << endl;
    }*/


}