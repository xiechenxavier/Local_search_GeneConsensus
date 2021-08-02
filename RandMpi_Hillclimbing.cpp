#include "SolveGene.hpp"
#include <iostream>
#include <time.h>
#include <fstream>
#include <cstring>
#include <random>
#include <mpi.h>
#include<math.h>
#include <chrono>

using std::cout;
using std::endl;
using std::fstream;
using std::ofstream;
using std::string;

using namespace std;

void merge_ranking(int tab1[],int taille1,int tab2[], int taille2, SolveGene genepb0){
    //create an array to contain tab1 and tab2 
    int* mergedTab = new int[taille1+taille2];
    //remplir les taille1 premiers éléments du mergedTab 
    for(int i=0;i<taille1;i++){
        mergedTab[i] = tab1[i];
    }
    //remplir le taille2 premiers éléments du mergedTab
    for(int j = taille1;j<taille1+taille2;j++){
        mergedTab[j] = tab2[j-taille1];
    }
    //Trier le tableau des gènes mergedTab par hill-climbing
    genepb0.setCurSol(mergedTab,taille1+taille2);
    genepb0.hillClimbing(true,true,true);
    vector<int> curSolMerged = genepb0.getCurSol(); // get the merged array
    //tab1 garde les premiers n/p éléments du mergedTab 
    for(int i=0;i<taille1;i++){
        tab1[i] = curSolMerged[i];
    }
    //tab2 garde les derniers n/p éléments du mergedTab
    for(int j=taille1;j<curSolMerged.size();j++){
        tab2[j-taille1] = curSolMerged[j];
    }
    //libérer le mémoire du mergedTab
    delete [] mergedTab;

}
void print_tab_distrib_gather(int tabLocal[],int tailleLocal,SolveGene genepb0) //这个函数的作用是把所有tabLocal在进程0处整合成tabGlobal然后打印
     {
	int id, p,j;
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Status status;
    //tableau fusionné et complète
	int* tabGlobal;
	int tailleGlobal = tailleLocal * p;
	tabGlobal = new int[tailleGlobal];

	MPI_Gather(tabLocal, tailleLocal, MPI_INT,tabGlobal, tailleLocal,MPI_INT, 0, MPI_COMM_WORLD);

	if (id==0)  
		{
	    	for(j=0;j<tailleGlobal;j++){
		    	std::cout  << " " << tabGlobal[j];  
            }
            	std::cout  << std::endl; 
            genepb0.setCurSol(tabGlobal,tailleGlobal);
            cout<< "best cout is:"<< genepb0.getCurrentProbability()<<endl; 
		}
	delete [] tabGlobal;
}

void init_tab_random(int tabLocal[],int tailleLocal,SolveGene genepb0)
     {
	
	int id, p,j;
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Status status;

	int* tabGlobal;
	int tailleGlobal = tailleLocal * p;

	srand(id+1);

	if (id==0){
        vector<int> vectorGlob;
        vectorGlob.resize(tailleGlobal);
        genepb0.randomizedConstr();
        vectorGlob = genepb0.getCurSol();
        tabGlobal = new int[tailleGlobal];
        memcpy(tabGlobal, &vectorGlob[0], tailleGlobal * sizeof(vectorGlob[0]));
	}
		//TODO: Scatter diffusion des morceaux du grand tableau
		MPI_Scatter(tabGlobal, tailleLocal, MPI_INT, tabLocal, tailleLocal,MPI_INT, 0, MPI_COMM_WORLD);
		// Libère la mémoire du grand tableau
	if (id==0){
		delete [] tabGlobal;
	}
	
	if (id==0) std::cout  << "Impression tableau initial" << std::endl;
	print_tab_distrib_gather(tabLocal,tailleLocal,genepb0);
	//print_tab_distrib(tabLocal,tailleLocal);

}


int main(int argc, char **argv)
{
    if (argc < 2)
    {
        //cerr << "Usage: tsp inputFile " << endl;
        return 1;
    }
    const char *instanceFile = argv[1];
    SolveGene genepb0;
    int nbMaxRestart = 2;
    genepb0.setVerboseMode(true);
    genepb0.setLSparameters(100, nbMaxRestart);

    genepb0.readInstanceFromFile(instanceFile);
    //display the instances
    genepb0.displayInstance();

    //MPI part
    int charge_unit; // unit taille

    double totaltime = 0;

    int tailleGlobalGenes = genepb0.getTailleGenes();
    MPI_Init(NULL, NULL);

    int rankID,j; //current process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &rankID);
    int nbprocs; // nb of process
    MPI_Comm_size(MPI_COMM_WORLD, &nbprocs);

    MPI_Status status;
    
    charge_unit = ceil((double)tailleGlobalGenes /
                       (double)nbprocs);
    
    auto start = std::chrono::steady_clock::now();
    srand(rankID+1);

    int *tabLocal;
    int *tabRecv;
    int tailleLocal = charge_unit;
    tabLocal = new int[tailleLocal];
    tabRecv = new int[tailleLocal];
    init_tab_random(tabLocal,tailleLocal,genepb0);

    genepb0.setCurSol(tabLocal, tailleLocal);
    genepb0.hillClimbing(true,true,true);
    
    vector<int> curSolLocal = genepb0.getCurSol();
    memcpy(tabLocal,&curSolLocal[0],tailleLocal*sizeof(curSolLocal[0]));

    for(j=0;j<nbprocs;j++){
        if(j%2 == 0){
            if(rankID%2 == 0 && rankID+1 <nbprocs){
                MPI_Recv(tabRecv,charge_unit,MPI_INT,rankID+1,0,MPI_COMM_WORLD,&status);
                merge_ranking(tabLocal,tailleLocal,tabRecv,charge_unit,genepb0);
                MPI_Send(tabRecv,tailleLocal,MPI_INT,rankID+1,0,MPI_COMM_WORLD);
            }
            if(rankID%2==1){
                MPI_Send(tabLocal,tailleLocal,MPI_INT,rankID-1,0,MPI_COMM_WORLD);
                MPI_Recv(tabLocal,tailleLocal,MPI_INT,rankID-1,0,MPI_COMM_WORLD,&status);
            }
            
        }

        if (j%2==1){
			//TODO: etape impaire
			if (rankID%2==1 && rankID+1 < nbprocs){ // un processus impair recoit du suivant et lui reenvoie son nouveau tableau
				MPI_Recv(tabRecv, charge_unit, MPI_INT, rankID+1, 0, MPI_COMM_WORLD,&status);
				merge_ranking(tabLocal,tailleLocal,tabRecv,charge_unit,genepb0);
				MPI_Send(tabRecv, tailleLocal, MPI_INT, rankID+1, 0, MPI_COMM_WORLD);
			}
			if (rankID%2==0 && rankID>0){
				MPI_Send(tabLocal, tailleLocal, MPI_INT, rankID-1, 0, MPI_COMM_WORLD);
				MPI_Recv(tabLocal, charge_unit, MPI_INT, rankID-1, 0, MPI_COMM_WORLD,&status);
			}
		}
        
    }

          //compter le temp
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    double count = elapsed_seconds.count();
    //compter le temp en final
    MPI_Allreduce(&count, &totaltime, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    //termination
    if(rankID==0) {
       std::cout << "\n Impression tableau final" << std::endl;
       cout<< "it takes " << totaltime << " s"<< endl;
    }
    print_tab_distrib_gather(tabLocal, tailleLocal,genepb0); 


    delete [] tabLocal;
    delete [] tabRecv;


    // delete [] tabLocal;
    MPI_Finalize();
    return 0;


}