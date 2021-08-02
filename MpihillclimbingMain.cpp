#include "SolveGene.hpp"
#include <iostream>
#include <time.h>
#include <fstream>
#include <cstring>
#include <mpi.h>
#include<math.h>
#include <chrono>
#include "omp.h"
#include <algorithm>

using std::cout;
using std::endl;
using std::fstream;
using std::ofstream;
using std::string;

using namespace std;

void merge_ranking(int tab1[],int taille1,int tab2[], int taille2, SolveGene genepb0){
    //create an array to contain tab1 and tab2 
    int* mergedTab = new int[taille1+taille2];
    #pragma omp parallel 
    {
        // cout<<omp_get_thread_num()<<endl;
        #pragma omp for
        for(int i=0;i<taille1+taille2;i++){
            if(i<taille1)
                mergedTab[i] = tab1[i];
            else
                mergedTab[i] = tab2[i-taille1];
        }
    }
    genepb0.setCurSol(mergedTab,taille1+taille2);
    genepb0.hillClimbing(true,true,true);
    vector<int> curSolMerged = genepb0.getCurSol(); // get the merged array
    #pragma omp parallel
    {
        #pragma omp for
        for (int i = 0; i < taille1+taille2; i++)
        {
            if(i<taille1) tab1[i] = curSolMerged[i];
                
            else tab2[i-taille1] = curSolMerged[i];      
        }
    }
    delete [] mergedTab;

}

void print_tab_distrib_gather(int tabLocal[],int tailleLocal,SolveGene genepb0) //这个函数的作用是把所有tabLocal在进程0处整合成tabGlobal然后打印
     {

	int id, p,j;
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Status status;

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


int main(int argc, char **argv)
{
    if (argc < 2)
    {
        //cerr << "Usage: tsp inputFile " << endl;
        return 1;
    }
    //l'instance importée
    const char *instanceFile = argv[1];
    //creer un objet 
    SolveGene genepb0;
    int nbMaxRestart = 2;
    genepb0.setVerboseMode(true);
    genepb0.setLSparameters(1000, nbMaxRestart);
    // selon le chimin de l'instance, lire l'instance
    genepb0.readInstanceFromFile(instanceFile);
    //display the instances
    genepb0.displayInstance();
    //pour initialier le classement des genes par l'algorithme glouton
    genepb0.naiveInit(); 
    // obtenir la solution du classement 
    vector<int> firstSol = genepb0.getCurSol();

    int tailleGlobalGenes = firstSol.size();
    // préparer un tableau contient la solution pour la parallélisation après
    int* tabfirstSol = new int[tailleGlobalGenes];

    //vector -> tableau
    memcpy(tabfirstSol, &firstSol[0], tailleGlobalGenes * sizeof(firstSol[0]));

    //MPI part
    int charge_unit; // unit de taille
    // temps pris par chaque processus
    double* recvtemp;

    MPI_Init(NULL, NULL);

    int rankID,j; //current process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &rankID);
    int nbprocs; // nb of process
    MPI_Comm_size(MPI_COMM_WORLD, &nbprocs);

    MPI_Status status;

    charge_unit = ceil((double)tailleGlobalGenes /
                       (double)nbprocs);
    
    if(rankID==0){
        recvtemp = new double[nbprocs];
    }
    
    auto start = std::chrono::steady_clock::now();
    //TODO hillclimbing distribue
    // ParallelHillClimbing(charge_unit, rankID, tabfirstSol,tailleGlobalGenes,genepb0);
    //charge_unit and rankID will indicate which part of curSol will be treated for the current process
    int *tabLocal;//tableau local pour chaque processus
    int *tabRecv; 
    int tailleLocal; //taille du tableau local
    tabRecv = new int[charge_unit];
    //pour diviser le tableau complete (ici, j'ai pensé si le dernier indice > taille du tableau complète)
    if ((rankID + 1) * charge_unit < tailleGlobalGenes)
    {
        tailleLocal = charge_unit;
        tabLocal = new int[charge_unit]; //initializer le tableau local
        // remplir le tableau local
        for (int i = rankID * charge_unit; i < (rankID + 1) * charge_unit; i++)
        {
            tabLocal[i - rankID * charge_unit] = tabfirstSol[i];
        }

    }else{
        tailleLocal = tailleGlobalGenes - rankID*charge_unit ;
        tabLocal = new int[tailleLocal];

        for(int i= rankID*charge_unit;i<tailleGlobalGenes;i++){
            tabLocal[i - rankID * charge_unit] = tabfirstSol[i];
        }
    }
    // après la distribution du tableau local, le tableau local se trie par l'hill-climbing
    genepb0.setCurSol(tabLocal, tailleLocal);
    genepb0.hillClimbing(true,true,true);
    
    vector<int> curSolLocal = genepb0.getCurSol();
    memcpy(tabLocal,&curSolLocal[0],tailleLocal*sizeof(curSolLocal[0]));
    
    for(j=0;j<nbprocs;j++){
        if(j%2 == 0){// si c'est l'étape en paire
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
            // un processus impair recoit du suivant et lui reenvoie son nouveau tableau
			if (rankID%2==1 && rankID+1 < nbprocs){ 
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
    //compter le temp en final, le vrai temp est le longueur temps pris par un process parmi tous
    // MPI_Allreduce(&count, &totaltime, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Gather(&count,1,MPI_DOUBLE,recvtemp,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    //termination
    if(rankID==0){
        std::cout << "\n Impression tableau final" << std::endl;
        // cout<< "it takes " << totaltime/nbprocs << " s"<< endl;
        double* vrai_temp;
        vrai_temp = max_element(recvtemp,recvtemp+nbprocs);
        cout<< "Maximum time: "<< *vrai_temp<< " s"<<endl;
        delete [] recvtemp;
    } 
    print_tab_distrib_gather(tabLocal, tailleLocal,genepb0); 


    delete [] tabLocal;
    delete [] tabRecv;

    // delete [] tabLocal;
    MPI_Finalize();
    delete [] tabfirstSol;
    return 0;
}