#include <iostream>
#include <algorithm>
#include <random>
#include "omp.h"
#include <vector>
#include <time.h>
#include <fstream>
#include <chrono>

using namespace std;

vector<int> GenerateRandomRanking(int nbGenes){
    vector<int> curSol;
    curSol.resize(nbGenes);
    for (auto i = 0; i < nbGenes; i++)
    {
        /* code */
        curSol[i]=i;
    }
    random_device rd;
    mt19937 g(rd());
    shuffle(begin(curSol), end(curSol), g);
    return curSol;
}
// Générer une matrice contenant le nombre donné des séquences  
vector<int> generateMatrix(int taille,int nbGenes){
    vector<int> Matrix;
    Matrix.resize(taille*nbGenes);
    #pragma omp parallel
    {
        #pragma omp for
        for(int i=0;i<taille;i++){
            vector<int> tempor = GenerateRandomRanking(nbGenes);
            for(int j = 0;j<nbGenes;j++){
                Matrix[i*nbGenes+j] = tempor[j];
            }
        }
    }
    
    return Matrix;
}
// la fonction pour générer la table des probabilités
vector<double> compterAllPairs(vector<int> Matrix,int nbGenes,int nbrankings){
    vector<double> wij_table; // la table des probabilites
    wij_table.resize(nbGenes*nbGenes); //resize la table wij
    // double start = omp_get_wtime();
    #pragma omp parallel
    {
        #pragma omp for
        for(int i = 0;i < nbrankings;i++){ 
            for(int j=0; j< nbGenes-1;j++){ //j th gene in i th ranking
                for(int z = j+1; z<nbGenes;z++){ // z th gene in i th ranking
                    int JthGene = Matrix[nbGenes*i+j];
                    int ZthGene = Matrix[nbGenes*i+z];
                    wij_table[JthGene*nbGenes+ZthGene]+=1;
                }
            }
        }
    }
     #pragma omp parallel
    {
            #pragma omp for
        for (int i = 0; i < nbGenes; i++)
        {
            for (int j = 0; j < nbGenes; j++)
            {
                wij_table[i*nbGenes+j] /= nbrankings; 
            }
            
        }
    }
    // double end = omp_get_wtime();
    // cout<<end-start<<endl;
    return wij_table;
    
}

void sauveGarderInstance(vector<double> Matrix,int nbGenes,string filename){
    //    char data[100];
 
   // 以写模式打开文件
   ofstream outfile;
   outfile.open(filename);

//    #pragma omp parallel for
    for (int i = 0; i < nbGenes; i++)
    {
        for (int j = 0; j < nbGenes; j++)
        {
            if (j!=nbGenes-1)
            {
                // 向文件写入用户输入的数据
                outfile << Matrix[i*nbGenes+j] << "\t";
            }else{
                outfile << Matrix[i*nbGenes+j];
            }
        }
        outfile<<endl;
    }
    
   // 关闭打开的文件
   outfile.close();

}

int main(){
    int nbGenes = 1000;
    int nbrankings = 20000;
    vector<int> Matrix = generateMatrix(nbrankings,nbGenes); //2000 classements de taille 1000

    // for(int i=0;i< nbrankings;i++){ // i th ranking
    //     for(int j=0;j<nbGenes;j++){ // i th ranking's j th gene
    //         cout<<Matrix[i*nbGenes+j]<<",";
    //     }
    //     cout<<endl;
    // }
    // for (int i = 0; i < 5; i++)
    // {
    //     cout<<endl;

    // }
    auto start = std::chrono::steady_clock::now();
    //compter table
    vector<double> wij_table = compterAllPairs(Matrix,nbGenes,nbrankings);
    sauveGarderInstance(wij_table,nbGenes,"afile.dat");
    //compter le temp
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    double count = elapsed_seconds.count();
    cout<<count<<endl;

    // for (int i = 0; i < nbGenes; i++)
    // {
    //     for (int j = 0; j < nbGenes; j++)
    //     {   
    //         if (j!=nbGenes-1)
    //         {
    //             cout<<wij_table[i*nbGenes+j]<<",";
    //         }else{
    //             cout<<wij_table[i*nbGenes+j];
    //         }
            
            
    //     }
    //     cout<<endl;
        
    // }
    
}
