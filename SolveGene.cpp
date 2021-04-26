#include <iostream>
#include <fstream>
#include <string.h>
#include <algorithm>
#include <random>
#include <map>
#include<time.h>
#include "SolveGene.hpp"


using namespace std;


/* Reads instance data. */
void SolveGene::readInstanceFromFile(const string& fileName) {
	ifstream infile;
	infile.exceptions(ifstream::failbit | ifstream::badbit);
	infile.open(fileName.c_str());

	// The input files follow the TSPLib "explicit" format.
	string str;
	char * pch;
	char* line;

	while (true) {
		getline(infile, str);
		line = strdup(str.c_str());
		pch = strtok(line, " :");
		if  (strcmp(pch, "DIMENSION") == 0){ //find the right line to get the nbGenes
			getline(infile, str);
			line = strdup(str.c_str());
			pch = strtok(NULL, " :");
			nbGenes = atoi(pch);
		} else if (strcmp(pch, "EDGE_WEIGHT_SECTION") == 0){
			break;
		}
	}
	// Distances from i to j
	wijMatrix.resize(nbGenes*nbGenes);
	for(int i = 0; i < nbGenes; i++) {
		for (int j = 0; j < nbGenes; j++) {
			infile >> wijMatrix[nbGenes * i+j];
		}
	}

	currentPossibility=-1;
	bestPossibility=-1;

	nbIterLS=0;

}

//To verify whether a solution is feasible, we need to verify whether it contains all elements
bool SolveGene::checkFeasibility() {
	vector<int> copyCurrentSol = curSol; //current solution list
	// copyCurrentSol.pop_back();//delete the last element
	std::sort (copyCurrentSol.begin(), copyCurrentSol.end());//sort soluiton and verify whether this solution is feasible
	bool b=true;
	int it=0;
	while (b && it< nbGenes){
		b = ((copyCurrentSol[it]-1)==it);
		// cout<<copyCurrentSol[it]<<endl;
		it++;
	}
	// cout<<endl;
	if (!b) cout << "Warning: infeasible solution found" << endl;
	return b;
}

void SolveGene::printStatus() {
	cout << "current solution with objective value " << currentPossibility << endl;

	for(int i=0; i <nbGenes; i++){
		cout << " " << curSol[i];
	}
	cout << endl;
}

void SolveGene::testGene(){
    curSol.clear();
    curSol.push_back(9);
    curSol.push_back(4);
    curSol.push_back(5);
    curSol.push_back(6);
    curSol.push_back(8);
    curSol.push_back(2);
    curSol.push_back(10);
    curSol.push_back(1);
    curSol.push_back(3);
    curSol.push_back(7);
    computeProbability();
    printStatus();
}

bool SolveGene::hillClimbingIter(bool swapMoves, bool revMoves, bool insertMoves){ //to do

    int modeImpr = -1 ; // pour dire quel mode de move à appliquer
	int param1 = -1; // pour dire que le paramètre 1 est
	int param2 = -1; // pour dire que le paramètre 2 est
	int i1, i2; // indice1 et indice2
	float delta = 0; // pour dire du progression maximale
	float temp = 0; // la progression temporaire
    int temp_var = 0; // variable temporaire pour finir l'échange entre deux éléments dans curSol
  
	if (swapMoves){
        //second case: between neighbours 
		for(i1=0; i1 < nbGenes-1 ; i1++){
			i2=i1+1;
			temp = sumProba(i2,i1) - sumProba(i1,i2);
			if (delta < temp ){
				delta=temp;
				param1 = i1; // i1 is going to exchange with i2
				param2 = i2;
				modeImpr=1; //the way number to improve
			}
		}
        // second case: on both sides of at least an element, exchange 
        for(i1 = 0; i1< nbGenes-2;i1++){
            for (i2 = i1+2; i2 < nbGenes; i2++)
            {   
                temp = 0; // reinitialiser temp
                for(int x = i1+1;x<i2+1;x++){
                    temp += sumProba(x,i1) - sumProba(i1,x);
                }
                for(int y = i1+1; y<i2;y++){
                    temp += sumProba(i2,y) - sumProba(y,i2);
                }
                if (delta < temp ){
                    delta=temp;
                    param1 = i1;
                    param2 = i2;
                    modeImpr=1;
                }
            }
        }
       
	}

    // the second mode: reverse
    if (revMoves) {
		// i2 == i1 + 1 same as insert and swap
		// i2 == i1 + 2 same as swap
		// i2 == i1 + 3 ...
		vector<int>::const_iterator first, last;
		for(i1 = 0; i1 < nbGenes - 3; i1++) {
			for (i2 = i1 + 2; i2 < nbGenes; i2++) {
				first = curSol.begin() + i1;
				last = curSol.begin() + i2 + 1;
				revSol.resize(i2 - i1 + 1);
				revSol.assign(first, last);
				reverseGenes();// firstly reverse and then to compute probability and sub
				temp = revSolutionProb - subSolutionProb;
				if (delta < temp ) {
					delta = temp;
					param1 = i1;
					param2 = i2;
					modeImpr = 2;
				}
			}
		}
	}

    if(insertMoves){
        // when i2 < i1, //i1 insert in head when i2 == 0
        for(i2 =0; i2 < nbGenes-1;i2++){
            for(i1=i2+1;i1<nbGenes;i1++){
                temp = 0; //initialize in each round
                for(int x = i2;x<i1;x++){
                    temp+= sumProba(i1,x) - sumProba(x,i1);
                }
                if(delta < temp){
                    delta = temp;
                    param1 = i1; //insert gene
                    param2 = i2; //insert position
                    modeImpr = 3;
                }
            }
        }

        // when i1 < i2, //i1 insert
        for(i2 = nbGenes-1;i2>0;i2--){
            for(i1 = i2-1;i1>=0;i1--){
                temp = 0;
                for(int x=i1+1;x<=i2;x++){
                    temp+= sumProba(x,i1) - sumProba(i1,x);
                }
                if(delta < temp){
                    delta = temp;
                    param1 = i1; //insert gene
                    param2 = i2; //insert position
                    modeImpr = 3;
                }
            }
        }

    }

    if (verboseMode) {
		cout << "best move : " << modeImpr << " param1 : " << param1 << " param2 : " << param2 << endl;
	}

    currentPossibility += delta; // update current possibility
	// update curSol de facon correspondante    
	if (modeImpr==1){
		temp_var = curSol[param2];
		curSol[param2]=curSol[param1];
		curSol[param1]=temp_var;
        countSwap++;
	}

    if (modeImpr==2){
		reverse(curSol.begin() + param1, curSol.begin() + param2 + 1);
        countRev++;
	}
    
    if( modeImpr == 3){
        if(param1>param2){
            curSol.insert(curSol.begin() + param2, curSol[param1]);
            curSol.erase(curSol.begin() + param1 + 1);
        }
        if(param1<param2){ //ex 1 -> 5
            curSol.insert(curSol.begin() + param2 + 1, curSol[param1]);
			curSol.erase(curSol.begin() + param1);
        }
        countInsert++;
    }
	//afficher le status
	if (verboseMode){ printStatus();}


	return delta > 0;
}

void SolveGene::hillClimbing(bool swapMoves, bool revMoves, bool insertMoves) {
    countSwap = 0;
    countInsert = 0;
    countRev =0;
	nbIterLS=0; // nbMaxIterLS = 20
	bool improving = true;
	while(improving && nbIterLS< nbMaxIterLS){
		improving = hillClimbingIter(swapMoves, revMoves, insertMoves);
		nbIterLS++;
	}
    //to print all results infos
    cout<<"After "<< nbIterLS << " iterations, optimization is finished"<<endl;
    cout<<"swapMoves: "<<countSwap<<endl;
    cout<<"revMoves: "<< countRev<<endl;
    cout<<"InsertMoves: "<<countInsert<<endl;
}

//Because we will only use this function after reversed vector,so don't worry about the order mistake
void SolveGene::computeReversedPossibility(){
    revSolutionProb = 0; // the reversed solution probability
	subSolutionProb = 0; // I think it's the normal order solution probability
	for(int i = 0; i < revSol.size() - 1; i++){
        for(int j=i+1;j<revSol.size();j++){
            revSolutionProb += ProbaIbeforeJ(revSol[i]-1,revSol[j]-1);
        }
    }
     
	for(int i = revSol.size() - 1; i >= 1; i--) {
        for(int j=i-1;j>=0;j--){
            subSolutionProb += ProbaIbeforeJ(revSol[i]-1,revSol[j]-1);
        }
	}
}

//reverse genes in ranking
void SolveGene::reverseGenes(){
    reverse(revSol.begin(), revSol.end()); //revSol is used to stock partially curSol and then reverse it
	computeReversedPossibility();
}

// public part
void SolveGene::updateBestSolution() {
	if (bestPossibility<0 || bestPossibility < currentPossibility){ //when the currentProba > bestProba, replace it
			bestSolution = curSol;
			bestPossibility=currentPossibility;
	}
}

//the correct way to calculate: wij*xij
void SolveGene::computeProbability() {
	currentPossibility = 0;
	for(int i=0; i < nbGenes-1; i++){
        for(int j=i+1;j<nbGenes;j++){
            currentPossibility += ProbaIbeforeJ(curSol[i]-1,curSol[j]-1);
        }
    }    
        // cout<< ProbaIbeforeJ(curSol[i]-1,curSol[i+1]-1)<<" ";
    // cout<<endl;
}

//display instance
void SolveGene::displayInstance() {
	cout << "nb of cities : " << nbGenes << endl;
	cout << "probability matrix : " << endl;

	for(int i1=0; i1 < nbGenes; i1++){
		for(int i2=0; i2 < nbGenes; i2++) cout << " " << ProbaIbeforeJ(i1,i2);
		cout << endl;
	}
}

//naive init
void SolveGene::naiveInit() {
	curSol.resize(nbGenes);
	for(int i=0; i < nbGenes; i++) curSol[i]=i+1;
	computeProbability();
    // printStatus();
}

//greedy init
void SolveGene::greedyInit() {
	curSol.resize(nbGenes); //firstly clear curSol
	curSol[0] = 2; //we set first gene 1
	int cur_index = 0;//current index starts from 0
	double cur_dist; //cur_dist is wi,j
	double max_proba; // we always wanna get the max local proba
	for (int i = 1; i < nbGenes; i++) {
		max_proba = INT_MIN + 0.01;
		// find the biggest proba gene
		for (int j = 0; j < nbGenes; j++) {
			// check if gene already in vector, find function allows to get (index+1) of j+1 in curSol
            if (find(curSol.begin(), curSol.begin() + i, j+1) != curSol.begin() + i) continue;
                cur_dist = ProbaIbeforeJ(curSol[i-1]-1, j);
			    if (max_proba < cur_dist) {
				    cur_index = j;
				    max_proba = cur_dist;
			    }
		        curSol[i] = cur_index+1;   
			
		}
	}
	computeProbability();
    printStatus();
}

//randomizedConstr
void SolveGene::randomizedConstr() {
	curSol.resize(nbGenes);
	random_device rd;
    mt19937 g(rd());
	// for (int i = 0; i < nbMaxRestart + 1; i++) {
	    curSol.clear();
		curSol.resize(nbGenes);
		for(int i=0; i < nbGenes; i++) curSol[i]=i+1;
		shuffle(begin(curSol), end(curSol) - 1, g);
		computeProbability();
		updateBestSolution();
    // computeProbability();
    printStatus();
}

void SolveGene::graspHC(bool swapMoves, bool revMoves, bool insertMoves) {
	for (int i = 0; i < nbMaxRestart + 1; i++) {
        clock_t start,end; // calculate the time by two slots
		greedyInit();
        start = clock();
		hillClimbing(swapMoves, revMoves, insertMoves);
        end = clock();
		updateBestSolution();
        cout<<"This turn takes "<<(double)(end-start)/CLOCKS_PER_SEC<<" s"<<endl;
	}
}

void SolveGene::ilsHC(bool swapMoves, bool revMoves, bool insertMoves) {
	for (int i = 0; i < nbMaxRestart + 1; i++) {
        ilsHCtime = 0;
		clock_t start,end; // calculate the time by two slots
		randomizedConstr();
        start = clock();
		hillClimbing(swapMoves, revMoves, insertMoves);
        end = clock();
		updateBestSolution();
        ilsHCtime = (double)(end-start)/CLOCKS_PER_SEC;
        cout<<"This turn takes "<<ilsHCtime<<" s"<<endl;
	}
}

void SolveGene::simulatedAnnealing() {
    nbIterSimulated = 0;
	auto t_start = 5000.0;
	auto t_end = 1e-8;//10^-8
	auto q = 0.8;
	random_device rd; //random_device 本身是均匀分布整数随机数生成器，通常仅用于播种; 如果是正态分布，那么有两个参数：均值和标准差
	mt19937 g(rd());//mt19937 是用来生成高性能随机数的算法，mt是Mernne Twister，19937表示最大位数2^19937 -1种可能
	uniform_int_distribution<> dis_int(0, 2);//在范围 0 ，2 中出现所有的整数都是等可能性的
	uniform_real_distribution<double> dis_real(0.0,1.0);//这也是一样，所有的小数出现可能性相同,0-1随机数
	int mode;
	bool swapMoves = false, insertMoves = false, revMoves = false;
	while (t_start > t_end) {
        swapMoves = false; insertMoves = false; revMoves = false;
		for (int i = 0; i < nbMaxRestart; i++) {
			vector<int> tmpSol = curSol;
			int tmpSolutionCost = 0;
			for(int i = 0; i < nbGenes-1; i++) {
                for(int j=i+1;j<nbGenes;j++)
                    tmpSolutionCost += ProbaIbeforeJ(tmpSol[i]-1,tmpSol[j]-1);
            }
            //get any mode
			mode = dis_int(g); // get a int between (0,2) 0,1,2
			switch (mode) {
				case 0:
					swapMoves = true; 
					break;
				case 1:
					revMoves = true;
					break;
				case 2:
					insertMoves = true;
					break;
				default:
					break;
			}
            cout<<"which move mode applying:"<<swapMoves<<revMoves<<insertMoves<<endl;
            //if there is no improvement for this move
			if (!hillClimbingIter(swapMoves, revMoves, insertMoves)) { //if there's an improvement
				if (exp((currentPossibility-tmpSolutionCost) / t_start) <= dis_real(g)) {
					curSol = tmpSol;
					computeProbability();
				} else updateBestSolution();
			} else updateBestSolution(); 
		}
		t_start *= q;
        //nombre de l'iterations
        nbIterSimulated++;
	}
}


