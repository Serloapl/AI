#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>

#ifndef MARKOV_H
#define MARKOV_H

using namespace std;

class Markov {

public:

	int NumFlightPattern = 5;
	int NumMovementDirections = 9;

	vector<vector<double>> A;
	vector<vector<double>> B;
	vector<vector<double>> pi;
	vector<int> oseq;
	
	vector<vector<double>> alpha;
	vector<vector<double>> beta;
	vector<vector<double>> gamma;
	vector<vector<vector<double>>> gamma_di;

	int NextDirection;
	double NextDirectionProb;
	int species;

	
	void InitialiseModel();
	void TrainModel(vector<int > oseq);
	void VisualizeModel();
	void PredictNextMov(vector<int> oseq, int*, double*);
	void Predict_Trial(vector<int> oseq, int*, double*);
	void PredictViterbi(vector<int > oseq, int*, double*);
	double CalcProbSeq(vector<int> oseq);
	



};

#endif