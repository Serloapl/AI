#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include "Markov.h"
#include <iomanip>
#include <algorithm>


using namespace std;

vector<vector< double> > multiplication(vector<vector< double> > A, vector<vector< double> > B)
{

	//function for matrix multiplication
	int N_rows_A = A.size();
	int N_cols_A = A[0].size();
	int N_rows_B = B.size();
	int N_cols_B = B[0].size();

	vector<vector< double> > C(N_rows_A, vector<double>(N_cols_B)); //Matrix Output

	if (N_cols_A != N_rows_B)
		cout << "Error: dimentions do not agree";
	else
	{
		for (int i = 0; i < N_rows_A; i++)
		{
			for (int j = 0; j < N_cols_B; j++)
			{
				C[i][j] = 0;
				for (int k = 0; k < N_cols_A; k++)
				{
					C[i][j] += A[i][k] * B[k][j];
				}
			}
		}
		return C;
	}
}

void Markov::VisualizeModel()
{	// Function for matrix visualization

	int N_rows_A = this->A.size();
	int N_cols_A = this->A[0].size();

	cerr << "A = [";

	for (int i = 0; i < N_rows_A; i++)
	{
		cerr << std::fixed << std::setprecision(3) << this->A[i][0];
		for (int j = 1; j < N_cols_A; j++)
		{
			cerr << setw(6) << " " << std::fixed << std::setprecision(3) << this->A[i][j];
		}
		cerr << ";" << endl << "     ";
	}

	cerr << endl;
	int N_rows_B = this->B.size();
	int N_cols_B = this->B[0].size();
	cerr << "B = [";

	for (int i = 0; i < N_rows_B; i++)
	{
		cerr << this->B[i][0];
		for (int j = 1; j < N_cols_B; j++)
		{
			cerr << " " << this->B[i][j];
		}
		cerr << ";" << endl << "     ";
	}
	cerr << endl;

	cerr << endl;
	int N_rows_pi = this->pi.size();
	int N_cols_pi = this->pi[0].size();
	cerr << "Pi =[";

	for (int i = 0; i < N_rows_pi; i++)
	{
		cerr << this->pi[i][0];
		for (int j = 1; j < N_cols_pi; j++)
		{
			cerr << " " << this->pi[i][j];
		}
		cerr << ";" << endl << "     ";
	}
	cerr << endl;
}

void Markov::InitialiseModel()
{
	// ---------------------------------------------------------Data reading (beginning)
	string line;
	vector<double> transition_aux;	//auxiliar vectors for reading
	vector<double> emission_aux;	//auxiliar vectors for reading
	vector<double> pi_aux;			//auxiliar vectors for reading


	////Vector reordenation into matrices
	int N_states = this->NumFlightPattern;
	int N_observ = this->NumMovementDirections;

	vector<vector< double> > A(N_states, vector<double>(N_states));	
	A[0] = {0.201, 0.200, 0.200, 0.200, 0.198};
	A[1] = {0.198, 0.198, 0.202, 0.199, 0.204};
	A[2] = {0.195, 0.199, 0.203, 0.198, 0.204};
	A[3] = {0.196, 0.204, 0.196, 0.204, 0.200};
	A[4] = {0.203, 0.198, 0.201, 0.200, 0.198};
	
	/*for (size_t i = 0; i < N_states; i++)
	{
		double total = 0;
		for (size_t j = 0; j < N_states; j++)
		{
			A[i][j] = (1 / N_states) +(((double)rand() / (RAND_MAX)) / 100);
			total += A[i][j];
		}

		for (size_t j = 0; j < N_states; j++)
		{
			A[i][j] = A[i][j] / total;
		}
	}*/

	vector<vector< double> > B(N_states, vector<double>(N_observ));

	B[0] = {0.106, 0.114, 0.107, 0.114, 0.110, 0.110, 0.112, 0.115, 0.115};
	B[1] = {0.109, 0.110, 0.112, 0.114, 0.114, 0.108, 0.111, 0.113, 0.109};
	B[2] = {0.112, 0.113, 0.105, 0.112, 0.114, 0.113, 0.110, 0.113, 0.109};
	B[3] = {0.114, 0.115, 0.115, 0.109, 0.108, 0.113, 0.109, 0.111, 0.107};
	B[4] = {0.111, 0.109, 0.114, 0.107, 0.114, 0.108, 0.113, 0.111, 0.112 };

	/*for (size_t i = 0; i < N_states; i++)
	{
		double total = 0;
		for (size_t j = 0; j < N_observ; j++)
		{
			B[i][j] = (1 / N_observ) +(((double)rand() / (RAND_MAX)) / 100);
			total += B[i][j];
		}

		for (size_t j = 0; j < N_observ; j++)
		{
			B[i][j] = B[i][j] / total;
		}
	}*/

	vector<vector< double> > pi(1, vector<double>(N_states));
	pi[0] = { 0.197, 0.198, 0.201, 0.203, 0.201 };
	/*double total = 0;
		for (size_t j = 0; j < N_states; j++)
		{
			pi[0][j] = (1 / N_states) +(((double)rand() / (RAND_MAX)) / 100);
			total += pi[0][j];
		}

		for (size_t j = 0; j < N_states; j++)
		{
			pi[0][j] = pi[0][j] / total;
		}*/

	
	this->A = A;
	this->B = B;
	this->pi = pi;
}

void Markov::TrainModel(vector<int> oseq)
{
	//------------------------------------------------------------- Computing alpha, Beta, Gamma parameters (beginning)

	int N_rows_A = this->A.size();
	int N_cols_A = this->A[0].size();
	int N_cols_B = B[0].size();
	int N_rows_B = this->B.size();
	int N_cols_oseq = oseq.size();
	int N_rows_pi = this->pi.size();
	int N_cols_pi = this->pi[0].size();

	vector<vector<double> > output_aux; //Result matrix
	vector<vector<double>> alpha(N_cols_A, vector<double>(N_cols_oseq));
	vector<vector<double>> beta(N_cols_A, vector<double>(N_cols_oseq));
	vector<vector<double>> gamma(N_cols_A, vector<double>(N_cols_oseq));

	vector<vector<vector<double>>> gamma_di(N_cols_oseq, vector<vector<double>>(N_cols_A, vector<double>(N_cols_A)));


	int N = N_cols_A; // number of states
	int M = N_cols_B; // number os possible different observations
	int T = N_cols_oseq;
	int iters = 1;
	int maxIters = 200;
	double logProb = -1 * 9999999998;
	double oldLogProb = -1 * 9999999999;
	vector<double> c(T);

	/*cerr << "A1 = [";

	for (int i = 0; i < N_rows_A; i++)
	{
		cerr << std::fixed << std::setprecision(3) << this->A[i][0];
		for (int j = 1; j < N_cols_A; j++)
		{
			cerr << setw(6) << " " << std::fixed << std::setprecision(3) << this->A[i][j];
		}
		cerr << ";" << endl << "     ";
	}

	cerr << endl;
	
	cerr << "B1 = [";

	for (int i = 0; i < N_rows_B; i++)
	{
		cerr << this->B[i][0];
		for (int j = 1; j < N_cols_B; j++)
		{
			cerr << " " << this->B[i][j];
		}
		cerr << ";" << endl << "     ";
	}
	cerr << endl;

	cerr << endl;

	cerr << "Pi1 =[";

	for (int i = 0; i < N_rows_pi; i++)
	{
		cerr << this->pi[i][0];
		for (int j = 1; j < N_cols_pi; j++)
		{
			cerr << " " << this->pi[i][j];
		}
		cerr << ";" << endl << "     ";
	}
	cerr << endl;*/



	while (iters < maxIters && logProb > oldLogProb)
	//while(iters < maxIters)
	{
		//---------------------------------- Alpha computation (beginning)
		// compute alpha0(i) - alpha[0][i]
		for (size_t i = 0; i < N; i++)
		{
			alpha[i][0] = pi[0][i] * B[i][oseq[0]];
			c[0] += alpha[i][0];
		}

		//scale the alpha0(i) - alpha[0][i]
		c[0] = 1 / c[0];
		for (size_t i = 0; i < N; i++)
		{
			alpha[i][0] = c[0] * alpha[i][0];
		}

		//compute alphat(i) - alpha[t][i]
		for (size_t t = 1; t < T; t++)
		{
			c[t] = 0;
			for (size_t i = 0; i < N; i++)
			{
				alpha[i][t] = 0;
				for (size_t j = 0; j < N; j++)
				{
					alpha[i][t] += alpha[j][t - 1] * A[j][i];
				}
				alpha[i][t] = alpha[i][t] * B[i][oseq[t]];
				c[t] += alpha[i][t];
			}

			//scale alphat(i)-alpha[t][i]
			c[t] = 1 / c[t];
			for (size_t i = 0; i < N; i++)
			{
				alpha[i][t] = c[t] * alpha[i][t];
			}
		}

		//visualization(alpha);
		//---------------------------------- Alpha computation (end)

		//---------------------------------- Beta computation (beginning)

		//Let betaT-1(i) - beta[i][T-1]
		for (size_t i = 0; i < N; i++)
		{
			beta[i][T - 1] = c[T - 1];
		}

		//beta-pass
		for (int t = T - 2; t >= 0; t--)
		{
			for (size_t i = 0; i < N; i++)
			{
				beta[i][t] = 0;
				for (size_t j = 0; j < N; j++)
				{
					beta[i][t] += A[i][j] * B[j][oseq[t + 1]] * beta[j][t + 1];
				}
				//scale betat(i)- beta[i][t] with same scale factor as alphat(i) - alpha[i][t]
				beta[i][t] = c[t] * beta[i][t];
			}
		}
		//---------------------------------- beta computation (end)

		//---------------------------------- gamma computation (beginning)
		double denom;
		for (size_t t = 0; t < T - 1; t++)
		{
			denom = 0;
			for (size_t i = 0; i < N; i++)
			{
				for (size_t j = 0; j < N; j++)
				{
					denom += alpha[i][t] * A[i][j] * B[j][oseq[t + 1]] * beta[j][t + 1];
				}
			}
			for (size_t i = 0; i < N; i++)
			{
				gamma[i][t] = 0;
				for (size_t j = 0; j < N; j++)
				{
					gamma_di[t][i][j] = (alpha[i][t] * A[i][j] * B[j][oseq[t + 1]] * beta[j][t + 1]) / denom;
					gamma[i][t] += gamma_di[t][i][j];
				}
			}
		}
		//---------------------------------- gamma computation (end)

		//------------------------------------------------------------- Computing alpha, Beta, Gamma parameters (end)

		//------------------------------------------------------------- reestimation of A,B,pi parameters (beginning)
		//---------------------------------- reestimation of pi parameters (beginning)

		for (size_t i = 0; i < N; i++)
		{
				pi[0][i] = gamma[i][0] + (((double)rand() / (RAND_MAX)) / 1000000000);
		}

		//---------------------------------- reestimation of pi parameters (end)
		//---------------------------------- reestimation of A parameters (beginning)

		double numer;
		for (size_t i = 0; i < N; i++)
		{
			for (size_t j = 0; j < N; j++)
			{
				numer = denom = 0;
				for (size_t t = 0; t < T - 1; t++)
				{
					numer += gamma_di[t][i][j];
					denom += gamma[i][t];
				}
				A[i][j] = numer / denom + (((double)rand() / (RAND_MAX)) / 1000000000);
			}
		}

		//---------------------------------- reestimation of A parameters (end)
		//---------------------------------- reestimation of B parameters (beginning)

		for (size_t i = 0; i < N; i++)
		{
			for (size_t j = 0; j < M; j++)
			{
				numer = denom = 0;
				for (size_t t = 0; t < T - 1; t++)
				{
					if (oseq[t] == j)
					{
						numer += gamma[i][t];
					}
					denom += gamma[i][t];
				}
				B[i][j] = numer / denom + (((double)rand() / (RAND_MAX)) / 1000000000);
			}
		}

		//---------------------------------- reestimation of B parameters (end)
		//------------------------------------------------------------- reestimation of A,B,pi parameters (end)

		//------------------------------------------------------------- log() computation and iteration decision (beginning)
		//---------------------------------- log() (beginning)
		oldLogProb = logProb;
		logProb = 0;
		for (size_t i = 0; i < T; i++)
		{
			logProb += log(c[i]);
		}
		logProb = -1 * logProb;
		//---------------------------------- log() (end)
		//---------------------------------- iteration decision (beginning)
		iters++;
		//cout << "oldLogprob = " << oldLogProb << endl;
		//cout << "Logprob = " << logProb << endl;
		//---------------------------------- iteration decision (end)
		//------------------------------------------------------------- log() computation and iteration decision (end)

		//cout << iters << endl;
		/*std::cerr << "BW Iter = " << iters << std::endl;
		cerr << "A_BW" << iters << "=[";

		for (int i = 0; i < N_rows_A; i++)
		{
			cerr << std::fixed << std::setprecision(3) << this->A[i][0];
			for (int j = 1; j < N_cols_A; j++)
			{
				cerr << setw(6) << " " << std::fixed << std::setprecision(3) << this->A[i][j];
			}
			cerr << ";" << endl << "     ";
		}

		cerr << endl;

		cerr << "B_BW" << iters << "=[";

		for (int i = 0; i < N_rows_B; i++)
		{
			cerr << this->B[i][0];
			for (int j = 1; j < N_cols_B; j++)
			{
				cerr << " " << this->B[i][j];
			}
			cerr << ";" << endl << "     ";
		}
		cerr << endl;

		cerr << endl;

		cerr << "Pi_BW" << iters << "=[";

		for (int i = 0; i < N_rows_pi; i++)
		{
			cerr << this->pi[i][0];
			for (int j = 1; j < N_cols_pi; j++)
			{
				cerr << " " << this->pi[i][j];
			}
			cerr << ";" << endl << "     ";
		}
		cerr << endl;*/

	}


	this -> alpha = alpha;
	this -> beta = beta;
	this -> gamma = gamma;
	this -> gamma_di = gamma_di;



}

void Markov::PredictNextMov(vector <int> oseq, int* IndexMovement, double* ProbDirection )
{
	int M = this->NumMovementDirections;
	int N = this->NumFlightPattern;
	int T = oseq.size();
	vector<vector<double>> alpha_aux(N, vector<double>(M));
	vector<vector<double>> alpha(N, vector<double>(T));
	vector<double> c_aux(M);
	vector<double> c(T);

	for (size_t i = 0; i < N; i++)
	{
		alpha[i][0] = pi[0][i] * B[i][oseq[0]];
		c[0] += alpha[i][0];
	}

	//scale the alpha0(i) - alpha[0][i]
	c[0] = 1 / c[0];
	for (size_t i = 0; i < N; i++)
	{
		alpha[i][0] = c[0] * alpha[i][0];
	}

	//compute alphat(i) - alpha[t][i]
	for (size_t t = 1; t < T; t++)
	{
		c[t] = 0;
		for (size_t i = 0; i < N; i++)
		{
			alpha[i][t] = 0;
			for (size_t j = 0; j < N; j++)
			{
				alpha[i][t] += alpha[j][t - 1] * A[j][i];
			}
			alpha[i][t] = alpha[i][t] * B[i][oseq[t]];
			c[t] += alpha[i][t];
		}

		//scale alphat(i)-alpha[t][i]
		c[t] = 1 / c[t];
		for (size_t i = 0; i < N; i++)
		{
			alpha[i][t] = c[t] * alpha[i][t];
		}
	}





	/*std::cerr << "Alpha = :";
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < T; j++)
		{
			std::cerr << alpha[i][j] << " ";
		}
		std::cerr << std::endl;
	}*/
	double denom = 0;

	// For every future direction
	for (size_t k = 0; k < M; k++)
	{
		//compute alphat(i) - alpha[t][i]
		c_aux[k] = 0;
		for (size_t i = 0; i < N; i++)
		{
			alpha_aux[i][k] = 0;
			for (size_t j = 0; j < N; j++)
			{
				alpha_aux[i][k] += alpha[j][T - 1] * A[j][i];
			}
			alpha_aux[i][k] = alpha_aux[i][k] * B[i][k];
			c_aux[k] += alpha_aux[i][k];
		}
		//std::cerr << c[k] << " ";
		denom += c_aux[k];
	}
	//std::cerr << std::endl;


	/*for (size_t i = 0; i < N; i++)
	{
		denom += alpha_aux[i][T - 1];
	}*/

	//std::cerr << "Denom =" << denom << std::endl;

	for (size_t k = 0; k < M; k++)
	{
		c_aux[k] = c_aux[k]/denom;
		//std::cerr << c[k] << " ";
	}
	//std:cerr << std::endl;


	////extract maximum and its direction

	auto max = max_element(begin(c_aux), end(c_aux));
	//this->NextDirectionProb = double(*max);
	*ProbDirection = double(*max);
	//this->NextDirection = distance(begin(c), max);
	*IndexMovement = distance(begin(c_aux), max);

}

void Markov::Predict_Trial(vector<int> oseq, int* IndexMovement, double* ProbDirection)
{
	int M = this->NumMovementDirections;
	int N = this->NumFlightPattern;
	int T = oseq.size();

	vector<vector<double>> curr_state;
	vector<vector<double>> curr_emiss;

	//vector<vector<double>> alpha_aux(N, vector<double>(M));
	//vector<double> c(M);

	curr_state = multiplication(pi, A);
	
	for (size_t t = 0; t < T - 1; t++)
	{
		curr_state = multiplication(curr_state, A);
	}
	
	curr_emiss = multiplication(curr_state, B);

	/*for (size_t i = 0; i < M; i++)
	{
		std::cerr << curr_emiss[0][i] << " ";
	}
	std::cerr << std::endl;*/

	auto max = max_element(begin(curr_emiss[0]), end(curr_emiss[0]));
	*ProbDirection = double(*max);
	*IndexMovement = distance(begin(curr_emiss[0]), max);

}

void Markov::PredictViterbi(vector<int> oseq, int* IndexMovement, double* ProbDirection)
{
	int M = this->NumMovementDirections;
	int N = this->NumFlightPattern;
	int T = oseq.size();

	//vector<vector<double> > output_aux; //Result matrix
	vector<vector<double>> delta(N, vector<double>(T));
	//vector<vector<int>> delta_idx(N_cols_A, vector<int>(N_cols_oseq));

	vector<double> aux;

	// compute delta0(i) - delta[0][i]
	for (size_t i = 0; i < N; i++)
	{
		delta[i][0] = pi[0][i] * B[i][oseq[0]];
	}

	//compute deltat(i) - delta[t][i]
	for (size_t t = 1; t < T; t++)
	{
		for (size_t i = 0; i < N; i++)
		{
			delta[i][t] = 0;
			for (size_t j = 0; j < N; j++)
			{
				aux.push_back(A[j][i] * B[i][oseq[t]] * delta[j][t - 1]);
			}
			auto max = max_element(begin(aux), end(aux));
			delta[i][t] = double(*max);
			//delta_idx[i][t] = distance(begin(aux), max);
			aux.clear();
		}
	}

	//Finding index of larger delta[i][T - 1]
	int index = 0;
	double acum = -1;

	for (size_t i = 0; i < N; i++)
	{
		if (delta[i][T - 1] > acum)
		{
			index = i;
			acum = delta[i][T - 1];
		}
	}

	auto max = max_element(begin(B[index]), end(B[index]));
	*ProbDirection = double(*max);
	*IndexMovement = distance(begin(B[index]), max);

	//this->NextDirectionProb = double(*max);

	//this->NextDirection = distance(begin(B[index]), max);

}

double Markov::CalcProbSeq(vector<int> oseq)
{
	int M = this->NumMovementDirections;
	int N = this->NumFlightPattern;
	int T = oseq.size();

	std::vector<std::vector<double>> alpha(N, std::vector<double>(T));

	vector<double> c(T);

	//this->VisualizeModel();

	// compute alpha0(i) - alpha[0][i]
	for (size_t i = 0; i < N; i++)
	{
		alpha[i][0] = pi[0][i] * B[i][oseq[0]];
		c[0] += alpha[i][0];
	}

	//scale the alpha0(i) - alpha[0][i]
	c[0] = 1 / c[0];
	for (size_t i = 0; i < N; i++)
	{
		alpha[i][0] = c[0] * alpha[i][0];
	}

	
	//compute alphat(i) - alpha[t][i]
	for (size_t t = 1; t < T; t++)
	{
		//std::cerr << "alpha col = " << t << std::endl;
		c[t] = 0;
		for (size_t i = 0; i < N; i++)
		{
			//std::cerr << "alpha row = " << i << std::endl;
			alpha[i][t] = 0;
			for (size_t j = 0; j < N; j++)
			{
				//std::cerr << "alpha update" << j;
				alpha[i][t] += alpha[j][t - 1] * A[j][i];
			}
			alpha[i][t] = alpha[i][t] * B[i][oseq[t]];
			c[t] += alpha[i][t];
		}
		//std::cerr << "He entrado a CalProb(4)" << std::endl;

		//std:cerr << "LAS CS VALEN:  " << c[t] << std::endl;

		//scale alphat(i)-alpha[t][i]
		//std::cerr << "c[" << t << "] = " << c[t] << std::endl;
		if (c[t] != 0)
		{
			c[t] = 1 / c[t];
		}
		else
		{
			c[t] = 100000000000;
		}
		for (size_t i = 0; i < N; i++)
		{
			alpha[i][t] = c[t] * alpha[i][t];
		}
	}
	//std::cerr << "He entrado a CalProb(5)" << std::endl;

	double LogProb = 0;
	for (size_t t = 0; t < T; t++)
	{
		LogProb += log(c[t]);
	}
	return -1 * LogProb;

	//std::cerr << "He entrado a CalProb(6)" << std::endl;

	/*double prob = 0;
	for (size_t i = 0; i < N; i++)
	{
		prob += alpha[i][T - 1];
	}*/


}