//HMM2 SERGIO LOPEZ AND DIEGO YUS
//ARTIFICIAL INTELLIGENCE DD2380


#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>

using namespace std;
vector<vector< double> > multiplication(vector<vector< double> > A, vector<vector< double> > B);
void visualization(vector<vector< double> > A);

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

void visualization(vector<vector< double> > A)
{	// Function for matrix visualization

	int N_rows_A = A.size();
	int N_cols_A = A[0].size();

	cout << "A = (";

	for (int i = 0; i < N_rows_A; i++)
	{
		cout << A[i][0];
		for (int j = 1; j < N_cols_A; j++)
		{
			cout << ", " << A[i][j];
		}
		cout << ";" << endl << "     ";
	}
	//cout << ")";
}

int main()
{
	// ---------------------------------------------------------Data reading (beginning)
	string line;
	vector<double> transition_aux;	//auxiliar vectors for reading
	vector<double> emission_aux;	//auxiliar vectors for reading
	vector<double> pi_aux;			//auxiliar vectors for reading
	vector<double> oseq_aux;			//auxiliar vectors for reading

	double d;


	//Getting Matrix A Transition
	getline(cin, line);
	stringstream lineStream(line);

	while (lineStream >> d)
	{
		transition_aux.push_back(d);
		//cout << d << endl;
	}

	//Getting Matrix B Emission
	getline(cin, line);
	lineStream = stringstream(line);
	while (lineStream >> d)
	{
		emission_aux.push_back(d);
	}

	//Getting Matrix Pi (Initial)
	getline(cin, line);
	lineStream = stringstream(line);
	while (lineStream >> d)
	{
		pi_aux.push_back(d);
	}

	//Getting the observation sequence
	getline(cin, line);
	lineStream = stringstream(line);
	while (lineStream >> d)
	{
		oseq_aux.push_back(d);
	}



	//Vector reordenation into matrices
	//Matrix A Transition
	int N_rows_A = int(transition_aux[0]);
	int N_cols_A = int(transition_aux[1]);
	int k = 2;

	vector<vector< double> > A(N_rows_A, vector<double>(N_cols_A));
	for (size_t i = 0; i < N_rows_A; i++)
	{
		for (size_t j = 0; j < N_cols_A; j++)
		{
			A[i][j] = transition_aux[k];
			k++;
		}
	}
	//visualization(A);

	//Matrix B Emission
	int N_rows_B = int(emission_aux[0]);
	int N_cols_B = int(emission_aux[1]);
	k = 2;

	vector<vector<double>> B(N_rows_B, vector<double>(N_cols_B));
	for (size_t i = 0; i < N_rows_B; i++)
	{
		for (size_t j = 0; j < N_cols_B; j++)
		{
			B[i][j] = emission_aux[k];
			k++;
		}
	}
	//visualization(B);

	//Matrix(vector) Pi (Initial)
	int N_cols_pi = int(pi_aux[1]);
	vector<vector<double > > pi(1, vector<double>(N_cols_pi));

	for (size_t i = 0; i < N_cols_pi; i++)
	{
		pi[0][i] = pi_aux[i + 2];
	}
	//visualization(pi);

	//Matrix(vector) oseq (Sequence of observation)
	int N_cols_oseq = int(oseq_aux[0]);
	vector<int > oseq;

	for (size_t i = 1; i < N_cols_oseq + 1; i++)
	{
		oseq.push_back(oseq_aux[i]);
	}

	

	//--------------------------------------------------------------Data reading (End)

	//------------------------------------------------------------- Computing alpha, Beta, Gamma parameters (beginning)

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
	double oldLogProb = -1*9999999999;
	vector<double> c(T);



	while (iters < maxIters && logProb > oldLogProb)
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
			pi[0][i] = gamma[i][0];
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
				A[i][j] = numer / denom;
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
				B[i][j] = numer / denom;
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
	}

	

	//visualization(A);

	cout << A.size() << " " << A[0].size();
	for (size_t i = 0; i < A.size(); i++)
	{
		for (size_t j = 0; j < A[0].size(); j++)
		{
			cout << " " << A[i][j];
		}
	}

	cout << endl;

	cout << B.size() << " " << B[0].size();
	for (size_t i = 0; i < B.size(); i++)
	{
		for (size_t j = 0; j < B[0].size(); j++)
		{
			cout << " " << B[i][j];
		}
	}

	//visualization(B);


	//-----------------------------------------------------------------Output (beginning)


	//----------------------------------------------------------------------Output (end)


	system("pause");
	return 0;
}