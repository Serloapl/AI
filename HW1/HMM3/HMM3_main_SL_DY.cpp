//HMM3 SERGIO LOPEZ AND DIEGO YUS
//ARTIFICIAL INTELLIGENCE DD2380


#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iomanip>


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


	for (int i = 0; i < N_rows_A; i++)
	{
		for (int j = 0; j < N_cols_A; j++)
		{
			cout << ", " << A[i][j];
		}
		cout << ";" << endl;
	}
	//cout << ")";
}

int main()
{
	vector<int> aux12 = {1,5,6,23};
	


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

	//------------------------------------------------------------- Multiplication (beginning)

	vector<vector<double> > output_aux; //Result matrix
	vector<vector<double>> delta(N_cols_A, vector<double>(N_cols_oseq));
	vector<vector<int>> delta_idx(N_cols_A, vector<int>(N_cols_oseq));

	int N = N_cols_A; // number of states
	int T = N_cols_oseq;
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
			delta_idx[i][t] = distance(begin(aux), max);
			aux.clear();
		}
	}

	//visualization(delta);

	/*for (size_t i = 0; i < N; i++)
	{
		for (size_t t = 0; t < T; t++)
		{
			cout << delta_idx[i][t] << " ";
		}
		cout << endl;
	}
	system("pause");
	*/

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

	int last_state = index; //Last state
	vector<int> output(T); //Output vector
	output[T - 1] = last_state;

	//Explanation of the recursivity for the states.
	//IMPORTANT TO REMEMBER THAT T = 4 BUT OUR LAST STATE IS STATE_3!!

	/*int state_2 = delta_idx[last_state][2]; 

	int state_1 = delta_idx[state_2][1];

	int state_0 = delta_idx[state_1][0];
	*/
	
	for (size_t t = T - 1; t > 0; t--) //First we compute the states backwards recursively
	{
		index =	delta_idx[index][t];
		output[t - 1] = index;
	}
	
	//----------------------------------------------------------------- Multiplication (end)

	//-----------------------------------------------------------------Output (beginning)

	for (size_t t = 0; t < T; t++) //And then we print them onwards
	{
		cout << output[t] << " ";
	}

	//----------------------------------------------------------------------Output (end)


	system("pause");
	return 0;
}









