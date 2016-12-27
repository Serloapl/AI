#include "Player.hpp"
#include <cstdlib>
#include <iostream>
#include <vector>
#include "Bird.hpp"
#include "Markov.h"



namespace ducks
{

	Player::Player()
	{
		NumberSpecies = 6;
		ModelsGlobal.resize(NumberSpecies);
		ModelSpecies.resize(NumberSpecies);
	}

	Markov Player::FuseModels(std::vector<Markov> BirdSameSpecies, int Ronda)
	{

		int NumBird = BirdSameSpecies.size();

		//int N = BirdSameSpecies[0].gamma.size();
		int N = BirdSameSpecies[0].NumFlightPattern;

		//int M = BirdSameSpecies[0].gamma_di[0][0].size();
		int M = BirdSameSpecies[0].NumMovementDirections;

		//int T = BirdSameSpecies[0].gamma_di.size();

		Markov Acum;
		//std::cerr << "Estoy dentro de Fuse" << std::endl;

		/*for (size_t k = 0; k < NumBird; k++)
		{
		BirdSameSpecies[k].VisualizeModel();
		}*/

		double numer;
		double denom;


		Acum.A.resize(N);
		for (size_t i = 0; i < N; i++)
		{
			Acum.A[i].resize(N);
			for (size_t j = 0; j < N; j++)
			{
				numer = denom = 0;
				for (size_t k = 0; k < NumBird; k++)
				{
					for (size_t t = 0; t < BirdSameSpecies[k].gamma_di.size(); t++)
					{
						numer += BirdSameSpecies[k].gamma_di[t][i][j];
						denom += BirdSameSpecies[k].gamma[i][t];
					}
				}
				Acum.A[i][j] = numer / denom;
			}
		}



		Acum.B.resize(N);
		for (size_t i = 0; i < N; i++)
		{
			Acum.B[i].resize(M);
			for (size_t j = 0; j < M; j++)
			{
				numer = denom = 0;
				for (size_t k = 0; k < NumBird; k++)
				{
					for (size_t t = 0; t < BirdSameSpecies[k].gamma_di.size(); t++)
					{
						if (BirdSameSpecies[k].oseq[t] == j)
						{
							numer += BirdSameSpecies[k].gamma[i][t];
						}
						denom += BirdSameSpecies[k].gamma[i][t];
					}
				}
				Acum.B[i][j] = numer / denom;
			}
		}



		double total;
		Acum.pi.resize(1);
		Acum.pi[0].resize(N);
		for (size_t i = 0; i < N; i++)
		{
			total = 0;
			for (size_t k = 0; k < NumBird; k++)
			{
				total += BirdSameSpecies[k].gamma[i][0];
			}
			Acum.pi[0][i] = total / NumBird;
		}


		return Acum;
	}

	Action Player::shoot(const GameState &pState, const Deadline &pDue)
	{


		int Ronda = pState.getRound();
		int Birdnumber = pState.getNumBirds();
		std::vector<Bird> Birds(Birdnumber);
		std::vector<std::vector<int>> BirdObservations;
		BirdObservations.resize(Birdnumber);


		//calculo el numero de pajaros aun vivos
		int BirdsAlive = 0;
		for (size_t i = 0; i < Birdnumber; i++)
		{
			//Birds.push_back(pState.getBird(i));
			Birds[i] = pState.getBird(i);
			if (Birds[i].isAlive())
			{
				BirdsAlive++;
			}
		}

		this->Models.resize(Birdnumber);
		double MaxAbsolut = -99999999999;
		double MaxInterno;
		int ind_dir = 1;
		int ind_bird = 1;
		int birdSpecies;
		int indexMovement;
		int Direction;
		double ProbDirection;
		double Prob;
		double Old_Prob;
		vector<int> seqAux;



		//////////////////////////////////////////////////////////
		//--------------------TRAINING--------------------------//
		//////////////////////////////////////////////////////////

		//entreno al final de la ronda o todo el rato si solo queda un pajaro vivo
		if ((BirdsAlive == 1) || (Birds[0].getSeqLength() > 60))
		{

			this->Models.clear();
			for (size_t i = 0; i < Birdnumber; i++)
			{
				this->Models[i].InitialiseModel();
				for (size_t j = 0; j < Birds[i].getSeqLength(); j++)
				{
					if (Birds[i].getObservation(j) == -1)
					{
						break;
					}
					BirdObservations[i].push_back(static_cast<int>(Birds[i].getObservation(j)));

				}

				Models[i].TrainModel(BirdObservations[i]);
				Models[i].oseq = BirdObservations[i];
			}
		}
		//////////////////////////////////////////////////////////
		//--------------------SHOOTING--------------------------//
		//////////////////////////////////////////////////////////

		//no disparo si es ronda 0
		if (Ronda > 5 && ModelsGlobal[5].size() == 0)
		{
			return cDontShoot;
		}
		//si no es ronda 0, disparo a partir de 60
		if (Birds[0].getSeqLength() > 60)
		{
			indexMovement = -1;
			ProbDirection = 0;
			for (size_t i = 0; i < Birdnumber; i++)
			{
				if (Birds[i].isDead())
				{
					continue;
				}

				//std:cerr << "---------------->  " << "HAY " << ModelsGlobal[5].size() << "Blackstorks y la prob es "<< ModelSpecies[5].CalcProbSeq(BirdObservations[i]) << std::endl;
				if (ModelsGlobal[5].size() != 0)
				{
					if (ModelSpecies[5].CalcProbSeq(BirdObservations[i]) > -300)
					{
						continue;
					}
				}
				

				Models[i].PredictNextMov(BirdObservations[i], &indexMovement, &ProbDirection);
				MaxInterno = ProbDirection;
				Direction = indexMovement;

				if (MaxInterno > MaxAbsolut)
				{
					MaxAbsolut = MaxInterno;
					ind_dir = Direction;
					ind_bird = i;
				}
			}

			//std::cerr << " DISPARA en la direccion " << Direction << " la prob es " << MaxInterno << "con secuencia de longitud " << seqAux.size() << std::endl;

			if (MaxAbsolut > 0.65)
			{
				std::cerr << "DIRECCION DE DISPARO: " << ind_dir << " con prob " << MaxAbsolut << " a pajaro " << ind_bird << std::endl;
			   	return Action(ind_bird, (EMovement)ind_dir);
			}
		}


		return cDontShoot;
	}

	std::vector<ESpecies> Player::guess(const GameState &pState, const Deadline &pDue)
	{
		std::vector<ESpecies> lGuesses(pState.getNumBirds(), SPECIES_UNKNOWN);



		int Round = pState.getRound();
		int Birdnumber = pState.getNumBirds();
		int NDirections = Models[0].NumMovementDirections;
		int NStates = Models[0].NumFlightPattern;
		std::vector<Bird> Birds(Birdnumber);
		double Prob;
		double Old_Prob;


		if (Round == 0)
		{
			for (size_t i = 0; i < Birdnumber; i++)
			{
				lGuesses[i] = SPECIES_RAVEN;
			}
		}
		else
		{
			for (size_t i = 0; i < Birdnumber; i++)
			{
				Old_Prob = -99999999999999;
				Prob = -99999999999998;
				for (size_t k = 0; k < NumberSpecies; k++)
				{
					if (ModelsGlobal[k].size() != 0)
					{
						//std::cerr << "ModelsGlobal no esta vacio" << std::endl;
						Prob = ModelSpecies[k].CalcProbSeq(Models[i].oseq);
						//std::cerr << "prob de seq" << Prob << std::endl;
						if (Prob > Old_Prob)
						{
							lGuesses[i] = (ESpecies)k;
							Old_Prob = Prob;
						}
					}
				}
			}

		}


		return lGuesses;

	}

	void Player::hit(const GameState &pState, int pBird, const Deadline &pDue)
	{
		std::cerr << "HIT BIRD!!!" << std::endl;
	}

	void Player::reveal(const GameState &pState, const std::vector<ESpecies> &pSpecies, const Deadline &pDue)
	{
		int Birdnumber = pState.getNumBirds();
		int especie;

		std::cerr << "Solucion ronda " << pState.getRound() << " es ";
		for (size_t i = 0; i < Birdnumber; i++)
		{
			especie = pSpecies[i];
			this->ModelsGlobal[especie].push_back(this->Models[i]);

			std::cerr << pSpecies[i];
		}
		std::cerr << std::endl;


		for (size_t k = 0; k < NumberSpecies; k++)
		{
			if (ModelsGlobal[k].size() != 0)
			{
				ModelSpecies[k] = FuseModels(ModelsGlobal[k], pState.getRound());
			}
		}



		/*
		* If you made any guesses, you will find out the true species of those birds in this function.
		*/



	}


} /*namespace ducks*/
