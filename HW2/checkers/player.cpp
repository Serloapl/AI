#include "player.hpp"
#include <cstdlib>
#include <algorithm>
#include <math.h>
#include <time.h>


#define A 1
#define B 2
#define d 8
#define dim 8
#define NUMPIECESTYPES 4
#define RANDOMRANGE 999999
#define HASHLENGTH  9999999
#define MODULE 9999999


namespace checkers
{
	struct HashState {
		GameState Board;
		double heuristic;
		int depth;
		int player;
	};



	std::vector<std::vector<HashState>> hash_vector(HASHLENGTH);

	bool BoardMatch(const GameState &Board1, const GameState &Board2)
	{
		for (size_t i = 0; i < dim; i++)
		{
			for (size_t j = 0; j < dim; j++)
			{
				if (Board1.at(i, j) != Board2.at(i, j))
				{
					return false;
				}
			}
		}
		return true;
	}

	int zobrist_hash(const GameState &lBoard)
	{
		//std::cerr << "LLEGA1" << std::endl;
		int h = 0;
		int q;
		int PrimeNumber = 31;

		for (size_t i = 0; i < dim; i++)
		{
			for (size_t j = 0; j < dim; j++)
			{
				h = (h*PrimeNumber + lBoard.at(i, j)) % MODULE;
			}

		}
		//std::cerr << "LLEGA2" << std::endl;
		//std::cerr << "--->" << h << std::endl;

		return h;
	}

	bool PlayerCanEat(const GameState &lBoard)
	{
		std::vector<GameState> lNextStates;
		lBoard.findPossibleMoves(lNextStates);
		int numPieces=0, numPieces_original=0;



		for (size_t i = 0; i < dim; i++)
		{
			for (size_t j = 0; j < dim; j++)
			{
				if (lBoard.at(i,j)!=0)
				{
					numPieces_original++;
				}
			}
		}

		for (size_t q = 0; q < lNextStates.size(); q++)
		{
			numPieces = 0;
			for (size_t i = 0; i < dim; i++)
			{
				for (size_t j = 0; j < dim; j++)
				{
					if (lNextStates[q].at(i, j) != 0)
					{
						numPieces++;
					}
				}
			}
			if (numPieces < numPieces_original)
			{
				return true;
			}
		}

		return false;

	}

	double heuristic(const GameState &lBoard)
	{
		
		double RedCont = 0;
		double WhiteCont = 0;
		double WhiteKingCont = 0;
		double RedKingCont = 0;

		std::vector<GameState> lNextStates;
		lBoard.findPossibleMoves(lNextStates);

		if (lBoard.isRedWin())
		{
			return 10000000;
		}
		if (lBoard.isWhiteWin())
		{
			return -100000000;
		}

		
		for (size_t i = 0; i < dim; i++)
		{
			for (size_t j = 0; j < dim; j++)
			{
				if (lBoard.at(i, j)&CELL_RED)
				{
					if (lBoard.at(i, j)&CELL_KING)
					{
						RedKingCont++;
					}
					else RedCont++;
				}
				if (lBoard.at(i, j)&CELL_WHITE)
				{
					if (lBoard.at(i, j)&CELL_KING)
					{
						WhiteKingCont++;
					}
					else WhiteCont++;
				}
			}
		}

		int CanEat = 0;
		if (lBoard.getNextPlayer() == A && PlayerCanEat(lBoard))
		{
			//std::cerr << "LLEGA3" << std::endl;
			CanEat = 100;
		}
		if (lBoard.getNextPlayer() == B && PlayerCanEat(lBoard))
		{
			//std::cerr << "LLEGA4" << std::endl;
			CanEat = -100;
		}

		
		return ((RedCont + 1.2 * RedKingCont) - (WhiteCont + 1.2 * WhiteKingCont) + CanEat);
	}

	int minimax(const GameState &Board, int player, int depth)
	{
		std::vector<GameState> AllNextStates;
		Board.findPossibleMoves(AllNextStates);
		int bestPossible, indexbestPossible;
		int v;
		int aux;


		//termination
		if (depth == d || AllNextStates.size() == 0) // 0-> GAME OVER
		{
			return heuristic(Board);
		}


		if (player == A) // player is 1/A
		{
			bestPossible = -100000000;
			for (size_t i = 0; i < AllNextStates.size(); i++)
			{
				v = minimax(AllNextStates[i], B, depth + 1);
				if (v > bestPossible)
				{
					bestPossible = v;
					indexbestPossible = i;
				}
			}
			return bestPossible;

		}
		else // player is 0/B
		{
			bestPossible = 100000000;
			for (size_t i = 0; i < AllNextStates.size(); i++)
			{
				v = minimax(AllNextStates[i], A, depth + 1);
				//std::cerr << AllNextStates[i].toString(A) << std::endl;
				//std::cerr << "Heuristica = " << v << std::endl;
				if (v < bestPossible)
				{
					bestPossible = v;
				}
			}
			return bestPossible;
		}

	}

	std::pair<double,int> alphabeta(const GameState &Board, int player, int depth, double alpha, double beta)
	{

		std::vector<GameState> AllNextStates;
		Board.findPossibleMoves(AllNextStates);
		double bestPossible;
		std::pair<double, int> v;
		int hash_visited,index, indexbestPossible;
		HashState aux;
		int vreal;

		//termination
		if (depth >= d || AllNextStates.size() == 0) // 0-> GAME OVER
		{
			return std::make_pair(heuristic(Board),-1);
		}

		if (player == A) // player is 1/A
		{
			//std::cerr << "LLEGA3" << std::endl;
			bestPossible = -100000000;
			for (size_t i = 0; i < AllNextStates.size(); i++)
			{
				v = alphabeta(AllNextStates[i], B, depth + 1, alpha, beta);

				//hash_visited = zobrist_hash(AllNextStates[i]);
				//if (hash_vector[hash_visited].size() != 0) // Boards at some hash index
				//{
				//	index = -1;
				//	for (size_t j = 0; j < hash_vector[hash_visited].size(); j++)
				//	{
				//		if ( BoardMatch(hash_vector[hash_visited][j].Board, AllNextStates[i]) && hash_vector[hash_visited][j].depth<=depth && hash_vector[hash_visited][j].player==A)
				//		{
				//			index = j;
				//			v.first = hash_vector[hash_visited][index].heuristic;
				//			//vreal = alphabeta(AllNextStates[i], B, depth + 1, alpha, beta);
				//			//if (v != vreal) contador_fallo++; else contador_correcto++;
				//			break;
				//		}
				//	}
				//	if (index == -1) // no Board is the one studied
				//	{
				//		v = alphabeta(AllNextStates[i], B, depth + 1, alpha, beta);
				//		aux.Board = AllNextStates[i]; aux.heuristic = v.first; aux.depth = depth; aux.player = A;
				//		hash_vector[hash_visited].push_back(aux);
				//	}
				//}
				//else // no Boards at that hash index
				//{
				//	v = alphabeta(AllNextStates[i], B, depth + 1, alpha, beta);
				//	aux.Board = AllNextStates[i]; aux.heuristic = v.first; aux.depth = depth; aux.player = A;
				//	hash_vector[hash_visited].push_back(aux);
				//}

				if (v.first > bestPossible)
				{
					bestPossible = v.first;
					indexbestPossible = i;
				}
				if (bestPossible > alpha)
				{
					alpha = bestPossible;
				}
				if (beta <= alpha) break;
			}
			return std::make_pair(bestPossible,indexbestPossible);

		}
		else // player is 0/B
		{
			bestPossible = 100000000;
			for (size_t i = 0; i < AllNextStates.size(); i++)
			{
				v = alphabeta(AllNextStates[i], A, depth + 1, alpha, beta);

				//hash_visited = zobrist_hash(AllNextStates[i]);
				////std::cerr << "-------> " << hash_visited;
				//if (hash_vector[hash_visited].size() != 0) // Boards at some hash index
				//{
				//	index = -1;
				//	for (size_t j = 0; j < hash_vector[hash_visited].size(); j++)
				//	{
				//		if (BoardMatch(hash_vector[hash_visited][j].Board, AllNextStates[i]) && hash_vector[hash_visited][j].depth<=depth && hash_vector[hash_visited][j].player == B)
				//		{
				//			index = j;
				//			v.first = hash_vector[hash_visited][index].heuristic;
				//			//vreal = alphabeta(AllNextStates[i], A, depth + 1, alpha, beta);
				//			//if (v != vreal) contador_fallo++; else contador_correcto++;
				//			break;
				//		}
				//	}
				//	if (index == -1) // no Board is the one studied
				//	{
				//		v = alphabeta(AllNextStates[i], A, depth + 1, alpha, beta);
				//		aux.Board = AllNextStates[i]; aux.heuristic = v.first; aux.depth = depth; aux.player = B;
				//		hash_vector[hash_visited].push_back(aux);
				//	}
				//}
				//else // no Boards at that hash index
				//{
				//	v = alphabeta(AllNextStates[i], A, depth + 1, alpha, beta);
				//	aux.Board = AllNextStates[i]; aux.heuristic = v.first; aux.depth = depth; aux.player = B;
				//	hash_vector[hash_visited].push_back(aux);
				//}


				if (v.first < bestPossible)
				{
					bestPossible = v.first;
					indexbestPossible = i;
				}
				if (bestPossible < beta)
				{
					beta = bestPossible;
				}
				if (beta <= alpha) break;
			}
			return std::make_pair(bestPossible, indexbestPossible);
		}
	}

	GameState Player::play(const GameState &pState, const Deadline &pDue)
	{
		std::vector<GameState> lNextStates;
		pState.findPossibleMoves(lNextStates);
		double bestPossible, indexbestPossible, alpha,beta;
		//init_zobrist();


		if (lNextStates.size() == 0) return GameState(pState, Move());

		alpha = -INFINITY;
		beta = INFINITY;
		std::pair<double,int> results;	


		results = alphabeta(pState, pState.getNextPlayer(), 0, alpha, beta);
		return lNextStates[results.second];
	}

	/*namespace checkers*/
}