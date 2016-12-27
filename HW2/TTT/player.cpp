#include "player.hpp"
#include <cstdlib>
#include <algorithm>
#include <math.h>

#define A 1
#define B 0
#define d 1




namespace TICTACTOE
{

int heuristic(GameState lBoard)
{
	int dims = 4;
	int contX1=0,contX2=0, contX3=0, contX4=0;
	int contO1=0, contO2=0, contO3=0, contO4=0;
	int total_cont=0;

	std::vector<std::vector<int>> TableHeuristic(dims + 1, std::vector<int>(dims + 1));
	TableHeuristic[0] = { 0    ,   -1, -10, -100, -1000 };
	TableHeuristic[1] = { 1    ,    0,   0,	   0,	  0 };
	TableHeuristic[2] = { 10   ,    0,   0,    0,	  0 };
	TableHeuristic[3] = { 100  ,    0,   0,    0,     0 };
	TableHeuristic[4] = { 1000 ,    0,   0,    0,	  0 };
	
		


	//computing sum of all X symbols on O free rows
	for (size_t i = 0; i < dims; i++)
	{
		contX1 = 0;
		contO1 = 0;
		contX2 = 0;
		contO2 = 0;
		for (size_t j = 0; j < dims; j++)
		{
			if (lBoard.at(i,j)&CELL_X)
			{
				contX1 += 1;
			}
			if (lBoard.at(i, j)&CELL_O)
			{
				contO1 += 1;
			}
			if (lBoard.at(j, i)&CELL_X)
			{
				contX2 += 1;
			}
			if (lBoard.at(j, i)&CELL_O)
			{
				contO2 += 1;
			}
		}
		total_cont += TableHeuristic[contX1][contO1] + TableHeuristic[contX2][contO2];
		if (lBoard.at(i, i)&CELL_X)
		{
			contX3 += 1;
		}
		if (lBoard.at(i, i)&CELL_O)
		{
			contO3 += 1;
		}
		if (lBoard.at(i, dims - 1 - i)&CELL_X)
		{
			contX4 += 1;
		}
		if (lBoard.at(i, dims - 1 - i)&CELL_O)
		{
			contO4 += 1;
		}
	}
	total_cont += TableHeuristic[contX3][contO3] + TableHeuristic[contX4][contO4];
	return total_cont;
}

int minimax(GameState Board, int player, int depth)
{
	std::vector<GameState> AllNextStates;
	Board.findPossibleMoves(AllNextStates);
	int bestPossible, indexbestPossible;
	int v;	
	int aux;


	//termination
	if ( depth==d || AllNextStates.size() == 0) // 0-> GAME OVER
	{
		return heuristic(Board);
	}


	if (player == A) // player is 1/A
	{
		bestPossible = -100000000;
		for (size_t i = 0; i < AllNextStates.size(); i++)
		{
			v = minimax(AllNextStates[i], B, depth + 1);
			if (v>bestPossible)
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
			if (v < bestPossible)
			{
				bestPossible = v;
			}
		}
		return bestPossible;
	}
	
}

int alphabeta(GameState Board, int player, int depth, int alpha, int beta)
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
			v = alphabeta(AllNextStates[i], B, depth + 1, alpha, beta);
			if (v > bestPossible)
			{
				bestPossible = v;
				//indexbestPossible = i;
			}
			if (bestPossible > alpha)
			{
				alpha = bestPossible;
			}
			if (beta <= alpha) break;
		}
		return bestPossible;

	}
	else // player is 0/B
	{
		bestPossible = 100000000;
		for (size_t i = 0; i < AllNextStates.size(); i++)
		{
			v = alphabeta(AllNextStates[i], A, depth + 1, alpha, beta);
			//std::cerr << AllNextStates[i].toString(A) << std::endl;
			//std::cerr << "Heuristica = " << v << std::endl;
			if (v < bestPossible)
			{
				bestPossible = v;
			}
			if (bestPossible < beta)
			{
				beta = bestPossible;
			}
			if (beta <= alpha) break;
		}
		return bestPossible;
	}




}


GameState Player::play(const GameState &pState,const Deadline &pDue)
{
	std::vector<GameState> lNextStates;
	pState.findPossibleMoves(lNextStates);
	int bestPossible, indexbestPossible, v;
	int alpha, beta;

	if (lNextStates.size() == 0) return GameState(pState, Move());

	bestPossible = -100000000;
	alpha = -400000000;
	beta = 400000000;

	for (size_t i = 0; i < lNextStates.size(); i++)
	{
		//v = minimax(lNextStates[i], B, 1);
		v = alphabeta(lNextStates[i], B, 1, alpha, beta);
		//std::cerr << lNextStates[i].toString(B) << std::endl;
		//std::cerr << "Heuristica = " << v << std::endl;
		if (v > bestPossible)
		{
			bestPossible = v;
			indexbestPossible = i;
		}
		if (bestPossible > alpha)
		{
			alpha = bestPossible;
		}
		if (beta <= alpha) break;
	}

	/*std::cerr << "El índice es:" << indexbestPossible << std::endl;
	std::cerr << lNextStates[indexbestPossible].toString(1) << std::endl;

	std::cerr << "el número de segundos es = " << pDue.getSeconds() << std::endl;*/

	return lNextStates[indexbestPossible];



    /*
     * Here you should write your clever algorithms to get the best next move, ie the best
     * next state. This skeleton returns a random move instead.
     */


    //return lNextStates[rand() % lNextStates.size()];
}

/*namespace TICTACTOE*/ }
