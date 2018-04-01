/********************************************************
*							*
*	CSCI426 Assignment 2				*
*	Student Name: Kuan Wen Ng			*
*	Subject Code: CSCI464				*
*	Student Number: 5078052				*
*	Email ID: kwn961				*
*	Filename: tsp.cpp (ass2)			*
*	Description: Genetic Algorithm for TSP		*
*							*
********************************************************/

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <cstring>
#include <fstream>
#include <cmath>
#include <iterator>
#include <sstream>
#include <cstdlib>
#include <algorithm>
using namespace std;

struct City;
struct Tour;

//Alias
typedef vector<City> Cities;	//Vector of city
typedef pair<City, City> CityPair;	//City pair
typedef map<CityPair, double> Distances;	//Map for CityPair distances
typedef vector<Tour> Population;	//Population of tours
typedef vector<vector<double> > CityType;	//City type
typedef vector<double> NormFitness;	//Normalized fitness

//Prototypes
bool operator == (const City & , const City & );
bool operator != (const City & , const City & );
bool operator < (const CityPair & , const CityPair & );
bool operator == (const Tour & , const Tour & );
ostream &operator << (ostream & , const City & );
ostream &operator << (ostream & , const CityPair & );
ostream &operator << (ostream & , const Tour & );
string pairString(CityPair );
double randDouble();
int randInt(int );

//Struct for city
struct City
{
	int x;
	int y;
	int type;

	bool operator < (const City &right)
	{
		bool status = false;

		if (x == right.x)
			return y < right.y;

		else
			return x < right.x;
	}
};

//Tour of cities
struct Tour
{
	Cities cities;
	double cost;
	double fitness;
};

//City comparator
bool operator == (const City &left, const City &right)
{
	return (left.x == right.x) && (left.y == right.y);
}

//City comparator
bool operator != (const City &left, const City &right)
{
	return (left.x != right.x) || (left.y != right.y);
}

//CityPair comparator
bool operator < (const CityPair &left, const CityPair &right)
{
	return pairString(left) < pairString(right);
}

//Tour comparator
bool operator == (const Tour &left, const Tour &right)
{
	bool same = true;

	for (int index = 0; index < right.cities.size() && same; index++)
	{
		if (left.cities[index] != right.cities[index])
			same= false;
	}

	return same;
}

//City ostream
ostream &operator << (ostream &strm, const City &city)
{
	strm << city.x << " " << city.y << " ";

	return strm;
}

//CityPair ostream
ostream &operator << (ostream &strm, const CityPair &cityPair)
{
	strm << cityPair.first << " " << cityPair.second << " ";

	return strm;
}

//Tour ostream
ostream &operator << (ostream &strm, const Tour &tour)
{
	for (int index = 0; index < tour.cities.size(); index++)
		strm << tour.cities[index] << endl;

	strm << tour.cities[0];

	return strm;
}

//String format for CityPair
string pairString(CityPair cityPair)
{
	string output;
	stringstream stream;

	stream << cityPair.first.x << "," << cityPair.first.y << ":" << cityPair.second.x << "," << cityPair.second.y; 

	output = stream.str();
	return output;
}

//Random number 0 to 1
double randDouble()
{
	return(rand() / (double)(RAND_MAX));
}

//Random number 0 to num
int randInt(int num)
{
	return int(rand() / (double(RAND_MAX) + 1) * num);
}

//Class for GA
class GA
{
private:
	static CityType cityType;	//Cost factor between cities
	int popSize;	//Population size
	int tournamentSize;	//Tournament size
	int generation;	//Evolution cycle
	double crossoverRate;	//Crossover rate
	double mutationRate;	//Mutation rate
	bool selection;	//Selection type
	Cities cities;	//Cities
	Distances distances;	//Map for CityPair distances
	Population population;	//Current population
	NormFitness normFitness;	//Normalized fitness


public:
	GA();
	~GA();
	static void initCost();
	void readData(char *);
	double cityDistance(CityPair );
	double tourCost(Tour );
	void initDistances();
	void init();
	Tour crossover(Tour , Tour );
	Tour mutate(Tour );
	Tour tournament();
	Tour rouletteWheel();
	Tour fittest();
	void normalizeFitness();
	void evolve();
	void printCities();
};

CityType GA::cityType;

//GA constructor
GA::GA()
{
	popSize = 0;
	tournamentSize = 0;
	generation = 0;
	crossoverRate = 1.0;
	mutationRate = 0.05;
	selection = false;	//true = tournament, false = roulette
}

//GA destructor
GA::~GA()
{
}

//Initialize cost between city types
void GA::initCost()
{
	srand(time(0));
	vector<double> tempCost;
	tempCost.push_back(10);
	tempCost.push_back(7.5);
	tempCost.push_back(5);
	cityType.push_back(tempCost);
	tempCost.clear();
	tempCost.push_back(7.5);
	tempCost.push_back(5);
	tempCost.push_back(2.5);
	cityType.push_back(tempCost);
	tempCost.clear();
	tempCost.push_back(5);
	tempCost.push_back(2.5);
	tempCost.push_back(1);
	cityType.push_back(tempCost);
}

//Read file input
void GA::readData(char *filename)
{
	string tempString;
	ifstream inputFile;
	City tempCity;

	inputFile.open(filename);

	if (inputFile.good())
	{
		inputFile >> tempString;

		//Read data and store cities
		while (inputFile.good() && !inputFile.eof())
		{
			if (inputFile >> tempCity.x >> tempCity.y >> tempCity.type)
				cities.push_back(tempCity);
		}
		inputFile.close();
	}

	else
		cout << "Invalid filename." << endl;
}

//Distance between CityPair
double GA::cityDistance(CityPair cityPair)
{
	double cost = 0;

	//Euclidean distance
	cost = sqrt(pow(static_cast<double>(cityPair.first.x - cityPair.second.x), 2.0) + pow(static_cast<double>(cityPair.first.y - cityPair.second.y), 2.0));
	cost *= cityType[cityPair.first.type - 1][cityPair.second.type - 1] * 0.01;	//Cost between city type in cents

	return cost;
}

//Tour cost
double GA::tourCost(Tour tour)
{
	double cost = 0;
	CityPair cityPair;

	//Compute total costs between all cities in tour
	for (int index = 0; index < tour.cities.size(); index++)
	{
		if (index == tour.cities.size() - 1)
		{
			cityPair.first = tour.cities[index];
			cityPair.second = tour.cities[0];
			cost += distances[cityPair];
		}

		else
		{
			cityPair.first = tour.cities[index];
			cityPair.second = tour.cities[index + 1];
			cost += distances[cityPair];
		}
	}
	return cost;
}

//Initialize distances
void GA::initDistances()
{
	CityPair tempPair;

	//Compute non repeatable possible permutatations of CityPair distance 
	for (int index1 = 0; index1 < cities.size(); index1++)
	{
		for (int index2 = 0; index2 < cities.size(); index2++)
		{
			if (cities[index1] != cities[index2])
			{
				tempPair.first = cities[index1];
				tempPair.second = cities[index2];
				distances[tempPair] = cityDistance(tempPair);
			}
		}
	}
}

void GA::printCities()
{
	for (int index = 0; index < cities.size(); index++)
	{
		cout << cities[index] << endl;
	}
}

//Initialize population
void GA::init()
{
	Cities tempCities;
	Tour tempTour;
	int tempIndex;

	//Initialize dynamic parameters
	popSize = cities.size();
	generation = popSize;
	tournamentSize = 0.05 * popSize;

	//Generate random tours for population
	for (int popIndex = 0; popIndex < popSize; popIndex++)
	{
		tempCities = cities;
		for (int cityIndex = 0; cityIndex < cities.size(); cityIndex++)
		{
			tempIndex = randInt(tempCities.size());
			tempTour.cities.push_back(tempCities[tempIndex]);
			tempCities.erase(tempCities.begin() + tempIndex);
		}

		tempTour.cost = tourCost(tempTour);
		tempTour.fitness = 1 / tempTour.cost;
		population.push_back(tempTour);
		tempTour.cities.clear();
	}
}

//Two point crossover
Tour GA::crossover(Tour parent1, Tour parent2)
{
	bool clear = false;
	Tour tour;
	Cities tempTour, clearTour;
	int left = randInt(parent1.cities.size()), right = randInt(parent1.cities.size() - (left + 1)) + left;	//Select left and right for cutoff

	for (int index = left; index < right + 1; index++)
		tempTour.push_back(parent1.cities[index]);

	clearTour = tempTour;

	for (int index = 0; index < parent2.cities.size() && !clear; index++)
	{
		for (int clearIndex = 0; clearIndex < clearTour.size(); clearIndex++)
		{
			if (parent2.cities[index] == clearTour[clearIndex])
			{
				parent2.cities.erase(parent2.cities.begin() + index);
				clearTour.erase(clearTour.begin() + clearIndex);
				index--;
			}
		}
		if (clearTour.empty())
			clear = true;
	}

	for (int index = 0; index < left; index++)
	{
		tour.cities.push_back(parent2.cities[0]);
		parent2.cities.erase(parent2.cities.begin() + 0);
	}

	for (int index = 0; index < tempTour.size(); index++)
		tour.cities.push_back(tempTour[index]);

	for (int index = 0; index < parent2.cities.size(); index++)
		tour.cities.push_back(parent2.cities[index]);

	tour.cost = tourCost(tour);
	tour.fitness = 1 / tour.cost;

	return tour;
}

//Mutate with city swap
Tour GA::mutate(Tour tour)
{
	Tour mutatedTour = tour;
	City tempCity;

	int left = randInt(mutatedTour.cities.size()), right;	//Select two different cities to swap

	do
	{
		right = randInt(mutatedTour.cities.size());
	} while (left == right);

	tempCity = mutatedTour.cities[left];
	mutatedTour.cities[left] = mutatedTour.cities[right];
	mutatedTour.cities[right] = tempCity;

	mutatedTour.cost = tourCost(mutatedTour);
	mutatedTour.fitness = 1 / mutatedTour.cost;

	return mutatedTour;
}

//Tournament selection
Tour GA::tournament()
{
	double fitness = -1;
	int winnerIndex = -1, candidate = 0;;

	//Select fittest from tournament
	for (int index = 0; index < tournamentSize; index++)
	{
		candidate = randInt(popSize);

		if (population[candidate].fitness > fitness)
		{
			winnerIndex = candidate;
			fitness = population[candidate].fitness;
		}
	}
	return population[winnerIndex];
}

//Roulette selection
Tour GA::rouletteWheel()
{
	double random, sumFitness = 0;

	random = randDouble() * sumFitness;

	//Select fittest from roulette wheel
	for (int index = 0; index < popSize; index++)
	{
		sumFitness += normFitness[index];

		if (sumFitness > random)
			return population[index];
	}
	return population[popSize - 1];
}

//Normalize fitness
void GA::normalizeFitness()
{
	double maxFitness = population[0].fitness, minFitness = population[0].fitness, tempFitness;

	normFitness.clear();
	normFitness.resize(popSize, 0);

	//Get min and max fitness
	for (int index = 1; index < popSize; index++)
	{
		if (population[index].fitness > maxFitness)
			maxFitness = population[index].fitness;

		if (population[index].fitness < minFitness)
			minFitness = population[index].fitness;
	}

	tempFitness = maxFitness - minFitness;

	//Normalize fitness 0 to 1
	for (int index = 0; index < popSize; index++)
	{
		normFitness[index] = population[index].fitness - minFitness;
		normFitness[index] /= tempFitness;
	}
}

//Select fittest
Tour GA::fittest()
{
	int winnerIndex = 0;
	double fitness = population[0].fitness;

	for (int index = 1; index < popSize; index++)
	{
		if (population[index].fitness > fitness)
		{
			fitness = population[index].fitness;
			winnerIndex = index;
		}
	}
	return population[winnerIndex];
}

//Evolve population
void GA::evolve()
{
	Population newPopulation;
	Tour parent1, parent2, child;
	int round = 0, tempIndex = 0, interval = generation / 20;
	double tempFitness = fittest().fitness;

	//Evolve population N times
	for (int generationIndex = 0; generationIndex < generation; generationIndex++)
	{
		//Preserve fittest tour
		newPopulation.push_back(fittest());

		if (tempFitness != newPopulation[0].fitness)
		{
			tempFitness = newPopulation[0].fitness;
			tempIndex = generationIndex;
		}

		//Reduce crossover and mutation rate slightly if there are no changes
		if ((generationIndex - tempIndex) == (generation * 0.01))
		{
			crossoverRate *= 0.99;
			mutationRate *= 0.99;
		}

		if ((generationIndex - tempIndex) == (generation * 0.1))
		{
			//if ((generationIndex - 1) % interval != interval - 1)
				//cout << "Generation: " << generationIndex << " Cost: " << newPopulation[0].cost << endl;
			break;
		}

		//Normalize fitness for roulette wheel
		if (!selection)
			normalizeFitness();

		//Generate new tours from crossover and mutation for new population
		for (int populationIndex = 0; populationIndex < popSize - 1; populationIndex++)
		{
			round = 0;

			//Select parent from tournament
			if (selection)
			{
				parent1 = tournament();
				parent2 = tournament();
			}			

			//Select parent from roulette wheel
			else
			{
				parent1 = rouletteWheel();
				parent2 = rouletteWheel();
			}

			if (parent1 == parent2)
				parent2 = newPopulation[0];

			//Crossover
			if (randDouble() <= crossoverRate)
			{
				child = crossover(parent1, parent2);

				if (randDouble() <= mutationRate)
					child = mutate(child);	//Mutate

				newPopulation.push_back(child);
			}
			else
			{
				if (parent1.fitness > parent2.fitness)
					child = parent1;

				else
					child = parent2;

				if (randDouble() <= mutationRate)
					child = mutate(child);	//Mutate

				newPopulation.push_back(child);
			}
		}
		population = newPopulation;	//Store new population
		newPopulation.clear();

		//if (generationIndex % interval == interval - 1)
			cout << "Generation: " << generationIndex + 1 << " Cost: " << population[0].cost << endl;
	}
	cout << "Tour: \n" << fittest() << endl;
}

int main(int argc, char *argv[])
{
	char filename[30];
	ifstream inputFile;

	if (argc > 2)
	{
		cout << "Invalid number of arguments." << endl;
		return 1;
	}

	else if (argc == 2)
	{
		strcpy(filename, argv[1]);
	}

	else
	{
		cout << "Please enter filename: ";
		cin.getline(filename, 20);
	}

	inputFile.open(filename);

	if (!inputFile.good())
	{
		cout << "Invalid filename." << endl;
		inputFile.close();
		return 1;
	}

	inputFile.close();

	GA ga;
	ga.initCost();
	ga.readData(filename);
	ga.initDistances();
	ga.init();
	ga.evolve();

	return 0;
}

/************************************************************************************************************************************************
*																		*
*	The population size and number of generations is dynamic and will depend on the number of cities that will be travelled			*
*	and is N * 2.5 where N is number of cities. Parents that are selected for crossover will be different to promote diversity.		*
*	It should also be noted that to achieve an optimal result with the increase in population size where the number of cities is		*
*	large, the compute time is long and could take hours on tsp500.txt. This optimal settings is used to demonstrate what it could		*
*	on tsp100.txt where compute time is not too long. The complexity of the algorithm is O(5/2 N^2) which is huge. The generation		*
*	number is twice the population size but the algorithm will terminate if there are no changes in 10% of the max generation number.	*
*	The crossover and mutation rate will also decrease slowly when there are no changes to fitness every 1% of generation.			*
*	Parameters used are shown below:													*
*	Crossover Rate = 1.00															*
*	Mutation Rate = 0.05															*
*	Tournament Size = 5% of population size													*
*	Selection Type = Tournament														*
*	Generations = Population size * 2													*
*																		*
************************************************************************************************************************************************/
