#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <vector>
#include <cmath>
#include <cstdlib>

using namespace std;

double massimo(double array[], int dim);

void CheckCoincidence(vector<vector<double>> &, int , Random &);

void WriteBestPath(vector<vector<vector<double>>> , int , ofstream &, const std::string& name, int n);

void CalculateL2(vector<vector<double>> &, int );

void mescola(vector<vector<double>> &, int, Random & );

void ordina(vector<vector<vector<double>>> &, int , int);

int seleziona(int, Random &);

void permCoppie(vector<vector<double>> &, int , Random &);

void permContig(vector<vector<double>> &, int, Random & );

void inversione(vector<vector<double>> &, int, Random & );

void shiftContig(vector<vector<double>> &, int, Random & );

void crossover(vector<vector<double>> &, vector<vector<double>> &, int, Random & );

vector<vector<vector<double>>> GeneraPopStart(vector<vector<double>> , int , int, Random & );

void NewPop(vector<vector<vector<double>>> &, int , int , double , double, Random & );
