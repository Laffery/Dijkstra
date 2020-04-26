#pragma once

#ifndef DIJKSTRA_H
#define DIJKSTRA_H

#include <iostream>
#include <fstream>
#include <cstdint>
#include <string>
#include <sstream>

using namespace std;

/**
 * This file declare the main class of Project2:DijkstraProject2.
 * You should implement the methods:`readFromFile`,`run1`and`run2` in Dijkstra.cpp.
 * You can add anything in DijkstraProject2 class to implement Project2.
 */
class DijkstraProject2 {
private:
	//You can declare your graph structure here.
	int nodeCount = 0;
	int edgeCount = 0;
	int **matrix; // neighbor matrix

	class Edge{
		private:
			int beg;
			int end;
			int weight;

		public:
			Edge(int b, int e, int w):beg(b), end(e), weight(w){}

			~Edge(){}

			int getBegin(){return this->beg;}

			int getEnd(){return this->end;}

			int getWeight(){return this->weight;}

			void display(){cout << beg << "," << end << "," << weight << endl;}
	};

	Edge *edgs;

public:
	
	DijkstraProject2(){}

	~DijkstraProject2(){
		freeMatrix();
		freeEdges();
	}

	/*
	 * Read graph from Param:`inputfile`.
	 */
	void readFromFile(const char* inputfile = "input.txt");
	
	/*
	 * Part 1, implement Dijstra algorithm to finish Part 1
	 * and save the result to Param:`outputFile`.
	 * Save the path as: node_1,node_2...node_n. (seperate nodes with comma)
	 */
	void run1(const char* outputFile = "output.txt");

	/*
	 * Part 2, find the monotonically increasing path to finish Part 2
	 * and save the result to Param:`outputFile`.
	 * Save the path as: node_1,node_2...node_n. (seperate nodes with comma)
	 */
	void run2(const char* outputFile = "output.txt");

	/* 
	 * To free matrix of neighbor
	 */
	void freeMatrix();

	/*
	 * To free array of all edges 
 	 */
	void freeEdges();

};

#endif // Dijkstra_h


/* 
 * Read data from input file
 */
void DijkstraProject2::readFromFile(const char* inputfile)
{
	cout << "readFromFile: " << inputfile << endl;

	//TODO
	int lines = 0;
	ifstream infile(inputfile, ios::in);

	if(!infile){
		cout << "Can not open" << inputfile << endl;
		return;
	}

	string tmp;
	char comma = ',';
	while (getline(infile, tmp)){
		if(lines == 0){
			stringstream ss(tmp);
			ss >> nodeCount >> comma >> edgeCount;
			
			matrix = new int*[nodeCount];
			for(int i = 0; i < nodeCount; ++i){
				matrix[i] = new int[nodeCount];

				for(int j = 0; j < nodeCount; ++j)
					matrix[i][j] = 0;
			}

			++lines;
		}

		else{
			int beg, end, weight;
			stringstream ss(tmp);
			ss >> beg >> comma >> end >> comma >> weight;

			matrix[beg][end] = weight;

			cout << beg << comma << end << comma << weight << endl;
		}
	}

	for(int i = 0; i < nodeCount; i++){
		for(int j = 0; j < nodeCount; j++)
			cout << matrix[i][j] << ' ';

		cout << endl;
	}
	infile.close();
}

void DijkstraProject2::run1(const char* outputFile)
{
	cout << "Save result to file:" << outputFile << endl;
	
	//TODO
}

void DijkstraProject2::run2(const char* outputFile)
{
	cout << "Save result to file:" << outputFile << endl;

	//TODO
}

/* 
 * To free matrix of neighbor
 */
void DijkstraProject2::freeMatrix()
{
	for(int i = 0; i < nodeCount; ++i){
		delete []matrix[i];
	}

	delete []matrix;
}

/*
 * To free array of all edges 
 */
void DijkstraProject2::freeEdges()
{
	delete []edgs;
}