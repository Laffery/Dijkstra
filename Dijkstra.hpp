#pragma once

#ifndef DIJKSTRA_H
#define DIJKSTRA_H

#include <iostream>
#include <fstream>
#include <cstdint>
#include <string>
#include <sstream>
#include <vector>
#include <stack>

#define INF INT_MAX

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

	class EdgeNode{
		private:
			int key;
			int weight = 0;
			EdgeNode *next = nullptr;

		public:
			EdgeNode(int k, int w):key(k), weight(w){}

			~EdgeNode(){}

			int getKey(){return this->key;}

			int getWeight(){return this->weight;}

			EdgeNode *getNext(){return this->next;}

			void setNext(EdgeNode *n){this->next = n;}
	};

	vector<EdgeNode *> graph;

public:
	
	DijkstraProject2(){}

	~DijkstraProject2(){
		freeGraph();
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
	 * It's a helper function od run1, to find minimum node whose flag is false 
	 */
	int run1Helper(bool flag[], int dist[]){
		int min = INF, minNode = -1;
		for(int i = 0; i < nodeCount; ++i){
			if(!flag[i] && dist[i] < min){
				min = dist[i];
				minNode = i;
			}
		}

		return minNode;
	}

	/* 
	 * Find all shortest paths
	 */
	void findAllPath(vector<stack<int>> &paths, vector<int> *prev, int node, stack<int> path){
		path.push(node);

		if(prev[node].at(0) == -1){
			paths.push_back(path);
			return;
		}

		for(int i = 0; i < prev[node].size(); ++i){
			findAllPath(paths, prev, prev[node].at(i), path);
		}
	}

	/*
	 * Part 2, find the monotonically increasing path to finish Part 2
	 * and save the result to Param:`outputFile`.
	 * Save the path as: node_1,node_2...node_n. (seperate nodes with comma)
	 */
	void run2(const char* outputFile = "output.txt");

	/* 
	 * To free the graph
	 */
	void freeGraph();

};

#endif // Dijkstra_h


/* 
 * Read data from input file
 */
void DijkstraProject2::readFromFile(const char* inputfile)
{
	cout << "readFromFile: " << inputfile << endl;

	//TODO
	bool firstLine = true;
	ifstream infile(inputfile, ios::in);

	if(!infile){
		cout << "Can not open" << inputfile << endl;
		return;
	}

	string tmp;
	char comma = ',';
	while (getline(infile, tmp)){
		if(firstLine){
			stringstream ss(tmp);
			ss >> nodeCount >> comma >> edgeCount;
			
			for(int i = 0; i < nodeCount; ++i){
				EdgeNode *helper;
				helper = new EdgeNode(i, 0);
				graph.push_back(helper);
			}

			firstLine = false;
		}

		else{
			int beg, end, weight;
			stringstream ss(tmp);
			ss >> beg >> comma >> end >> comma >> weight;

			EdgeNode *helper;
			helper = new EdgeNode(end, weight);

			EdgeNode *curr = graph.at(beg);
			while(curr->getNext() != nullptr)
			{
				curr = curr->getNext();
			}

			curr->setNext(helper);
		}
	}
	infile.close();
}

void DijkstraProject2::run1(const char* outputFile)
{
	cout << "Save result to file:" << outputFile << endl;
	
	//TODO
	fstream outfile(outputFile, ios::out | ios::app);
	if(!outfile){
		cout << "Can not open file " << outputFile << endl;
		outfile.close();
		return;
	}


	/* first node has no out degree, so no path to n */
	if(graph.at(0)->getNext() == nullptr){
		outfile << "Infinite" << endl << 0 << endl <<endl;
		outfile.close();
		return;
	}


	/* Parsing graph */
	int pathNum = 1;
	int dist[nodeCount];
	vector<int> *prev;// path[nodeCount];
	prev = new vector<int>[nodeCount];
	bool flag[nodeCount];
	for(int i = 0; i < nodeCount; ++i){
		dist[i] = INF;
		prev[i].push_back(-1);
		flag[i] = false;
	}
	dist[0] = 0;


	/* Core step */
	int currNode = run1Helper(flag, dist), endNode = nodeCount - 1;
	while(currNode != -1){
		flag[currNode] = true;
		
		while(currNode != endNode){
			EdgeNode *curr = graph.at(currNode);
			int nearNode = currNode;
			int nearDist = INF;
			int tmpk, tmpw;

			/* Traverse all next nodes of current node */
			while(curr->getNext() != nullptr){
				curr = curr->getNext();
				tmpk = curr->getKey();
				tmpw = curr->getWeight();

				// cout << currNode << "(" << dist[currNode] << "): " << tmpk << endl;

				/* If path through this node is much nearer, update distance */
				if(dist[tmpk] > tmpw + dist[currNode])
				{
					dist[tmpk] = tmpw + dist[currNode];
					// cout << tmpk << " :>" << dist[tmpk] << endl;
					pathNum -= (prev[tmpk].size() - 1);
					prev[tmpk].clear();
					prev[tmpk].push_back(currNode);
				}
				else if (dist[tmpk] == tmpw + dist[currNode]){
					// cout << tmpk << " :=\n";
					bool helper = false;
					for(int i = 0; i < prev[tmpk].size(); ++i){
						if(prev[tmpk].at(i) == currNode){
							helper = true;
							break;
						}
					}

					if(!helper){
						prev[tmpk].push_back(currNode);
						pathNum++;
					}
				}

				/* Find nearest node */
				for(int i = 0; i < nodeCount; ++i){
					// cout << i <<"(" << flag[i] << ") ";
					// if(dist[i] == INF)
					// 	cout << "INF,";
					// else
					// 	cout << dist[i] << ",";
					

					if(!flag[i] && dist[i] < nearDist){
						nearDist = dist[i];
						nearNode = i;
					}
				}
				// cout << endl << currNode << "(" << dist[currNode] << "), " 
				// 	<< nearNode << "(" << dist[nearNode] << ")\n";
			}

			/* out-degree is 0 */
			if(currNode != nearNode)
			{
				currNode = nearNode;
				flag[currNode] = true;
			}
		}

		currNode = run1Helper(flag, dist);
	}


	/* output result */
	outfile << dist[endNode] << endl << pathNum << endl;
	vector<stack<int>> Spath;
	stack<int> path;
	findAllPath(Spath, prev, endNode, path);

	for(int i = 0; i < pathNum; ++i){
		outfile << Spath.at(i).top();
		Spath.at(i).pop();

		while(!Spath.at(i).empty()){
			outfile << "," << Spath.at(i).top();
			Spath.at(i).pop();
		}
		outfile << endl;
	}

	outfile << endl;

	outfile.close();
}

void DijkstraProject2::run2(const char* outputFile)
{
	cout << "Save result to file:" << outputFile << endl;

	//TODO
	fstream outfile(outputFile, ios::out | ios::app);
	if(!outfile){
		cout << "Can not open file " << outputFile << endl;
		outfile.close();
		return;
	}

	int pathNum = 0;

	/* first node has no out degree, so no path to n */
	if(graph.at(0)->getNext() == nullptr){
		outfile << "Infinite" << endl << pathNum << endl << "end" << endl <<endl;
		outfile.close();
		return;
	}

	outfile.close();
}

void DijkstraProject2::freeGraph()
{
	for(int i = 0; i < nodeCount; ++i){
		EdgeNode *curr = graph.at(i);
		EdgeNode *helper;

		while(curr->getNext() != nullptr){
			helper = curr;
			curr = curr->getNext();
			delete helper;	
		}
	}

	graph.clear();
}