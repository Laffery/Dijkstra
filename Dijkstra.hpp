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

	class Answer{
		private:
			int SDist = -1;
			int pathNum = 0;
			vector<stack<int>> Spath;
		
		public:
			Answer(){}
			Answer(int s, int p, vector<stack<int>> sp){
				this->SDist = s;
				this->pathNum = p;
				for(int i = 0; i < sp.size(); ++i)
					Spath.push_back(sp.at(i));
			}

			~Answer(){};

			void outputAnswer(fstream &outfile){
				outfile << SDist << endl << pathNum << endl;
				for(int i = 0; i < pathNum; ++i){
					outfile << Spath.at(i).top();
					Spath.at(i).pop();

					while(!Spath.at(i).empty()){
						outfile << "," << Spath.at(i).top();
						Spath.at(i).pop();
					}
					outfile << endl;
				}
			}
	};

	class Graph{
		private:
			int nodeCount = 0;
			int edgeCount = 0;
			vector<EdgeNode *> matrix;
			Answer answer[2];

		public:
			Graph(){}

			Graph(int n, int e):nodeCount(n), edgeCount(e){}

			~Graph(){}

			int getNodeCount(){return this->nodeCount;}

			void setNodeCount(int n){this->nodeCount = n;}

			int getEdgeCount(){return this->edgeCount;}

			void setEdgeCount(int e){this->edgeCount = e;}

			vector<EdgeNode *> getMatrix(){return this->matrix;}

			void setMatrix(vector<EdgeNode *> m){
				for(int i = 0; i < m.size(); ++i)
					this->matrix.push_back(m.at(i));
			}

			void showGraph(){
				char comma = ',';
				cout << nodeCount << comma << edgeCount << endl;

				for(int i = 0; i < matrix.size(); ++i){
					EdgeNode *curr = matrix.at(i);
					int beg = curr->getKey();

					while(curr->getNext() != nullptr){
						curr = curr->getNext();
						cout << beg << comma << curr->getKey() << comma << curr->getWeight() << endl;
					}
				}
			}

			void setAnswer(int index, Answer ans){this->answer[index] = ans;}

			void getAnswer(int index, fstream &outfile){this->answer[index].outputAnswer(outfile);}
	};


	vector<Graph> graphs;

public:
	
	DijkstraProject2(){}

	~DijkstraProject2(){}

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

	void run1MultiGraphHelper(int index, vector<EdgeNode *> graph, int nodeCount);

	/*
	 * It's a helper function od run1, to find minimum node whose flag is false 
	 */
	int run1Helper(bool flag[], int dist[], int nodeCount);

	/* 
	 * Find all shortest paths
	 */
	void findAllPath(vector<stack<int>> &paths, vector<int> *prev, int node, stack<int> path);

	/*
	 * Part 2, find the monotonically increasing path to finish Part 2
	 * and save the result to Param:`outputFile`.
	 * Save the path as: node_1,node_2...node_n. (seperate nodes with comma)
	 */
	void run2(const char* outputFile = "output.txt");

	void run2MultiGraphHelper(int index, vector<EdgeNode *> graph, int nodeCount);

	void ShowAnswer(const char* outputFile = "output.txt");

};

#endif // Dijkstra_h


/* 
 * Read data from input file
 */
void DijkstraProject2::readFromFile(const char* inputfile)
{
	cout << "readFromFile: " << inputfile << endl;

	//TODO
	ifstream infile(inputfile, ios::in);

	if(!infile){
		cout << "Can not open" << inputfile << endl;
		return;
	}
	
	/* Read */
	string tmp;
	char comma = ',';

	while (!infile.eof())
	{	
		int nodeCount = 0;
		int edgeCount = 0;
		
		getline(infile, tmp);
		stringstream ss(tmp);
		ss >> nodeCount >> comma >> edgeCount;
		// cout << nodeCount << comma << edgeCount <<endl;

		Graph graph(nodeCount, edgeCount);
		vector<EdgeNode *> matrix;
			
		for(int i = 0; i < nodeCount; ++i){
			EdgeNode *helper;
			helper = new EdgeNode(i, 0);
			matrix.push_back(helper);
		}

		/* Read edges */
		while(edgeCount){
			getline(infile, tmp);
			int beg, end, weight;
			stringstream ss(tmp);
			ss >> beg >> comma >> end >> comma >> weight;
			// cout << beg << comma << end << comma << weight <<endl;

			EdgeNode *helper;
			helper = new EdgeNode(end, weight);

			EdgeNode *curr = matrix.at(beg);
			while(curr->getNext() != nullptr)
			{
				curr = curr->getNext();
			}

			curr->setNext(helper);

			edgeCount--;
		}

		graph.setMatrix(matrix);
		graphs.push_back(graph);

		if(!infile.eof())
			getline(infile, tmp);
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

	for(int i = 0; i < graphs.size(); ++i){
		Graph graph = graphs.at(i);
		int nodeCount = graph.getNodeCount();
		vector<EdgeNode *> matrix = graph.getMatrix();

		run1MultiGraphHelper(i, matrix, nodeCount);
	}

	outfile.close();
}

void DijkstraProject2::run1MultiGraphHelper(int index, vector<EdgeNode *> graph, int nodeCount)
{
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
	int currNode = run1Helper(flag, dist, nodeCount), endNode = nodeCount - 1;
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


				/* If path through this node is much nearer, update distance */
				if(dist[tmpk] > tmpw + dist[currNode])
				{
					dist[tmpk] = tmpw + dist[currNode];
					pathNum -= (prev[tmpk].size() - 1);
					prev[tmpk].clear();
					prev[tmpk].push_back(currNode);
				}
				else if (dist[tmpk] == tmpw + dist[currNode]){
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
					if(!flag[i] && dist[i] < nearDist){
						nearDist = dist[i];
						nearNode = i;
					}
				}
			}

			/* out-degree is 0 */
			if(currNode != nearNode)
			{
				currNode = nearNode;
				flag[currNode] = true;
			}
		}

		currNode = run1Helper(flag, dist, nodeCount);
	}


	/* output result */
	vector<stack<int>> Spath;
	stack<int> path;
	findAllPath(Spath, prev, endNode, path);

	Answer ans(dist[endNode], pathNum, Spath);
	graphs.at(index).setAnswer(0, ans);
}

int DijkstraProject2::run1Helper(bool flag[], int dist[], int nodeCount)
{
	int min = INF, minNode = -1;
	for(int i = 0; i < nodeCount; ++i){
		if(!flag[i] && dist[i] < min){
			min = dist[i];
			minNode = i;
		}
	}

	return minNode;
}

void DijkstraProject2::findAllPath(vector<stack<int>> &paths, vector<int> *prev, int node, stack<int> path)
{
	path.push(node);

	if(prev[node].at(0) == -1){
		paths.push_back(path);
		return;
	}

	for(int i = 0; i < prev[node].size(); ++i){
		findAllPath(paths, prev, prev[node].at(i), path);
	}
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

	for(int i = 0; i < graphs.size(); ++i){
		Graph graph = graphs.at(i);
		int nodeCount = graph.getNodeCount();
		vector<EdgeNode *> matrix = graph.getMatrix();

		run2MultiGraphHelper(i, matrix, nodeCount);
	}

	outfile.close();
}

void DijkstraProject2::run2MultiGraphHelper(int index, vector<EdgeNode *> graph, int nodeCount)
{
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
	int currNode = run1Helper(flag, dist, nodeCount), endNode = nodeCount - 1;
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


				/* If path through this node is much nearer, update distance */
				if(dist[tmpk] > tmpw + dist[currNode])
				{
					dist[tmpk] = tmpw + dist[currNode];
					pathNum -= (prev[tmpk].size() - 1);
					prev[tmpk].clear();
					prev[tmpk].push_back(currNode);
				}
				else if (dist[tmpk] == tmpw + dist[currNode]){
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
					if(!flag[i] && dist[i] < nearDist){
						nearDist = dist[i];
						nearNode = i;
					}
				}
			}

			/* out-degree is 0 */
			if(currNode != nearNode)
			{
				currNode = nearNode;
				flag[currNode] = true;
			}
		}

		currNode = run1Helper(flag, dist, nodeCount);
	}


	/* output result */
	vector<stack<int>> Spath;
	stack<int> path;
	findAllPath(Spath, prev, endNode, path);

	Answer ans(dist[endNode], pathNum, Spath);
	graphs.at(index).setAnswer(1, ans);
}

void DijkstraProject2::ShowAnswer(const char* outputFile)
{

	fstream outfile(outputFile, ios::out | ios::app);
	if(!outfile){
		cout << "Can not open file " << outputFile << endl;
		outfile.close();
		return;
	}

	int gsize = graphs.size();

	for(int i = 0; i < gsize; ++i){
		graphs.at(i).getAnswer(0, outfile);
		outfile << endl;
		graphs.at(i).getAnswer(1, outfile);
		
		if(i < gsize - 1){
			outfile << "end" << endl << endl;
		}
		else{
			outfile << "End";
		}
	}

	outfile.close();
}