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

#define INF INT16_MAX

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
			int Sdist = -1;
			int pathNum = 0;
			vector<stack<int>> Spath;
		
		public:
			Answer(){}
			Answer(int s, int p, vector<stack<int>> sp){
				this->Sdist = s;
				this->pathNum = p;
				for(int i = 0; i < sp.size(); ++i)
					Spath.push_back(sp.at(i));
			}

			~Answer(){};

			int getSdist(){return this->Sdist;}

			void setSdist(int s){this->Sdist = s;}

			int getPathNum(){return this->pathNum;}

			void setPathNum(int p){this->pathNum = p;}

			void addPathNum(int p){this->pathNum += p;}

			vector<stack<int>> getSpath(){return this->Spath;}

			void setSpath(vector<stack<int>> p){
				if(!Spath.empty())
					Spath.clear();

				for(int i = 0; i < p.size(); ++i)
					this->Spath.push_back(p.at(i));
			}

			void addSpath(vector<stack<int>> p){
				for(int i = 0; i < p.size(); ++i)
					this->Spath.push_back(p.at(i));
			}

			void outputAnswer(fstream &outfile){
				outfile << Sdist << endl << pathNum << endl;
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

			int getAnswerSdist(int index){return this->answer[index].getSdist();}

			int getAnswerPathNum(int index){return this->answer[index].getPathNum();}

			vector<stack<int>> getAnswerSpath(int index){return this->answer[index].getSpath();}

			void showAnswer(int index, fstream &outfile){this->answer[index].outputAnswer(outfile);}
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
	int run1NearHelper(bool flag[], int dist[], int nodeCount);

	/* 
	 * Find all shortest paths
	 */
	void run1findAllPath(vector<stack<int>> &paths, vector<int> *prev, int node, stack<int> path);

	/*
	 * Part 2, find the monotonically increasing path to finish Part 2
	 * and save the result to Param:`outputFile`.
	 * Save the path as: node_1,node_2...node_n. (seperate nodes with comma)
	 */
	void run2(const char* outputFile = "output.txt");

	void run2MultiGraphHelper(int index, vector<EdgeNode *> graph, int nodeCount, bool asc);

	int run2NearHelper(vector<EdgeNode *> matrix, bool flag[], int dist[], int nodeCount, vector<int> *prev, bool asc);

	int run2WeightHelper(vector<EdgeNode *> matrix, int beg, int end);

	bool run2AscHelper(vector<EdgeNode *> matrix, vector<int> *prev, int beg, int end, bool asc);

	void run2findAllPath(vector<EdgeNode *> matrix, vector<stack<int>> &paths, vector<int> *prev, int node, stack<int> path, bool asc);

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
	for(int i = 0; i < graphs.size(); ++i){
		Graph graph = graphs.at(i);
		int nodeCount = graph.getNodeCount();
		vector<EdgeNode *> matrix = graph.getMatrix();

		run1MultiGraphHelper(i, matrix, nodeCount);
	}
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
	int currNode = run1NearHelper(flag, dist, nodeCount), endNode = nodeCount - 1;
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
			else
				break;
		}

		currNode = run1NearHelper(flag, dist, nodeCount);
	}


	/* store result */
	if(dist[endNode] == INF)
		return;
	
	vector<stack<int>> Spath;
	stack<int> path;
	run1findAllPath(Spath, prev, endNode, path);

	Answer ans(dist[endNode], pathNum, Spath);
	graphs.at(index).setAnswer(0, ans);
	delete []prev;
}

int DijkstraProject2::run1NearHelper(bool flag[], int dist[], int nodeCount)
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

void DijkstraProject2::run1findAllPath(vector<stack<int>> &paths, vector<int> *prev, int node, stack<int> path)
{
	path.push(node);

	if(prev[node].at(0) == -1){
		paths.push_back(path);
		return;
	}

	for(int i = 0; i < prev[node].size(); ++i){
		run1findAllPath(paths, prev, prev[node].at(i), path);
	}
}

void DijkstraProject2::run2(const char* outputFile)
{
	cout << "Save result to file:" << outputFile << endl;

	//TODO
	for(int i = 0; i < graphs.size(); ++i){
		Graph graph = graphs.at(i);
		int nodeCount = graph.getNodeCount();
		vector<EdgeNode *> matrix = graph.getMatrix();

		/* Asc */
		run2MultiGraphHelper(i, matrix, nodeCount, true);

		/* Desc */
		run2MultiGraphHelper(i, matrix, nodeCount, false);
	}
}

void DijkstraProject2::run2MultiGraphHelper(int index, vector<EdgeNode *> graph, int nodeCount, bool asc)
{
	/* Parsing graph */
	// cout << index << asc << endl;
	int pathNum = 1;
	int dist[nodeCount];
	vector<int> *prev;
	prev = new vector<int>[nodeCount];
	bool flag[nodeCount];
	for(int i = 0; i < nodeCount; ++i){
		dist[i] = INF;
		prev[i].push_back(-1);
		flag[i] = false;
	}
	dist[0] = 0;


	/* Core step */
	int currNode = 0, endNode = nodeCount - 1;
	while(currNode != -1){
		flag[currNode] = true;
		int nearNode = currNode;
		int nearDist = INF;
		int tmpk, tmpw;

		while(currNode != endNode){
			/* Traverse all next nodes of current node */
			EdgeNode *curr = graph.at(currNode);
			while(curr->getNext() != nullptr){
				curr = curr->getNext();
				tmpk = curr->getKey();
				tmpw = curr->getWeight();

				/* If path through this node is much nearer, update distance */
				bool ascHelper = run2AscHelper(graph, prev, currNode, tmpk, asc);

				// cout << "(" << currNode << "," << dist[currNode] 
				// 		<< "),(" << tmpk << "," << dist[tmpk] 
				// 		<< ")," << tmpw << ":" << ascHelper << endl;
				

				if(dist[tmpk] > tmpw + dist[currNode] && ascHelper)
				{
					dist[tmpk] = tmpw + dist[currNode];
					pathNum -= (prev[tmpk].size() - 1);
					prev[tmpk].clear();
					prev[tmpk].push_back(currNode);
				}
				else if (dist[tmpk] == tmpw + dist[currNode] && ascHelper){
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
			}

			/* Find nearest node */
			for(int i = 0; i < nodeCount; ++i){
				if(!flag[i] && dist[i] < nearDist){
					nearDist = dist[i];
					nearNode = i;
				}
			}
			
			// cout << currNode << "," << nearNode << ")\n";
			if(currNode != nearNode)
			{
				currNode = nearNode;
				flag[currNode] = true;
			}
			else
				break;
		}

		// cout << "curr (" << currNode << "," << dist[currNode] 
		//  << "),(" << nearNode << "," << dist[nearNode] << ")\n";

		currNode = run2NearHelper(graph, flag, dist, nodeCount, prev, asc);

		// cout << "next (" << currNode << "," << dist[currNode] 
		//  << "),(" << nearNode << "," << dist[nearNode] << ")\n";
	}


	/* store result */
	if(dist[endNode] == INF)
		return;
	
	vector<stack<int>> Spath;
	stack<int> path;
	run2findAllPath(graph, Spath, prev, endNode, path, asc);

	int currSdist = graphs.at(index).getAnswerSdist(1);
	int currPathNum = graphs.at(index).getAnswerPathNum(1);
	vector<stack<int>> currSpath = graphs.at(index).getAnswerSpath(1);
	Answer ans(currSdist, currPathNum, currSpath);

	if(dist[endNode] < currSdist || (dist[endNode] > currSdist && currSdist == -1)){
		ans.setSdist(dist[endNode]);
		ans.setPathNum(pathNum);
		ans.setSpath(Spath);
	}
	else if(dist[endNode] == currSdist)
	{
		ans.addPathNum(pathNum);
		ans.addSpath(Spath);
	}
	graphs.at(index).setAnswer(1, ans);
	delete []prev;
}

int DijkstraProject2::run2WeightHelper(vector<EdgeNode *> matrix, int beg, int end)
{
	EdgeNode *curr = matrix.at(beg);
	while (curr->getKey() != end){
		curr = curr->getNext();
	}

	return curr->getWeight();
} 

bool DijkstraProject2::run2AscHelper(vector<EdgeNode *> matrix, vector<int> *prev, int beg, int end, bool asc)
{
	if(prev[beg].at(0) == -1)
		return true;
	
	// cout << currNode << " " << prev[currNode].size() << endl;
	int weight = run2WeightHelper(matrix, beg, end);
	for(int i = 0; i < prev[beg].size(); ++i){
		int preNode = prev[beg].at(i);
		int prew = run2WeightHelper(matrix, preNode, beg);
		// cout << currNode << ' ' << tmpk << ' ' << prew << endl;
		if(asc && weight > prew)
			return true;

		if(!asc && weight < prew)
			return true;
	}

	return false;
}

int DijkstraProject2::run2NearHelper(vector<EdgeNode *> matrix, bool flag[], int dist[], int nodeCount, 
									vector<int> *prev, bool asc)
{
	int min = INF, minNode = -1;
	for(int i = 0; i < nodeCount; ++i){
		if(!flag[i] && dist[i] < min){
			bool ascHelper = false;
			for(int j = 0; j < matrix.size(); ++j){
				EdgeNode *curr = matrix.at(j);
				while(curr->getNext() != nullptr){
					curr = curr->getNext();
					if(curr->getKey() == i){
						ascHelper = run2AscHelper(matrix, prev, j, i, asc);
						break;
					}
				}

				if(ascHelper)
					break;
			}

			if(ascHelper){
				min = dist[i];
				minNode = i;
			}
		}
	}

	return minNode;
}

void DijkstraProject2::run2findAllPath(vector<EdgeNode *> matrix, vector<stack<int>> &paths, vector<int> *prev,
									 int node, stack<int> path, bool asc)
{
	path.push(node);

	if(prev[node].at(0) == -1){
		paths.push_back(path);
		return;
	}

	for(int i = 0; i < prev[node].size(); ++i){
		int beg = prev[node].at(i);

		if(run2AscHelper(matrix, prev, beg, node, asc))
			run2findAllPath(matrix, paths, prev, beg, path, asc);
	}
}

void DijkstraProject2::ShowAnswer(const char* outputFile)
{

	fstream outfile(outputFile, ios::out);
	if(!outfile){
		cout << "Can not open file " << outputFile << endl;
		outfile.close();
		return;
	}

	int gsize = graphs.size();

	for(int i = 0; i < gsize; ++i){
		graphs.at(i).showAnswer(0, outfile);
		outfile << endl;
		graphs.at(i).showAnswer(1, outfile);
		
		if(i < gsize - 1){
			outfile << "end" << endl << endl;
		}
		else{
			outfile << "End";
		}
	}

	outfile.close();
}