#pragma once

#include <cstdint>
#include <string>


/**
 * This file declare the main class of Project2:DijstraProject2.
 * You should implement the methods:`readFromFile`,`run1`and`run2` in Dijstra.cpp.
 * You can add anything in DijstraProject2 class to implement Project2.
 */
class DijstraProject2 {
private:
	//You can declare your graph structure here.


public:

	/**
	 * Read graph from Param:`inputfile`.
	 * 
	 */
	void readFromFile(const char* inputfile="input.txt");
	
	/**
	 * Part 1, implement Dijstra algorithm to finish Part 1
	 * and save the result to Param:`outputFile`.
	 * Save the path as: node_1,node_2...node_n. (seperate nodes with comma)
	 *
	 */
	void run1(const char* outputFile = "output.txt");

	/**
	 * Part 2, find the monotonically increasing path to finish Part 2
	 * and save the result to Param:`outputFile`.
	 * Save the path as: node_1,node_2...node_n. (seperate nodes with comma)
	 *
	 */
	void run2(const char* outputFile = "output.txt");

};