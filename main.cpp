#include <iostream>
#include "Dijkstra.hpp"

/*
 * You can use this file to test your code.
 */
int main()
{
    DijkstraProject2 pro;
    pro.readFromFile("input.txt");
    pro.run1("output.txt");
    pro.run2("output.txt");
    pro.ShowAnswer("output.txt");
    return 0;
}