README: CoreExact and Exact Algorithms
Project Overview
This project implements two algorithms from the paper "Efficient Algorithms for Densest Subgraph Discovery":

- CoreExact Algorithm: Finds the maximum core number of an undirected graph.
- Exact Algorithm: Finds the densest subgraph based on average degree or pattern-based densities (e.g., 4-cliques).

Both are implemented in C++ and work efficiently on large graphs.
Requirements
C++ compiler (supporting at least C++11).
Files Included
- coreExact.cpp : Implementation of CoreExact algorithm.
- exact.cpp : Implementation of Exact densest subgraph algorithm.
Compilation Procedures
While implementing and compiling the code, one must write the names of dataset parameters such as edge, triangle and 3clique etc. Any other form while compiling will show error as it won't match the pattern, as we wanted to maintain the actual names from the assignment paper.
CoreExact
g++ -std=c++17 -O2 -o coreExact coreExact.cpp
Exact
g++ -std=c++11 -O2 -o exact exact.cpp
Running the Programs
CoreExact
./coreExact <graph_file> <h>

Example:
./coreExact as19991210.txt 100
Exact
./exact <graph_file> [pattern]

Example:
./exact as19991210.txt
./exact output_graph.txt 4clique
Input Format
Both programs expect the input file to contain edges:

- Each line: two integers representing an undirected edge.
- Comment lines start with # and are ignored.

Example Input:
1 2
2 3
3 4
4 1

Indexing Notes:
- CoreExact expects 1-indexed files.
- Exact expects 0-indexed files for pattern-based searches.
Output Format
CoreExact
Outputs the maximum core number found in the graph.

Example Output:
Maximum Core Number: 36
Exact
Outputs maximum average degree or maximum pattern density.
Lists nodes included in the densest subgraph found.

Example Output (Average Degree):
Maximum Average Degree: 7.24
Densest Subgraph Size: 120 nodes
Nodes: 1, 5, 10, 12, ...

Example Output (4-Clique Pattern):
Maximum 4-Clique Density: 242.25
Subgraph Size: 20 nodes
Nodes: 645, 1429, 1447, ...
Notes
- Ensure correct indexing depending on the algorithm.
- Larger datasets may take longer but are handled efficiently.
- Pattern searches (like 4-cliques) may require more computational time.
References
Yixiang Fang, Kaiqiang Yu, Reynold Cheng, Laks V.S. Lakshmanan, Xuemin Lin,
"Efficient Algorithms for Densest Subgraph Discovery", PVLDB 2019.
