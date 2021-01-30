#include <stdlib.h>
#include <time.h>
#include <random>
#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <limits>
#include <memory>
#include <fstream>
#include <sstream>
#include <string>
//compilation: g++ -O3 -std=c++11 *.cpp


uint32_t loadGraph(uint32_t& nvertices, uint32_t*& edges, uint32_t*& rowpointers, std::string filename) {
  std::ifstream infile(filename);

  int currentvertex = -1;
  uint32_t ii = 0; //current row pointer
  uint32_t i = 0; //current edge
  std::string line;
  uint32_t a, b, nedges;

  std::getline(infile, line);
  std::istringstream iss(line);
  if (iss >> a >> b) {
    nvertices = a;
    nedges = b;
  }
  std::cout << nvertices << " " << nedges << std::endl;

  edges = new uint32_t[nedges];
  rowpointers = new uint32_t[nvertices];

  while(std::getline(infile, line)) {
    std::istringstream iss(line);
    if (iss >> a >> b) {
      if (a != currentvertex) { //this comparison only works because unsigned value of -1 is really large and we know we start at 0.
        int missed = a - currentvertex; //this also only works because we know we start at 0, you should actually never mix signed and unsigned!
        for (int j = 0; j < missed; j++) {
          rowpointers[ii] = i;
          ii += 1;
        }
        currentvertex = a;
      }
      edges[i] = b;
      i += 1;
    }
  }
  return nedges;
}

/* Sequential graph coloring based on adjacency list traversals.

    This code
    * Loads a graph from a text file and stores it in two arrays: one of continuous
      vertex numbers, and one array that points to where the next vertex starts.
      From this start index to the next start index, all the vertices it shares
      an edge with is listed.

    * Then uses the first fit approach to color the graph. in particular, nodes are
      colored sequentially by stepping through the vertices by their index.
      Colors are stored in a separate data structure, which is in this case a vector of size n.
      The element at index n is the current color of the node n.

    * This code uses unsigned 32bit integers for the biggest part to handle huge
      graphs. Maximum number of edges allowed is 4'294'967'295.

 */


int main(int argc, char* argv[]) {

  if (argc < 2) {
    std::cout << "specify input file..";
  }
  std::cout << "processing file: " << argv[1] << std::endl;

  std::string filename = argv[1];

  uint32_t nglob, nedges;
  uint32_t* rowpointers;
  uint32_t* edges;
  uint32_t* nvertices;

  nedges = loadGraph(nglob,edges,rowpointers,filename);


  int maxcolors = 500; //what's a good heuristic here?


  std::vector<int> colors(nglob,0); //vector of size n, colors[i] is the color of node i.
  std::vector<uint32_t> forbidden(maxcolors,0); //vector of size maxcolors, initialized to -1.

  auto wcts = std::chrono::system_clock::now();

  uint32_t start, end;

  for (uint32_t i = 0; i < nglob; i++) { // for all nodes in adjacency matrix

    start = rowpointers[i];
    end = (i != nglob-1) ? rowpointers[i+1] : nedges;
    for (uint32_t j = start; j < end; j++) {
      uint32_t neighbor = edges[j];
      int nodecolor = colors[neighbor];
      forbidden[nodecolor] = i;
    }

    //first fit assignment
    for (int k = 0; k < maxcolors; k++) { //iterate over forbidden array
      if (forbidden[k] != i) { //assign the lowest color that wasn't forbidden in previous step
        colors[i] = k;
        break;
      }
    }
  }
  //parallel version needs conflict resolution code at this point.

  //stop timer
  std::chrono::duration<double> wctduration = (std::chrono::system_clock::now() - wcts);
  std::cout << "Finished in " << wctduration.count() << " seconds [Wall Clock]" << std::endl;

  //check correctness (trivial in sequential, but maybe usable for parallel)
  for (uint32_t i = 0; i < nglob; i++) {
    int owncolor = colors[i];
    start = rowpointers[i];
    end = (i != nglob-1) ? rowpointers[i+1] : nedges;
    for (uint32_t j = start; j < end; j++) {
      uint32_t edge = edges[j];
      if (owncolor == colors[edge]) {
        std::cout << "mistake at vertex " << i << std::endl;
      }
    }
  }

  //print number of colors
  auto max = *std::max_element(colors.begin(),colors.end());
  std::cout << max << " colors used." << std::endl;

  std::ofstream results("results_sequential.csv", std::ios_base::app | std::ios_base::out);
  results << filename << "," << 1 << "," << nglob << "," << nedges << "," << wctduration.count() << "," << 1 << "," << max << "," << 1 << "," << 1 << "\n";

  //print colors
  /*
  for (int i = 0; i < n; i++) {
    std::cout << colors[i];
  }
  */

  return 0;
}
