#include <stdlib.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <boost/sort/spreadsort/spreadsort.hpp>
#include <boost/bind.hpp>
#include <limits>
#include <memory>
#include <fstream>
#include <sstream>
#include <string>
#include <iterator>

/* Transform the graph from an edgelist into a form readable by the mpi implementation

    -> instead of only one pair per edge, have a pair in each direction
    -> sort by first index, to simplify building the csr representation
    -> store back to a file, because the sorting takes quite some time and memory.

    Some assumptions about the final graph to work with the mpi:
    - graph is sorted
    - contains both directions of each edge
    - start with index 0

*/

int main(int argc, char* argv[]) {

  if (argc < 2) {
    std::cout << "specify input file..";
  }
  std::cout << "sorting file: " << argv[1] << std::endl;

  std::string filename = argv[1];
  std::ifstream infile(argv[1]);

  uint32_t nedges = 0;
  uint32_t nvertices = 0;
  uint32_t minimum = 10; //maybe not the best init here..
  std::string nstring, estring, h;
  std::string line;
  std::vector<std::pair<uint32_t,uint32_t>> edgepairs;
  while(std::getline(infile, line)) {
    std::istringstream iss(line);
    uint32_t a, b;
    if (iss >> a >> b) {
      if (a != b) {
        minimum = (a < minimum) ? a : minimum; //check for base index
        minimum = (b < minimum) ? b : minimum;
        edgepairs.push_back(std::make_pair(a,b));
        edgepairs.push_back(std::make_pair(b,a));
        nedges += 2;
      }
    } else { // get the number of vertices in the graph
      std::istringstream iss2(line);
      if (iss2 >> h >> nstring >> a >> estring >> b) {
        if (nstring == "Nodes:") {
          nvertices = a;
        }
      }
    }
  }

  std::sort(edgepairs.begin(), edgepairs.end(), [](auto &left, auto &right) {
    return left.first < right.first || (!(right.first < left.first) && left.second < right.second);
  });

  std::cout << "size before: " << edgepairs.size() << std::endl;

  edgepairs.erase(std::unique(edgepairs.begin(), edgepairs.end()), edgepairs.end());

  std::cout << "size after: " << edgepairs.size() << std::endl;


  if (minimum != 0) {
    std::cout << "base is: " << minimum << std::endl;
    for (auto &i: edgepairs) {
      i.first -= minimum;
      i.second -= minimum;
    }
  }

  filename.append("sorted.txt");
  std::cout << "saving to: " << filename << std::endl;
  std::ofstream outfile(filename);

  nvertices = edgepairs.back().first + 1;

  //encode number of vertices and number of edges in the first line
  outfile << nvertices << " " << edgepairs.size() << std::endl;

  for(auto &x:edgepairs) {
    outfile << x.first << " " << x.second << std::endl;
  }

}
