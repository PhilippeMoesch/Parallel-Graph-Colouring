// Copyright (C) 2005, 2006 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Douglas Gregor
//           Andrew Lumsdaine
//
//  compilation: mpic++ -O3 -std=c++14 -ldl -lboost_mpi -lboost_serialization -lboost_graph_parallel -lboost_system
#define PBGL_ACCOUNTING

#include <boost/graph/use_mpi.hpp>
#include <boost/config.hpp>
#include <boost/throw_exception.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/graph/distributed/boman_et_al_graph_coloring.hpp>
#include <boost/graph/distributed/mpi_process_group.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/graph/parallel/distribution.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/graph/distributed/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <iostream>
#include <boost/random.hpp>
#include <boost/test/minimal.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/accounting.hpp>
#include <fstream>
#include <time.h>
#include <chrono>

#ifdef BOOST_NO_EXCEPTIONS
void
boost::throw_exception(std::exception const& ex)
{
    std::cout << ex.what() << std::endl;
    abort();
}
#endif

using namespace boost;
using boost::graph::distributed::mpi_process_group;

void
test_distributed_graph_coloring(std::string filename, int superstep)
{
  typedef adjacency_list<vecS,
                         distributedS<mpi_process_group, vecS>,
                         undirectedS> Graph;

  typedef property_map<Graph, vertex_index_t>::type vertex_index_map;

  std::ifstream infile(filename);
  if (!infile) {
    std::cout << "Error loading file" << std::endl;
  }

  std::string line;
  uint32_t nvertices;
  uint32_t nedges;
  uint32_t a;
  uint32_t b;

  std::getline(infile, line);
  std::istringstream iss(line);
  if (iss >> a >> b) {
    nvertices = a;
    nedges = b;
  }

  Graph g(nvertices);

  if (process_id(g.process_group()) == 0) {
    while (std::getline(infile,line)) {
    std::istringstream iss(line);
      if (iss >> a >> b) {
        add_edge(vertex(a,g),vertex(b,g),g);
      }
    }
  }

  // Set up color map
  std::vector<int> colors_vec(num_vertices(g));
  iterator_property_map<int*, vertex_index_map> color(&colors_vec[0],
                                                      get(vertex_index, g));

  auto wtcs = std::chrono::system_clock::now();
  // Run the graph coloring algorithm
  graph::boman_et_al_graph_coloring(g, color, superstep);

  std::chrono::duration<double> wtcduration = std::chrono::system_clock::now() - wtcs;
  if (process_id(g.process_group()) == 0) {
    std::cout << "Finished in " << wtcduration.count() << " seconds [Wall time]" << std::endl;
  }
  if (process_id(g.process_group()) == 0) {
    graph::distributed::boman_et_al_graph_coloring_stats.print(std::cout);
  }

  graph::distributed::boman_et_al_graph_coloring_stats_t stats = graph::distributed::boman_et_al_graph_coloring_stats;

  if (process_id(g.process_group()) == 0) {
    std::ofstream results("results_baseline.csv", std::ios_base::app | std::ios_base::out);
    results << filename <<  "," << num_processes(g.process_group()) << "," << nvertices << "," << nedges << "," << /*boost::graph::accounting::print_time(stats.execution_time)*/ stats.execution_time << "," << -1 << "," << stats.num_colors << "," << superstep << "," << 1 << std::endl; 
  }
}

int test_main(int argc, char* argv[])
{
  std::string filename;
  int superstep;

  mpi::environment env(argc, argv);
  mpi::communicator comm;

  if (argc < 3) {
    if (comm.rank() == 0) {
      std::cout << "invalid number of arguments";
    }
    boost::mpi::environment::abort(-1);
  }

  if (comm.rank() == 0) {
    std::cout << "Processing file: " << argv[1] << std::endl;
  }

  filename = argv[1];

  std::istringstream bla(argv[2]);
  if (!(bla >> superstep)) {
    if (comm.rank() == 0) {
      std::cout << "invalid superstep" << std::endl;
    }
    boost::mpi::environment::abort(-1);
  }

  test_distributed_graph_coloring(filename, superstep);

  return 0;
}
