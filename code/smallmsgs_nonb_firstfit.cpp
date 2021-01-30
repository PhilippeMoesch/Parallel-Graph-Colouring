#include <stdlib.h>
#include <time.h>
#include <random>
#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <list>
#include <limits>
#include <memory>
#include <mpi.h>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>

// For the memory usage logging
#include <stdio.h>
#include <string.h>
#include <memory.h>

#define ERROR -1

//compilation: mpic++ -O3 -std=c++14 *.cpp

/*
 * Source: https://hpcf.umbc.edu/general-productivity/checking-memory-usage/
 * Gets the VmRSS and VmSize from the /proc/self/status
 * Input:
 *      long* vmrss_kb: the VmRSS size taken from /proc/self/status
 *      long* vmsize_kb: the VmSize size taken from /proc/self/status
 * Output:
 *      int: 0 for successful termination or something else if an error occurred
 */
int get_memory_usage_kb(long* vmrss_kb, long* vmsize_kb) {
    // Get the the current process' status file from the proc filesystem
    FILE* procfile = fopen("/proc/self/status", "r");

    long to_read = 8192;
    char buffer[to_read];
    int read = fread(buffer, sizeof (char), to_read, procfile);
    fclose(procfile);

    short found_vmrss = 0;
    short found_vmsize = 0;
    char* search_result;

    // Look through proc status contents line by line
    char delims[] = "\n";
    char* line = strtok(buffer, delims);

    while (line != NULL && (found_vmrss == 0 || found_vmsize == 0)) {
        search_result = strstr(line, "VmRSS:");
        if (search_result != NULL) {
            sscanf(line, "%*s %ld", vmrss_kb);
            found_vmrss = 1;
        }

        search_result = strstr(line, "VmSize:");
        if (search_result != NULL) {
            sscanf(line, "%*s %ld", vmsize_kb);
            found_vmsize = 1;
        }

        line = strtok(NULL, delims);
    }

    return (found_vmrss == 1 && found_vmsize == 1) ? 0 : 1;
}

/*
 * Source: https://hpcf.umbc.edu/general-productivity/checking-memory-usage/
 * Gets the local cluster memory usage.
 * Input:
 *      long* vmrss_per_process: the array with the VmRSS size for all clusters, taken from /proc/self/status
 *      long* vmsize_per_process: the array with the VmSize size for all clusters, taken from /proc/self/status
 *      int np: number of processors
 * Output:
 *      int: 0 for successful termination or something else if an error occurred
 */
int get_cluster_memory_usage_kb(long* vmrss_per_process, long* vmsize_per_process, int root, int np) {
    long vmrss_kb;
    long vmsize_kb;
    int ret_code = get_memory_usage_kb(&vmrss_kb, &vmsize_kb);

    if (ret_code != 0) {
        printf("Could not gather memory usage!\n");
        return ret_code;
    }

    MPI_Gather(&vmrss_kb, 1, MPI_UNSIGNED_LONG,
            vmrss_per_process, 1, MPI_UNSIGNED_LONG,
            root, MPI_COMM_WORLD);

    MPI_Gather(&vmsize_kb, 1, MPI_UNSIGNED_LONG,
            vmsize_per_process, 1, MPI_UNSIGNED_LONG,
            root, MPI_COMM_WORLD);

    return 0;
}

/*
 * Source: https://hpcf.umbc.edu/general-productivity/checking-memory-usage/
 * Gets the global memory usage.
 * Input:
 *      long* global_vmrss: the global VmRSS size, taken from /proc/self/status
 *      long* global_vmsize: the global VmSize size, taken from /proc/self/status
 *      int np: number of processors
 * Output:
 *      int: 0 for successful termination or something else if an error occurred
 */
int get_global_memory_usage_kb(long* global_vmrss, long* global_vmsize, int np) {
    long vmrss_per_process[np];
    long vmsize_per_process[np];
    int ret_code = get_cluster_memory_usage_kb(vmrss_per_process, vmsize_per_process, 0, np);

    if (ret_code != 0) {
        return ret_code;
    }

    *global_vmrss = 0;
    *global_vmsize = 0;
    for (int i = 0; i < np; i++) {
        *global_vmrss += vmrss_per_process[i];
        *global_vmsize += vmsize_per_process[i];
    }

    return 0;
}

/*
 * Loads a graph from the file (needs to be preprocessed by sort_graph.cpp)
 * Input:
 *      uint32_t& nvertices: the number of vertices the graph has
 *      uint32_t*& edges: the number of edges the graph has
 *      uint32_t*& rowpointers: the position in the edge array that the next vertex starts
 *      std::string filename: the name of the file containing the graph
 * Output:
 *      uint32_t: the number of edges the graph has or ERROR if there was an error while opening the file
 */
uint32_t loadGraph(uint32_t& nvertices, uint32_t*& edges, uint32_t*& rowpointers, std::string filename) {

    // Opening the file and checking for an error
    std::ifstream infile(filename);
    if (!infile) {
        std::cout << "An error occurred while opening the graph file" << std::endl;
        return ERROR;
    }

    int currentvertex = -1;
    uint32_t crow = 0; // Current row pointer
    uint32_t cedge = 0; // Current edge
    std::string line;
    uint32_t a, b, nedges;

    // Reading the first line that contains the number of edges and vertices of the graph
    std::getline(infile, line);
    std::istringstream iss(line);
    if (iss >> a >> b) {
        nvertices = a;
        nedges = b;
    }
    std::cout << nvertices << " " << nedges << std::endl;

    // Creating the graph arrays (vertices and edges)
    edges = new uint32_t[nedges];
    rowpointers = new uint32_t[nvertices];

    // Reading all edges from the file and saving them to the arrays
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        if (iss >> a >> b) {
            if (a != currentvertex) {
                // if we have gaps in the vertex indices, set the rowpointers for the vertices before too
                int missed = a - currentvertex; //this only works because we know the first vertex id is 0.
                for (int j = 0; j < missed; j++) {
                    rowpointers[crow] = cedge;
                    crow += 1;
                }
                currentvertex = a;
            }
            edges[cedge] = b;
            cedge += 1;
        }
    }

    infile.clear();
    return nedges;
}

/* MPI Graph coloring

  Smallmessages: Only the changes are sent after each superstep.

  This happens as follows:
   - in the first iteration, after each superstep all the colors of the superstep are broadcast everyone.
     This is a contiguous block of memory and works nicely.
   - in subsequent iterations, two arrays are broadcast after each superstep:
      - one array with the indices of the vertices that where recolored
      - another array with the actual colors
      - due to the unknown size, the first element of the index array contains the number vertices that where recolored.
         Thus, color i belongs to index i+1, and is stored in colors[i].

  The graph is laid out like this:

  The basic data structures are 2 arrays: one contiguous array of edge connections, and a rowpointer that
  specifies at which position in the edge array the next vertex starts. These 2 data structures
  are distributes to the MPI ranks and each rank works on their portion of these arrays.

  Colors are stored in local array of size n, and at each communication point this array is updated
  with the new values from the other ranks.

  This code for the greatest part works with unsigned 32bit integers. This means that
  the maximum number of edges in the graphs permitted is 4'294'967'295.

  The color array that is sent over the network uses unsigned 16bit integers. That should
  be enough to handle most colorings.

  Communication is done asynchronously, which complicates things heavily. The code is annotated as well
  as possible to make sure why we need everything.

 */



int main(int argc, char* argv[]) {
    std::string filename; // The graph filename
    int superstep; // The number of supersteps
    int stats; // If we monitor stats or not

    int rank, size;
    uint32_t nedges;

    // MPI Initialization
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Initialization from input

    // Getting the input and printing relevant errors
    if (argc < 4) {
        if (rank == 0) {
            std::cout << "Error!" << std::endl;
            std::cout << "Please run the program in the following format:" << std::endl;
            std::cout << "\t" << argv[0] << " input_file[string] superstep_size[integer] monitor_stats[0 or 1]" << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        return ERROR;
    }
    if (rank == 0) {
        std::cout << "Processing file: " << argv[1] << std::endl;
    }

    filename = argv[1];

    std::istringstream bla(argv[2]);
    if (!(bla >> superstep) || superstep < 1) {
        if (rank == 0) {
            std::cout << "Please enter a valid number for the superstep size" << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        return ERROR;
    }

    std::istringstream mem(argv[3]);
    if (!(mem >> stats) || (stats < 0 && stats > 1)) {
        if (rank == 0) {
            std::cout << "Please enter 0 or 1 for the stats monitoring" << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        return ERROR;
    }

    uint32_t* rowpointers;
    uint32_t* edges;
    uint32_t n;

    // Loading graph from file and checking for errors
    if (rank == 0) {
        nedges = loadGraph(n, edges, rowpointers, filename);
    }

    if (nedges < 0) {
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        return ERROR;
    }

    // Making sure all ranks wait until the graph is loaded
    MPI_Barrier(MPI_COMM_WORLD);

    // Sending "n" to all processes to get ready
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Number of vertices per rank
    uint32_t chunksize = (n + size - 1) / size; // Ceiling division

    // Number of vertices per rank, with the last rank being shrunk to its actual size.
    uint32_t localchunksize;
    if (rank != (size - 1)) {
        localchunksize = chunksize;
    } else {
        localchunksize = n - ((size - 1) * chunksize); // Last process needs to deal with less work
    }

    // Debugging: show the graph
    /*
    if (rank == 0) {

      for (int i = 0; i < n; i++) {
        int row_start = rowpointers[i];
        int row_end = (i == (n-1)) ? nedges : rowpointers[i+1];
        std::cout << "vertex " << i << ": ";
        for (int j = row_start; j < row_end; j++) {
          std::cout << edges[j] << ", ";
        }
        std::cout << "\n";
      }

      std::cout << "rowpointers: ";
      for (int i = 0; i < n; i++) {
        std::cout << i << " -> " << rowpointers[i] << ", ";
      }
      std::cout << "\n";

    }
     */

    uint32_t offset = rank * chunksize;
    uint32_t* localrowpointers = new uint32_t[chunksize];

    // Sending pieces to each rank: allocate receiving buffer and scatter
    MPI_Scatter(rowpointers, chunksize, MPI_INT, localrowpointers, chunksize, MPI_INT, 0, MPI_COMM_WORLD);

    uint32_t localnedges, blockend, blockstart, blocksize;

    // Sending size of adjacency list to each rank
    if (rank == 0) {
        for (int i = 0; i < size; i++) {
            blockstart = rowpointers[i * chunksize];
            blockend = (i < (size - 1)) ? rowpointers[(i + 1) * chunksize] : nedges;
            blocksize = blockend - blockstart;
            if (i == 0) {
                localnedges = blocksize;
            } else {
                MPI_Send(&blocksize, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
            }
        }
    } else {
        // recv for every rank
        MPI_Recv(&localnedges, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    uint32_t* localadjlist = new uint32_t[localnedges];

    // Sending local adjacency array. use scatterv?
    if (rank == 0) {
        for (int i = 0; i < size; i++) {
            blockstart = rowpointers[i * chunksize];
            blockend = (i < (size - 1)) ? rowpointers[(i + 1) * chunksize] : nedges;
            blocksize = blockend - blockstart;
            if (i == 0) {
                std::copy(edges, edges + blocksize, localadjlist);
            } else {
                MPI_Send(&edges[blockstart], blocksize, MPI_INT, i, 2, MPI_COMM_WORLD);
            }
        }
    } else {
        MPI_Recv(localadjlist, localnedges, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    //should be maxdegree + 1?
    int maxcolors = 500; //is this a good heuristic?

    uint16_t* colors = new uint16_t[size * chunksize](); // Array of size number of vertices, colors[i] is color of vertex i. make it too big to have everything nice when sending stuff...
    std::vector<int> forbidden(maxcolors, -1); // Forbidden data structure

    uint32_t pointeroffset = localrowpointers[0]; // Shifting the local row pointers to have base 0
    for (uint32_t i = 0; i < localchunksize; i++) {
        localrowpointers[i] -= pointeroffset;
    }

    // Finding the boundary points before starting to color
    std::vector<uint32_t> boundarypoints;
    for (uint32_t i = 0; i < localchunksize; i++) {
        uint32_t start = localrowpointers[i];
        uint32_t end = (localchunksize - 1 == i) ? localnedges : localrowpointers[i + 1];
        for (uint32_t j = start; j < end; j++) {
            int neighbor = localadjlist[j];
            if ((neighbor < offset) || (neighbor >= (offset + localchunksize))) {
                boundarypoints.push_back(i);
                break;
            }
        }
    }

    std::vector<uint32_t> conflicts; // Container for the vertices to work on
    std::vector<uint32_t> oldconflicts; //second container for conflict resolution
    for (uint32_t i = 0; i < localchunksize; i++) {
        conflicts.push_back(i);
    }

    // Status variable
    bool finished = false;
    int counter = 0;
    uint32_t max_conflicts = chunksize; // Initialized to do enough supersteps

    //we can not have bigger supersteps than chunksizes
    if (superstep > chunksize) {
      superstep = chunksize;
    }

    //these two arrays are the working set, the updates during each superstep are saved here
    uint32_t* reassigned = new uint32_t[superstep+1]; //reassigned vertices, first index holds number of recolored vertices
    uint16_t* reassigned_colors = new uint16_t[superstep]; //reassigned colors

    //those 3 arrays are needed as containers for every process and serve as send/receive buffers
    uint32_t** changed_vertices = new uint32_t*[size]; //changed vertices, array of pointers of size n. of procs
    uint16_t** changed_colors = new uint16_t*[size]; //changed colors, array of pointers of size n. of procs
    uint16_t** first_round_colors = new uint16_t*[size]; //colors in the first round, pointers of size n. of procs

    //init: one array per process
    for (int i = 0; i < size; i++) {
      changed_vertices[i] = new uint32_t[superstep+1];
      changed_colors[i] = new uint16_t[superstep];
      first_round_colors[i] = new uint16_t[superstep];
    }

    //mpi request array to make sure we wait for everyone
    MPI_Request* reqcolors = new MPI_Request[size];
    MPI_Request* reqvertices = new MPI_Request[size];
    MPI_Request* req = new MPI_Request[size];

    //flag is 0 now, set to 1 after first superstep. Needed because we need to start waiting for messages only in the second superstep.
    int it_flag = 0;

    // Memory usage stats
    long vmrss_per_process[size];
    long vmsize_per_process[size];
    long global_vmrss, global_vmsize;
    double average_global_vmrss = 0.0, average_global_vmsize = 0.0;

    MPI_Barrier(MPI_COMM_WORLD); // All threads should reach this place before we start the real work
    double starttime = MPI_Wtime();

    /*
     *
     * Starting the parallel execution
     *
     */

    while (finished == false) {

        int nit = 0; //iterations in this round
        int n_colored = 0; //vertices colored in this superstep
        counter += 1; //current round
        int n_superstep = 0; //current superstep
        int n_supersteps = (max_conflicts + (superstep - 1)) / superstep; //ceiling division - get max number of supersteps for this round
        //for all vertices to be colored

        size_t conf_size = conflicts.size();
        for (uint32_t ii = 0; ii < max_conflicts; ii++) { //iterate for as long as the process with the most conflicts
          if (ii < conf_size) { //if we still have conflicts to resolve
            size_t i = conflicts[ii];
            uint32_t start = localrowpointers[i];
            uint32_t end = (localchunksize-1 == i) ? localnedges : localrowpointers[i+1]; //make sure not to go out of bounds
            for (uint32_t j = start; j < end; j++) { //for each neighbor
              uint32_t neighbor = localadjlist[j];
              uint32_t vertexcolor = colors[neighbor];
              forbidden[vertexcolor] = i;
            }
            //this step makes sure starting from round 2 we store the vertices we changed. This was we can only send the changes.
            if (counter > 1) {
              n_colored += 1;
              reassigned[n_colored] = offset + i;
            }
            //first fit assignment
            for (int k = 0; k < maxcolors; k++) {
              if (forbidden[k] != i) {
                colors[offset + i] = k;
                if (counter > 1) {
                  reassigned_colors[n_colored-1] = k;
                }
                break;
              }
            }
          }

          //coloring is done here, communication starts
          if ((((ii+1) % superstep) == 0) || (ii == max_conflicts-1)) { //superstep done, communicate

            // This block waits for communcation from earlier loop iterations to complete,
            // and then places the received colors in the correct positions in the color array.
            if (it_flag == 1) { //only from 2nd superstep
              if (counter == 1) { //if we are in the first round -> take advantage of contiguous memory layout
                MPI_Waitall(size, req, MPI_STATUS_IGNORE);
                for (int p = 0; p < size; p++) {
                  if (p != rank) {
                    uint32_t startchunk = (p * chunksize) + ((n_superstep-1) * superstep); //current starting index
                    uint32_t sendsize = superstep;
                    if ((startchunk + superstep) > (size * chunksize)) { //make sure we don't overflow at the end of color array
                      sendsize = (size * chunksize) - startchunk;
                    }
                    std::memcpy(&colors[startchunk], first_round_colors[p], sendsize*2); //size is in BYTES!
                  }
                }
              } else { //from second round, it's complicated
                MPI_Waitall(size, reqcolors, MPI_STATUS_IGNORE);
                MPI_Waitall(size, reqvertices, MPI_STATUS_IGNORE);
                for (int p = 0; p < size; p ++) {
                  if (p != rank) {
                    uint32_t* localchanges = changed_vertices[p];
                    uint16_t* localcolors = changed_colors[p];
                    for (int t = 0; t < localchanges[0]; t++) {
                        colors[localchanges[t+1]] = localcolors[t]; //place at the right locations
                    }
                  }
                }
              }
            }
            it_flag = 1;
            //receiving is done here, we can start sending new stuff

            if (counter == 1) { //if first round, take advantage of contiguous layout
              for (int p = 0; p < size; p++) {
                if (p == rank) {
                  uint32_t startchunk = (p * chunksize) + (n_superstep * superstep); //current starting index
                  uint32_t sendsize = superstep;
                  if ((startchunk + superstep) > (size * chunksize)) { //make sure we don't overflow at the end of color array
                    sendsize = (size * chunksize) - startchunk;
                  }
                  std::memcpy(first_round_colors[p], &colors[startchunk], sendsize*2); //copy into a buffer such that we can immediately reuse the color array
                }
                MPI_Ibcast(first_round_colors[p], superstep, MPI_SHORT, p, MPI_COMM_WORLD, &req[p]);

                n_colored = 0;
              }
            } else { //from second round, it's compicated..
              reassigned[0] = n_colored; //store the number of elements first
              for (int p = 0; p < size; p++) {
                  if (rank == p) { //at the owning process, copy into send/recv buffer
                    std::memcpy(changed_colors[p], reassigned_colors, 2*superstep);
                    std::memcpy(changed_vertices[p], reassigned, 4*(superstep+1));
                  }
                  MPI_Ibcast(changed_colors[p], superstep, MPI_SHORT, p, MPI_COMM_WORLD, &reqcolors[p]);
                  MPI_Ibcast(changed_vertices[p], superstep+1, MPI_INT, p, MPI_COMM_WORLD, &reqvertices[p]);
              }
              n_colored = 0;
            }
            n_superstep += 1;
          }
          nit++;
        }

        //do all the receiving again because at the end of the round we have outstanding messages.
        if (counter == 1) {
          MPI_Waitall(size, req, MPI_STATUS_IGNORE);
          for (int p = 0; p < size; p++) {
            if (p != rank) {
              uint32_t startchunk = (p * chunksize) + ((n_superstep-1) * superstep); //current starting index
              uint32_t sendsize = superstep;
              if ((startchunk + superstep) > (size * chunksize)) { //make sure we don't overflow at the end of color array
                sendsize = (size * chunksize) - startchunk;
              }
              std::memcpy(&colors[startchunk], first_round_colors[p], sendsize*2);
            }
          }
        } else {
          MPI_Waitall(size, reqcolors, MPI_STATUS_IGNORE);
          MPI_Waitall(size, reqvertices, MPI_STATUS_IGNORE);
          for (int p = 0; p < size; p ++) {
            if (p != rank) {
              uint32_t* localchanges = changed_vertices[p];
              uint16_t* localcolors = changed_colors[p];
              for (int t = 0; t < localchanges[0]; t++) {
                  colors[localchanges[t+1]] = localcolors[t];
              }
            }
          }
        }
        it_flag = 0;

        // Conflict detection
        oldconflicts.swap(conflicts); //swap
        if (counter == 1) {
          oldconflicts = boundarypoints;
        }
        conflicts.clear();
        for (auto &i : oldconflicts) { // All boundary vertices
            uint32_t start = localrowpointers[i];
            uint32_t end = (localchunksize - 1 == i) ? localnedges : localrowpointers[i + 1];
            uint32_t owncolor = colors[i + offset];
            for (uint32_t j = start; j < end; j++) {
                uint32_t neighbor = localadjlist[j];
                uint32_t vertexcolor = colors[neighbor];
                if (vertexcolor == owncolor && (i + offset) <= neighbor) { // The smaller rank recolors
                    conflicts.push_back(i);
                    break;
                }
            }
        }

        // Get status from all processes.
        max_conflicts = conflicts.size();
        MPI_Allreduce(MPI_IN_PLACE, &max_conflicts, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if (max_conflicts == 0) {
            finished = true;
        }
        if (rank == 0) { //status update
            std::cout << "Round " << counter << " done.\n";
            std::cout << max_conflicts << " conflicts found\n";
        }

        // Getting memory stats
        if (stats == 1) {
            get_cluster_memory_usage_kb(vmrss_per_process, vmsize_per_process, 0, size);

            if (rank == 0) {
                for (int k = 0; k < size; k++) {
                    printf("Process %03d: VmRSS = %6ld KB, VmSize = %6ld KB\n", k, vmrss_per_process[k], vmsize_per_process[k]);
                }
            }

            long global_vmrss, global_vmsize;
            get_global_memory_usage_kb(&global_vmrss, &global_vmsize, size);
            average_global_vmrss += global_vmrss;
            average_global_vmsize += global_vmsize;
            if (rank == 0) {
                printf("Global memory usage: VmRSS = %6ld KB, VmSize = %6ld KB\n", global_vmrss, global_vmsize);
            }
        }
    }


    /*
     *
     *  Parallel execution ended
     *
     */

    MPI_Barrier(MPI_COMM_WORLD);
    double endtime = MPI_Wtime();

    // Output statistics
    if (rank == 0) {
        /*
        std::cout << "colors: ";
        for (int i = 0; i < n; i ++) {
          std::cout << colors[i] << ", ";
        }
         */
        std::cout << "\nfinished in " << counter << " rounds.\n";
        std::cout << endtime - starttime << " seconds [Wall time]\n";
        auto max = *std::max_element(colors, colors + n);
        std::cout << max << " colors used." << std::endl;

        // Checking for the validity of our results
        bool valid = true;

        for (uint32_t i = 0; i < n; i++) {
            int owncolor = colors[i];
            uint32_t start = rowpointers[i];
            uint32_t end = (i != n - 1) ? rowpointers[i + 1] : nedges;
            for (int j = start; j < end; j++) {
                int edge = edges[j];
                if (owncolor == colors[edge]) {
                    valid = false;
                    std::cout << "invalid results" << std::endl;
                    break;
                }
            }
        }

        std::ofstream results("results_smallmsg_nonb_firstfit.csv", std::ios_base::app | std::ios_base::out);
        if (stats == 0) {
            results << filename << "," << size << "," << n << "," << nedges << "," << (endtime - starttime) <<
                    "," << counter << "," << max << "," << superstep << "," << valid << "\n";
        } else {
            results << filename << "," << size << "," << n << "," << nedges << "," << (endtime - starttime) <<
                    "," << counter << "," << max << "," << superstep << "," << average_global_vmrss / (double) counter << "," <<
                    average_global_vmsize / (double) counter << "," << valid << "\n";
        }
    }

    MPI_Finalize();

    //Free memory
    if (rank == 0) {
      delete[] edges;
      edges = 0;
      delete[] rowpointers;
      rowpointers = 0;
    }
    delete[] reqcolors;
    reqcolors = 0;
    delete[] reqvertices;
    reqvertices =0;
    delete[] req;
    req = 0;
    delete[] localrowpointers;
    localrowpointers = 0;
    delete[] localadjlist;
    localadjlist = 0;
    delete[] colors;
    colors = 0;
    delete [] reassigned;
    reassigned = 0;
    delete [] reassigned_colors;
    reassigned_colors = 0;
    for (int i = 0; i < size; i++) {
      delete [] changed_colors[i];
      delete [] changed_vertices[i];
      delete [] first_round_colors[i];
    }
    delete [] changed_colors;
    changed_colors = 0;
    delete [] changed_vertices;
    changed_vertices = 0;
    delete [] first_round_colors;

    return 0;
}
