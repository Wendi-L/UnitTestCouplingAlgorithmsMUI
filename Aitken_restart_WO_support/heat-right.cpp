/*
 * heat-right.cpp
 *
 *  Created on: Jul 16, 2019
 *      Author: W. Liu
 */

#include "mui.h"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

/*                             Left Domain                       Right Domain             
 * Coarse : +-------+-------+-------+-------o=======+=======o-------+-------+-------+-------+
 *          0       1       2       3       4       5       6       7       8       9      10
 * +: grid points
 * o: interface points
 * -: single domain zone
 * =: overlapping zone
 */

// USAGE: mpirun -np 1 ./heat-left : -np 1 ./heat-right

int main( int argc, char ** argv ) {
    using namespace mui;

    const static int N = 11;
    double u1[N], u2[N]; 

    /// Initialise values from file
    std::string inoutFilenameL = "right_AITKEN.csv";
	std::fstream inFile;
	std::vector<std::vector<std::string>> content;
	std::vector<std::string> row;
    std::string line, word;
    inFile.open(inoutFilenameL, std::ios::in);
    if(inFile.is_open())
    {
		while(getline(inFile, line)) {
			row.clear();
			std::stringstream str(line);
			while (std::getline(str, word, ',')) {
				row.push_back(word);
		   }
		   content.push_back(row);
		}

		for ( int i = 4; i <  11; i++ ) u1[i] = stod(content[i-3][1]);

	} else {
		std::cerr<<"right_AITKEN.csv missing" << std::endl;
		for ( int i = 4; i <  11; i++ ) u1[i] = 0.0;
	}

    inFile.close();

    uniface1d interface( "mpi://right/ifs" );

	MPI_Comm  world = mui::mpi_split_by_app();
	MPI_Comm*  Cppworld = &world;    
    int rankLocal = MPI::COMM_WORLD.Get_rank();
    int sizeLocal = MPI::COMM_WORLD.Get_size();
    
    int rank, size;
    MPI_Comm_rank( world, &rank );
    MPI_Comm_size( world, &size );

    /// Create rbf matrix folder
    std::string makedirMString = "results_right" + std::to_string(rank);
    mkdir(makedirMString.c_str(), 0777);
    std::string fileAddress(makedirMString);

    double        k = 0.515, H = 1;
    double *      u = u1, *v = u2;

    std::vector<std::pair<mui::point1d, double>> ptsVluInit;

    for ( int i = 4; i <  11; i++ ) {
		mui::point1d pt(i);
		ptsVluInit.push_back(std::make_pair(pt,u1[i]));
	}

    // fetch data from the other solver
    sampler_pseudo_nearest_neighbor1d<double> s1(0.1);
    chrono_sampler_exact1d  s2;
	algo_aitken1d aitken(0.01,1.0);
	//algo_aitken1d aitken(0.514535,1.0,world,ptsVluInit);

     // Print off a hello world message
    printf("Hello world from Right rank %d out of %d MUI processors\n",
           rank, size);
           
     // Print off a hello world message
    printf("Hello world from Right rank %d out of %d local processors\n",
           rankLocal, sizeLocal);

    /// Output
    std::ofstream outputFileRight;
    std::string filenameR = "results_right" + std::to_string(rank) + "/solution-right_AITKEN_0.csv";
    outputFileRight.open(filenameR);
    outputFileRight << "\"X\",\"u\"\n";
    for ( int i = 4; i <  11; i++ ) outputFileRight << i * H << "," << u[i] << ", \n";
    outputFileRight.close();

    for ( int t = 1; t <= 1000; ++t ) {
        printf( "Right grid step %d\n", t );

            u[4] = interface.fetch( "u", 4 * H, t, s1, s2, aitken );
            printf( "Right under relaxation factor at t= %d is %f\n", t, aitken.get_under_relaxation_factor(t));
            printf( "Right residual L2 Norm at t= %d is %f\n", t, aitken.get_residual_L2_Norm(t));

            // calculate 'interior' points
            for ( int i = 5; i <  11; i++ ) v[i] = u[i] + k / ( H * H ) * ( u[i - 1] + u[i + 1] - 2 * u[i] );
            // calculate 'boundary' points
            v[N - 1] = 0.0;
            v[4]     = u[4    ];

            // push data to the other solver
            interface.push( "u0", 6 * H, u[6] );
            interface.commit( t );
        // I/O
        std::swap( u, v );

        /// Output
        std::ofstream outputFileRight;
        std::string filenameR = "results_right" + std::to_string(rank) + "/solution-right_AITKEN_"
          + std::to_string(t) + ".csv";
        outputFileRight.open(filenameR);
        outputFileRight << "\"X\",\"u\"\n";
        for ( int i = 4; i <  11; i++ ) outputFileRight << i * H << "," << u[i] << ", \n";
        outputFileRight.close();
    }

    return 0;
}
