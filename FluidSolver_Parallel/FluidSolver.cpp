#include "precice/SolverInterface.hpp"
#include "FluidSolver.h"

using namespace precice;
using namespace precice::constants;

int main(int argc, char** argv){

	std::cout << "Starting Fluid Solver..." << std::endl;
	if ( argc != 5 )
	{
		std::cout << std::endl;
		std::cout << "Fluid: Usage: mpiexec -np <#procs> " << argv[0] << " <configurationFileName> <N> <tau> <kappa>" << std::endl;
		std::cout << std::endl;
		std::cout << "N:     Number of mesh elements, needs to be equal for fluid and structure solver." << std::endl;
		std::cout << "tau:   Dimensionless time step size." << std::endl;
		std::cout << "kappa: Dimensionless structural stiffness." << std::endl;

		return -1;
	}

    MPI_Init(&argc, &argv);

    int domainSize, gridOffset, rank, size, chunkLength, *vertexIDs;
    double *pressure, *pressure_n, *crossSectionLength, *crossSectionLength_n, *velocity, *velocity_n, *grid, tau, kappa;   // Declare dataset

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    domainSize = atoi(argv[2]);
    tau = 		 atof(argv[3]);
    kappa = 	 atof(argv[4]);

    std::cout << "tau: " << tau << ", kappa: " << kappa << std::endl;

    if ((domainSize+1) % size == 0) {
    	chunkLength = (domainSize+1)/size;
    	gridOffset = rank * chunkLength;
    } else if (rank < (domainSize+1) % size) {
    	chunkLength = (domainSize+1)/size + 1;
    	gridOffset = rank * chunkLength;
    } else {
    	chunkLength = (domainSize+1)/size;
    	gridOffset = ((domainSize+1) % size) * ((domainSize+1)/size + 1) + (rank - ((domainSize+1) % size)) * (domainSize+1)/size;
    }

    pressure = 			new double[chunkLength];
    pressure_n = 		new double[chunkLength];
    crossSectionLength = 	new double[chunkLength];
    crossSectionLength_n = 	new double[chunkLength];
    velocity = 			new double[chunkLength];
    velocity_n = 		new double[chunkLength];

    fluidDataDisplay(pressure, chunkLength);

    std::string configFileName(argv[1]);
    std::string solverName = "FLUID";

    SolverInterface interface(solverName, rank, size);
    interface.configure(configFileName);

    int meshID = interface.getMeshID("Fluid_Nodes");
    int pressureID = interface.getDataID("Stresses", meshID);
    int crossSectionLengthID = interface.getDataID("Displacements", meshID);

    int dimensions = interface.getDimensions();
    grid = new double[dimensions * chunkLength];
    vertexIDs = new int[chunkLength];

	for (int i = 0; i < chunkLength; i++) {
		pressure[i] = 0.0;
		pressure_n[i] = 0.0;
		crossSectionLength[i] = 1.0;
		crossSectionLength_n[i] = 1.0;
		velocity[i] = 1.0 / kappa;
		velocity_n[i] = 1.0 / kappa;
		std::cout << "[  ";
		for (int j = 0; j < dimensions; j++) {
			grid[i*dimensions + j] = j == 0 ? gridOffset + (double)i : 0.0;
			std::cout << grid[i*dimensions + j] << "  ";
		}
		std::cout << "]" << std::endl;
	}

    interface.setMeshVertices(meshID, chunkLength, grid, vertexIDs);

    interface.initialize();

    double t = 0;
    double dt = 1.0;

	if (interface.isActionRequired(actionWriteInitialData())) {
		std::cout << "Fluid: Writing initial data.." << std::endl;
		interface.writeBlockScalarData(pressureID, chunkLength, vertexIDs, pressure);
		interface.fulfilledAction(actionWriteInitialData());
	}

    interface.initializeData();

	if (interface.isReadDataAvailable()) {
		std::cout << "Fluid: Reading initial data.." << std::endl;
		interface.readBlockScalarData(crossSectionLengthID, chunkLength, vertexIDs, crossSectionLength);
	}

    while (interface.isCouplingOngoing()){
    	int convergenceCounter = 0;
    	if (interface.isActionRequired(actionWriteIterationCheckpoint())){
    		interface.fulfilledAction(actionWriteIterationCheckpoint());
		}

        fluidComputeSolution(rank, size, domainSize, chunkLength, kappa, tau, 1, t, pressure, pressure_n, pressure, crossSectionLength, crossSectionLength_n, velocity, velocity_n);// Call "Solver"

        fluidDataDisplay(pressure, chunkLength);
        fluidDataDisplay(crossSectionLength, chunkLength);

		interface.writeBlockScalarData(pressureID, chunkLength, vertexIDs, pressure);

        interface.advance(dt);

		interface.readBlockScalarData(crossSectionLengthID, chunkLength, vertexIDs, crossSectionLength);

		if (interface.isActionRequired(actionReadIterationCheckpoint())){ // i.e. not yet converged
			interface.fulfilledAction(actionReadIterationCheckpoint());
			convergenceCounter++;
		}
		else {
			std::cout << "Fluid" << rank << ": CONVERGED in " << convergenceCounter << " iterations." << std::endl;
			t += dt;
			for (int i = 0; i < chunkLength; i++) {
				pressure_n[i] = pressure[i];
				velocity_n[i] = velocity[i];
				crossSectionLength_n[i] = crossSectionLength[i];
			}
		}
    }

    std::cout << "Fluid" << rank << ": Time-loop finished." << std::endl;

    delete(pressure);
    delete(pressure_n);
    delete(crossSectionLength);
    delete(crossSectionLength_n);
    delete(velocity);
    delete(velocity_n);
    delete(grid);

    std::cout << "Fluid" << rank << ": Memory deallocation done." << std::endl;

    MPI_Finalize();

    std::cout << "Fluid: MPIFinalize() done. Exiting.." << std::endl;

    return 0;
}
