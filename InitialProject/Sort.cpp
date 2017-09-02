#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct Point
{
	int x;
	int y;
	double distance;
};

enum eDirections
{
	ROWS,
	COLS
};

// Methods
int readPointFromFile(Point*& arrPoints, char* fileName)
{	
	int i = 0;
	int numOfPoints = 0;
	FILE *file;

	// Open the points' file
	file = fopen(fileName, "r");

	if (file) 
	{
		// Reading number of points from file
		fscanf (file, "%d", &numOfPoints);    

		// Allocation the array of points
		arrPoints = new Point[numOfPoints];

		// Reading points from file
		for(i = 0; i < numOfPoints; i++)
		{
			fscanf (file, "%d", &arrPoints[i].x);
			fscanf (file, "%d", &arrPoints[i].y);
		}

		fclose(file);
	}

	return numOfPoints;
}

void createPointMpiDataType(MPI_Datatype& PointMPIType)
{
	struct Point point;	
	MPI_Datatype type[3] = { MPI_INT, MPI_INT, MPI_DOUBLE };
	int blocklen[3] = { 1, 1, 1 };
	MPI_Aint disp[3];

	disp[0] = (char *) &point.x - (char *) &point;
	disp[1] = (char *) &point.y - (char *) &point;
	disp[2] = (char *) &point.distance - (char *) &point;

	MPI_Type_create_struct(3, blocklen, disp, type, &PointMPIType);
	MPI_Type_commit(&PointMPIType);	
}

void createCart(int numOfProcs, MPI_Comm& comm)
{
	int dim[2], period[2], reorder;		

	dim[0] = dim[1] = (int)(sqrt((float)numOfProcs));
	period[0] = period[1] = 0;		
	reorder = 1;		
	MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &comm); 	
}

double calcDistance(int x, int y)
{
	return sqrt((double)(x * x + y * y));
}

void oddEvenStep(int* myCoords, int* neighbourCoords, Point& myPoint, char compareOpetor, MPI_Comm cartComm,const MPI_Datatype& PointMPIType, MPI_Status& status)
{
	int neighbourRank;
	Point neighbourPoint;

	// Getting the neighbour's rank
	MPI_Cart_rank(cartComm, neighbourCoords, &neighbourRank);

	// Sending my point to the neighbour and receiving its point back
	MPI_Sendrecv(&myPoint, 1, PointMPIType, neighbourRank, 0, &neighbourPoint, 1, PointMPIType, neighbourRank, 0, MPI_COMM_WORLD, &status);

	// Checking if a switch between points is needed
	if ((compareOpetor == '>' && myPoint.distance > neighbourPoint.distance) || 
		(compareOpetor == '<' && myPoint.distance < neighbourPoint.distance))
	{
		myPoint.x = neighbourPoint.x;
		myPoint.y = neighbourPoint.y;
		myPoint.distance = neighbourPoint.distance;
	}
}

void oddEvenSort(int numOfElements, eDirections direction, int rank, MPI_Comm& cartComm, Point& myPoint, const MPI_Datatype& PointMPIType, MPI_Status& status)
{
	int coords[2], neighbourCoords[2], neighbourRank;
	Point neighbourPoint;

	// Getting the coords of the process in the cart
	MPI_Cart_coords(cartComm, rank, 2, coords);

	for (int i = 0; i < numOfElements; i++)
	{
		// Even step of odd-even sort
		if (i % 2 == 0)
		{
			if (direction == COLS)
			{
				if (coords[0] % 2 == 0) // Even row
				{
					if (coords[0] != numOfElements - 1)
					{
						// Setting the neighbour's coords (the top one)
						neighbourCoords[1] = coords[1];
						neighbourCoords[0] = coords[0] + 1;

						oddEvenStep(coords, neighbourCoords, myPoint, '>', cartComm, PointMPIType, status);
					}
				}
				else // Odd row
				{
					// Setting the neighbour's coords (the bottom one)
					neighbourCoords[1] = coords[1];
					neighbourCoords[0] = coords[0] - 1;

					oddEvenStep(coords, neighbourCoords, myPoint, '<', cartComm, PointMPIType, status);
				}
			}
			else // ROWS
			{
				if (coords[1] % 2 == 0) // Even col
				{
					if (coords[1] != numOfElements - 1)
					{
						// Setting the neighbour's coords (the right one)
						neighbourCoords[1] = coords[1] + 1;
						neighbourCoords[0] = coords[0];

						if(coords[0] % 2 == 0) // Even row
							oddEvenStep(coords, neighbourCoords, myPoint, '>', cartComm, PointMPIType, status);
						else // Odd row
							oddEvenStep(coords, neighbourCoords, myPoint, '<', cartComm, PointMPIType, status);
					}
				}
				else // // Odd col
				{
					// Setting the neighbour's coords (the left one)
					neighbourCoords[1] = coords[1] - 1;
					neighbourCoords[0] = coords[0];

					if(coords[0] % 2 == 0) // Even row
						oddEvenStep(coords, neighbourCoords, myPoint, '<', cartComm, PointMPIType, status);
					else // Odd row
						oddEvenStep(coords, neighbourCoords, myPoint, '>', cartComm, PointMPIType, status);
				}
			}
		}
		else // Odd step of odd-even sort
		{
			if (direction == COLS)
			{
				if (coords[0] % 2 == 0) // Even row
				{
					if (coords[0] != 0)
					{
						// Setting the neighbour's coords (the bottom one)
						neighbourCoords[1] = coords[1];
						neighbourCoords[0] = coords[0] - 1;

						oddEvenStep(coords, neighbourCoords, myPoint, '<', cartComm, PointMPIType, status);
					}
				}
				else // Odd row
				{
					if (coords[0] != numOfElements - 1)
					{
						// Setting the neighbour's coords (the top one)
						neighbourCoords[1] = coords[1];
						neighbourCoords[0] = coords[0] + 1;

						oddEvenStep(coords, neighbourCoords, myPoint, '>', cartComm, PointMPIType, status);
					}
				}
			}
			else // ROWS
			{
				if (coords[1] % 2 == 0) // Even col
				{
					if (coords[1] != 0)
					{
						// Setting the neighbour's coords (the left one)
						neighbourCoords[1] = coords[1] - 1;
						neighbourCoords[0] = coords[0];

						if(coords[0] % 2 == 0) // Even row
							oddEvenStep(coords, neighbourCoords, myPoint, '<', cartComm, PointMPIType, status);
						else // Odd row
							oddEvenStep(coords, neighbourCoords, myPoint, '>', cartComm, PointMPIType, status);
					}
				}
				else // // Odd col
				{
					if (coords[1] != numOfElements - 1)
					{
						// Setting the neighbour's coords (the right one)
						neighbourCoords[1] = coords[1] + 1;
						neighbourCoords[0] = coords[0];

						if(coords[0] % 2 == 0) // Even row
							oddEvenStep(coords, neighbourCoords, myPoint, '>', cartComm, PointMPIType, status);
						else // Odd row
							oddEvenStep(coords, neighbourCoords, myPoint, '<', cartComm, PointMPIType, status);
					}
				}
			}
		}

	}
}

int main(int argc,char *argv[])
{
	int  namelen, numprocs, myid;
	MPI_Status status;
	MPI_Datatype PointMPIType;
	Point* arrPoints = nullptr;
	Point myPoint;	
	MPI_Comm cartComm;
	int coord[2];
	int numOfPoints;
	int shearSortCounter = 0;
	int maxShearSortIterations;
	eDirections direction;
	int currCoords[2];
	int currRank;
	int index = 0;
	Point* arrSortedPoints;

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);			

	if (myid == 0) 
	{		
		char FILE_NAME[] = "C:\\Users\\Zuzu\\Desktop\\Points.txt";
		numOfPoints = readPointFromFile(arrPoints, FILE_NAME);

		// Calculating the distance for each point
		for(int i = 0; i < numOfPoints; i++)
			arrPoints[i].distance = calcDistance(arrPoints[i].x, arrPoints[i].y);

		// Sending the number of points to all processes
		for(int i = 1; i < numprocs; i++)
			MPI_Send(&numOfPoints, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
	}
	else
	{
		// Receiving the number of points from process 0
		MPI_Recv(&numOfPoints, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	}

	// Creating MPI data type of point
	createPointMpiDataType(PointMPIType);

	// Scatter the points to all processes
	MPI_Scatter(arrPoints, 1, PointMPIType, &myPoint, 1, PointMPIType, 0, MPI_COMM_WORLD);
	delete []arrPoints;

	// Create cart
	createCart(numprocs, cartComm);
	MPI_Cart_coords(cartComm, myid, 2, coord); // Coord[0] - Y axis, Coord[1] - X axis

	// Calculating the number of iterations that need to be executed in share sort
	maxShearSortIterations = (log((double)numOfPoints) / log(2.0)) + 1;
	
	// Shear sort
	while(shearSortCounter < maxShearSortIterations)
	{
		// Even iteration
		if(shearSortCounter % 2 == 0)
		{
			// Sorting by rows
			direction = ROWS;			
		}
		else // Odd iteration
		{
			// Sorting by columns
			direction = COLS;
		}

		// Odd even sort
		oddEvenSort(sqrt((float)numOfPoints),direction, myid, cartComm, myPoint, PointMPIType, status);

		shearSortCounter++;
	}

	if (myid == 0)
	{
		arrSortedPoints = new Point[numOfPoints];
		for (int row = (int)sqrt((float)numOfPoints) - 1; row >= 0; row--)
		{			
			if (row % 2 != 0)
			{
				for (int col = 0; col < sqrt((float)numOfPoints); col++)
				{
					currCoords[1] = col;
					currCoords[0] = row;

					MPI_Cart_rank(cartComm, currCoords, &currRank);
					if (currRank != 0)
						MPI_Recv(arrSortedPoints + index, 1, PointMPIType, currRank, 0, MPI_COMM_WORLD, &status);
					else
						arrSortedPoints[index] = myPoint;
					index++;
				}
			}
			else
			{
				for (int col = sqrt((float)numOfPoints) - 1; col >= 0; col--)
				{
					currCoords[1] = col;
					currCoords[0] = row;

					MPI_Cart_rank(cartComm, currCoords, &currRank);
					if (currRank != 0)
						MPI_Recv(arrSortedPoints + index, 1, PointMPIType, currRank, 0, MPI_COMM_WORLD, &status);
					else
						arrSortedPoints[index] = myPoint;
					index++;
				}
			}
		}

		for (int i = 0; i < numOfPoints; i++)
		{
			printf("X = %d, Y = %d, Dist = %.2f\n", arrSortedPoints[i].x, arrSortedPoints[i].y, arrSortedPoints[i].distance);fflush(stdout);
		}
		delete []arrSortedPoints;
	}
	else
	{
		// All the processes except 0 send their point to process 0
		MPI_Send(&myPoint, 1, PointMPIType, 0, 0, MPI_COMM_WORLD);
	}

	MPI_Finalize();	
	return 0;
}
