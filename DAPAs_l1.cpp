/* 
 * File:   main.cpp
 * Author: cse10ssn
 *
 * Created on May 24, 2015, 1:43 PM
 */

#include <stdio.h>
#include<iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sysexits.h>
#include <math.h>
#include <sys/time.h>
#include <mpi/mpi.h>
#include <stddef.h>
using namespace std;
#define debug 0
/*
 * 
 */
// Define Infinite (Using INT_MAX caused overflow problems)
#define INF 10000

static double gettime(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + 1e-6*tv.tv_usec;
}

typedef struct p {
    float x;
    float y;
} Point;


typedef enum {
	TURN_RIGHT = 1,
	TURN_LEFT = 2
} turn_t;

// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are colinear
// 1 --> Clockwise
// 2 --> Counterclockwise

int orientation(Point p, Point q, Point r) {
	float val = ((q.x - p.x) * (r.y - p.y)) - ((r.x - p.x) * (q.y - p.y));

	if (val <= 0)   // right turn  or collinear
		return 1;
	else         //left turn
		return 2;
}
/*
Point getPoint(float *points, int k)
{
	return Point
}
*/
// Prints convex hull of a set of n points.
void convexHull(Point points[], int n, int *Hull) {
	// There must be at least 3 points (6 floats)
	if (n < 3)
		return;

	// Initialize Result
//	int *Hull = new int[n];
	//int next[n];
	for (int k = 0; k < n; k++)
		Hull[k] = -1;

	// Find the leftmost point
	int l = 0;
	for (int k = 0; k < n; k++){
		if (points[k].x < points[l].x)
			l = k;
	}
	// Start from leftmost point, keep moving counterclockwise
	// until reach the start point again
	int q = l;
	int i = 0;
	do {
		Hull[i] = q;
		i = i + 1;
		int p = q;
		q = 0; // index of the first point in the points array
		// Search for a point 'q' such that orientation(p, i, q) is
		// counterclockwise for all points 'i'
		
		if (p == q)	// if P and Q become same point move the Q to the next point in the array(Assumption 1)
		{
			q = q + 1;
		}
		for (int j = 1; j < n; j++){
			if (j == p)		// If P and R(j) become same point then skip that point.
				continue;

			if (orientation(points[p], points[q], points[j]) == TURN_LEFT)
			{
				q = j;
			}
		}
	} while (q != Hull[0]);

/*	// Print Result
	for (int i = 0; i < n; i++) {
		if (Hull[i] != -1)
			cout << "(" << points[Hull[i]].x << ", " << points[Hull[i]].y << ")\n";
	}
*/
}

bool comparePoints(Point M, Point N)
{
	if(M.x == N.x && M.y == N.y) return true;
	else return false;
}	

void leftMost(void *M, void *N, int *len, MPI_Datatype *dptr)
{
	Point *m = (Point *)M;
	Point *n = (Point *)N;
	if(m->x < n->x )
	{ 
		n->x = m->x;
		n->y = m->y;
	}
}



Point FindNextHullPoint(Point p, Point points[], int n) {
	int q = 0; // index of the first point in the points array
	// Search for a point 'q' such that orientation(p, q, r(j)) is
	// clockwise for all points 'j'

	if (comparePoints(p, points[q]))	// if P and Q become same point move the Q to the next point in the array(Assumption 1)
	{
		q = q + 1;
	}
	for (int j = 1; j < n; j++){
		if (comparePoints(p, points[j]))		// If P and R(j) become same point then skip that point.
			continue;

		if (orientation(p, points[q], points[j]) == TURN_LEFT)
		{
			q = j;
		}
	}
	return points[q];
}

// Driver program to test above functions

float rand_coord(float limit) {
  float safe_limit = limit == 0 ? 1l : limit;
  return (fmod((float) rand() , (float) safe_limit)); 
}

Point random_point(float limit) {
  return (Point) { rand_coord(limit), rand_coord(limit) };
}


int main(int argc, char *argv[]) {

	int myid, numprocs, i = 0, m = 40;

 if (argc >= 2) {
	  
  m = atoi(argv[1]);
    if (m < 40)
      m = 40;
  }else{
		printf("you forgot total number of points!");
		return 0;
	}

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	/* create a type for struct car */
	const int nitems = 2;
	int blocklengths[2] = {1, 1};
	MPI_Datatype types[2] = {MPI_FLOAT, MPI_FLOAT};
	MPI_Datatype mpi_point_type;
	MPI_Aint offsets[2];
	offsets[0] = offsetof(Point, x);
	offsets[1] = offsetof(Point, y);
	MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_point_type);
	MPI_Type_commit(&mpi_point_type);
	MPI_Op MPI_Leftmost;
	MPI_Op_create( leftMost, true, &MPI_Leftmost );

	Point *send = NULL;
	Point *Hull = NULL;
	Point *sendQ = NULL;
	Point *recvQ = NULL;

	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	double t1 = gettime();
	if (myid == 0) {
//			int 
			send=new Point[m];
			for(i = 0; i<m;++i)
			{
				send[i] = random_point(100);
			}
		if(debug){
			printf("Root has elements");
			for (i=0; i<m; i++) {
				printf(" (%.0f, %.0f)", send[i].x, send[i].y);
				if(i%9 == 0)	printf("\n");
			}
			printf("\n");
		}
		Hull = new Point[m];
	}
	int n = 10;
	n = m/numprocs;
	Point *recv=new Point[n];

	MPI_Scatter(send, n, mpi_point_type, recv, n, mpi_point_type, 0,MPI_COMM_WORLD);

	if(debug){
		printf("Process %d has elements:", myid);
		for (i=0; i<n; i++) {
			printf(" (%.0f, %.0f)", recv[i].x, recv[i].y);
		}
		printf("\n");
	}

//	int *Hull = new int[n];
//	if (n < 3)
//		return;

	// Initialize Result
//	int *Hull = new int[n];
	//int next[n];
//	for (int k = 0; k < n; k++)
//		Hull[k] = -1;

	// Find the local leftmost point
	int l = 0;
	for (int k = 0; k < n; k++){
		if (recv[k].x < recv[l].x)
			l = k;
	}
	printf("I'm %d and i have (%f, %f) \n", myid, recv[l].x, recv[l].y);	
	MPI_Barrier( MPI_COMM_WORLD ) ;
	Point llm = recv[l], glm;
	glm.x = 2.4;
	glm.y = 5.4;
	MPI_Allreduce(&llm, &glm, 1, mpi_point_type, MPI_Leftmost, MPI_COMM_WORLD);
	printf("I'm %d and i have (%f, %f) \n", myid, glm.x, glm.y);
	MPI_Barrier( MPI_COMM_WORLD ) ;

	// Start from leftmost point, keep moving counterclockwise
	// until reach the start point again
	Point q = glm;
	i = 0;
	sendQ = new Point[numprocs];
	recvQ = new Point[numprocs];

	do {
		if(myid == 0){
			Hull[i] = q;
		}
		Point p = q;
		q = FindNextHullPoint(p, recv, n); 
		for (int k = 0; k < numprocs; k++) sendQ[k] = q;
		printf("\n I'm %d and i've got (%f, %f)\n", myid, q.x, q.y);	
		MPI_Barrier( MPI_COMM_WORLD ) ;
		MPI_Alltoall(sendQ, 1, mpi_point_type, recvQ, 1, mpi_point_type, MPI_COMM_WORLD);
		if(myid == 0){
			for (int k = 0; k < numprocs; k++) 
				printf("\n I'm %d and i've got (%f, %f)\n", myid, recvQ[k].x, recvQ[k].y);	
		}
		q = FindNextHullPoint(p, recvQ, numprocs); 
		i = i + 1;
	//printf("\n I'm %d and i'm running %d time\n", myid, i);	
	} while((i<n) && !comparePoints(q, glm));
//	printf("\n I'm %d and i've got (%f, %f)\n", myid, q.x, q.y);	
	if(myid == 0){
		// Print Result
		printf("Process %d found Hull( elements):", myid);
		for (int i = 0; i < m; i++) {
				    cout << "(" << Hull[i].x << ", " << Hull[i].y << ") ";		      
			printf("\n");
		}
	}


	MPI_Barrier( MPI_COMM_WORLD ) ;
	double t2 = gettime();
	double t3 = t2 - t1;
	printf("\n Run Time = %.6lf\n", t3);
	MPI_Type_free(&mpi_point_type);
	MPI_Finalize();
	return 0;
}
