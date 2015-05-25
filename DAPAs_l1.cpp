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
using namespace std;

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

struct Point {
	float x;
	float y;
};

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
// Driver program to test above functions

float rand_coord(float limit) {
  float safe_limit = limit == 0 ? 1l : limit;
  return (fmod((float) rand() , (float) safe_limit)); 
}

Point random_point(float limit) {
  return (Point) { rand_coord(limit), rand_coord(limit) };
}

int main(int argc, char *argv[]) {
	      
	int m = 100;

	// generate random inputs points.
	m = 80;
	float *points1 = new float[m];
	int q = 20;
	float *points2 = new float[q];
		 
	int myid, numprocs, i = 0;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	printf("Who am i = %d\n", myid);

	if(myid == 0){
		for(i = 0; i<m;++i)
		{
			points1[i] = rand_coord(100);
		}

/* Print out own portion of the scattered array */
    printf("Root has elements");
    for (i=0; i<m; i++) {
      printf(" %.0f", points1[i]);
			if(i%19 == 0)	printf("\n");
    }
    printf("\n");

  }

	MPI_Scatter(points1, q, MPI_FLOAT, points2, q, MPI_FLOAT, 0,MPI_COMM_WORLD);

	int pointsLength = 10;
	Point *myPoints = new Point[pointsLength];
	int *Hull = new int[pointsLength];
  for (i=0; i<pointsLength; i++) {
		myPoints[i].x = points2[i*2];
		myPoints[i].y = points2[i*2+1];
  }

  printf("Process %d has elements:", myid);
    for (i=0; i<pointsLength; i++) {
//      printf(" %.0f", points2[i]);
			cout << "(" << myPoints[i].x << ", " << myPoints[i].y << ") ";
    }
    printf("\n");

	double t1 = gettime();
	convexHull(myPoints, pointsLength, Hull);
	double t2 = gettime();
	double t3 = t2-t1;

// Print Result
	printf("Process %d found Hull( elements):", myid);
	for (int i = 0; i < pointsLength; i++) {
		if (Hull[i] != -1)
			cout << "(" << myPoints[Hull[i]].x << ", " << myPoints[Hull[i]].y << ") ";
	}
    printf("\n");
	printf("\n Run Time = %.6lf\n", t3);

        
        
  MPI_Finalize();
	return 0;
}
