#pragma once

typedef struct Descripter {
	int x;
	int y;
	int direc;
	double vector[128];
} descripter;

typedef struct Point {
	int x;
	int y;
} point;


typedef struct PointPair {
	struct Point p1;
	struct Point p2;
} pointpair;

typedef struct affinematrix {
	//[x']   [ A1, B1, e ][x]
	//[y'] = [ A2, B2, f ][y]
	//[1 ]   [ 0,  0,  1 ][1]
	double A1;
	double B1;
	double A2;
	double B2;
	double e;
	double f;
} affinematrix;
