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

typedef struct HomographicMatrix {
	//[x']   [ H11, H12, H13 ][x]
	//[y'] = [ H21, H22, H23 ][y]
	//[1 ]   [ H31, H32,  1  ][p]
	double H11;
	double H12;
	double H13;
	double H21;
	double H22;
	double H23;
	double H31;
	double H32;
} homographicmatrix;