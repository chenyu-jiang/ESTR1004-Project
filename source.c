#include <stdio.h>
#include "qdbmp.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define DEBUG
// REGARDS TO http://qdbmp.sourceforge.net/

#define PI 3.14159
#define LOGE 2.71828

void RGBtoBW(BMP* bmp);
void invert(BMP* bmp);
BMP* zoom(BMP* bmp, int k);
void reduceLevel(BMP* bmp, int level);
void BritandCntr(BMP* bmp, int brightness, double contrast);
BMP* AMRotation(BMP* bmp, double theta);
BMP* GaussianBlur(BMP* bmp,double theta,int radius);

int main()
{
	BMP*    bmp;
	char filename[] = "E:\\Lesson units\\ESTR 1004\\Project\\Area rotation\\lena.bmp";
	bmp = BMP_ReadFile(filename);
	BMP_CHECK_ERROR(stderr, -1);
	/////////////////////////////////////////////////////////////////////////
	//Your code in between
	BMP* newbmp = AMRotation(bmp, PI/6);
	/////////////////////////////////////////////////////////////////////////
	/* Save result */
	BMP_WriteFile(newbmp, "myout.bmp");
	BMP_CHECK_ERROR(stderr, -2);
	/* Free all memory allocated for the image */
	BMP_Free(bmp);
	BMP_Free(newbmp);
	//BMP_Free(output);
	return 0;
}

void RGBtoBW(BMP* bmp)
{
	UCHAR   r, g, b;
	UINT    width, height;
	UINT    x, y;
	width = BMP_GetWidth(bmp);
	height = BMP_GetHeight(bmp);
	for (x = 0; x < width; x++)
	{
		for (y = 0; y < height; y++)
		{
			BMP_GetPixelRGB(bmp, x, y, &r, &g, &b);
			BMP_SetPixelRGB(bmp, x, y, 0.3*r + 0.59*g + 0.11*b, 0.3*r + 0.59*g + 0.11*b, 0.3*r + 0.59*g + 0.11*b);
		}
	}
}

void invert(BMP* bmp)
{
	UCHAR   r, g, b;
	UINT    width, height;
	UINT    x, y;
	width = BMP_GetWidth(bmp);
	height = BMP_GetHeight(bmp);
	for (x = 0; x < width; x++)
	{
		for (y = 0; y < height; y++)
		{
			BMP_GetPixelRGB(bmp, x, y, &r, &g, &b);
			BMP_SetPixelRGB(bmp, x, y, 255 - r, 255 - g, 255 - b);
		}
	}
}

BMP* zoom(BMP* bmp, int k)
{
	UINT    width, height;
	width = BMP_GetWidth(bmp);
	height = BMP_GetHeight(bmp);
	BMP* newim = BMP_Create(width*k, height*k, BMP_GetDepth(bmp));
	for (UINT iX = 0; iX < width; iX++)
	{
		for (UINT iY = 0; iY < height; iY++)
		{
			UCHAR difr, difg, difb;
			UCHAR tmpr, tmpg, tmpb;
			BMP_GetPixelRGB(bmp, iX, iY, &tmpr, &tmpg, &tmpb);
			BMP_GetPixelRGB(bmp, iX, iY + 1, &difr, &difg, &difb);
			difr = -(tmpr - difr) / k;
			difg = -(tmpg - difg) / k;
			difb = -(tmpb - difb) / k;
			for (int i = 0; i < k; i++)
			{
				BMP_SetPixelRGB(newim, iX*k, iY*k + i, tmpr + difr*i, tmpg + difg*i, tmpb + difb*i);
			}
		}
	}
	for (UINT iY = 0; iY < height*k; iY++)
	{
		for (UINT iX = 0; iX < width; iX++)
		{
			UCHAR difr, difg, difb;
			UCHAR tmpr, tmpg, tmpb;
			BMP_GetPixelRGB(newim, iX*k, iY, &tmpr, &tmpg, &tmpb);
			BMP_GetPixelRGB(newim, iX*k + k, iY, &difr, &difg, &difb);
			difr = -(tmpr - difr) / k;
			difg = -(tmpg - difg) / k;
			difb = -(tmpb - difb) / k;
			for (int i = 0; i < k; i++)
			{
				BMP_SetPixelRGB(newim, iX*k + i, iY, tmpr + difr*i, tmpg + difg*i, tmpb + difb*i);
			}
		}
	}
	return newim;
}

void reduceLevel(BMP* bmp, int level)
{
	UINT    width, height;
	USHORT depth;
	width = BMP_GetWidth(bmp);
	height = BMP_GetHeight(bmp);
	depth = BMP_GetDepth(bmp);
	for (UINT i = 0; i < width; i++)
	{
		for (UINT j = 0; j < height; j++)
		{
			UCHAR r,g,b;
			BMP_GetPixelRGB(bmp, i, j, &r, &g, &b);
			BMP_SetPixelRGB(bmp, i, j, (UCHAR)(r /(255/level)/(1.0*level-1)*255), (UCHAR)(g / (255 / level) / (1.0*level - 1) * 255), (UCHAR)(b / (255 / level) / (1.0*level - 1) * 255));
			BMP_CHECK_ERROR(stdout, -1);
		}
	}
}

void BritandCntr(BMP* bmp, int brightness,double contrast) //Brightness: amount of RGB value you want to add
//Contrast: | 1 means original | <1 less contrast | >1 more contrast
{
	UINT    width, height;
	width = BMP_GetWidth(bmp);
	height = BMP_GetHeight(bmp);
	for (UINT i = 0; i < width; i++)
	{
		for (UINT j = 0; j < height; j++)
		{
			UCHAR Ur, Ug, Ub;
			BMP_GetPixelRGB(bmp, i, j, &Ur, &Ug, &Ub);
			int r, g, b;
			r = (int)(contrast*(Ur - 128)) + 128 + brightness;
			g = (int)(contrast*(Ug - 128)) + 128 + brightness;
			b = (int)(contrast*(Ub - 128)) + 128 + brightness;
			if (r < 0) r = 0;
			if (g < 0) g = 0;
			if (b < 0) b = 0;
			if (r > 255) r = 255;
			if (g > 255) g = 255;
			if (b > 255) b = 255;
			BMP_SetPixelRGB(bmp, i, j, (UCHAR)r, (UCHAR)g, (UCHAR)b);
			BMP_CHECK_ERROR(stdout, -1);
		}
	}
}

double GaussFunction(double theta, int dimension, double* gaussmatrix)
{
	double sum = 0;
	for (int i = 0; i < dimension; i++)
	{
		for (int j = 0; j < dimension; j++)
		{
			int x = i - (dimension - 1) / 2;
			int y = j - (dimension - 1) / 2;
			double gauss = 1.0 / (2 * PI*theta*theta)*pow(LOGE, 1.0*(-x*x - y*y) / (2 * theta*theta));
			gaussmatrix[i*dimension + j] = gauss;
			sum += gauss;
		}
	}
	return sum;
}

double maximum(double a, double b, double c, double d)
{
	double par = a;
	double var[3] = { b,c,d };
	for (int i = 0; i < 3; i++)
		par = (par > var[i]) ? a : var[i];
	return par;
}

double minimum(double a, double b, double c, double d)
{
	double par = a;
	double var[3] = { b,c,d };
	for (int i = 0; i < 3; i++)
		par = (par < var[i]) ? a : var[i];
	return par;
}

double r_compute(double var_x, double var_y, BMP* bmp, double cost, double sint)
{
	UCHAR R[4] = { 0,0,0,0 };
	UCHAR G[4] = { 0,0,0,0 };
	UCHAR B[4] = { 0,0,0,0 };

	BMP_GetPixelRGB(bmp, (unsigned long)var_x, (unsigned long)var_y, R, G, B);
	BMP_GetPixelRGB(bmp, (unsigned long)var_x+1, (unsigned long)var_y, R+1, G+1, B+1);
	BMP_GetPixelRGB(bmp, (unsigned long)var_x, (unsigned long)var_y+1, R+2, G+2, B+2);
	BMP_GetPixelRGB(bmp, (unsigned long)var_x+1, (unsigned long)var_y+1, R+3, G+3, B+3);
	
	//Area average of the colors
	return ((1-var_x+(unsigned long)var_x)*(1-var_y+(unsigned long)var_y)*R[0] + (var_x - (unsigned long)var_x)*(1 - var_y + (unsigned long)var_y)*R[1] + (1 - var_x + (unsigned long)var_x)*(var_y - (unsigned long)var_y)*R[2] + (var_x - (unsigned long)var_x)*(var_y - (unsigned long)var_y)*R[3]);
}

double g_compute(double var_x, double var_y, BMP* bmp, double cost, double sint)
{
	UCHAR R[4] = { 0,0,0,0 };
	UCHAR G[4] = { 0,0,0,0 };
	UCHAR B[4] = { 0,0,0,0 };

	BMP_GetPixelRGB(bmp, (unsigned long)var_x, (unsigned long)var_y, R, G, B);
	BMP_GetPixelRGB(bmp, (unsigned long)var_x + 1, (unsigned long)var_y, R + 1, G + 1, B + 1);
	BMP_GetPixelRGB(bmp, (unsigned long)var_x, (unsigned long)var_y + 1, R + 2, G + 2, B + 2);
	BMP_GetPixelRGB(bmp, (unsigned long)var_x + 1, (unsigned long)var_y + 1, R + 3, G + 3, B + 3);

	//Area average of the colors
	return ((1 - var_x + (unsigned long)var_x)*(1 - var_y + (unsigned long)var_y)*G[0] + (var_x - (unsigned long)var_x)*(1 - var_y + (unsigned long)var_y)*G[1] + (1 - var_x + (unsigned long)var_x)*(var_y - (unsigned long)var_y)*G[2] + (var_x - (unsigned long)var_x)*(var_y - (unsigned long)var_y)*G[3]);
}

double b_compute(double var_x, double var_y, BMP* bmp, double cost, double sint)
{
	UCHAR R[4] = { 0,0,0,0 };
	UCHAR G[4] = { 0,0,0,0 };
	UCHAR B[4] = { 0,0,0,0 };

	BMP_GetPixelRGB(bmp, (unsigned long)var_x, (unsigned long)var_y, R, G, B);
	BMP_GetPixelRGB(bmp, (unsigned long)var_x + 1, (unsigned long)var_y, R + 1, G + 1, B + 1);
	BMP_GetPixelRGB(bmp, (unsigned long)var_x, (unsigned long)var_y + 1, R + 2, G + 2, B + 2);
	BMP_GetPixelRGB(bmp, (unsigned long)var_x + 1, (unsigned long)var_y + 1, R + 3, G + 3, B + 3);

	//Area average of the colors
	return ((1 - var_x + (unsigned long)var_x)*(1 - var_y + (unsigned long)var_y)*B[0] + (var_x - (unsigned long)var_x)*(1 - var_y + (unsigned long)var_y)*B[1] + (1 - var_x + (unsigned long)var_x)*(var_y - (unsigned long)var_y)*B[2] + (var_x - (unsigned long)var_x)*(var_y - (unsigned long)var_y)*B[3]);
}

BMP* AMRotation(BMP* bmp, double theta)
{
	//Compute some important coefficient
	double cost = cos(theta);
	double sint = sin(theta);
	double pre_width = BMP_GetWidth(bmp);
	double pre_height = BMP_GetHeight(bmp);
	double depth = BMP_GetDepth(bmp);

	//Analysing size of new picture;
	double max_x = maximum((pre_width*cost - pre_height*sint), -pre_height*sint, 0.1 * 0, pre_width*cost);
	double min_x = minimum((pre_width*cost - pre_height*sint), -pre_height*sint, 0.1 * 0, pre_width*cost);
	double max_y = maximum((pre_width*sint + pre_height*cost), pre_height*cost, 0.1 * 0, pre_width*sint);
	double min_y = minimum((pre_width*sint + pre_height*cost), pre_height*cost, 0.1 * 0, pre_width*sint);
	double width = max_x - min_x;
	double height = max_y - min_y;

	BMP* newbmp = BMP_Create(width, height, depth);

	//x`=Ax + (-minx,-miny);
	for (UINT x = 0; x < width; x++)
	{
		for (UINT y = 0; y < height; y++)
		{
			double r = 0, g = 0, b = 0;

			//Compute original coordinate in newbmp system
			double pre_x = cost*(x + min_x) + sint*(y + min_y);
			double pre_y = -sint*(x + min_x) + cost*(y + min_y);

			if ((pre_x > 0 && pre_x < pre_width) && (pre_y > 0 && pre_y < pre_height))
			{
				r = r_compute(pre_x, pre_y, bmp, cost, sint);
				g = g_compute(pre_x, pre_y, bmp, cost, sint);
				b = b_compute(pre_x, pre_y, bmp, cost, sint);
			}
			BMP_SetPixelRGB(newbmp, x, y, (UCHAR)r, (UCHAR)g, (UCHAR)b);
		}
	}

	return newbmp;
}

BMP* GaussianBlur(BMP* bmp,double theta,int radius)
{
	double* gaussmatrix = malloc(sizeof(double)*radius*radius);
	double sum=GaussFunction(theta, radius, gaussmatrix);
#ifdef DEBUG
	for (int i = 0; i < radius; i++)
	{
		for (int j = 0; j < radius; j++)
		{
			printf("%10.4lf", gaussmatrix[i*radius + j]);
		}
		printf("\n");
	}
#endif // DEBUG

	UINT    width, height;
	USHORT depth;
	width = BMP_GetWidth(bmp);
	height = BMP_GetHeight(bmp);
	depth = BMP_GetDepth(bmp);
	BMP* newbmp=BMP_Create(width, height, depth);
	//Go through all pixels
	for (UINT x = 0; x < width; x++)
	{
		for (UINT y = 0; y < height; y++)
		{
			//Fill in colors
			double r = 0, g = 0, b = 0;
			for (int i = 0; i < radius; i++)
			{
				for (int j = 0; j < radius; j++)
				{
					//Reflect boundary
					int xi=x, yi=y;
					if (xi - (radius - 1) / 2 + i < 0)
					{
						xi = -(xi - (radius - 1) / 2 + i);
					}
					else
					{
						xi = (xi - (radius - 1) / 2 + i);
					}
					if (yi - (radius - 1) / 2 + j < 0)
					{
						yi = -(yi - (radius - 1) / 2 + j);
					}
					else
					{
						yi = (yi - (radius - 1) / 2 + j);
					}
					if (xi >= width)
					{
						xi = width - (xi - width) - 1;
					}
					if (yi >= height)
					{
						yi = height - (yi - height) - 1;
					}
					UCHAR R, G, B;
					BMP_GetPixelRGB(bmp, xi, yi, &R, &G, &B);
					r += R*gaussmatrix[i*radius + j];
					g += G*gaussmatrix[i*radius + j];
					b += B*gaussmatrix[i*radius + j];
				}
			}

			//Calculate total
			r /= sum;
			g /= sum;
			b /= sum;
			BMP_SetPixelRGB(newbmp, x, y, (UCHAR)r, (UCHAR)g, (UCHAR)b);
		}
	}
	free(gaussmatrix);
	return newbmp;
}