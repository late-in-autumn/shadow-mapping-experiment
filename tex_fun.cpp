/* Texture functions for cs580 GzLib	*/
#include    "stdafx.h" 
#include	"stdio.h"
#include	"Gz.h"
#include	"rend.h"
#include	<limits>

GzColor* image = nullptr;
GzColor randomColors[NUM_QUADRANTS]{};
int xs, ys;
int reset = 1;

int texture_value(int x, int y, GzColor out)
{
	if (image == nullptr) return GZ_FAILURE;
	else
	{
		if (x < 0) x = 0;
		if (x > xs - 1) x = xs - 1;
		if (y < 0) y = 0;
		if (y > ys - 1) y = ys - 1;
		memcpy(out, image + MATRIX_OFFSET(xs, x, y), sizeof(GzColor));
		return GZ_SUCCESS;
	}
}

int interpolate_color(float s, float t, GzColor corners[4], GzColor out)
{
	out[RED] = s * t * corners[3][RED] + (1 - s) * t * corners[2][RED]
		+ s * (1 - t) * corners[1][RED] + (1 - s) * (1 - t) * corners[0][RED];
	out[GREEN] = s * t * corners[3][GREEN] + (1 - s) * t * corners[2][GREEN]
		+ s * (1 - t) * corners[1][GREEN] + (1 - s) * (1 - t) * corners[0][GREEN];
	out[BLUE] = s * t * corners[3][BLUE] + (1 - s) * t * corners[2][BLUE]
		+ s * (1 - t) * corners[1][BLUE] + (1 - s) * (1 - t) * corners[0][BLUE];
	return GZ_SUCCESS;
}

/* Image texture function */
int tex_fun(float u, float v, GzColor color)
{
	unsigned char		pixel[3];
	unsigned char     dummy;
	char  		foo[8];
	int   		i;
	FILE* fd;

	if (reset)
	{          /* open and load texture file */
		fd = fopen("texture", "rb");
		if (fd == nullptr)
		{
#ifdef GZ_DEBUG_MODE
			OutputDebugStringA("texture file not found\n");
#endif
			exit(-1);
		}
		fscanf(fd, "%s %d %d %c", foo, &xs, &ys, &dummy);
		image = (GzColor*)malloc(sizeof(GzColor) * (xs + 1) * (ys + 1));
		if (image == nullptr)
		{
#ifdef GZ_DEBUG_MODE
			OutputDebugStringA("malloc for texture image failed\n");
#endif
			exit(-1);
		}

		for (i = 0; i < xs * ys; i++) {	/* create array of GzColor values */
			fread(pixel, sizeof(pixel), 1, fd);
			image[i][RED] = (float)((int)pixel[RED]) * (static_cast<float>(1) / static_cast<float>(255));
			image[i][GREEN] = (float)((int)pixel[GREEN]) * (static_cast<float>(1) / static_cast<float>(255));
			image[i][BLUE] = (float)((int)pixel[BLUE]) * (static_cast<float>(1) / static_cast<float>(255));
		}

		reset = 0;          /* init is done */
		fclose(fd);
	}

	/* bounds-test u,v to make sure nothing will overflow image array bounds */
	/* determine texture cell corner values and perform bilinear interpolation */
	/* set color to interpolated GzColor value and return */
	u = CLAMP(u, 0.0f, 1.0f);
	v = CLAMP(v, 0.0f, 1.0f);

	float scaledU = u * (xs - 1);
	float scaledV = v * (ys - 1);
	int floorU = (int)floorf(scaledU);
	int floorV = (int)floorf(scaledV);
	GzColor corners[NUM_QUADRANTS]{};

	texture_value(floorU, floorV, corners[0]);
	texture_value(floorU + 1, floorV, corners[1]);
	texture_value(floorU, floorV + 1, corners[2]);
	texture_value(floorU + 1, floorV + 1, corners[3]);
	interpolate_color(scaledU - floorU, scaledV - floorV, corners, color);

	return GZ_SUCCESS;
}

/* Procedural texture function */
int ptex_fun(float u, float v, GzColor color)
{
	if (reset) // generate four random colors for the four quadrants
	{
#ifdef GZ_DEBUG_MODE
		OutputDebugStringA("Initializing procedural texture with four random colors...\n");
#endif
		for (int i = 0; i < NUM_QUADRANTS; i++)
		{
			randomColors[i][RED] = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
			randomColors[i][GREEN] = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
			randomColors[i][BLUE] = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
		}
		reset = 0;
	}

	u = CLAMP(u, 0.0f, 1.0f);
	v = CLAMP(v, 0.0f, 1.0f);

	if (FP_LESS(v, 0.5f, std::numeric_limits<float>::epsilon())
		&& FP_LESS(u, 0.5f, std::numeric_limits<float>::epsilon())) memcpy(color, randomColors[0], sizeof(GzColor));
	else if (FP_LESS(v, 0.5f, std::numeric_limits<float>::epsilon())
		&& !FP_LESS(u, 0.5f, std::numeric_limits<float>::epsilon())) memcpy(color, randomColors[1], sizeof(GzColor));
	else if (!FP_LESS(v, 0.5f, std::numeric_limits<float>::epsilon())
		&& FP_LESS(u, 0.5f, std::numeric_limits<float>::epsilon())) memcpy(color, randomColors[2], sizeof(GzColor));
	else memcpy(color, randomColors[3], sizeof(GzColor));

	return GZ_SUCCESS;
}

/* Free texture memory */
int GzFreeTexture()
{
	if (image != nullptr)
		free(image);
	return GZ_SUCCESS;
}
