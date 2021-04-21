/* CS580 Homework 3 */

#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"
#include	<algorithm>
#include	<fstream>
#include	<limits>
#include	<new>
#include	<string>

int GzRender::GzRotXMat(float degree, GzMatrix mat)
{
	memset(mat, 0, sizeof(GzMatrix));
	mat[0][0] = 1;
	mat[1][1] = cosf(TO_RADIAN(degree));
	mat[1][2] = -sinf(TO_RADIAN(degree));
	mat[2][1] = sinf(TO_RADIAN(degree));
	mat[2][2] = cosf(TO_RADIAN(degree));
	mat[3][3] = 1;

	return GZ_SUCCESS;
}

int GzRender::GzRotYMat(float degree, GzMatrix mat)
{
	memset(mat, 0, sizeof(GzMatrix));
	mat[1][1] = 1;
	mat[0][0] = cosf(TO_RADIAN(degree));
	mat[0][2] = sinf(TO_RADIAN(degree));
	mat[2][0] = -sinf(TO_RADIAN(degree));
	mat[2][2] = cosf(TO_RADIAN(degree));
	mat[3][3] = 1;

	return GZ_SUCCESS;
}

int GzRender::GzRotZMat(float degree, GzMatrix mat)
{
	memset(mat, 0, sizeof(GzMatrix));
	mat[2][2] = 1;
	mat[0][0] = cosf(TO_RADIAN(degree));
	mat[0][1] = -sinf(TO_RADIAN(degree));
	mat[1][0] = sinf(TO_RADIAN(degree));
	mat[1][1] = cosf(TO_RADIAN(degree));
	mat[3][3] = 1;

	return GZ_SUCCESS;
}

int GzRender::GzTrxMat(GzCoord translate, GzMatrix mat)
{
	memset(mat, 0, sizeof(GzMatrix));
	mat[0][0] = 1;
	mat[1][1] = 1;
	mat[2][2] = 1;
	mat[3][3] = 1;
	mat[0][3] = translate[0];
	mat[1][3] = translate[1];
	mat[2][3] = translate[2];

	return GZ_SUCCESS;
}

int GzRender::GzScaleMat(GzCoord scale, GzMatrix mat)
{
	memset(mat, 0, sizeof(GzMatrix));
	mat[0][0] = scale[0];
	mat[1][1] = scale[1];
	mat[2][2] = scale[2];
	mat[3][3] = 1;

	return GZ_SUCCESS;
}

int GzRender::GzSetShadowMap(bool in)
{
	use_shadow_map = in;
	return GZ_SUCCESS;
}

GzRender::GzRender(int xRes, int yRes)
{
	// default to not using shadowmap
	use_shadow_map = false;

	// update x and y resolutions
	xres = xRes;
	yres = yRes;

	// both the pixel and frame buffers are stored linearly
	// allocate them with the correct size
	pixelbuffer = (GzPixel*)calloc(xres * yres, sizeof(GzPixel));
	framebuffer = (char*)calloc(3 * xres * yres, sizeof(char));

	// in case of memory allocation failures
	if (pixelbuffer == nullptr || framebuffer == nullptr) throw std::bad_alloc();

	// initialize number of lights
	numlights = 0;

	memset(Xsp, 0, sizeof(GzMatrix));
	Xsp[0][0] = xres / static_cast<float>(2);
	Xsp[1][1] = -yres / static_cast<float>(2);
	Xsp[2][2] = static_cast<float>((std::numeric_limits<int>::max)());
	Xsp[3][3] = 1;
	Xsp[0][3] = xres / static_cast<float>(2);
	Xsp[1][3] = yres / static_cast<float>(2);

	GzDefaultCamera();
}

GzRender::~GzRender()
{
	// clean up the pixel and frame buffers
	free(pixelbuffer);
	free(framebuffer);
}

int GzRender::GzDefault()
{
	// initialize current triangle and current offset
	memset(currentTriangle, 0, sizeof(GzTriangle));
	memset(currentOffset, 0, sizeof(GzOffsetsXY));

	// initialize everything to white and max z distance
	for (int offset = 0; offset < xres * yres; offset++)
	{
		(pixelbuffer + offset)->red = 4095;
		(pixelbuffer + offset)->green = 4095;
		(pixelbuffer + offset)->blue = 4095;
		(pixelbuffer + offset)->alpha = 4095;
		(pixelbuffer + offset)->z = (std::numeric_limits<int>::max)();
	}

	// update the framebuffer accordingly
	GzFlushDisplay2FrameBuffer();

	return GZ_SUCCESS;
}

int GzRender::GzBeginRender()
{
	GzMatrix Riw{};

	ExtractRotation(m_camera.Xiw, Riw);
	GzDefault();

	memcpy(Ximage[0], Xsp, sizeof(GzMatrix));
	MultiplyMatrices(Ximage[0], m_camera.Xpi, Ximage[1]);
	MultiplyMatrices(Ximage[1], m_camera.Xiw, Ximage[2]);

	memcpy(Xnorm[0], IDENTITY, sizeof(GzMatrix));
	memcpy(Xnorm[1], IDENTITY, sizeof(GzMatrix));
	memcpy(Xnorm[2], Riw, sizeof(GzMatrix));

	matlevel = 2;

	return GZ_SUCCESS;
}

int GzRender::GzDefaultCamera()
{
	m_camera.position[0] = DEFAULT_IM_X;
	m_camera.position[1] = DEFAULT_IM_Y;
	m_camera.position[2] = DEFAULT_IM_Z;
	m_camera.lookat[0] = 0;
	m_camera.lookat[1] = 0;
	m_camera.lookat[2] = 0;
	m_camera.worldup[0] = 0;
	m_camera.worldup[1] = 1.0;
	m_camera.worldup[2] = 0.0;
	m_camera.FOV = DEFAULT_FOV;
	UpdateCameraMatrices();

	return GZ_SUCCESS;
}

int GzRender::GzPutCamera(GzCamera camera)
{
	m_camera = camera;
	UpdateCameraMatrices();

	return GZ_SUCCESS;
}

int GzRender::GzPushMatrix(GzMatrix	matrix)
{
	if (matlevel + 1 >= MATLEVELS) return GZ_FAILURE;

	GzMatrix rotation;

	ExtractRotation(matrix, rotation);

	if (matlevel == -1)
	{
		matlevel++;
		memcpy(Ximage[matlevel], matrix, sizeof(GzMatrix));
		memcpy(Xnorm[matlevel], rotation, sizeof(GzMatrix));
	}
	else
	{
		matlevel++;
		MultiplyMatrices(Ximage[matlevel - 1], matrix, Ximage[matlevel]);
		MultiplyMatrices(Xnorm[matlevel - 1], rotation, Xnorm[matlevel]);
	}

	return GZ_SUCCESS;
}

int GzRender::GzPopMatrix()
{
	if (matlevel - 1 < -1) return GZ_FAILURE;

	matlevel--;

	return GZ_SUCCESS;
}

int GzRender::GzPut(int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
	// normalize x and y inputs
	int x, y;
	x = i > xres - 1 ? xres - 1 : (i < 0 ? 0 : i);
	y = j > yres - 1 ? yres - 1 : (j < 0 ? 0 : j);

	// set the pixel (with out-of-range values clipped)
	(pixelbuffer + MATRIX_OFFSET(xres, x, y))->red = r > 4095 ? 4095 : (r < 0 ? 0 : r);
	(pixelbuffer + MATRIX_OFFSET(xres, x, y))->green = g > 4095 ? 4095 : (g < 0 ? 0 : g);
	(pixelbuffer + MATRIX_OFFSET(xres, x, y))->blue = b > 4095 ? 4095 : (b < 0 ? 0 : b);
	(pixelbuffer + MATRIX_OFFSET(xres, x, y))->alpha = a > 4095 ? 4095 : (a < 0 ? 0 : a);
	(pixelbuffer + MATRIX_OFFSET(xres, x, y))->z = z;

	return GZ_SUCCESS;
}


int GzRender::GzGet(int i, int j, GzIntensity* r, GzIntensity* g, GzIntensity* b, GzIntensity* a, GzDepth* z)
{
	// normalize x and y inputs
	int x, y;
	x = i > xres - 1 ? xres - 1 : (i < 0 ? 0 : i);
	y = j > yres - 1 ? yres - 1 : (j < 0 ? 0 : j);

	// get the pixel values
	*r = (pixelbuffer + MATRIX_OFFSET(xres, x, y))->red;
	*g = (pixelbuffer + MATRIX_OFFSET(xres, x, y))->green;
	*b = (pixelbuffer + MATRIX_OFFSET(xres, x, y))->blue;
	*a = (pixelbuffer + MATRIX_OFFSET(xres, x, y))->alpha;
	*z = (pixelbuffer + MATRIX_OFFSET(xres, x, y))->z;

	return GZ_SUCCESS;
}


int GzRender::GzFlushDisplay2File(FILE* outfile)
{
	// create the output stream from file
	std::fstream output(outfile);

	// write PPM header (P3 for ASCII encoded PPM: http://paulbourke.net/dataformats/ppm/)
	output << "P3" << std::endl
		<< xres << " " << yres << " " << "255" << std::endl
		<< "# ASCII-encoded PPM output from the Gz library" << std::endl;

	// write PPM body (ASCII encoded PPM)
	for (int offset = 0; offset < xres * yres; offset++)
		output << ((pixelbuffer + offset)->red >> 4) << " "
		<< ((pixelbuffer + offset)->green >> 4) << " "
		<< ((pixelbuffer + offset)->blue >> 4) << std::endl;

	return GZ_SUCCESS;
}

int GzRender::GzFlushDisplay2FrameBuffer()
{
	// offsets for pixel and frame buffer elements
	int pixel_offset = 0;
	int frame_offset = 0;

	// copy contents of the pixel buffer to the frame buffer
	for (int y = 0; y < yres; y++)
		for (int x = 0; x < xres; x++)
		{
			*(framebuffer + frame_offset) = ((pixelbuffer + pixel_offset)->blue >> 4);
			*(framebuffer + frame_offset + 1) = ((pixelbuffer + pixel_offset)->green >> 4);
			*(framebuffer + frame_offset + 2) = ((pixelbuffer + pixel_offset)->red >> 4);

			pixel_offset++;
			frame_offset += 3;
		}

	return GZ_SUCCESS;
}

/***********************************************/
/* HW2 methods: implement from here */

int GzRender::GzPutAttribute(int numAttributes, GzToken* nameList, GzPointer* valueList)
{
	/* HW 2.1
	-- Set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
	-- In later homeworks set shaders, interpolaters, texture maps, and lights
	*/

	/*
	- GzPutAttribute() must accept the following tokens/values:

	- GZ_RGB_COLOR					GzColor		default flat-shade color
	- GZ_INTERPOLATE				int			shader interpolation mode
	- GZ_DIRECTIONAL_LIGHT			GzLight
	- GZ_AMBIENT_LIGHT            	GzLight		(ignore direction)
	- GZ_AMBIENT_COEFFICIENT		GzColor		Ka reflectance
	- GZ_DIFFUSE_COEFFICIENT		GzColor		Kd reflectance
	- GZ_SPECULAR_COEFFICIENT       GzColor		Ks coef's
	- GZ_DISTRIBUTION_COEFFICIENT   float		spec power
	*/
	for (int i = 0; i < numAttributes; i++)
	{
		switch (*(nameList + i))
		{
		case GZ_RGB_COLOR:
			memcpy(flatcolor, ((GzColor*)valueList[i]), sizeof(GzColor));
			break;
		case GZ_INTERPOLATE:
			interp_mode = *((int*)valueList[i]);
			break;
		case GZ_DIRECTIONAL_LIGHT:
			if (numlights < MAX_LIGHTS)
			{
				lights[numlights] = *((GzLight*)valueList[i]);
				numlights++;
			}
			break;
		case GZ_AMBIENT_LIGHT:
			ambientlight = *((GzLight*)valueList[i]);
			break;
		case GZ_AMBIENT_COEFFICIENT:
			memcpy(Ka, ((GzColor*)valueList[i]), sizeof(GzColor));
			break;
		case GZ_DIFFUSE_COEFFICIENT:
			memcpy(Kd, ((GzColor*)valueList[i]), sizeof(GzColor));
			break;
		case GZ_SPECULAR_COEFFICIENT:
			memcpy(Ks, ((GzColor*)valueList[i]), sizeof(GzColor));
			break;
		case GZ_DISTRIBUTION_COEFFICIENT:
			spec = *((float*)valueList[i]);
			break;
		case GZ_TEXTURE_MAP:
			tex_fun = (GzTexture)valueList[i];
			break;
		case GZ_AASHIFTX:
			currentOffset[0] = *((float*)valueList[i]);
			break;
		case GZ_AASHIFTY:
			currentOffset[1] = *((float*)valueList[i]);
			break;
		default:
			break;
		}
	}

	return GZ_SUCCESS;
}

int GzRender::GzPutTriangle(int numParts, GzToken* nameList, GzPointer* valueList)
/* numParts - how many names and values */
{
	/* HW 2.2
	-- Pass in a triangle description with tokens and values corresponding to
		  GZ_NULL_TOKEN:		do nothing - no values
		  GZ_POSITION:		3 vert positions in model space
	-- Return error code
	*/
	/*
	-- Xform positions of verts using matrix on top of stack
	-- Clip - just discard any triangle with any vert(s) behind view plane
			- optional: test for triangles with all three verts off-screen (trivial frustum cull)
	-- invoke triangle rasterizer
	*/
	for (int i = 0; i < numParts; i++)
	{
		switch (*(nameList + i))
		{
		case GZ_POSITION:
			memcpy(currentTriangle[0].coord, ((GzCoord*)valueList[i]), sizeof(GzCoord));
			memcpy(currentTriangle[1].coord, ((GzCoord*)valueList[i]) + 1, sizeof(GzCoord));
			memcpy(currentTriangle[2].coord, ((GzCoord*)valueList[i]) + 2, sizeof(GzCoord));
			break;
		case GZ_NORMAL:
			memcpy(currentTriangle[0].normal, ((GzCoord*)valueList[i]), sizeof(GzCoord));
			memcpy(currentTriangle[1].normal, ((GzCoord*)valueList[i]) + 1, sizeof(GzCoord));
			memcpy(currentTriangle[2].normal, ((GzCoord*)valueList[i]) + 2, sizeof(GzCoord));
			break;
		case GZ_TEXTURE_INDEX:
			memcpy(currentTriangle[0].uv, ((GzTextureIndex*)valueList[i]), sizeof(GzTextureIndex));
			memcpy(currentTriangle[1].uv, ((GzTextureIndex*)valueList[i]) + 1, sizeof(GzTextureIndex));
			memcpy(currentTriangle[2].uv, ((GzTextureIndex*)valueList[i]) + 2, sizeof(GzTextureIndex));
			RenderTriangle(&(currentTriangle[0]), &(currentTriangle[1]), &(currentTriangle[2]));
			break;
		case GZ_NULL_TOKEN:
			break;
		default:
			break;
		}
	}

	return GZ_SUCCESS;
}

void GzRender::UpdateCameraMatrices()
{
	GzCoord zAxis{};
	GzCoord yAxis{};
	GzCoord xAxis{};
	float divideFactor, upDot, xDotOrigin, yDotOrigin, zDotOrigin;

	divideFactor = tanf(TO_RADIAN(m_camera.FOV / 2));
	zAxis[0] = m_camera.lookat[0] - m_camera.position[0];
	zAxis[1] = m_camera.lookat[1] - m_camera.position[1];
	zAxis[2] = m_camera.lookat[2] - m_camera.position[2];
	NormalizeVector(3, zAxis);
	upDot = DotProduct(3, m_camera.worldup, zAxis);
	yAxis[0] = m_camera.worldup[0] - upDot * zAxis[0];
	yAxis[1] = m_camera.worldup[1] - upDot * zAxis[1];
	yAxis[2] = m_camera.worldup[2] - upDot * zAxis[2];
	NormalizeVector(3, yAxis);
	CrossProduct(yAxis, zAxis, xAxis);
	NormalizeVector(3, xAxis);
	xDotOrigin = DotProduct(3, xAxis, m_camera.position);
	yDotOrigin = DotProduct(3, yAxis, m_camera.position);
	zDotOrigin = DotProduct(3, zAxis, m_camera.position);

	memset(m_camera.Xpi, 0, sizeof(GzMatrix));
	m_camera.Xpi[0][0] = 1;
	m_camera.Xpi[1][1] = 1;
	m_camera.Xpi[2][2] = divideFactor;
	m_camera.Xpi[3][2] = divideFactor;
	m_camera.Xpi[3][3] = 1;

	memset(m_camera.Xiw, 0, sizeof(GzMatrix));
	m_camera.Xiw[0][0] = xAxis[0];
	m_camera.Xiw[0][1] = xAxis[1];
	m_camera.Xiw[0][2] = xAxis[2];
	m_camera.Xiw[1][0] = yAxis[0];
	m_camera.Xiw[1][1] = yAxis[1];
	m_camera.Xiw[1][2] = yAxis[2];
	m_camera.Xiw[2][0] = zAxis[0];
	m_camera.Xiw[2][1] = zAxis[1];
	m_camera.Xiw[2][2] = zAxis[2];
	m_camera.Xiw[0][3] = -xDotOrigin;
	m_camera.Xiw[1][3] = -yDotOrigin;
	m_camera.Xiw[2][3] = -zDotOrigin;
	m_camera.Xiw[3][3] = 1;

	memset(m_camera.Xwi, 0, sizeof(GzMatrix));
	m_camera.Xwi[0][0] = xAxis[0];
	m_camera.Xwi[1][0] = xAxis[1];
	m_camera.Xwi[2][0] = xAxis[2];
	m_camera.Xwi[0][1] = yAxis[0];
	m_camera.Xwi[1][1] = yAxis[1];
	m_camera.Xwi[2][1] = yAxis[2];
	m_camera.Xwi[0][2] = zAxis[0];
	m_camera.Xwi[1][2] = zAxis[1];
	m_camera.Xwi[2][2] = zAxis[2];
	m_camera.Xwi[0][3] = xDotOrigin;
	m_camera.Xwi[1][3] = yDotOrigin;
	m_camera.Xwi[2][3] = zDotOrigin;
	m_camera.Xwi[3][3] = 1;
}

void GzRender::RenderTriangle(GzVertex* v1, GzVertex* v2, GzVertex* v3)
{
	// transform the three coordinates here
	GzVertex t1{}, t2{}, t3{};

	t1 = *v1;
	t2 = *v2;
	t3 = *v3;

	TranslateVertex(&t1, Ximage[matlevel], Xnorm[matlevel]);
	TranslateVertex(&t2, Ximage[matlevel], Xnorm[matlevel]);
	TranslateVertex(&t3, Ximage[matlevel], Xnorm[matlevel]);

	if (FP_LESS(t1.coord[2], 0.0f, std::numeric_limits<float>::epsilon())
		|| FP_LESS(t2.coord[2], 0.0f, std::numeric_limits<float>::epsilon())
		|| FP_LESS(t3.coord[2], 0.0f, std::numeric_limits<float>::epsilon())) return;

	LinearEvaulator(&t1, &t2, &t3);
}

void GzRender::LinearEvaulator(GzVertex* v1, GzVertex* v2, GzVertex* v3)
{
	// step 0: set up the vertices, edges, and the planes needed for computation
	GzVertex va{}, vb{}, vc{};
	GzEdge ac{}, e1{}, e2{}, e3{};
	GzPlane plane{}, pr{}, pg{}, pb{};
	GzNormalPlane pn{}, pshadow{};
	GzUvPlane puv{};
	long start[2]{}, end[2]{};

	// step 1: populate the vertices
	va = *v1;
	vb = *v2;
	vc = *v3;

	// step 2: compute the bound box
	ComputeBoundBox(&va, &vb, &vc, start, end);

	// step 3: sort the three vertices by y
	SortVertices(&va, &vb, &vc);

	// step 4: compute the colors of the three vertices
	if (tex_fun == nullptr)
	{
		ComputeColor(&va, Ka, Kd, Ks, const_cast<GzMatrix&>(IDENTITY));
		ComputeColor(&vb, Ka, Kd, Ks, const_cast<GzMatrix&>(IDENTITY));
		ComputeColor(&vc, Ka, Kd, Ks, const_cast<GzMatrix&>(IDENTITY));
	}
	else
	{
		ComputeColor(&va, const_cast<GzColor&>(FULL_INTENSITY), const_cast<GzColor&>(FULL_INTENSITY),
			const_cast<GzColor&>(FULL_INTENSITY), const_cast<GzMatrix&>(IDENTITY));
		ComputeColor(&vb, const_cast<GzColor&>(FULL_INTENSITY), const_cast<GzColor&>(FULL_INTENSITY),
			const_cast<GzColor&>(FULL_INTENSITY), const_cast<GzMatrix&>(IDENTITY));
		ComputeColor(&vc, const_cast<GzColor&>(FULL_INTENSITY), const_cast<GzColor&>(FULL_INTENSITY),
			const_cast<GzColor&>(FULL_INTENSITY), const_cast<GzMatrix&>(IDENTITY));
	}

	if (FP_EQUALS(va.coord[1], vb.coord[1], std::numeric_limits<float>::epsilon()))
	{
		// ab
		e1.start = va;
		e1.end = vb;

		// bc
		e2.start = vb;
		e2.end = vc;

		// ca
		e3.start = vc;
		e3.end = va;
	}
	else if (FP_EQUALS(vb.coord[1], vc.coord[1], std::numeric_limits<float>::epsilon()))
	{
		// ac
		e1.start = va;
		e1.end = vc;

		// cb
		e2.start = vc;
		e2.end = vb;

		// ba
		e3.start = vb;
		e3.end = va;
	}
	else
	{
		// step 5: set up edge ac needed for left/right edge detection
		ac.start = va;
		ac.end = vc;
		ComputeEdgeEquasion(&ac);

		// step 6: set up the three edges (clockwise) according to the x coordinate of the point on ac with the y coordinate of b
		if (FP_LESS(-(ac.c + ac.b * vb.coord[1]) / ac.a, vb.coord[1], std::numeric_limits<float>::epsilon())) // ac is the left edge
		{
			// ab
			e1.start = va;
			e1.end = vb;

			// bc
			e2.start = vb;
			e2.end = vc;

			// ca
			e3.start = vc;
			e3.end = va;
		}
		else // ac is the right edge
		{
			// ac
			e1.start = va;
			e1.end = vc;

			// cb
			e2.start = vc;
			e2.end = vb;

			// ba
			e3.start = vb;
			e3.end = va;
		}
	}

	// step 7: compute the a, b, c for  e1, e2, and e3
	ComputeEdgeEquasion(&e1);
	ComputeEdgeEquasion(&e2);
	ComputeEdgeEquasion(&e3);

	// step 8: compute the planes the current triangle lies in (coordinates, normals, and color channels)
	plane.p1 = e1.start;
	plane.p2 = e2.start;
	plane.p3 = e3.start;

	pr.p1 = e1.start;
	pr.p2 = e2.start;
	pr.p3 = e3.start;

	pg.p1 = e1.start;
	pg.p2 = e2.start;
	pg.p3 = e3.start;

	pb.p1 = e1.start;
	pb.p2 = e2.start;
	pb.p3 = e3.start;

	pn.p1 = e1.start;
	pn.p2 = e2.start;
	pn.p3 = e3.start;

	puv.p1 = e1.start;
	puv.p2 = e2.start;
	puv.p3 = e3.start;

	ComputePlaneEquasion(&plane);
	ComputeColorPlaneEquasions(&pr, &pg, &pb);
	ComputeNormalPlaneEquasion(&pn);
	ComputeUvPlaneEquasion(&puv);

	// step 9: fill the bound box
	FillBoundBox(&e1, &e2, &e3, &plane, &pr, &pg, &pb, &pn, &puv, start, end);
}

void GzRender::SortVertices(GzVertex* v1, GzVertex* v2, GzVertex* v3)
{
	Point coords[3]{};

	memcpy(&(coords[0]), v1, sizeof(Point));
	memcpy(&(coords[1]), v2, sizeof(Point));
	memcpy(&(coords[2]), v3, sizeof(Point));

	std::qsort(coords, 3, sizeof(Point), CompareVertices);

	memcpy(v1, &(coords[0]), sizeof(Point));
	memcpy(v2, &(coords[1]), sizeof(Point));
	memcpy(v3, &(coords[2]), sizeof(Point));
}

void GzRender::ComputeBoundBox(GzVertex* v1, GzVertex* v2, GzVertex* v3, long start[2], long end[2])
{
	float x[3]{}, y[3]{};

	x[0] = v1->coord[0];
	x[1] = v2->coord[0];
	x[2] = v3->coord[0];

	y[0] = v1->coord[1];
	y[1] = v2->coord[1];
	y[2] = v3->coord[1];

	std::sort(x, x + 3);
	std::sort(y, y + 3);

	start[0] = static_cast<long>(floorl(x[0]) < 0 ? 0 : floorl(x[0]));
	start[1] = static_cast<long>(floorl(y[0]) < 0 ? 0 : floorl(y[0]));

	end[0] = static_cast<long>(ceill(x[2]) > xres ? xres : ceill(x[2]));
	end[1] = static_cast<long>(ceill(y[2]) > yres ? yres : ceill(y[2]));
}

void GzRender::ComputeEdgeEquasion(GzEdge* e)
{
	e->a = e->end.coord[1] - e->start.coord[1];
	e->b = -(e->end.coord[0] - e->start.coord[0]);
	e->c = (-e->b * e->start.coord[1]) - (e->a * e->start.coord[0]);
}

void GzRender::ComputePlaneEquasion(GzPlane* p)
{
	GzCoord v1{};
	GzCoord v2{};
	GzCoord product{};

	v1[0] = p->p2.coord[0] - p->p1.coord[0];
	v1[1] = p->p2.coord[1] - p->p1.coord[1];
	v1[2] = p->p2.coord[2] - p->p1.coord[2];

	v2[0] = p->p3.coord[0] - p->p2.coord[0];
	v2[1] = p->p3.coord[1] - p->p2.coord[1];
	v2[2] = p->p3.coord[2] - p->p2.coord[2];

	CrossProduct(v1, v2, product);

	p->a = product[0];
	p->b = product[1];
	p->c = product[2];
	p->d = -(p->c * p->p2.coord[2]) - (p->b * p->p2.coord[1]) - (p->a * p->p2.coord[0]);
}

void GzRender::ComputeNormalPlaneEquasion(GzNormalPlane* p)
{
	GzCoord a{}, b{}, product{};

	// plane equasion for Nx
	a[0] = p->p2.coord[0] - p->p1.coord[0];
	a[1] = p->p2.coord[1] - p->p1.coord[1];
	a[2] = p->p2.normal[0] - p->p1.normal[0];

	b[0] = p->p3.coord[0] - p->p2.coord[0];
	b[1] = p->p3.coord[1] - p->p2.coord[1];
	b[2] = p->p3.normal[0] - p->p2.normal[0];

	CrossProduct(a, b, product);

	p->ax = product[0];
	p->bx = product[1];
	p->cx = product[2];
	p->dx = -(p->cx * p->p2.normal[0]) - (p->bx * p->p2.coord[1]) - (p->ax * p->p2.coord[0]);

	// plane equasion for Ny
	a[0] = p->p2.coord[0] - p->p1.coord[0];
	a[1] = p->p2.coord[1] - p->p1.coord[1];
	a[2] = p->p2.normal[1] - p->p1.normal[1];

	b[0] = p->p3.coord[0] - p->p2.coord[0];
	b[1] = p->p3.coord[1] - p->p2.coord[1];
	b[2] = p->p3.normal[1] - p->p2.normal[1];

	CrossProduct(a, b, product);

	p->ay = product[0];
	p->by = product[1];
	p->cy = product[2];
	p->dy = -(p->cy * p->p2.normal[1]) - (p->by * p->p2.coord[1]) - (p->ay * p->p2.coord[0]);

	// plane equasion for Nz
	a[0] = p->p2.coord[0] - p->p1.coord[0];
	a[1] = p->p2.coord[1] - p->p1.coord[1];
	a[2] = p->p2.normal[2] - p->p1.normal[2];

	b[0] = p->p3.coord[0] - p->p2.coord[0];
	b[1] = p->p3.coord[1] - p->p2.coord[1];
	b[2] = p->p3.normal[2] - p->p2.normal[2];

	CrossProduct(a, b, product);

	p->az = product[0];
	p->bz = product[1];
	p->cz = product[2];
	p->dz = -(p->cz * p->p2.normal[2]) - (p->bz * p->p2.coord[1]) - (p->az * p->p2.coord[0]);
}

void GzRender::ComputeColorPlaneEquasions(GzPlane* pr, GzPlane* pg, GzPlane* pb)
{
	GzCoord a{}, b{}, product{};

	// plane equasion for red
	a[0] = pr->p2.coord[0] - pr->p1.coord[0];
	a[1] = pr->p2.coord[1] - pr->p1.coord[1];
	a[2] = pr->p2.color[0] - pr->p1.color[0];

	b[0] = pr->p3.coord[0] - pr->p2.coord[0];
	b[1] = pr->p3.coord[1] - pr->p2.coord[1];
	b[2] = pr->p3.color[0] - pr->p2.color[0];

	CrossProduct(a, b, product);

	pr->a = product[0];
	pr->b = product[1];
	pr->c = product[2];
	pr->d = -(pr->c * pr->p2.color[0]) - (pr->b * pr->p2.coord[1]) - (pr->a * pr->p2.coord[0]);

	// plane equasion for green
	a[0] = pg->p2.coord[0] - pg->p1.coord[0];
	a[1] = pg->p2.coord[1] - pg->p1.coord[1];
	a[2] = pg->p2.color[1] - pg->p1.color[1];

	b[0] = pg->p3.coord[0] - pg->p2.coord[0];
	b[1] = pg->p3.coord[1] - pg->p2.coord[1];
	b[2] = pg->p3.color[1] - pg->p2.color[1];

	CrossProduct(a, b, product);

	pg->a = product[0];
	pg->b = product[1];
	pg->c = product[2];
	pg->d = -(pg->c * pg->p2.color[1]) - (pg->b * pg->p2.coord[1]) - (pg->a * pg->p2.coord[0]);

	// plane equasion for blue
	a[0] = pb->p2.coord[0] - pb->p1.coord[0];
	a[1] = pb->p2.coord[1] - pb->p1.coord[1];
	a[2] = pb->p2.color[2] - pb->p1.color[2];

	b[0] = pb->p3.coord[0] - pb->p2.coord[0];
	b[1] = pb->p3.coord[1] - pb->p2.coord[1];
	b[2] = pb->p3.color[2] - pb->p2.color[2];

	CrossProduct(a, b, product);

	pb->a = product[0];
	pb->b = product[1];
	pb->c = product[2];
	pb->d = -(pb->c * pb->p2.color[2]) - (pb->b * pb->p2.coord[1]) - (pb->a * pb->p2.coord[0]);
}

void GzRender::ComputeUvPlaneEquasion(GzUvPlane* p)
{
	GzCoord v1{};
	GzCoord v2{};
	GzCoord product{};

	v1[0] = p->p2.coord[0] - p->p1.coord[0];
	v1[1] = p->p2.coord[1] - p->p1.coord[1];
	v1[2] = p->p2.uv[0] - p->p1.uv[0];

	v2[0] = p->p3.coord[0] - p->p2.coord[0];
	v2[1] = p->p3.coord[1] - p->p2.coord[1];
	v2[2] = p->p3.uv[0] - p->p2.uv[0];

	CrossProduct(v1, v2, product);

	p->au = product[0];
	p->bu = product[1];
	p->cu = product[2];
	p->du = -(p->cu * p->p2.uv[0]) - (p->bu * p->p2.coord[1]) - (p->au * p->p2.coord[0]);

	v1[0] = p->p2.coord[0] - p->p1.coord[0];
	v1[1] = p->p2.coord[1] - p->p1.coord[1];
	v1[2] = p->p2.uv[1] - p->p1.uv[1];

	v2[0] = p->p3.coord[0] - p->p2.coord[0];
	v2[1] = p->p3.coord[1] - p->p2.coord[1];
	v2[2] = p->p3.uv[1] - p->p2.uv[1];

	CrossProduct(v1, v2, product);

	p->av = product[0];
	p->bv = product[1];
	p->cv = product[2];
	p->dv = -(p->cv * p->p2.uv[1]) - (p->bv * p->p2.coord[1]) - (p->av * p->p2.coord[0]);
}

void GzRender::FillBoundBox(GzEdge* e1, GzEdge* e2, GzEdge* e3,
	GzPlane* p, GzPlane* r, GzPlane* g, GzPlane* b,
	GzNormalPlane* n, GzUvPlane* uv, long start[2], long end[2])
{
	float edgeEval[3]{};
	GzIntensity rValue, gValue, bValue, aValue;
	GzDepth zValue, zCurrent;
	GzVertex vCurrent{};
	GzColor cCurrent{};
	GzMatrix Xinvert;

	if (use_shadow_map)
	{
		// Generate Xinvert matrix (screen space to world space)
		int result = InvertMatrix(Ximage[matlevel], Xinvert);
#ifdef GZ_DEBUG_MODE
		if (result == GZ_FAILURE) OutputDebugStringA("Invert didn't work!\n");
#endif
	}

	for (int x = start[0]; x < end[0]; x++)
		for (int y = start[1]; y < end[1]; y++)
		{
			edgeEval[0] = e1->a * x + e1->b * y + e1->c;
			edgeEval[1] = e2->a * x + e2->b * y + e2->c;
			edgeEval[2] = e3->a * x + e3->b * y + e3->c;
			zCurrent = static_cast<GzDepth>(floorl(InterpolateParameter(p, static_cast<float>(x), static_cast<float>(y))));

			if ((edgeEval[0] < 0 && edgeEval[1] < 0 && edgeEval[2] < 0)
				|| (edgeEval[0] > 0 && edgeEval[1] > 0 && edgeEval[2] > 0)
				|| (FP_EQUALS(edgeEval[0], 0.0f, std::numeric_limits<float>::epsilon()) && edgeEval[1] > 0 && edgeEval[2] > 0)
				|| (FP_EQUALS(edgeEval[0], 0.0f, std::numeric_limits<float>::epsilon()) && edgeEval[1] < 0 && edgeEval[2] < 0))
			{
				GzGet(x, y, &rValue, &gValue, &bValue, &aValue, &zValue);
				if (zCurrent < zValue)
				{
					vCurrent.coord[0] = static_cast<float>(x);
					vCurrent.coord[1] = static_cast<float>(y);
					vCurrent.coord[2] = InterpolateParameter(p, vCurrent.coord[0], vCurrent.coord[1]);

					switch (interp_mode)
					{
					case GZ_NORMALS:
						InterpolateNormal(n, vCurrent.coord[0], vCurrent.coord[1], vCurrent.normal);
						if (tex_fun == nullptr) ComputeColor(&vCurrent, Ka, Kd, Ks, Xinvert);
						else
						{
							InterpolateUv(uv, vCurrent.coord[0], vCurrent.coord[1], vCurrent.uv);
							vCurrent.uv[0] *= ComputeWrapFactor(vCurrent.coord[2]);
							vCurrent.uv[1] *= ComputeWrapFactor(vCurrent.coord[2]);
							tex_fun(vCurrent.uv[0], vCurrent.uv[1], cCurrent);
							ComputeColor(&vCurrent, cCurrent, cCurrent, Ks, Xinvert);
						}
						break;
					case GZ_COLOR:
						vCurrent.color[0] = InterpolateParameter(r, vCurrent.coord[0], vCurrent.coord[1]);
						vCurrent.color[1] = InterpolateParameter(g, vCurrent.coord[0], vCurrent.coord[1]);
						vCurrent.color[2] = InterpolateParameter(b, vCurrent.coord[0], vCurrent.coord[1]);
						if (tex_fun != nullptr)
						{
							InterpolateUv(uv, vCurrent.coord[0], vCurrent.coord[1], vCurrent.uv);
							vCurrent.uv[0] *= ComputeWrapFactor(vCurrent.coord[2]);
							vCurrent.uv[1] *= ComputeWrapFactor(vCurrent.coord[2]);
							tex_fun(vCurrent.uv[0], vCurrent.uv[1], cCurrent);
							vCurrent.color[0] *= cCurrent[0];
							vCurrent.color[1] *= cCurrent[1];
							vCurrent.color[2] *= cCurrent[2];
						}
						break;
					case GZ_FLAT:
						vCurrent.color[0] = (e1->start.color[0] + e2->start.color[0] + e3->start.color[0]) / 3;
						vCurrent.color[1] = (e1->start.color[1] + e2->start.color[1] + e3->start.color[1]) / 3;
						vCurrent.color[2] = (e1->start.color[2] + e2->start.color[2] + e3->start.color[2]) / 3;
						break;
					case GZ_SHADOWMAP:
						vCurrent.color[0] = SHADOW_MAP_IMAGE_INTENSITY_FACTOR<float> *
							vCurrent.coord[2] / static_cast<float>((std::numeric_limits<int>::max)());
						vCurrent.color[1] = SHADOW_MAP_IMAGE_INTENSITY_FACTOR<float> *
							vCurrent.coord[2] / static_cast<float>((std::numeric_limits<int>::max)());
						vCurrent.color[2] = SHADOW_MAP_IMAGE_INTENSITY_FACTOR<float> *
							vCurrent.coord[2] / static_cast<float>((std::numeric_limits<int>::max)());
						break;
					default:
						memcpy(vCurrent.color, flatcolor, sizeof(GzColor));
						break;
					}
					GzPut(x, y, ctoi(vCurrent.color[0]), ctoi(vCurrent.color[1]), ctoi(vCurrent.color[2]), aValue, zCurrent);
				}
			}
		}
}

void GzRender::TranslateVertex(GzVertex* v, GzMatrix ximage, GzMatrix xnorm)
{
	float input[4]{}, output[4]{};

	// translate vertex
	memcpy(input, v->coord, sizeof(GzCoord));
	input[3] = static_cast<float>(1);

	TranslateCoord(ximage, input, output);

	v->coord[0] = (output[0] / output[3]) + currentOffset[0];
	v->coord[1] = (output[1] / output[3]) + currentOffset[1];
	v->coord[2] = output[2] / output[3];

	// translate normal
	memcpy(input, v->normal, sizeof(GzCoord));
	input[3] = static_cast<float>(1);

	TranslateCoord(xnorm, input, output);

	v->normal[0] = output[0] / output[3];
	v->normal[1] = output[1] / output[3];
	v->normal[2] = output[2] / output[3];

	NormalizeVector(3, v->normal);

	// translate UV
	v->uv[0] /= ComputeWrapFactor(v->coord[2]);
	v->uv[1] /= ComputeWrapFactor(v->coord[2]);
}

void GzRender::CrossProduct(GzCoord a, GzCoord b, GzCoord result)
{
	result[0] = a[1] * b[2] - a[2] * b[1];
	result[1] = a[2] * b[0] - a[0] * b[2];
	result[2] = a[0] * b[1] - a[1] * b[0];
}

void GzRender::NormalizeVector(long size, float v[])
{
	float norm = 0;

	for (int i = 0; i < size; i++)
		norm += powf(v[i], 2);

	norm = sqrtf(norm);

	for (int i = 0; i < size; i++)
		v[i] /= norm;
}

void GzRender::MultiplyMatrices(GzMatrix a, GzMatrix b, GzMatrix result)
{
	memset(result, 0, sizeof(GzMatrix));
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			for (int k = 0; k < 4; k++)
				result[i][j] += a[i][k] * b[k][j];
}

void GzRender::TranslateCoord(GzMatrix a, float b[4], float result[4])
{
	memset(result, 0, sizeof(float[4]));
	result[0] = a[0][0] * b[0] + a[0][1] * b[1] + a[0][2] * b[2] + a[0][3] * b[3];
	result[1] = a[1][0] * b[0] + a[1][1] * b[1] + a[1][2] * b[2] + a[1][3] * b[3];
	result[2] = a[2][0] * b[0] + a[2][1] * b[1] + a[2][2] * b[2] + a[2][3] * b[3];
	result[3] = a[3][0] * b[0] + a[3][1] * b[1] + a[3][2] * b[2] + a[3][3] * b[3];

}

void GzRender::ExtractRotation(GzMatrix input, GzMatrix output)
{
	GzCoord row{};
	memset(output, 0, sizeof(GzMatrix));
	for (int i = 0; i < 3; i++)
	{
		memcpy(row, input[i], sizeof(GzCoord));
		NormalizeVector(3, row);
		memcpy(output[i], row, sizeof(GzCoord));
	}
	output[3][3] = 1;
}

void GzRender::ComputeColor(GzVertex* v, GzColor ka, GzColor kd, GzColor ks, GzMatrix Xinvert)
{
	bool inShadow = false;
	int lightSpaceX, lightSpaceY;
	float eDotR, lDotN, eDotN, screenSpaceCoord[4], lightSpaceCoordHomogenous[4];
	GzIntensity rShadow, gShadow, bShadow, aShadow;
	GzDepth zShadow;
	GzCoord vEye{}, vReflection{}, vLight{}, vNormal{};

	memcpy(vEye, EYE, sizeof(GzCoord));
	memcpy(vNormal, v->normal, sizeof(GzCoord));

	// ambient light
	v->color[0] = ka[0] * ambientlight.color[0];
	v->color[1] = ka[1] * ambientlight.color[1];
	v->color[2] = ka[2] * ambientlight.color[2];

	// specular and diffuse for each light
	for (int i = 0; i < numlights; i++)
	{
		if (use_shadow_map && interp_mode == GZ_NORMALS && (GzRender*)(lights[i].shadow_map_renderer) != nullptr)
		{
			// Create screen space x,y,z coordinates to homogenous coordinates (used for Xls matrix multiplication) (needs to be a 4x1 matrix)
			memcpy(screenSpaceCoord, v->coord, sizeof(GzCoord));
			screenSpaceCoord[3] = static_cast<float>(1);

			// Transfrom screen space coordinates to light space coordinates
			TranslateCoord(Xinvert, screenSpaceCoord, lightSpaceCoordHomogenous);

			// Convet from homogenous coordinates to cartesian coordinates
			screenSpaceCoord[0] = lightSpaceCoordHomogenous[0] / lightSpaceCoordHomogenous[3];
			screenSpaceCoord[1] = lightSpaceCoordHomogenous[1] / lightSpaceCoordHomogenous[3];
			screenSpaceCoord[2] = lightSpaceCoordHomogenous[2] / lightSpaceCoordHomogenous[3];
			screenSpaceCoord[3] = static_cast<float>(1);

			TranslateCoord(((GzRender*)(lights[i].shadow_map_renderer))->Ximage[matlevel], screenSpaceCoord, lightSpaceCoordHomogenous);
			v->shadow[0] = lightSpaceCoordHomogenous[0] / lightSpaceCoordHomogenous[3];
			v->shadow[1] = lightSpaceCoordHomogenous[1] / lightSpaceCoordHomogenous[3];
			v->shadow[2] = lightSpaceCoordHomogenous[2] / lightSpaceCoordHomogenous[3];

			NearestNeighbor(v->shadow, &lightSpaceX, &lightSpaceY);
			((GzRender*)(lights[i].shadow_map_renderer))->GzGet(lightSpaceX, lightSpaceY, &rShadow, &gShadow, &bShadow, &aShadow, &zShadow);

			if (abs(v->shadow[2] - zShadow) > SHADOW_BIAS<float>)
			{
				inShadow = true;
#ifdef GZ_DEBUG_MODE
				OutputDebugStringA("In shadow!\n");
#endif
			}
			else
			{
				inShadow = false;
#ifdef GZ_DEBUG_MODE
				OutputDebugStringA("Not in shadow!\n");
#endif
			}
		}

		memcpy(vLight, lights[i].direction, sizeof(GzCoord));
		NormalizeVector(3, vLight);

		eDotN = DotProduct(3, vEye, vNormal);
		lDotN = DotProduct(3, vLight, vNormal);

		if ((eDotN > 0 && lDotN < 0)
			|| (eDotN < 0 && lDotN > 0))
			continue;
		else if (eDotN < 0 && lDotN < 0)
		{
			vNormal[0] = -vNormal[0];
			vNormal[1] = -vNormal[1];
			vNormal[2] = -vNormal[2];
			lDotN = DotProduct(3, vLight, vNormal);
		}

		vReflection[0] = (2 * lDotN) * vNormal[0] - vLight[0];
		vReflection[1] = (2 * lDotN) * vNormal[1] - vLight[1];
		vReflection[2] = (2 * lDotN) * vNormal[2] - vLight[2];
		NormalizeVector(3, vReflection);
		eDotR = DotProduct(3, vEye, vReflection);

		if (inShadow)
		{
			v->color[0] += SHADOW_KS[0] * ks[0] * lights[i].color[0] * powf(eDotR, spec);
			v->color[1] += SHADOW_KS[1] * ks[1] * lights[i].color[1] * powf(eDotR, spec);
			v->color[2] += SHADOW_KS[2] * ks[2] * lights[i].color[2] * powf(eDotR, spec);

			v->color[0] += SHADOW_KD[0] * kd[0] * lights[i].color[0] * lDotN;
			v->color[1] += SHADOW_KD[1] * kd[1] * lights[i].color[1] * lDotN;
			v->color[2] += SHADOW_KD[2] * kd[2] * lights[i].color[2] * lDotN;
		}
		else
		{
			v->color[0] += ks[0] * lights[i].color[0] * powf(eDotR, spec);
			v->color[1] += ks[1] * lights[i].color[1] * powf(eDotR, spec);
			v->color[2] += ks[2] * lights[i].color[2] * powf(eDotR, spec);

			v->color[0] += kd[0] * lights[i].color[0] * lDotN;
			v->color[1] += kd[1] * lights[i].color[1] * lDotN;
			v->color[2] += kd[2] * lights[i].color[2] * lDotN;
		}
	}

	v->color[0] = CLAMP(v->color[0], 0.0f, 1.0f);
	v->color[1] = CLAMP(v->color[1], 0.0f, 1.0f);
	v->color[2] = CLAMP(v->color[2], 0.0f, 1.0f);
}

void GzRender::InterpolateNormal(GzNormalPlane* p, float x, float y, GzCoord n)
{
	n[0] = -(p->dx + p->bx * y + p->ax * x) / p->cx;
	n[1] = -(p->dy + p->by * y + p->ay * x) / p->cy;
	n[2] = -(p->dz + p->bz * y + p->az * x) / p->cz;

	NormalizeVector(3, n);
}

void GzRender::InterpolateUv(GzUvPlane* p, float x, float y, GzTextureIndex uv)
{
	uv[0] = -(p->du + p->bu * y + p->au * x) / p->cu;
	uv[1] = -(p->dv + p->bv * y + p->av * x) / p->cv;
}

void GzRender::NearestNeighbor(GzCoord in, int* outX, int* outY)
{
	GzTextureIndex neighbors[4] = {
		{floorf(in[0]), floorf(in[1])},
		{ceilf(in[0]), floorf(in[1])},
		{floorf(in[0]), ceilf(in[1])},
		{ceilf(in[0]), ceilf(in[1])}
	};
	float distances[4]{};
	float minDistance = (std::numeric_limits<float>::max)();
	int minDistanceIndex = 0;

	for (int i = 0; i < 3; i++)
		distances[i] = sqrtf(powf(in[0] - neighbors[i][0], 2) + powf(in[1] - neighbors[1][1], 2));

	for (int i = 0; i < 3; i++)
		if (FP_LESS(distances[i], minDistance, std::numeric_limits<float>::epsilon()))
		{
			minDistanceIndex = i;
			minDistance = distances[i];
		}

	*outX = static_cast<int>(neighbors[minDistanceIndex][0]);
	*outY = static_cast<int>(neighbors[minDistanceIndex][1]);
}

float GzRender::InterpolateParameter(GzPlane* p, float x, float y)
{
	return (-(p->d + p->b * y + p->a * x) / p->c);
}

float GzRender::DotProduct(long size, float a[], float b[])
{
	float result = 0;

	for (int i = 0; i < size; i++)
		result += a[i] * b[i];

	return result;
}

float GzRender::ComputeWrapFactor(float Vzs)
{
	return (1 + (Vzs / ((std::numeric_limits<int>::max)() - Vzs)));
}

// generic matrix inversion, adapted from the Mesa 3D's GLUT implementation
int GzRender::InvertMatrix(GzMatrix in, GzMatrix out)
{
	GzMatrix invert;
	float det;
	int i;

	invert[0][0] = in[0][5] * in[0][10] * in[0][15] -
		in[0][5] * in[0][11] * in[0][14] -
		in[0][9] * in[0][6] * in[0][15] +
		in[0][9] * in[0][7] * in[0][14] +
		in[0][13] * in[0][6] * in[0][11] -
		in[0][13] * in[0][7] * in[0][10];

	invert[0][4] = -in[0][4] * in[0][10] * in[0][15] +
		in[0][4] * in[0][11] * in[0][14] +
		in[0][8] * in[0][6] * in[0][15] -
		in[0][8] * in[0][7] * in[0][14] -
		in[0][12] * in[0][6] * in[0][11] +
		in[0][12] * in[0][7] * in[0][10];

	invert[0][8] = in[0][4] * in[0][9] * in[0][15] -
		in[0][4] * in[0][11] * in[0][13] -
		in[0][8] * in[0][5] * in[0][15] +
		in[0][8] * in[0][7] * in[0][13] +
		in[0][12] * in[0][5] * in[0][11] -
		in[0][12] * in[0][7] * in[0][9];

	invert[0][12] = -in[0][4] * in[0][9] * in[0][14] +
		in[0][4] * in[0][10] * in[0][13] +
		in[0][8] * in[0][5] * in[0][14] -
		in[0][8] * in[0][6] * in[0][13] -
		in[0][12] * in[0][5] * in[0][10] +
		in[0][12] * in[0][6] * in[0][9];

	invert[0][1] = -in[0][1] * in[0][10] * in[0][15] +
		in[0][1] * in[0][11] * in[0][14] +
		in[0][9] * in[0][2] * in[0][15] -
		in[0][9] * in[0][3] * in[0][14] -
		in[0][13] * in[0][2] * in[0][11] +
		in[0][13] * in[0][3] * in[0][10];

	invert[0][5] = in[0][0] * in[0][10] * in[0][15] -
		in[0][0] * in[0][11] * in[0][14] -
		in[0][8] * in[0][2] * in[0][15] +
		in[0][8] * in[0][3] * in[0][14] +
		in[0][12] * in[0][2] * in[0][11] -
		in[0][12] * in[0][3] * in[0][10];

	invert[0][9] = -in[0][0] * in[0][9] * in[0][15] +
		in[0][0] * in[0][11] * in[0][13] +
		in[0][8] * in[0][1] * in[0][15] -
		in[0][8] * in[0][3] * in[0][13] -
		in[0][12] * in[0][1] * in[0][11] +
		in[0][12] * in[0][3] * in[0][9];

	invert[0][13] = in[0][0] * in[0][9] * in[0][14] -
		in[0][0] * in[0][10] * in[0][13] -
		in[0][8] * in[0][1] * in[0][14] +
		in[0][8] * in[0][2] * in[0][13] +
		in[0][12] * in[0][1] * in[0][10] -
		in[0][12] * in[0][2] * in[0][9];

	invert[0][2] = in[0][1] * in[0][6] * in[0][15] -
		in[0][1] * in[0][7] * in[0][14] -
		in[0][5] * in[0][2] * in[0][15] +
		in[0][5] * in[0][3] * in[0][14] +
		in[0][13] * in[0][2] * in[0][7] -
		in[0][13] * in[0][3] * in[0][6];

	invert[0][6] = -in[0][0] * in[0][6] * in[0][15] +
		in[0][0] * in[0][7] * in[0][14] +
		in[0][4] * in[0][2] * in[0][15] -
		in[0][4] * in[0][3] * in[0][14] -
		in[0][12] * in[0][2] * in[0][7] +
		in[0][12] * in[0][3] * in[0][6];

	invert[0][10] = in[0][0] * in[0][5] * in[0][15] -
		in[0][0] * in[0][7] * in[0][13] -
		in[0][4] * in[0][1] * in[0][15] +
		in[0][4] * in[0][3] * in[0][13] +
		in[0][12] * in[0][1] * in[0][7] -
		in[0][12] * in[0][3] * in[0][5];

	invert[0][14] = -in[0][0] * in[0][5] * in[0][14] +
		in[0][0] * in[0][6] * in[0][13] +
		in[0][4] * in[0][1] * in[0][14] -
		in[0][4] * in[0][2] * in[0][13] -
		in[0][12] * in[0][1] * in[0][6] +
		in[0][12] * in[0][2] * in[0][5];

	invert[0][3] = -in[0][1] * in[0][6] * in[0][11] +
		in[0][1] * in[0][7] * in[0][10] +
		in[0][5] * in[0][2] * in[0][11] -
		in[0][5] * in[0][3] * in[0][10] -
		in[0][9] * in[0][2] * in[0][7] +
		in[0][9] * in[0][3] * in[0][6];

	invert[0][7] = in[0][0] * in[0][6] * in[0][11] -
		in[0][0] * in[0][7] * in[0][10] -
		in[0][4] * in[0][2] * in[0][11] +
		in[0][4] * in[0][3] * in[0][10] +
		in[0][8] * in[0][2] * in[0][7] -
		in[0][8] * in[0][3] * in[0][6];

	invert[0][11] = -in[0][0] * in[0][5] * in[0][11] +
		in[0][0] * in[0][7] * in[0][9] +
		in[0][4] * in[0][1] * in[0][11] -
		in[0][4] * in[0][3] * in[0][9] -
		in[0][8] * in[0][1] * in[0][7] +
		in[0][8] * in[0][3] * in[0][5];

	invert[0][15] = in[0][0] * in[0][5] * in[0][10] -
		in[0][0] * in[0][6] * in[0][9] -
		in[0][4] * in[0][1] * in[0][10] +
		in[0][4] * in[0][2] * in[0][9] +
		in[0][8] * in[0][1] * in[0][6] -
		in[0][8] * in[0][2] * in[0][5];

	det = in[0][0] * invert[0][0] + in[0][1] * invert[0][4] + in[0][2] * invert[0][8] + in[0][3] * invert[0][12];

	if (det == 0)
		return GZ_FAILURE;

	det = 1 / det;

	for (i = 0; i < 16; i++)
		out[0][i] = invert[0][i] * det;

	return GZ_SUCCESS;
}

int CompareVertices(const void* left, const void* right)
{
	Point* l = (Point*)left;
	Point* r = (Point*)right;
	int result = 0;

	if (FP_LESS(l->coord[1], r->coord[1], std::numeric_limits<float>::epsilon())) result = -1;
	else if (FP_MORE(l->coord[1], r->coord[1], std::numeric_limits<float>::epsilon())) result = 1;
	else
	{
		if (FP_LESS(l->coord[0], r->coord[0], std::numeric_limits<float>::epsilon())) result = -1;
		else if (FP_MORE(l->coord[0], r->coord[0], std::numeric_limits<float>::epsilon())) result = 1;
		else result = 0;
	}

	return result;
}
