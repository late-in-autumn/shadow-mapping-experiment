#include	"gz.h"
#ifndef GZRENDER_
#define GZRENDER_


/* Camera defaults */
constexpr auto DEFAULT_FOV = 35.0;
constexpr auto DEFAULT_IM_Z = (-10.0);  /* world coords for image plane origin */
constexpr auto DEFAULT_IM_Y = (5.0);  /* default look-at point = 0,0,0 */
constexpr auto DEFAULT_IM_X = (-10.0);

constexpr auto DEFAULT_AMBIENT = { 0.1, 0.1, 0.1 };
constexpr auto DEFAULT_DIFFUSE = { 0.7, 0.6, 0.5 };
constexpr auto DEFAULT_SPECULAR = { 0.2, 0.3, 0.4 };
constexpr auto DEFAULT_SPEC = 32;

constexpr auto MATLEVELS = 100;	/* how many matrix pushes allowed */
constexpr auto MAX_LIGHTS = 10;	/* how many lights allowed */

// constants
// some of them, like PI and IDENTITY, should never be changed
// others, can be changed, read the symbol name for exact meaning
template<typename T> constexpr T PI = T(3.14159265358979323846);
template<typename T> constexpr T SHADOW_MAP_IMAGE_INTENSITY_FACTOR = T(0.75);
template<typename T> constexpr T SHADOW_BIAS = T(0.5);
template<typename T> constexpr T SHADOW_DIFFUSE_FACTOR = T(0.5);
template<typename T> constexpr T SHADOW_SPECULAR_FACTOR = T(0.7);
constexpr GzMatrix IDENTITY = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 };
constexpr GzCoord EYE = { 0, 0, -1 };
constexpr GzColor FULL_INTENSITY = { 1.0f, 1.0f, 1.0f };
constexpr GzColor SHADOW_KD = { SHADOW_DIFFUSE_FACTOR<float>, SHADOW_DIFFUSE_FACTOR<float>, SHADOW_DIFFUSE_FACTOR<float> };
constexpr GzColor SHADOW_KS = { SHADOW_SPECULAR_FACTOR<float>, SHADOW_SPECULAR_FACTOR<float>, SHADOW_SPECULAR_FACTOR<float> };
constexpr int NUM_QUADRANTS = 4;

// a more generic version of the array offset finder, assumes a square matrix
inline int MATRIX_OFFSET(int size, int x, int y) {
	return (x + y * size);
}

// floating point comparison helper functions based on Don Knuth's TAOCP
// left ~= right
// Don Knuth calls this "essentially equal" in TAOCP
template<typename T> bool FP_EQUALS(T left, T right, T epsilon)
{
	return std::abs(left - right) <= ((std::abs(left) > std::abs(right) ? std::abs(right) : std::abs(left)) * epsilon);
}

// left < right
// Don Knuth calls this "definitely less than" in TAOCP
template<typename T> bool FP_LESS(T left, T right, T epsilon)
{
	return (right - left) > ((std::abs(left) < std::abs(right) ? std::abs(right) : std::abs(left)) * epsilon);
}

// left > right
// Don Knuth calls this "definitely more than" in TAOCP
template<typename T> bool FP_MORE(float left, float right, T epsilon)
{
	return (left - right) > ((std::abs(left) < std::abs(right) ? std::abs(right) : std::abs(left)) * epsilon);
}

// a simple clamp function
template <typename T> T CLAMP(T input, T floor, T ceil)
{
	if (input < floor) input = floor;
	if (input > ceil) input = ceil;
	return input;
}

// convert degree to radian
template <typename T> T TO_RADIAN(T degree) {
	return (degree * (PI<T> / 180));
}

typedef struct {
	GzCoord coord;
	GzCoord normal;
	GzColor color;
	GzTextureIndex uv;
	GzCoord shadow;
} Point;

int CompareVertices(const void* left, const void* right);

class GzRender {			/* define a renderer */
public:
	unsigned short	xres;
	unsigned short	yres;
	GzPixel* pixelbuffer;		/* frame buffer array */
	char* framebuffer;
	bool  use_shadow_map;

	GzCamera	m_camera;
	short		matlevel;	        /* top of stack - current xform */
	GzMatrix	Ximage[MATLEVELS];	/* stack of xforms (Xsm) */
	GzMatrix	Xnorm[MATLEVELS];	/* xforms for norms (Xim) */
	GzMatrix	Xsp;		        /* NDC to screen (pers-to-screen) */
	GzColor		flatcolor;          /* color state for flat shaded triangles */
	int			interp_mode;
	int			numlights;
	GzLight		lights[MAX_LIGHTS];
	GzLight		ambientlight;
	GzColor		Ka, Kd, Ks;
	float		spec;		/* specular power */
	GzTexture	tex_fun;    /* tex_fun(float u, float v, GzColor color) */

	// Constructors
	GzRender(int xRes, int yRes);
	~GzRender();

	// Function declaration

	// HW1: Display methods
	int GzDefault();
	int GzBeginRender();
	int GzPut(int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z);
	int GzGet(int i, int j, GzIntensity* r, GzIntensity* g, GzIntensity* b, GzIntensity* a, GzDepth* z);

	int GzFlushDisplay2File(FILE* outfile);
	int GzFlushDisplay2FrameBuffer();

	// HW2: Render methods
	int GzPutAttribute(int numAttributes, GzToken* nameList, GzPointer* valueList);
	int GzPutTriangle(int numParts, GzToken* nameList, GzPointer* valueList);

	// HW3
	int GzDefaultCamera();
	int GzPutCamera(GzCamera camera);
	int GzPushMatrix(GzMatrix	matrix);
	int GzPopMatrix();

	// Extra methods: NOT part of API - just for general assistance */
	inline int ARRAY(int x, int y) { return (x + y * xres); }	/* simplify fbuf indexing */
	inline short	ctoi(float color) { return(short)((int)(color * ((1 << 12) - 1))); }		/* convert float color to GzIntensity short */


	// Object Translation
	int GzRotXMat(float degree, GzMatrix mat);
	int GzRotYMat(float degree, GzMatrix mat);
	int GzRotZMat(float degree, GzMatrix mat);
	int GzTrxMat(GzCoord translate, GzMatrix mat);
	int GzScaleMat(GzCoord scale, GzMatrix mat);

	// methods needed for shadow mapping
	int GzSetShadowMap(bool in);

private:
	typedef struct {
		GzCoord	coord;
		GzCoord normal;
		GzColor color;
		GzTextureIndex uv;
		GzCoord shadow;
	} GzVertex;

	typedef struct {
		GzVertex start;
		GzVertex end;
		float a;
		float b;
		float c;
	} GzEdge;

	typedef struct {
		GzVertex p1;
		GzVertex p2;
		GzVertex p3;
		float a;
		float b;
		float c;
		float d;
	} GzPlane;

	typedef struct {
		GzVertex p1;
		GzVertex p2;
		GzVertex p3;
		float ax;
		float ay;
		float az;
		float bx;
		float by;
		float bz;
		float cx;
		float cy;
		float cz;
		float dx;
		float dy;
		float dz;
	} GzNormalPlane;

	typedef struct {
		GzVertex p1;
		GzVertex p2;
		GzVertex p3;
		float au;
		float bu;
		float cu;
		float du;
		float av;
		float bv;
		float cv;
		float dv;
	} GzUvPlane;

	typedef GzVertex GzTriangle[3];

	typedef float GzOffsetsXY[3];

	GzTriangle currentTriangle;
	GzOffsetsXY currentOffset;

	// renderer subroutines
	void UpdateCameraMatrices();
	void RenderTriangle(GzVertex* v1, GzVertex* v2, GzVertex* v3);
	void LinearEvaulator(GzVertex* v1, GzVertex* v2, GzVertex* v3);
	void SortVertices(GzVertex* v1, GzVertex* v2, GzVertex* v3);
	void ComputeBoundBox(GzVertex* v1, GzVertex* v2, GzVertex* v3, long start[2], long end[2]);
	void ComputeEdgeEquasion(GzEdge* e);
	void ComputePlaneEquasion(GzPlane* p);
	void ComputeNormalPlaneEquasion(GzNormalPlane* p);
	void ComputeColorPlaneEquasions(GzPlane* pr, GzPlane* pg, GzPlane* pb);
	void ComputeUvPlaneEquasion(GzUvPlane* p);
	void FillBoundBox(GzEdge* e1, GzEdge* e2, GzEdge* e3,
		GzPlane* p, GzPlane* r, GzPlane* g, GzPlane* b,
		GzNormalPlane* n, GzUvPlane* uv, long start[2], long end[2]);
	void TranslateVertex(GzVertex* v, GzMatrix ximage, GzMatrix xnorm);
	void ExtractRotation(GzMatrix input, GzMatrix output);
	void CrossProduct(GzCoord a, GzCoord b, GzCoord result);
	void NormalizeVector(long size, float v[]);
	void MultiplyMatrices(GzMatrix a, GzMatrix b, GzMatrix result);
	void TranslateCoord(GzMatrix a, float b[4], float result[4]);
	void ComputeColor(GzVertex* v, GzColor ka, GzColor kd, GzColor ks, GzMatrix Xls);
	void InterpolateNormal(GzNormalPlane* p, float x, float y, GzCoord n);
	void InterpolateUv(GzUvPlane* p, float x, float y, GzTextureIndex uv);
	void NearestNeighbor(GzCoord in, int* outX, int* outY);
	float InterpolateParameter(GzPlane* p, float x, float y);
	float DotProduct(long size, float a[], float b[]);
	float ComputeWrapFactor(float Vzs);
	int InvertMatrix(GzMatrix in, GzMatrix out); // generic matrix inversion, adapted from the Mesa 3D's GLUT implementation
};
#endif
