#ifndef _stone2_h_
#define _stone2_h_

/*
	File di inclusione
*/
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include "mmath.h"
#include "random2.h"

/*
	Definizioni
*/
#define FAULT			99999999
#define PFAULT			99
#define STOP_CELL		-1
#define MAX_LEN_STRING	256


/* 
	Definizioni Constanti 
*/
#define G           9.805
#define INV_G       0.101988781
#define SQ2         1.414213562
#define PI          3.141592654
#define PI4         0.785398164
#define RAD_TO_GRAD 57.29577951


typedef enum { FLY, ROLL } StoneStatus;

typedef struct
{
	P3d		pos, v;
	double	dz;
	char	deleted;

} typePath;

typedef struct
{
    long    id;

	double  x, y, v, dz;

} typeInfoStat;

typedef struct 
{
	char	type;
	double	Tx, Ty, Tz, Rx, Ry, Rz;
	P3d		p0, p1, p2;
	double	rToPlane[16], rFromPlane[16], toPlane[16], fromPlane[16];

} typePlane;

typedef struct
{
	double	sw_x, sw_y, ne_x, ne_y;
	short	cell;
	short	rows, cols;

} typeGeometry; 

typedef struct
{
	double	v0, min_v, min_v2, tab;
	double	fly_step, roll_step, short_bounce2, fly_roll_thresh2;
	int		max_path, gen_3d_vect, gen_2d_vect;
	char	elev_f[MAX_LEN_STRING]; 
	char	stst_f[MAX_LEN_STRING]; 
	char	v_elas_f[MAX_LEN_STRING]; 
	char	h_elas_f[MAX_LEN_STRING]; 
	char	frict_f[MAX_LEN_STRING];
	char	OUT_STONES_FILE[MAX_LEN_STRING];
	char	OUT_DTM_FILE[MAX_LEN_STRING];
	char	OUT_POINT_2D_FILE[MAX_LEN_STRING];
	char	OUT_ATTRIBUTES_2D_FILE[MAX_LEN_STRING];
	char	OUT_COUNTERS_FILE[MAX_LEN_STRING];
	char	OUT_MAX_VEL_FILE[MAX_LEN_STRING];
	char	OUT_MAX_DZ_FILE[MAX_LEN_STRING];
    char    STAT_INPUT_FILE[MAX_LEN_STRING];
	int		stoc_flag, stoc_angle, stoc_vel, stoc_hel, stoc_frict;
	int		EnabledVideoOutput;
	int		EnabledLogFile;
	int     SwitchVelType;
	char	AccMtrxFile[MAX_LEN_STRING]; 
	double  FromAccToVel;

	double gdQuotaUMFactor;
    double gdVeloUMFactor;
    double gdTab2;
    double gFlagInfoStat;
    int giRockType;
	int randomGenerator;
	int mANG_STOCH_FUNC;
	int mVREST_STOCH_FUNC;

} typeParams;

typedef struct
{
	typeGeometry* gGeometry;
	int giRows;
	int giCols;
	long glRowsxCols;
	double gdCell;
	double gdInvCell;
	double gdInvCell0001;
	double gdInvCellSq20001;
	double gdOffset;

	double gdInvMaxRandPlusOne;
} globalParams;


typedef struct
{
	double gdStocAngleR;
	double gdStocVelR;
	double gdStocHelR;
	double gdStocFrictR;
	double gdStocAngle2R;
	double gdStocVel2R;
	double gdStocHel2R;
	double gdStocFrict2R;
	
	long *gplQuota;
	long *gplStartStop;
	long *gplFrict;
	long *gplCountStones;
	long *gplVelo;
	long *gplMaxQuota;
	long *gplStartVel;
	short *gpcVElas;
	short *gpcHElas;

	long glConter2d;
	int giCountRuns;
	long glIdMasso;
	long glLastPiv;
	int giUstSlice;

	typePlane gPlane;
	typePath *gpPathRoot;
    typePath *gpPathCur;
    P3d*    gP3dZero;
    short  gsKernel9[9];
} runtimeParams;

int NewPlane(typeParams* stoneRunParams, runtimeParams* stoneRunTimeParams, globalParams* stoneGlobalParams, P3d *cp, MathContext *ctx);
long Pivot(globalParams* gParams, double X, double Y);
void RoundP3d(P3d* p);
void Filter();
void MarkVelo();
void MarkPath();
void MarkQuota();
void WritePath(runtimeParams* rtParams, typeParams* tParams, globalParams* gParams, StoneStatus status);
void WriteFile3d(StoneStatus status);
void WriteFiles2d();
void PivToXY(long piv, double* X, double* Y);
double P3dDist2(P3d* p1, P3d* p2);
double RoundDouble(double d);
double GetRandAngle(typeParams* gParams, runtimeParams* rtParams, globalParams* globParams, UniSave* uniData, int dir);
double GetRandVElas(typeParams* gParams, runtimeParams* rtParams, globalParams* globParams, UniSave* uniData, long piv);
double GetRandHElas(typeParams* gParams, runtimeParams* rtParams, globalParams* globParams, UniSave* uniData, long piv);
double GetRandFrict(typeParams* gParams, runtimeParams* rtParams, globalParams* globParams, UniSave* uniData, long piv);

void runStone(typeParams* stoneRunParams, globalParams* stoneGlobalParams);
void WriteP3dToLog(const char *label, P3d* p);

// extern typePlane gPlane;

// extern typeGeometry gGeometry;

// extern typeParams gParams;

// extern char *gpcVElas;
// extern char *gpcHElas;

// extern long *gplQuota;
// extern long *gplStartStop;
// extern long *gplFrict;
// extern long *gplCountStones;
// extern long *gplVelo;
// extern long *gplMaxQuota;
// extern long *gplStartVel;

// extern typePath *gpPathRoot;
// extern typePath *gpPathCur;

// extern P3d    gP3dZero;
// extern int    giRows, giCols;
// extern long   glRowsxCols;
// extern int    giUstSlice;
// extern long   glLastPiv;
// extern long   glConter2d;
// extern short  gsKernel9[9]; /* offset alle celle attigue ad una cella data in una matrice */
// extern double gdCell;
// extern double gdInvCell;
// extern double gdInvCell0001;
// extern double gdInvCellSq20001;
// extern double gdTab2;
// extern double gdOffset;
// extern double gdInvMaxRandPlusOne; 
// extern double gdRadius; 

// extern double gdStocAngleR, gdStocVelR, gdStocHelR, gdStocFrictR;
// extern double gdStocAngle2R, gdStocVel2R, gdStocHel2R, gdStocFrict2R;

// extern FILE *gFileLog;
// extern FILE *gFileUst1, *gFileUst2;
// extern FILE *gFilePoints, *gFileAttrib;
// extern FILE *gFileInfoStat;

void TrackPoints(typeParams* stoneRunParams, globalParams* stoneGlobalParams, 
				runtimeParams* stoneRuntimeParams, UniSave* uniData, MathContext* ctx);
void TrackNoPoints();

// extern int     giRockType, giCountRuns; /* 0: point; 1: sphere; 2: cylinder; 3:disco  */
// extern int     gFlagInfoStat;
// extern double  gdRadius, gdRadiusxRadius, gdMinRadius, gdMaxRadius, gdDensity;
// extern double  gdQuotaUMFactor, gdVeloUMFactor;
// extern double  gdWidth, gdVolume, gdMassa, gdInvMassa, gdFrontalArea, gdDragCoefficient;

// extern long glIdMasso;

#endif