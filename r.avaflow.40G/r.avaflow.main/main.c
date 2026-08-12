/*
##############################################################################
#
# MODULE:       r.avaflow.main.c
#
# AUTHORS:      Martin Mergili and Shiva P. Pudasaini
# CONTRIBUTORS: Massimiliano Alvioli, Wolfgang Fellin, Jan-Thomas Fischer,
#               Sigridur S. Gylfadottir, Markus Metz, Markus Neteler, 
#               Alexander Ostermann, Matthias Rauter
#
# PURPOSE:      The simulation model for avalanche and debris flows
#               Flow routing script
#
# COPYRIGHT:    (c) 2008 - 2025 by the authors
#               (c) 2020 - 2025 by the University of Graz
#               (c) 2013 - 2021 by the BOKU University, Vienna
#               (c) 2015 - 2020 by the University of Vienna
#               (c) 2014 - 2020 by the University of Bonn
#               (c) 2020 - 2025 by the Technical University of Munich
#               (c) 2000 - 2025 by the GRASS Development Team
#
# VERSION:      20250902 (2 September 2025)
#
#               This program is free software under the GNU General Public
#               License (>=v2). Read the file COPYING that comes with GRASS
#               for details.
#
##############################################################################
*/


#define WITHGRASS // use of GRASS GIS (REMOVE THIS LINE IF GRASS IS NOT USED)


#include <fcntl.h> // libraries
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include <limits.h> 


#ifdef WITHGRASS


    #include <grass/gmath.h>
    #include <grass/gis.h>
    #include <grass/gprojects.h>
    #include <grass/glocale.h>
    #include <grass/raster.h>
    #include <grass/segment.h>


#endif


struct ico { // global structure for input constants and control variables

    char *MAINMAPSET; int MULT; int LIMITER; int MODEL; int M; int N; int IMAX; int ENTR; int ENTR2; int ENTR3; int ZONES; int CENTR; int CDEFORM;
    int ELEV; int RELM; int RELM2; int RELM3; int RELV; int RELV2; int RELV3; int PHI; int PHI2; int PHI3; int DELTA; int DELTA2; int DELTA3; 
    int ADDFRI; int ADDFRI2; int ADDFRI3; int COHESION; int COHESION2; int COHESION3; int NY; int NY2; int NY3; int AMBDRAG; int TRANS12; int TRANS13; int TRANS23; 
    int ZFRAG; int TRELEASE; int TRELSTOP; int STOPTIME; int TSLIDE; float PI; float UNDEF; float HFLOWMIN; float CSZ; float BDWEST; float BDNORTH; float BDSOUTH; 
    float BDEAST; float GRAVITY; float CFL[2]; float IMPTHR[4]; int CURVCTRL; int SURFACE; int ENTRAINMENT; int STOPPING; int NVECTMIN; int PMAX; int MESH; 
    float SLOMO; float FFLOW; int HYDADD; int ORIGINAL; int SEPFLUX; int COLLAPSE; int AFLAG; int LAYERS0; int LAYERS; int XREL; int YREL; float XDIST; float YDIST; 
    int GLACIER; int NONHYDRO; int PROFILE; int IMPACTAREA; int HDEPOSIT; int CTRLPOINTS; int PBG; int TSUNAMI; int TOPOFACT; int TEMP;


    #ifdef WITHGRASS


        float NSEGRC; float NSEGS;


    #endif


};

struct flow { // global structure for flow parameters

    float RHO1; float RHO2; float RHO3; float PHI0[3]; float PHI[3]; float DELTA0[3]; float DELTA[3]; float ADDFRI0[3]; float ADDFRI[3]; 
    float COHESION0[3]; float COHESION[3]; float NY0[3]; float NY[3]; float CDEFORM0; float CDEFORM[3]; float AMBDRAG0; float AMBDRAG; 
    float REP; int DRAG_je; float J; float TRANS120; float TRANS130; float TRANS230; float TRANS12; float TRANS13; float TRANS23; float FRAGSML; float FRAGML; 
    float ESP; float US; float CENTR0; float CENTR; float THETAS; float CSTOP; float VM_n0; float VM_l; float VM_n; float VM_fact; 
    float DRAG_k; float DRAG_m; int DRAG_n; float DRAG_rep; float DRAG_vterm; float VIS_chi[3]; float VIS_xi[3]; float VIS_a[3]; float VIS_ry; float NY_exp; 
    float KPMAX[3]; float sep_SHEAR; float disp_MULT; float CCONST; float EKINCOEF; float FRICMIN; float FRI_exp; float BETAPLAIN; float LAMBDA_FL; float LAMBDA_NU; 
    float GZMIN; float HREF; float TEMP0; float CR; float CI; float CF; float LAMBDAS; float LAMBDAB; float LAMBDA; float KAPPAR; float KAPPAI; float KAPPAF; 
    float TEMPB; float TEMPS; float L; float CL; float Q; float AS; float SLIDEEXP[3]; float SLIDEDEF[3];
};


// -- START -- Functions ----------------------------------------------------------------------------------------


float ffmax ( float gval1, float gval2 ) { // function for maximum value (float)

    float gmax;
    if ( gval2 > gval1 ) gmax = gval2;
    else gmax = gval1;
    return gmax;
}

float ffmin ( float gval1, float gval2 ) { // function for minimum value (float)

    float gmin;
    if ( gval2 < gval1 ) gmin = gval2;
    else gmin = gval1;
    return gmin;
}

float fabsmin ( float gval1, float gval2 ) { // function for minimum absolute value (float)

    float gmin;
    if ( fabs(gval2) < fabs(gval1) ) gmin = gval2;
    else gmin = gval1;
    return gmin;
}

float fsign ( float gval ) { // function for sign

    float gsign;
    if ( gval < 0 ) gsign = -1;
    else if ( gval == 0 ) gsign = 0;
    else gsign = 1;
    return gsign;
}

int fvalid ( float ghflow, float ghflow2, float ghflow3, struct ico sico ) { // function for minimum flow height

    int gvalid;

    gvalid = 0;
    if ( ghflow > sico.HFLOWMIN ) gvalid = 1;
    if ( sico.MODEL == 7 && ghflow2 > sico.HFLOWMIN ) gvalid = 1;
    if ( sico.MODEL == 7 && ghflow3 > sico.HFLOWMIN ) gvalid = 1;

    return gvalid;
}

int **alloc_imatrix(int mrows, int mcols) { // function for allocating 2D integer arrays
    int **mm;
    int mi;
    mm = (int **) calloc(mrows, sizeof(int *));
    mm[0] = (int *) calloc(mrows * mcols, sizeof(int));
    for (mi = 1; mi < mrows; mi++)
    mm[mi] = mm[mi - 1] + mcols;
    return mm;
}

float **alloc_dmatrix (int mrows, int mcols) { // function for allocating 2D float arrays
    float **mm;
    int mi;
    mm = (float **) calloc(mrows, sizeof(float *));
    mm[0] = (float *) calloc(mrows * mcols, sizeof(float));
    for (mi = 1; mi < mrows; mi++)
    mm[mi] = mm[mi - 1] + mcols;
    return mm;
}

char **alloc_cmatrix (int mrows, int mcols) { // function for allocating 2D char arrays
    char **mm;
    int mi;
    mm = (char **) calloc(mrows, sizeof(char *));
    mm[0] = (char *) calloc(mrows * mcols, sizeof(char));
    for (mi = 1; mi < mrows; mi++)
    mm[mi] = mm[mi - 1] + mcols;
    return mm;
}

float ***alloc_dmatrix3 (int mrows, int mcols, int mdepths) { // function for allocating 3D float arrays
    float ***mm;
    int mi, mj;
    mm = (float***)calloc(mrows, sizeof(float**));

    for (mi=0; mi<mrows; mi++) {
        mm[mi] = (float**)calloc(mcols, sizeof(float*));

        for (mj=0; mj<mcols; mj++) {
            mm[mi][mj] = (float*)calloc(mdepths, sizeof(float));
        }
    }
    return mm;
}

void free_dmatrix3 (float ***mm,int mi,int mj) { // function for freeing 3D float arrays
    int i,j;

    for(i=0;i<mi;i++)
    {
        for(j=0;j<mj;j++)
        {
            free(mm[i][j]);
        }
        free(mm[i]);
    }
    free(mm);
}

int fiparam ( FILE *gfparam ) { // function for reading integer parameters

    int gparam;
    char gsparam[50];
    char **strtodhlp=NULL;
    if ( fgets(gsparam, 50, gfparam) == NULL ) return -1;
    gparam=strtod(gsparam, strtodhlp);
    return gparam;
}

float fdparam ( FILE *gfparam ) { // function for reading float parameters

    float gparam;
    char gsparam[50];
    char **strtodhlp=NULL;
    if ( fgets(gsparam, 50, gfparam) == NULL ) return -1;
    gparam=strtod(gsparam, strtodhlp);
    return gparam;
}

char *fcparam ( FILE *gfparam ) { // function for reading character parameters

    char *gsparam = (char*) calloc ( 1010, sizeof(char));
    if ( fgets(gsparam, 1010, gfparam) == NULL ) return "-1";
    gsparam[strnlen(gsparam, 1010) - 1] = '\0';
    return gsparam;
}

char *flparam ( FILE *gfparam ) { // function for reading ascii raster line

    int asclinemax = 100000;
    char *gsparam2 = (char*) calloc (( asclinemax + 20 ), sizeof(char));
    if ( fgets(gsparam2, asclinemax, gfparam) == NULL ) return "-1";
    gsparam2[strnlen(gsparam2, asclinemax+10) - 1] = '\0';
    return gsparam2;
}

int *fin ( int gi, float *gpelev, struct ico sico ) { // function for identifying cell environment

    int *gin;

    gin = (int*) calloc( 9, sizeof(int));

    gin[0] = gi;
    if ( gi%sico.N == 0 || (gi+1)%sico.N == 0 || gi < sico.N || gi > sico.IMAX - sico.N || gpelev[gi] == sico.UNDEF ) { // edge cells (to be excluded from computation)
    
        gin[1] = sico.IMAX;
        gin[2] = sico.IMAX;
        gin[3] = sico.IMAX;
        gin[4] = -1;
        gin[5] = -1;
        gin[6] = -1;
        gin[7] = sico.IMAX;
        gin[8] = -1;
    }
    else { // cells to be included in computation
        gin[1] = gi + 1;
        gin[2] = gi + sico.N;
        gin[3] = gi + sico.N + 1;
        gin[4] = gi - 1;
        gin[5] = gi - sico.N;
        gin[6] = gi - sico.N - 1;
        gin[7] = gi + sico.N - 1;
        gin[8] = gi - sico.N + 1;
    }

    return gin;
}

int *f0noosc ( float *ghflow, struct ico sico ) { // function for controlling slope at flow boundaries (full cells)

    float gh[4];
    int *gcin;

    gcin = (int*) calloc( 6, sizeof(int));

    gh[0] = ghflow[2];
    gh[1] = ghflow[5];
    gh[2] = ghflow[1];
    gh[3] = ghflow[4];
    
    if ( gh[0] > sico.HFLOWMIN || sico.SURFACE == 0 || sico.SURFACE == 2 ) gcin[0] = 2; else gcin[0] = 0;
    if ( gh[1] > sico.HFLOWMIN || sico.SURFACE == 0 || sico.SURFACE == 2 ) gcin[1] = 5; else gcin[1] = 0;
    if ( gh[2] > sico.HFLOWMIN || sico.SURFACE == 0 || sico.SURFACE == 2 ) gcin[2] = 1; else gcin[2] = 0;
    if ( gh[3] > sico.HFLOWMIN || sico.SURFACE == 0 || sico.SURFACE == 2 ) gcin[3] = 4; else gcin[3] = 0;

    if ( gcin[0] == 0 || gcin[1] == 0 ) gcin[4] = 1; else gcin[4] = 2;

    if ( gcin[2] == 0 || gcin[3] == 0 ) gcin[5] = 1; else gcin[5] = 2;

    return gcin;
}

int *fnoosc ( float **ghflow, int **gin, int gi, struct ico sico ) { // function for controlling slope at flow boundaries (full cells)

    float gh[4], ghall;
    int *gcin;

    gcin = (int*) calloc( 6, sizeof(int));

    if ( sico.MODEL <= 3 ) {

        gh[0] = ghflow[gin[gi][5]][0];
        gh[1] = ghflow[gin[gi][2]][0];
        gh[2] = ghflow[gin[gi][4]][0];
        gh[3] = ghflow[gin[gi][1]][0];
 
    } else {
    
        gh[0] = ghflow[gin[gi][5]][0] + ghflow[gin[gi][5]][3] + ghflow[gin[gi][5]][6];
        gh[1] = ghflow[gin[gi][2]][0] + ghflow[gin[gi][2]][3] + ghflow[gin[gi][2]][6];
        gh[2] = ghflow[gin[gi][4]][0] + ghflow[gin[gi][4]][3] + ghflow[gin[gi][4]][6];
        gh[3] = ghflow[gin[gi][1]][0] + ghflow[gin[gi][1]][3] + ghflow[gin[gi][1]][6];
    }
    
    ghall =  gh[0] + gh[1] + gh[2] + gh[3];
    
    if ( gh[0] > sico.HFLOWMIN || ghall <= sico.HFLOWMIN || sico.SURFACE == 0 || sico.SURFACE == 2 ) gcin[0] = gin[gi][5]; else gcin[0] = gin[gi][0];
    if ( gh[1] > sico.HFLOWMIN || ghall <= sico.HFLOWMIN || sico.SURFACE == 0 || sico.SURFACE == 2 ) gcin[1] = gin[gi][2]; else gcin[1] = gin[gi][0];
    if ( gh[2] > sico.HFLOWMIN || ghall <= sico.HFLOWMIN || sico.SURFACE == 0 || sico.SURFACE == 2 ) gcin[2] = gin[gi][4]; else gcin[2] = gin[gi][0];
    if ( gh[3] > sico.HFLOWMIN || ghall <= sico.HFLOWMIN || sico.SURFACE == 0 || sico.SURFACE == 2 ) gcin[3] = gin[gi][1]; else gcin[3] = gin[gi][0];

    if ( gcin[0] == gin[gi][0] || gcin[1] == gin[gi][0] ) gcin[4] = 1; else gcin[4] = 2;
    if ( gcin[2] == gin[gi][0] || gcin[3] == gin[gi][0] ) gcin[5] = 1; else gcin[5] = 2;

    return gcin;
}

int *fwnoosc ( float ***ghflow, int **gin, int gi, int gj , struct ico sico ) { // function for controlling slope at flow boundaries (half cells)

    float gh[4];
    int *gcin;

    gcin = (int*) calloc( 6, sizeof(int));

    if ( sico.MODEL <= 3 ) {

        gh[0] = ghflow[gin[gi][5]][gj][0];
        gh[1] = ghflow[gin[gi][2]][gj][0];
        gh[2] = ghflow[gin[gi][4]][gj][0];
        gh[3] = ghflow[gin[gi][1]][gj][0];
 
    } else {
    
        gh[0] = ghflow[gin[gi][5]][gj][0] + ghflow[gin[gi][5]][gj][3] + ghflow[gin[gi][5]][gj][6];
        gh[1] = ghflow[gin[gi][2]][gj][0] + ghflow[gin[gi][2]][gj][3] + ghflow[gin[gi][2]][gj][6];
        gh[2] = ghflow[gin[gi][4]][gj][0] + ghflow[gin[gi][4]][gj][3] + ghflow[gin[gi][4]][gj][6];
        gh[3] = ghflow[gin[gi][1]][gj][0] + ghflow[gin[gi][1]][gj][3] + ghflow[gin[gi][1]][gj][6];
    }
    
    if ( gh[0] > sico.HFLOWMIN || sico.SURFACE == 0 || sico.SURFACE == 2 ) gcin[0] = gin[gi][5]; else gcin[0] = gin[gi][0];
    if ( gh[1] > sico.HFLOWMIN || sico.SURFACE == 0 || sico.SURFACE == 2 ) gcin[1] = gin[gi][2]; else gcin[1] = gin[gi][0];
    if ( gh[2] > sico.HFLOWMIN || sico.SURFACE == 0 || sico.SURFACE == 2 ) gcin[2] = gin[gi][4]; else gcin[2] = gin[gi][0];
    if ( gh[3] > sico.HFLOWMIN || sico.SURFACE == 0 || sico.SURFACE == 2 ) gcin[3] = gin[gi][1]; else gcin[3] = gin[gi][0];

    if ( gcin[0] == gin[gi][0] || gcin[1] == gin[gi][0] ) gcin[4] = 1; else gcin[4] = 2;
    if ( gcin[2] == gin[gi][0] || gcin[3] == gin[gi][0] ) gcin[5] = 1; else gcin[5] = 2;

    return gcin;
}

int *fw0noosc ( float *ghflow, struct ico sico ) { // function for controlling slope at flow boundaries (half cells)


    float gh[4];
    int *gcin;

    gcin = (int*) calloc( 6, sizeof(int));

    gh[0] = ghflow[5];
    gh[1] = ghflow[2];
    gh[2] = ghflow[4];
    gh[3] = ghflow[1];
 
    if ( gh[0] > sico.HFLOWMIN || sico.SURFACE == 0 || sico.SURFACE == 2 ) gcin[0] = 5; else gcin[0] = 0;
    if ( gh[1] > sico.HFLOWMIN || sico.SURFACE == 0 || sico.SURFACE == 2 ) gcin[1] = 2; else gcin[1] = 0;
    if ( gh[2] > sico.HFLOWMIN || sico.SURFACE == 0 || sico.SURFACE == 2 ) gcin[2] = 4; else gcin[2] = 0;
    if ( gh[3] > sico.HFLOWMIN || sico.SURFACE == 0 || sico.SURFACE == 2 ) gcin[3] = 1; else gcin[3] = 0;

    if ( gcin[0] == 0 || gcin[1] == 0 ) gcin[4] = 1; else gcin[4] = 2;
    if ( gcin[2] == 0 || gcin[3] == 0 ) gcin[5] = 1; else gcin[5] = 2;

    return gcin;
}

float fround ( float ginval, int gplaces ) { // function for rounding values to a given number of places

    float ground;
    ground = (float)(round( pow( 10, gplaces ) * ginval )) / (float)(pow( 10, gplaces ));
    return ground;
}

float fdiv ( float gabove, float gbelow, float gcriterion ) { // function to avoid division by zero

    float gdiv;
    if ( gbelow > gcriterion ) gdiv = gabove / gbelow; else gdiv = 0;
    return gdiv;
}

float fbeta ( float gelevm, float gelevp, float gspan, struct ico sico ) { // function for topographic slopes in x and y direction

    float gbeta;
    gbeta = atan(( gelevm - gelevp ) / ( gspan * sico.CSZ ));
    return gbeta;
}

float fbetaxy ( float gbetax, float gbetay ) { // function for topographic slope

    float gbetaxy;
    if ( gbetax != 0 || gbetay != 0 ) gbetaxy = atan( pow( tan ( gbetax) * tan ( gbetax) + tan( gbetay ) * tan( gbetay ), 0.5 ));
    else gbetaxy = 0;
    return gbetaxy;
}

float falpha ( float gbetax, float gbetay, struct ico sico ) { // function for topographic aspect

    float galpha = 0;

    if ( gbetax == 0 && gbetay == 0 ) galpha = 0;
    else if ( gbetax == 0 && gbetay > 0 ) galpha = sico.PI * 0.5;
    else if ( gbetax == 0 && gbetay < 0 ) galpha = 3 * sico.PI * 0.5;
    else if ( gbetax >= 0 ) galpha = atan( tan( gbetay ) / tan( gbetax ));
    else if ( gbetax < 0 ) galpha = atan( tan( gbetay ) / tan( gbetax )) + sico.PI;

    if ( galpha < 0 ) galpha += 2 * sico.PI;

    return galpha;
}

float falphav ( float gmomx, float gmomy, struct ico sico ) { // function for direction of movement

    float galphav = 0;

    if ( gmomx == 0 && gmomy == 0 ) galphav = 0;
    else if ( gmomx == 0 && gmomy > 0 ) galphav = sico.PI * 0.5;
    else if ( gmomx == 0 && gmomy < 0 ) galphav = 3 * sico.PI * 0.5;
    else if ( gmomx >= 0 ) galphav = atan( gmomy / gmomx );
    else if ( gmomx < 0 ) galphav = atan( gmomy / gmomx ) + sico.PI;

    if ( galphav < 0 ) galphav += 2 * sico.PI;

    return galphav;
}

float **finhyd ( char *gindir, char *gname, int ghydtmax, int ghydtmaxx ) { // function for input of hydrograph

    FILE *gfhyd;
    char gpath[1000], *gtesth, *gtesth2;
    int gi, gj;
    char *tok;
    char gdelim[] = " \t\n";
    char **strtodhlp=NULL;

    float **ghyd = alloc_dmatrix( ghydtmaxx+1, 14 );

    sprintf(gpath, "%s%s", gindir, gname );
    gfhyd=fopen(gpath, "r");

    gtesth = fcparam ( gfhyd );
    free( gtesth );
    for ( gi=0; gi<=ghydtmax; gi++ ) {
        gtesth = fcparam ( gfhyd );
        tok=strtok(gtesth, gdelim);
        gj = 0;
        while ( tok != NULL ) {
            gtesth2=tok;
            tok=strtok(NULL, gdelim);
            ghyd[gi][gj] = strtod (gtesth2, strtodhlp);
            gj += 1;
        }
        free(gtesth);
    }

    fclose(gfhyd);
    return ghyd;
}

int *fhydp ( int ghydx, int ghydy, float ghydl, int *gpx, int *gpy, float galpha, float ghydlmax, struct ico sico ) { // function for input hydrograph profile

    int gdistx, gdisty, gi, gj, gx0, gy0, gdprex, gdprey, gctrl, *gx, *gy;
    float gd, gdxtest, gdytest;

    int *ghydp = (int*) calloc( (int)(2*ghydlmax/sico.CSZ+2), sizeof(int) );

    gdistx = (int)( ghydl * fabs( cos( galpha )) / ( 2 * sico.CSZ ) + 0.5 );
    gdisty = (int)( ghydl * fabs( sin( galpha )) / ( 2 * sico.CSZ ) + 0.5 );

    gdxtest = fabs((float)( gdistx ) / cos( galpha ));
    gdytest = fabs((float)( gdisty ) / sin( galpha ));

    gx = (int*) calloc( 4 * (( gdistx + 1 ) * ( gdisty + 1 )), sizeof(int));
    gy = (int*) calloc( 4 * (( gdistx + 1 ) * ( gdisty + 1 )), sizeof(int));

    gctrl = 1;

    if ( galpha == 0 || galpha == sico.PI * 0.5 || gdytest == 0 ) {

        gx0 = ghydx - gdistx;
        gy0 = ghydy;

        gx[gctrl] = gx0;
        gy[gctrl] = gy0;

        while ( gx[gctrl] <= ghydx + gdistx ) {

            gctrl += 1;
            gx[gctrl] = gx0 + gctrl;
            gy[gctrl] = gy0;
        }
        
    } else if ( galpha == sico.PI || galpha == 3 * sico.PI * 0.5 || gdxtest == 0 ) {

        gx0 = ghydx;
        gy0 = ghydy - gdisty;

        gx[gctrl] = gx0;
        gy[gctrl] = gy0;

        while ( gy[gctrl] <= ghydy + gdisty ) {

            gctrl += 1;
            gx[gctrl] = gx0;
            gy[gctrl] = gy0 + gctrl;
        }

    } else {

        if ( gdxtest > gdytest ) {

            gx0 = ghydx - gdistx;

            if (( galpha > 0 && galpha < sico.PI * 0.5 ) || ( galpha > sico.PI && galpha < 3 * sico.PI * 0.5 ))
                gy0 = ghydy - gdisty;
            else gy0 = ghydy + gdisty;

            gdprex = 1;
            if (( galpha > 0 && galpha < sico.PI * 0.5 ) || ( galpha > sico.PI && galpha < 3 * sico.PI * 0.5 )) gdprey = 1;
            else gdprey = -1;
            
        } else {

            if (( galpha > 0 && galpha < sico.PI * 0.5 ) || ( galpha > sico.PI && galpha < 3 * sico.PI * 0.5 ))
                gx0 = ghydx - gdistx;
            else gx0 = ghydx + gdistx;
            gy0 = ghydy - gdisty;
            if (( galpha > 0 && galpha < sico.PI * 0.5 ) || ( galpha > sico.PI && galpha < 3 * sico.PI * 0.5 )) gdprex = 1;
            else gdprex = -1;
            gdprey = 1;
        }
        
        if ( fabs( 1 / cos( galpha )) < fabs(1 / sin ( galpha ))) gd = fabs( 1 / cos( galpha )); else gd = fabs(1 / sin ( galpha ));

        gx[gctrl] = gx0;
        gy[gctrl] = gy0;

        while ( gx[gctrl] <= ghydx + gdistx && gy[gctrl] <= ghydy + gdisty
            && gx[gctrl] >= ghydx - gdistx && gy[gctrl] >= ghydy - gdisty ) {

            gctrl += 1;
            gx[gctrl] = gx0 + gdprex * (int)( gctrl * gd * fabs(cos( galpha )) + 0.5 );
            gy[gctrl] = gy0 + gdprey * (int)( gctrl * gd * fabs(sin( galpha )) + 0.5 );
        }
    }
    
    gctrl -= 1;
    ghydp[0] = gctrl;

    for ( gi=1; gi<sico.IMAX; gi++ ) {
        for ( gj=1; gj<=gctrl; gj++ ) {
            if ( gpx[gi] == gx[gj] && gpy[gi] == gy[gj] ) {
                ghydp[gj] = gi;
            }
        }
    }

    free(gx); free(gy);

    return ghydp;
}

float fk ( float gh, float ghflow, float gvflowx, float gvflowy, float gdvflow, float gekin, int gp, struct ico sico, struct flow sflow ) {
    // function for earth pressure coefficients

    float ggk, gka, gkp, gdelta, gphi;
   
    gdelta = sflow.DELTA[gp];
    gphi = sflow.PHI[gp];

    if ( gphi < gdelta ) gphi = gdelta;
    if ( gphi == sico.UNDEF ) gphi = gdelta; // friction angles

    gka = 2 * ( 1 - pow( 1 - pow( cos( gphi ) / cos( gdelta ), 2 ), 0.5 )) / pow( cos( gphi ), 2 ) - 1; // active earth pressure coefficient
    gkp = 2 * ( 1 + pow( 1 - pow( cos( gphi ) / cos( gdelta ), 2 ), 0.5 )) / pow( cos( gphi ), 2 ) - 1; // passive eath pressure coefficient

    if ( gdvflow >= 0 ) ggk=gka;
    else ggk = ffmin(sflow.KPMAX[gp], gkp ); // earth pressure coefficient

    return ggk;
}

float *fvm ( float ghflow, float ghflow2, float ghflow3, float gvflowx, float gvflowx2, float gvflowx3, float gvflowy, float gvflowy2, float gvflowy3, 
    struct ico sico, struct flow sflow ) { // function for virtual mass

    int go;
    float *gvm, gh, galpha[3], ggamma[3], gu[3], gv[3];
    gvm = (float*) calloc( 21, sizeof(float));

    if ( sico.MODEL == 7 && sico.SLOMO <= 1.0 && sico.LAYERS == 0 ) {

        gu[0] = gvflowx; gu[1] = gvflowx2; gu[2] = gvflowx3;
        gv[0] = gvflowy; gv[1] = gvflowy2; gv[2] = gvflowy3;

        gh = ghflow + ghflow2 + ghflow3;
        if ( gh > sico.HFLOWMIN ) {

            if ( ghflow > sico.HFLOWMIN ) galpha[0] = ghflow / gh; else galpha[0] = 0;
            if ( ghflow2 > sico.HFLOWMIN ) galpha[1] = ghflow2 / gh; else galpha[1] = 0;
            if ( ghflow3 > sico.HFLOWMIN ) galpha[2] = ghflow3 / gh; else galpha[2] = 0;

        } else { galpha[0] = 0; galpha[1] = 0; galpha[2] = 0; }

        if ( galpha[0] > 0  && galpha[2] > 0 ) ggamma[0] = sflow.RHO3 / sflow.RHO1; else ggamma[0] = 0;
        if ( galpha[0] > 0  && galpha[1] > 0 ) ggamma[1] = sflow.RHO2 / sflow.RHO1; else ggamma[1] = 0;
        if ( galpha[1] > 0  && galpha[2] > 0 ) ggamma[2] = sflow.RHO3 / sflow.RHO2; else ggamma[2] = 0;

        if ( galpha[2] > 0 && galpha[0] > 0 )
            gvm[9] = ( sflow.VM_n0 * ( sflow.VM_l + pow( galpha[0], sflow.VM_n )) - 1 ) / ( galpha[0] / galpha[2] + ggamma[0] ); else gvm[9] = 0;
        if ( galpha[1] > 0 && galpha[0] > 0 )
            gvm[10] = ( sflow.VM_n0 * ( sflow.VM_l + pow( galpha[0], sflow.VM_n )) - 1 ) / ( galpha[0] / galpha[1] + ggamma[1] ); else gvm[10] = 0;
        if ( galpha[2] > 0 && galpha[1] > 0 )
            gvm[11] = ( sflow.VM_n0 * ( sflow.VM_l + pow( galpha[1], sflow.VM_n )) - 1 ) / ( galpha[1] / galpha[2] + ggamma[2] ); else gvm[11] = 0;

        if ( galpha[0] > 0 ) {

            gvm[0] = ggamma[0] * gvm[9] * ( gu[2] - gu[0] ) + ggamma[1] * gvm[10] * ( gu[1] - gu[0] );
            gvm[1] = ggamma[0] * gvm[9] * ( gu[2] * gu[2] - gu[0] * gu[0]) + ggamma[1] * gvm[10] * ( gu[1] * gu[1] - gu[0] * gu[0]);
            gvm[2] = ggamma[0] * gvm[9] * ( gu[2] * gv[2] - gu[0] * gv[0] ) + ggamma[1] * gvm[10] * ( gu[1] * gv[1] - gu[0] * gv[0] );

            if ( fabs( gvm[1] ) > sflow.VM_fact * fabs( gu[2] - gu[0] ) && fabs( gvm[1] ) > sflow.VM_fact * fabs( gu[1] - gu[0] ))
                gvm[1] = sflow.VM_fact * fsign( ffmax( fabs( gu[2] - gu[0] ), fabs( gu[1] - gu[0] )));
            if ( fabs( gvm[2] ) > sflow.VM_fact * fabs( gv[2] - gv[0] ) && fabs( gvm[2] ) > sflow.VM_fact * fabs( gv[1] - gv[0] ))
                gvm[2] = sflow.VM_fact * fsign( ffmax( fabs( gv[2] - gv[0] ), fabs( gv[1] - gv[0] )));

            gvm[12] = ggamma[0] * gvm[9] * ( gv[2] - gv[0] ) + ggamma[1] * gvm[10] * ( gv[1] - gv[0] );
            gvm[13] = ggamma[0] * gvm[9] * ( gv[2] * gv[2] - gv[0] * gv[0]) + ggamma[1] * gvm[10] * ( gv[1] * gv[1] - gv[0] * gv[0]);
            gvm[14] = ggamma[0] * gvm[9] * ( gu[2] * gv[2] - gu[0] * gv[0] ) + ggamma[1] * gvm[10] * ( gu[1] * gv[1] - gu[0] * gv[0] );

            if ( fabs( gvm[13] ) > sflow.VM_fact * fabs( gv[2] - gv[0] ) && fabs( gvm[13] ) > sflow.VM_fact * fabs( gv[1] - gv[0] ))
                gvm[13] = sflow.VM_fact * fsign( ffmax( fabs( gv[2] - gv[0] ), fabs( gv[1] - gv[0] )));
            if ( fabs( gvm[14] ) > sflow.VM_fact * fabs( gu[2] - gu[0] ) && fabs( gvm[14] ) > sflow.VM_fact * fabs( gu[1] - gu[0] ))
                gvm[14] = sflow.VM_fact * fsign( ffmax( fabs( gu[2] - gu[0] ), fabs( gu[1] - gu[0] )));

        }
        else { gvm[0] = 0; gvm[1] = 0; gvm[2] = 0; gvm[12] = 0; gvm[13] = 0; gvm[14] = 0; }

        if ( galpha[1] > 0 ) {

            gvm[3] = ggamma[2] * gvm[11] * ( gu[2] - gu[1] ) - galpha[0] / galpha[1] * gvm[10] * ( gu[1] - gu[0] );
            gvm[4] = ggamma[2] * gvm[11] * ( gu[2] * gu[2] - gu[1] * gu[1]) - galpha[0] / galpha[1] * gvm[10] * ( gu[1] * gu[1] - gu[0] * gu[0]);
            gvm[5] = ggamma[2] * gvm[11] * ( gu[2] * gv[2] - gu[1] * gv[1] ) - galpha[0] / galpha[1] * gvm[10] * ( gu[1] * gv[1] - gu[0] * gv[0] );

            if ( fabs( gvm[4] ) > sflow.VM_fact * fabs( gu[2] - gu[1] ) && fabs( gvm[4] ) > sflow.VM_fact * fabs( gu[0] - gu[1] ))
                gvm[4] = sflow.VM_fact * fsign( ffmax( fabs( gu[2] - gu[1] ), fabs( gu[0] - gu[1] )));
            if ( fabs( gvm[5] ) > sflow.VM_fact * fabs( gv[2] - gv[1] ) && fabs( gvm[5] ) > sflow.VM_fact * fabs( gv[0] - gv[1] ))
                gvm[5] = sflow.VM_fact * fsign( ffmax( fabs( gv[2] - gv[1] ), fabs( gv[0] - gv[1] )));

            gvm[15] = ggamma[2] * gvm[11] * ( gv[2] - gv[1] ) - galpha[0] / galpha[1] * gvm[10] * ( gv[1] - gv[0] );
            gvm[16] = ggamma[2] * gvm[11] * ( gv[2] * gv[2] - gv[1] * gv[1]) - galpha[0] / galpha[1] * gvm[10] * ( gv[1] * gv[1] - gv[0] * gv[0]);
            gvm[17] = ggamma[2] * gvm[11] * ( gu[2] * gv[2] - gu[1] * gv[1] ) - galpha[0] / galpha[1] * gvm[10] * ( gu[1] * gv[1] - gu[0] * gv[0] );

            if ( fabs( gvm[16] ) > sflow.VM_fact * fabs( gv[2] - gv[1] ) && fabs( gvm[16] ) > sflow.VM_fact * fabs( gv[0] - gv[1] ))
                gvm[16] = sflow.VM_fact * fsign( ffmax( fabs( gv[2] - gv[1] ), fabs( gv[0] - gv[1] )));
            if ( fabs( gvm[17] ) > sflow.VM_fact * fabs( gu[2] - gu[1] ) && fabs( gvm[17] ) > sflow.VM_fact * fabs( gu[0] - gu[1] ))
                gvm[17] = sflow.VM_fact * fsign( ffmax( fabs( gu[2] - gu[1] ), fabs( gu[0] - gu[1] )));

        } else { gvm[3] = 0; gvm[4] = 0; gvm[5] = 0; gvm[15] = 0; gvm[16] = 0; gvm[17] = 0; }

        if ( galpha[2] > 0 ) {

            gvm[6] = galpha[0] / galpha[2] * gvm[9] * ( gu[2] - gu[0] ) + galpha[1] / galpha[2] * gvm[11] * ( gu[2] - gu[1] );
            gvm[7] = galpha[0] / galpha[2] * gvm[9] * ( gu[2] * gu[2] - gu[0] * gu[0]) + galpha[1] / galpha[2] * gvm[11] * ( gu[2] * gu[2] - gu[1] * gu[1]);
            gvm[8] = galpha[0] / galpha[2] * gvm[9] * ( gu[2] * gv[2] - gu[0] * gv[0] ) + galpha[1] / galpha[2] * gvm[11] * ( gu[2] * gv[2] - gu[1] * gv[1] );

            if ( fabs( gvm[7] ) > sflow.VM_fact * fabs( gu[2] - gu[0] ) && fabs( gvm[7] ) > sflow.VM_fact * fabs( gu[2] - gu[1] ))
                gvm[7] = sflow.VM_fact * fsign( ffmax( fabs( gu[2] - gu[0] ), fabs( gu[2] - gu[1] )));
            if ( fabs( gvm[8] ) > sflow.VM_fact * fabs( gv[2] - gv[0] ) && fabs( gvm[8] ) > sflow.VM_fact * fabs( gv[2] - gv[1] ))
                gvm[8] = sflow.VM_fact * fsign( ffmax( fabs( gv[2] - gv[0] ), fabs( gv[2] - gv[1] )));

            gvm[18] = galpha[0] / galpha[2] * gvm[9] * ( gv[2] - gv[0] ) + galpha[1] / galpha[2] * gvm[11] * ( gv[2] - gv[1] );
            gvm[19] = galpha[0] / galpha[2] * gvm[9] * ( gv[2] * gv[2] - gv[0] * gv[0]) + galpha[1] / galpha[2] * gvm[11] * ( gv[2] * gv[2] - gv[1] * gv[1]);
            gvm[20] = galpha[0] / galpha[2] * gvm[9] * ( gu[2] * gv[2] - gu[0] * gv[0] ) + galpha[1] / galpha[2] * gvm[11] * ( gu[2] * gv[2] - gu[1] * gv[1] );

            if ( fabs( gvm[19] ) > sflow.VM_fact * fabs( gv[2] - gv[0] ) && fabs( gvm[19] ) > sflow.VM_fact * fabs( gv[2] - gv[1] ))
                gvm[19] = sflow.VM_fact * fsign( ffmax( fabs( gv[2] - gv[0] ), fabs( gv[2] - gv[1] )));
            if ( fabs( gvm[20] ) > sflow.VM_fact * fabs( gu[2] - gu[0] ) && fabs( gvm[20] ) > sflow.VM_fact * fabs( gu[2] - gu[1] ))
                gvm[20] = sflow.VM_fact * fsign( ffmax( fabs( gu[2] - gu[0] ), fabs( gu[2] - gu[1] )));
                
        } else { gvm[6] = 0; gvm[7] = 0; gvm[8] = 0; gvm[18] = 0; gvm[19] = 0; gvm[20] = 0; }
        
    } else { for ( go = 0; go < 21; go++ ) gvm[go] = 0; }

    return gvm;
}

float *fdrag ( float ghflow, float ghflow2, float ghflow3, struct ico sico, struct flow sflow ) { // function for drag coefficient in mixture

    float *gcdrag, gh, galpha[3], galphac[3], ggamma[3], gf, gg, gp, gsp;
    gcdrag = (float*) calloc( 3, sizeof(float));

    if ( sico.MODEL == 7 && sico.LAYERS == 0 && sflow.DRAG_k >= 0 ) {

        gh = ghflow + ghflow2 + ghflow3;
        if ( gh > sico.HFLOWMIN ) {

            if ( ghflow > sico.HFLOWMIN ) galpha[0] = ghflow / gh; else galpha[0] = 0;
            if ( ghflow2 > sico.HFLOWMIN ) galpha[1] = ghflow2 / gh; else galpha[1] = 0;
            if ( ghflow3 > sico.HFLOWMIN ) galpha[2] = ghflow3 / gh; else galpha[2] = 0;

        } else { galpha[0] = 0; galpha[1] = 0; galpha[2] = 0; }

        if ( galpha[0] > 0 && galpha[2] > 0 ) ggamma[0] = sflow.RHO3 / sflow.RHO1; else ggamma[0] = 0;
        if ( galpha[0] > 0 && galpha[1] > 0 ) ggamma[1] = sflow.RHO2 / sflow.RHO1; else ggamma[1] = 0;
        if ( galpha[1] > 0 && galpha[2] > 0 ) ggamma[2] = sflow.RHO3 / sflow.RHO2; else ggamma[2] = 0;

        galphac[0] = galpha[0] * galpha[2];
        galphac[1] = galpha[0] * galpha[1];
        galphac[2] = galpha[1] * galpha[2];

        if ( galpha[0] > 0 && galpha[2] > 0 ) {

            gf = ggamma[0] / 180 * pow( galpha[2] / galpha[0], 3 ) * sflow.DRAG_rep;
            gg = pow( galpha[2], sflow.DRAG_m - 1 );
            gp = pow( galpha[0], sflow.DRAG_n );
            gsp = ( gp / galpha[0] + ( 1 - gp ) / galpha[2] ) * sflow.DRAG_k;
            gcdrag[0] = pow( 1, sflow.DRAG_je )
                * galphac[0] * ( 1 - ggamma[0] ) * sico.GRAVITY / pow( sflow.DRAG_vterm * ( gp * gf + ( 1 - gp ) * gg ) + gsp, sflow.DRAG_je );

        } else gcdrag[0] = 0;

        if ( galpha[0] > 0 && galpha[1] > 0 ) {

            gf = ggamma[0] / 180 * pow( galpha[1] / galpha[0], 3 ) * sflow.DRAG_rep;
            gg = pow( galpha[1], sflow.DRAG_m - 1 );
            gp = pow( galpha[0], sflow.DRAG_n );
            gsp = ( gp / galpha[0] + ( 1 - gp ) / galpha[1] ) * sflow.DRAG_k;
            gcdrag[1] = pow( 1, sflow.DRAG_je )
                * galphac[1] * ( 1 - ggamma[1] ) * sico.GRAVITY / pow( sflow.DRAG_vterm * ( gp * gf + ( 1 - gp ) * gg ) + gsp, sflow.DRAG_je );

        } else gcdrag[1] = 0;

        if ( galpha[1] > 0 && galpha[2] > 0 ) {

            gf = ggamma[1] / 180 * pow( galpha[2] / galpha[1], 3 ) * sflow.DRAG_rep;
            gg = pow( galpha[2], sflow.DRAG_m - 1 );
            gp = pow( galpha[1], sflow.DRAG_n );
            gsp = ( gp / galpha[1] + ( 1 - gp ) / galpha[2] ) * sflow.DRAG_k;
            gcdrag[2] = pow( 1, sflow.DRAG_je )
                * galphac[2] * ( 1 - ggamma[2] ) * sico.GRAVITY / pow( sflow.DRAG_vterm * ( gp * gf + ( 1 - gp ) * gg ) + gsp, sflow.DRAG_je );

        } else gcdrag[2] = 0;

    } else { gcdrag[0] = 0; gcdrag[1] = 0; gcdrag[2] = 0; }

    return gcdrag;
}

float *fgze ( float ghflow, float ghflow2, float ghflow3, float gvflowx, float gvflowx2, float gvflowx3, float gvflowy, float gvflowy2, float gvflowy3, 
    float *gnh, float gghx, float gghy, float *ggux, float *ggvy, float *gguy, float *ggvx, float *gnuss, float *gnvss, float *gnufs, float *gnvfs, float *gnuff, float *gnvff, 
    float *ggalpha10, float *ggalpha20, float *ggalpha30, float gdx, float gdy, float *ggrav, float gbetax, float gbetay, float *gcdrag, float gtsum, 
    struct ico sico, struct flow sflow ) { // function for enhanced gravity

    int gb, gj, gl, gm, go;
    float *ggze, gh, gu[3], gv[3], gnu[3][9], gnv[3][9], ggalpha[3][9], gw[3], gww[3], guvw[3], galpha[3], ggamma[3], ggamma2[3], gdrag[3], gght, ggtlength[3], 
        gtestx[3], gtesty[3], gtestxy[3];
        
    ggze = (float*) calloc( 11, sizeof(float));

    if ( sico.MODEL == 7 && sico.SLOMO <= 1.0 && sico.LAYERS == 0 ) { // non-layered multi-phase model for extremely rapid flows (gravity enhanced for buoyancy)

        for ( gl = 0; gl < sico.PMAX; gl++ ) {

            gtestx[gl] = ggrav[0]; gtesty[gl] = ggrav[1]; gtestxy[gl] = ggrav[4];
        }
   
        gh = ghflow + ghflow2 + ghflow3; // total flow height

        if ( gh > sico.HFLOWMIN ) {

            if ( ghflow > sico.HFLOWMIN ) galpha[0] = ghflow / gh; else galpha[0] = 0; // phase fractions
            if ( ghflow2 > sico.HFLOWMIN ) galpha[1] = ghflow2 / gh; else galpha[1] = 0;
            if ( ghflow3 > sico.HFLOWMIN ) galpha[2] = ghflow3 / gh; else galpha[2] = 0;

        } else { galpha[0] = 0; galpha[1] = 0; galpha[2] = 0; }

        if ( galpha[0] > 0 && galpha[2] > 0 && sflow.PHI[2] == 0 ) ggamma[0] = sflow.RHO3 / sflow.RHO1; else ggamma[0] = 0; // density ratios
        if ( galpha[0] > 0 && galpha[1] > 0 && sflow.PHI[1] == 0 ) ggamma[1] = sflow.RHO2 / sflow.RHO1; else ggamma[1] = 0;
        if ( galpha[1] > 0 && galpha[2] > 0 && sflow.PHI[2] == 0 ) ggamma[2] = sflow.RHO3 / sflow.RHO2; else ggamma[2] = 0;

        if ( galpha[0] > 0 && galpha[2] > 0 && sflow.PHI[2] == 0 ) ggamma2[0] = sflow.RHO3 / sflow.RHO1; else ggamma2[0] = 1;
        if ( galpha[0] > 0 && galpha[1] > 0 && sflow.PHI[1] == 0 ) ggamma2[1] = sflow.RHO2 / sflow.RHO1; else ggamma2[1] = 1;
        if ( galpha[1] > 0 && galpha[2] > 0 && sflow.PHI[2] == 0 ) ggamma2[2] = sflow.RHO3 / sflow.RHO2; else ggamma2[2] = 1;

        ggze[0] = gtestx[0] * ( 1 - ggamma[0] ); // gz in x direction, PHASE 1
        ggze[1] = gtestx[1] * ggamma2[2]; // gz in x direction, PHASE 2        
        ggze[2] = gtestx[2]; // gz in x direction, PHASE 3        

        ggze[4] = gtesty[0] * ( 1 - ggamma[0] ); // gz in y direction, PHASE 1        
        ggze[5] = gtesty[1] * ggamma2[2]; // gz in y direction, PHASE 2
        ggze[6] = gtesty[2]; // gz in y direction, PHASE 3
        
        ggze[3] = gtestx[0] * ggamma2[0]; // gz in x direction, solid PHASE 1 (for dh/dx)
        ggze[7] = gtesty[0] * ggamma2[0]; // gz in y direction, solid PHASE 1 (for dh/dy)

        ggze[8] = gtestxy[0] * ( 1 - ggamma[0] ); // gz, PHASE 1
        ggze[9] = gtestxy[1] * ggamma2[2]; // gz, PHASE 2        
        ggze[10] = gtestxy[2]; // gz, PHASE 3        
        
        if ( sico.NONHYDRO > 1 ) { // advanced enhanced gravity model !!! CHECK for topographic adjustment
        
            for( gj=0; gj<9; gj++ ) { 

                gnu[0][gj] = gnuss[gj]; gnu[1][gj] = gnufs[gj]; gnu[2][gj] = gnuff[gj]; gnv[0][gj] = gnvss[gj]; gnv[1][gj] = gnvfs[gj]; gnv[2][gj] = gnvff[gj];
                ggalpha[0][gj] = ggalpha10[gj]; ggalpha[1][gj] = ggalpha20[gj]; ggalpha[2][gj] = ggalpha30[gj];
            }

            gu[0] = gvflowx; gu[1] = gvflowx2; gu[2] = gvflowx3;
            gv[0] = gvflowy; gv[1] = gvflowy2; gv[2] = gvflowy3;

            gght = -0.5 / gdx * ( gnh[2] * ( ggalpha[0][2] * gnu[0][2] + ggalpha[1][2] * gnu[1][2] + ggalpha[2][2] * gnu[2][2] ) - gnh[5] * ( ggalpha[0][5] * gnu[0][5] 
                + ggalpha[1][5] * gnu[1][5] + ggalpha[2][5] * gnu[2][5] )) - 0.5 / gdy * ( gnh[1] * ( ggalpha[0][1] * gnv[0][1] + ggalpha[1][1] * gnv[1][1] 
                + ggalpha[2][1] * gnv[2][1] ) - gnh[4] * ( ggalpha[0][4] * gnv[0][4] + ggalpha[1][4] * gnv[1][4] + ggalpha[2][4] * gnv[2][4] ));

            for ( go = 0; go < 3; go++ ) {

                ggtlength[go] = -0.5 * ( -gh * ( pow(ggux[go], 2 ) + 2 * ggvx[go] * gguy[go] + pow( ggvy[go], 2 )) 
                    + gght * ( ggux[go] + ggvy[go] ) + ( gu[go] * gghx + gv[go] * gghy ) * ( ggux[go] + ggvy[go] ));

                gw[go] = gu[go] * tan(gbetax) + gv[go] * tan(gbetay);
                gww[go] = (gu[go] * tan(gbetax) + gv[go] * tan( gbetay )) - ( ggux[go] + ggvy[go] ) * gh * 0.5;
                guvw[go] = pow( pow( gu[go], 2 ) + pow( gv[go], 2 ) + pow( gw[go], 2 ), 0.5 );
            }

            gdrag[0] = gcdrag[0] * pow( fabs( guvw[2] - guvw[0] ), sflow.DRAG_je - 1 );
            gdrag[1] = gcdrag[1] * pow( fabs( guvw[1] - guvw[0] ), sflow.DRAG_je - 1 );
            gdrag[2] = gcdrag[2] * pow( fabs( guvw[2] - guvw[1] ), sflow.DRAG_je - 1 );
            
            for ( gm=0; gm<2; gm++ ) {

                gb = 4 * gm;

                if ( galpha[0] > 0 ) ggze[gb+0] = ggze[gb+0] + ggtlength[0] - 1 / galpha[0] * ( gdrag[0] * ( gww[2] - gww[0] ) + gdrag[1] * ( gww[1] - gww[0] ));
                if ( galpha[1] > 0 ) ggze[gb+1] = ggze[gb+1] + ggtlength[1] - 1 / galpha[1] * ( -1 / ggamma2[1] * gdrag[1] * ( gww[1] - gww[0] ) 
                    + gdrag[2] * ( gww[2] - gww[1] ));
                if ( galpha[2] > 0 ) ggze[gb+2] = ggze[gb+2] + ggtlength[2] + 1 / galpha[2] * ( 1 / ggamma2[0] * gdrag[0] * ( gww[2] - gww[0] ) 
                    + 1 / ggamma2[2] * gdrag[2] * ( gww[2] - gww[1] )); // !!! CHECK for x, y, xy directions
            }
        }
        
    } else if ( sico.LAYERS != 0 ) { // layered model:

        ggze[0] = ggrav[0]; ggze[1] = ggrav[15]; ggze[2] = ggrav[17]; ggze[3] = ggrav[0];
        ggze[4] = ggrav[1]; ggze[5] = ggrav[16]; ggze[6] = ggrav[18]; ggze[7] = ggrav[1];
        ggze[8] = ggrav[4]; ggze[9] = ggrav[19]; ggze[10] = ggrav[20];

    } else { // all other models:
        
        ggze[0] = ggrav[0]; ggze[1] = ggrav[0]; ggze[2] = ggrav[0]; ggze[3] = ggrav[0];
        ggze[4] = ggrav[1]; ggze[5] = ggrav[1]; ggze[6] = ggrav[1]; ggze[7] = ggrav[1];
        ggze[8] = ggrav[4]; ggze[9] = ggrav[4]; ggze[10] = ggrav[4];
    }

    return ggze;
}

float *fdisp ( float ghflow, float ghflow2, float ghflow3, float gvflowx, float gvflowx2, float gvflowx3, float gvflowy, float gvflowy2, float gvflowy3, 
    float *ggux, float *ggvy, float *gguy, float *ggvx, float gbetax, float gbetay, float gbetaxy, float *gkx, float *gky, float *gvm, float *gcdrag, float gtsum,
    struct ico sico, struct flow sflow ) { // function for dispersion

    int go, gp;
    float *gdisp, ggtlength[3], gh, galpha[3], guu[3], gu[3], gv[3], ghnorm;
    gdisp = (float*) calloc( 6, sizeof(float));

    if ( sico.NONHYDRO == 1 || sico.NONHYDRO == 3 ) {

        gh = ghflow + ghflow2 + ghflow3;
        ghnorm = gh * cos( gbetaxy ); ghnorm = gh * cos( gbetaxy );

        gu[0] = gvflowx / cos( gbetaxy ); gu[1] = gvflowx2 / cos( gbetaxy ); gu[2] = gvflowx3 / cos( gbetaxy );
        gv[0] = gvflowy / cos( gbetaxy ); gv[1] = gvflowy2 / cos( gbetaxy ); gv[2] = gvflowy3 / cos( gbetaxy );
       
        if ( gh > sico.HFLOWMIN ) {

            if ( ghflow > sico.HFLOWMIN ) galpha[0] = ghflow / gh; else galpha[0] = 0;
            if ( ghflow2 > sico.HFLOWMIN ) galpha[1] = ghflow2 / gh; else galpha[1] = 0;
            if ( ghflow3 > sico.HFLOWMIN ) galpha[2] = ghflow3 / gh; else galpha[2] = 0;

        } else { galpha[0] = 0; galpha[1] = 0; galpha[2] = 0; }

        for ( go = 0; go < 3; go++ ) {

            guu[go] = fabs( sflow.disp_MULT ) * ( ggux[go] + ggvy[go] );
            ggtlength[go] = 2 * pow( guu[go], 2 ) - 2 * ggux[go] * ggvy[go] + 2 * ggvx[go] * gguy[go] - 2 * 0.0019 * pow( pow( gu[go], 2 ) + pow( gv[go], 2 ), 0.5 ) * guu[go];
        }

        if ( galpha[0] > 0 ) gdisp[0] = gkx[0] * ghnorm * ghnorm / 12 * ggtlength[0] * cos( gbetax );
        else gdisp[0] = 0;

        if ( galpha[1] > 0 ) gdisp[1] = ghnorm * ghnorm / 12 * ggtlength[1] * cos( gbetax );
        else gdisp[1] = 0;

        if ( galpha[2] > 0 ) gdisp[2] = ghnorm * ghnorm / 12 * ggtlength[2] * cos( gbetax );
        else gdisp[2] = 0;

        if ( galpha[0] > 0 ) gdisp[3] = gky[0] * ghnorm * ghnorm / 12 * ggtlength[0] * cos( gbetay );
        else gdisp[3] = 0;

        if ( galpha[1] > 0 ) gdisp[4] = ghnorm * ghnorm / 12 * ggtlength[1] * cos( gbetay );
        else gdisp[4] = 0;

        if ( galpha[2] > 0 ) gdisp[5] = ghnorm * ghnorm / 12 * ggtlength[2] * cos( gbetay );
        else gdisp[5] = 0;

        for ( gp=0; gp<6; gp++ ) gdisp[gp] *= fabs( sflow.disp_MULT );
        
    } else { 
    
        for ( gp = 0; gp < 6; gp++ ) gdisp[gp] = 0;
    }

    return gdisp;
}

float *ff ( float *ggh, float ghflow, float ghflow2, float ghflow3, float gvflowx, float gvflowx2, float gvflowx3, float gvflowy,
    float gvflowy2, float gvflowy3, float *gkx, float *ggz, float *gvm, float gghflowx, float ggalphax, float *gcdrag, float *ggrav, float *gdisp, 
    float gbetax, float gbetay, float gbetaxy, float gzfrag, float gtemp, struct flow sflow, struct ico sico ) { // function for fluxes in x direction

    int gl, gt = 0;
    float *ggf, gh, gu[3], gv[3], galpha[3], ggamma, glambdas, glambdaf, gseprate_x, gsepflux_x[3], ggtest, ghnorm, ghnorm1, ghnorm3, gux;
    ggf = (float*) calloc( sico.NVECTMIN, sizeof(float));

    gh = ghflow + ghflow2 + ghflow3;
    gu[0] = gvflowx; gu[1] = gvflowx2; gu[2] = gvflowx3;
    gv[0] = gvflowy; gv[1] = gvflowy2; gv[2] = gvflowy3;

    if ( gh > sico.HFLOWMIN ) {

        if ( ghflow > sico.HFLOWMIN ) galpha[0] = ghflow / gh; else galpha[0] = 0;
        if ( ghflow2 > sico.HFLOWMIN ) galpha[1] = ghflow2 / gh; else galpha[1] = 0;
        if ( ghflow3 > sico.HFLOWMIN ) galpha[2] = ghflow3 / gh; else galpha[2] = 0;

    } else { galpha[0] = 0; galpha[1] = 0; galpha[2] = 0; }


    // Separation fluxes

    if ( sico.MODEL == 7 && sico.LAYERS == 0 ) {

        if ( sico.SEPFLUX == 1 && ghflow > sico.HFLOWMIN && ghflow3 > sico.HFLOWMIN ) {

            if ( ggh[2] > sico.HFLOWMIN && ggh[5] > sico.HFLOWMIN ) ggalphax = ggalphax; else ggalphax = 0;

            ghnorm = gh * cos( gbetaxy ); ghnorm1 = ghflow * cos( gbetaxy ); ghnorm3 = ghflow3 * cos( gbetaxy );

            ggamma = sflow.RHO3 / sflow.RHO1;

            glambdas = ggamma / ( galpha[0] + ggamma * ( 1 - galpha[0] ));
            glambdaf = 1 / ( galpha[0] + ggamma * ( 1 - galpha[0] ));
   
            if ( gcdrag[0] > 0 ) {

                gseprate_x = -1 / sflow.RHO1 / gcdrag[0] * ( 
                    -fsign( gvflowx ) * tan( sflow.DELTA[0] ) * ( 1 - ggamma ) * galpha[0] * sflow.RHO1 * ggrav[4]
                    + 0.5 * gkx[0] * ( 1 - ggamma ) * sflow.RHO1 * ggrav[4] * ghnorm * ggalphax 
                    + gkx[0] * ( 1 - ggamma ) * sflow.RHO1 * ggrav[4] * galpha[0] * gghflowx * cos( gbetaxy )
                    + ggamma * sflow.RHO1 * ggrav[4] * galpha[0] * gghflowx - sflow.RHO1 * ggrav[4] * galpha[0] * tan( gbetax )
                    + sflow.sep_SHEAR * galpha[0] * sflow.NY[0] * sflow.RHO1 * gvflowx / cos( gbetaxy ) / ( ghnorm1 * ghnorm1));
                        // separation rate in x direction !!! CHECK for topographic adjustment
                
            } else gseprate_x = 0;

            gsepflux_x[0] = glambdas * gseprate_x * galpha[0] * ( 1 - galpha[0] ) * ( ghnorm1 + ghnorm3 ); // solid separation flux in x direction
            gsepflux_x[1] = 0; 
            gsepflux_x[2] = glambdaf * gseprate_x * galpha[0] * ( 1 - galpha[0] ) * ( ghnorm1 + ghnorm3 ) * -1; // fluid separation flux in x direction

        } else { gsepflux_x[0] = 0; gsepflux_x[1] = 0; gsepflux_x[2] = 0; }

    } else { gsepflux_x[0] = 0; gsepflux_x[1] = 0; gsepflux_x[2] = 0; }
    
    for ( gl = 0; gl < sico.PMAX; gl++ ) {


        // Flux terms

        if ( galpha[gl] > 0 ) {

            ggf[3*gl] = galpha[gl] * gh * gu[gl] + gsepflux_x[gl];
            ggf[3*gl+1] = galpha[gl] * gh * ( gu[gl] * gu[gl] - gvm[3*gl+1] + gdisp[gl] );

            if ( gdisp[0] != 0 ) {
            
                ggtest = galpha[gl] * gh * ( gu[gl] * gu[gl] - gvm[3*gl+1] );
                if ( fsign( ggtest) != fsign( ggf[3*gl+1] )) ggf[3*gl+1] = 0;
                else if ( fabs( ggf[3*gl+1] ) > sflow.CCONST * fabs( ggtest )) ggf[3*gl+1] = sflow.CCONST * ggtest;
            }

            ggf[3*gl+2] = galpha[gl] * gh * ( gu[gl] * gv[gl] - gvm[3*gl+2] );

        } else { ggf[3*gl] = 0; ggf[3*gl+1] = 0; ggf[3*gl+2] = 0; }
    }

    if ( sico.MODEL <= 3 ) gt=3; else gt=9;

    if ( sico.TEMP > 0 ) { // temperature term

        gtemp = 0;
        gux = gu[0] * galpha[0] + gu[1] * galpha[1] + gu[2] * galpha[2];

        ggf[gt] = gh * gtemp * gux;

    } else ggf[gt] = 0; 

    return ggf;
}

float *fg ( float *ggh, float ghflow, float ghflow2, float ghflow3, float gvflowx, float gvflowx2, float gvflowx3, float gvflowy, 
    float gvflowy2, float gvflowy3, float *gky, float *ggz, float *gvm, float gghflowy, float ggalphay, float *gcdrag, float *ggrav, float *gdisp, 
    float gbetax, float gbetay, float gbetaxy, float gzfrag, float gtemp, struct flow sflow, struct ico sico ) { // function for fluxes in y direction

    int gl, gt = 0;
    float *ggg, gh, gu[3], gv[3], galpha[3], ggamma, glambdas, glambdaf, gseprate_y, gsepflux_y[3], ggtest, ghnorm, ghnorm1, ghnorm3, gvx;
    ggg = (float*) calloc( sico.NVECTMIN, sizeof(float));

    gh = ghflow + ghflow2 + ghflow3;
    gu[0] = gvflowx; gu[1] = gvflowx2; gu[2] = gvflowx3;
    gv[0] = gvflowy; gv[1] = gvflowy2; gv[2] = gvflowy3;

    if ( gh > sico.HFLOWMIN ) {

        if ( ghflow > sico.HFLOWMIN ) galpha[0] = ghflow / gh; else galpha[0] = 0;
        if ( ghflow2 > sico.HFLOWMIN ) galpha[1] = ghflow2 / gh; else galpha[1] = 0;
        if ( ghflow3 > sico.HFLOWMIN ) galpha[2] = ghflow3 / gh; else galpha[2] = 0;

    } else { galpha[0] = 0; galpha[1] = 0; galpha[2] = 0; }


    // Separation fluxes

    if ( sico.MODEL == 7 && sico.LAYERS == 0 ) {

        if ( sico.SEPFLUX == 1 && ghflow > sico.HFLOWMIN && ghflow3 > sico.HFLOWMIN ) {

            ghnorm = gh * cos( gbetaxy ); ghnorm1 = ghflow * cos( gbetaxy ); ghnorm3 = ghflow3 * cos( gbetaxy );

            if ( ggh[1] > sico.HFLOWMIN && ggh[4] > sico.HFLOWMIN ) ggalphay = ggalphay; else ggalphay = 0;

            ggamma = sflow.RHO3 / sflow.RHO1;

            glambdas = ggamma / ( galpha[0] + ggamma * ( 1 - galpha[0] ));
            glambdaf = 1 / ( galpha[0] + ggamma * ( 1 - galpha[0] ));
      
            if ( gcdrag[0] > 0 ) {
         
                gseprate_y = -1 / sflow.RHO1 / gcdrag[0] * ( 
                    -fsign( gvflowy ) * tan( sflow.DELTA[0] ) * ( 1 - ggamma ) * galpha[0] * sflow.RHO1 * ggrav[4]
                    + 0.5 * gky[0] * ( 1 - ggamma ) * sflow.RHO1 * ggrav[4] * ghnorm * ggalphay 
                    + gky[0] * ( 1 - ggamma ) * sflow.RHO1 * ggrav[4] * galpha[0] * gghflowy * cos( gbetaxy )
                    + ggamma * sflow.RHO1 * ggrav[4] * galpha[0] * gghflowy - sflow.RHO1 * ggrav[4] * galpha[0] * tan( gbetay )
                    + sflow.sep_SHEAR * galpha[0] * sflow.NY[0] * sflow.RHO1 * gvflowy / cos( gbetaxy ) / ( ghnorm1 * ghnorm1 ));
                        // separation rate in y direction !!! CHECK for topographic adjustment
                
            } else  gseprate_y = 0;

            gsepflux_y[0] = glambdas * gseprate_y * galpha[0] * ( 1 - galpha[0] ) * ( ghnorm + ghnorm3 ); // solid separation flux in y direction
            gsepflux_y[1] = 0; 
            gsepflux_y[2] = glambdaf * gseprate_y * galpha[0] * ( 1 - galpha[0] ) * ( ghnorm + ghnorm3 ) * -1; // fluid separation flux in y direction

        } else { gsepflux_y[0] = 0; gsepflux_y[1] = 0; gsepflux_y[2] = 0; }

    } else { gsepflux_y[0] = 0; gsepflux_y[1] = 0; gsepflux_y[2] = 0; }
    
    for ( gl = 0; gl < sico.PMAX; gl++ ) {


        // Flux terms

        if ( galpha[gl] > 0 ) {

            ggg[3*gl] = galpha[gl] * gh * gv[gl] + gsepflux_y[gl];
            ggg[3*gl+1] = galpha[gl] * gh * ( gu[gl] * gv[gl] - gvm[3*gl+14] );
            ggg[3*gl+2] = galpha[gl] * gh * ( gv[gl] * gv[gl] - gvm[3*gl+13] + gdisp[3] );
            
            if ( gdisp[3] != 0 ) {
            
                ggtest = galpha[gl] * gh * ( gv[gl] * gv[gl] - gvm[3*gl+13] );
                if ( fsign( ggtest) != fsign( ggg[3*gl+2] )) ggg[3*gl+2] = 0;
                else if ( fabs( ggg[3*gl+2] ) > sflow.CCONST * fabs( ggtest )) ggg[3*gl+2] = sflow.CCONST * ggtest;
            }
            
        } else { ggg[3*gl] = 0; ggg[3*gl+1] = 0; ggg[3*gl+2] = 0; }
    }

    if ( sico.MODEL <= 3 ) gt=3; else gt=9;

    if ( sico.TEMP > 0 ) { // temperature term

        gtemp = 0;
        gvx = gv[0] * galpha[0] + gv[1] * galpha[1] + gv[2] * galpha[2];

        ggg[gt] = gh * gtemp * gvx;

    } else ggg[gt] = 0; 

    return ggg;
}

float *fs ( float *gh, float ghflow, float ghflow2, float ghflow3, float gvflowx, float gvflowx2, float gvflowx3, float gvflowy, float gvflowy2, float gvflowy3, 
    float gghx, float gghy, float gghx1, float gghy1, float gghx2, float gghy2, float ggalphassx, float ggalphassy, float ggalphafsx, float ggalphafsy, 
    float ggalphaffx, float ggalphaffy, float *ggrav, float *ggz, float gbetax, float gbetay, float gbetaxy, float *gcdrag, 
    float *gkx, float *gky, int gcslide, float *gangx, float *gangy, float gtsum, float gtslide, float gzfrag, float gtemp, struct ico sico, struct flow sflow ) {
    // function for source terms (accelerating components and slow flow)

    int gl, go, gctrl_deform, gt = 0;
    float *ggs, gu[3], gv[3], gw[3], galpha[3], ggalphax[3], ggalphay[3], ggamma0, ggamma1, ggamma2, ggamma3, ggamma4, ggamma5, ghnorm,
        ghnorm1s, ghnorm2s, ghnorm3s, guvw[3], gcdragx[3], gcdragy[3], gdragx[3], gdragy[3], gcompx[3], gcompy[3], gslideterm, gslidefractb, gghxx[3], gghyy[3], 
        glocx[3], glocy[3], glocz[3], gpbetax, gpbetay, gfact_slidex, gfact_slidey, gpx, gpy, galphari, gmyri, grhori, grho, guvwz, gub, gvb, gkappa, gc, gtau, 
        ggtempx2, ggtempy2, ggkappax, ggkappay, ggtempx, ggtempy, rmelt;
        
    ggs = (float*) calloc( sico.NVECTMIN, sizeof(float));

    if ( sico.MODEL <= 3 ) gt=3; else gt=9;
    if ( sico.MODEL <= 3 ) { gghx1 = gghx; gghy1 = gghy; }

    gu[0] = gvflowx / cos( gbetaxy ); gu[1] = gvflowx2 / cos( gbetaxy ); gu[2] = gvflowx3 / cos( gbetaxy );
    gv[0] = gvflowy / cos( gbetaxy ); gv[1] = gvflowy2 / cos( gbetaxy ); gv[2] = gvflowy3 / cos( gbetaxy );

    if ( sico.NONHYDRO > 0 ) { // non-hydrostatic effects

        gw[0] = gu[0] * tan(gbetax) + gv[0] * tan(gbetay);
        gw[1] = gu[1] * tan(gbetax) + gv[1] * tan(gbetay);
        gw[2] = gu[2] * tan(gbetax) + gv[2] * tan(gbetay);

    } else { gw[0] = 0; gw[1] = 0; gw[2] = 0; }

    if ( gh[0] > sico.HFLOWMIN ) {

        if ( ghflow > sico.HFLOWMIN ) galpha[0] = ghflow / gh[0]; else galpha[0] = 0;
        if ( ghflow2 > sico.HFLOWMIN ) galpha[1] = ghflow2 / gh[0]; else galpha[1] = 0;
        if ( ghflow3 > sico.HFLOWMIN ) galpha[2] = ghflow3 / gh[0]; else galpha[2] = 0;

    } else { galpha[0] = 0; galpha[1] = 0; galpha[2] = 0; }

    ghnorm = gh[0] * cos( gbetaxy );
    ghnorm1s = ghflow * cos( gbetaxy );
    ghnorm2s = ghflow2 * cos( gbetaxy );
    ghnorm3s = ghflow3 * cos( gbetaxy );
    
    gghx *= cos( gbetaxy ); gghy *= cos( gbetaxy );
    gghx1 *= cos( gbetaxy ); gghy1 *= cos( gbetaxy );
    gghx2 *= cos( gbetaxy ); gghy2 *= cos( gbetaxy );

    if ( galpha[0] > 0 && galpha[2] > 0 ) ggamma0 = sflow.RHO3 / sflow.RHO1; else ggamma0 = 1;
    if ( galpha[0] > 0 && galpha[1] > 0 ) ggamma1 = sflow.RHO2 / sflow.RHO1; else ggamma1 = 1;
    if ( galpha[1] > 0 && galpha[2] > 0 ) ggamma2 = sflow.RHO3 / sflow.RHO2; else ggamma2 = 1;

    if ( galpha[0] > 0 && galpha[2] > 0 && sflow.PHI[2] == 0 ) ggamma3 = sflow.RHO3 / sflow.RHO1; else ggamma3 = 0;
    if ( galpha[1] > 0 && galpha[2] > 0 && sflow.PHI[2] == 0 ) ggamma4 = sflow.RHO3 / sflow.RHO2; else ggamma4 = 0;
    ggamma5 = 1;

    for ( go = 0; go < 3; go++ ) guvw[go] = pow( pow( gu[go], 2 ) + pow( gv[go], 2 ) + pow( gw[go], 2 ), 0.5 );

    if ( gh[2] > sico.HFLOWMIN && gh[5] > sico.HFLOWMIN ) ggalphax[0] = ggalphassx; else ggalphax[0] = 0;
    if ( gh[2] > sico.HFLOWMIN && gh[5] > sico.HFLOWMIN ) ggalphax[1] = ggalphafsx; else ggalphax[1] = 0;
    if ( gh[2] > sico.HFLOWMIN && gh[5] > sico.HFLOWMIN ) ggalphax[2] = ggalphaffx; else ggalphax[2] = 0;
    if ( gh[1] > sico.HFLOWMIN && gh[4] > sico.HFLOWMIN ) ggalphay[0] = ggalphassy; else ggalphay[0] = 0;
    if ( gh[1] > sico.HFLOWMIN && gh[4] > sico.HFLOWMIN ) ggalphay[1] = ggalphafsy; else ggalphay[1] = 0;
    if ( gh[1] > sico.HFLOWMIN && gh[4] > sico.HFLOWMIN ) ggalphay[2] = ggalphaffy; else ggalphay[2] = 0;

    if ( sico.LAYERS > 0 ) {

        gghxx[0] = gghx1; gghxx[1] = gghx2 - gghx1; gghxx[2] = gghx - gghx2;
        gghyy[0] = gghy1; gghyy[1] = gghy2 - gghy1; gghyy[2] = gghy - gghy2;

        glocx[0] = ggrav[2]; glocx[1] = ggrav[7]; glocx[2] = ggrav[11];
        glocy[0] = ggrav[3]; glocy[1] = ggrav[8]; glocy[2] = ggrav[12];

    } else {

        gghxx[0] = gghx; gghxx[1] = gghx; gghxx[2] = gghx;
        gghyy[0] = gghy; gghyy[1] = gghy; gghyy[2] = gghy;
                
        glocx[0] = ggrav[2]; glocx[1] = ggrav[2]; glocx[2] = ggrav[2];
        glocy[0] = ggrav[3]; glocy[1] = ggrav[3]; glocy[2] = ggrav[3];
    }

    glocz[0] = ggz[8]; glocz[1] = ggz[9]; glocz[2] = ggz[10];


    if ( sico.SLOMO <= 1.0 ) {

        for ( gl = 0; gl < sico.PMAX; gl++ ) {

            if ( galpha[gl] > 0 ) { 

                if ( fabs( sflow.CDEFORM[gl] ) < 1 ) gctrl_deform = 1; else gctrl_deform = 0;


                // Hydraulic pressure coefficients

                if ( gctrl_deform == 1 ) {
                
                    if ( sflow.CDEFORM[gl] >= 0 ) {
                
                        gfact_slidex = fmin(9999, fmax(0, 1 - (( 1 - sflow.CDEFORM[gl] ) * ( 1 - gzfrag )) * sflow.GZMIN * 0.95 / ( ggz[gl] ) + sflow.FRAGSML * gzfrag ));
                        gfact_slidey = fmin(9999, fmax(0, 1 - (( 1 - sflow.CDEFORM[gl] ) * ( 1 - gzfrag )) * sflow.GZMIN * 0.95 / ( ggz[gl+4] ) + sflow.FRAGSML * gzfrag ));
                        
                    } else {

                        gfact_slidex = fmin(9999, fmax(0, 1 - ( 1 + sflow.CDEFORM[gl] ) * ( 1 - gzfrag ) + sflow.FRAGSML * gzfrag ));
                        gfact_slidey = fmin(9999, fmax(0, 1 - ( 1 + sflow.CDEFORM[gl] ) * ( 1 - gzfrag ) + sflow.FRAGSML * gzfrag ));
                    }
                    
                } else { 
                
                    gfact_slidey = 1;
                    gfact_slidex = 1; // factors for deformation control and fragmentation (x)
                }
        
                if ( sflow.DELTA[gl] > 0 ) gpbetax = gfact_slidex * gkx[gl] * ggz[gl] * gghxx[gl] /  cos( gbetaxy ); 
                else gpbetax = ggz[gl] * gghxx[gl] /  cos( gbetaxy );


                if ( sflow.DELTA[gl] > 0 ) gpbetay = gfact_slidey * gky[gl] * ggz[gl+4] * gghyy[gl] /  cos( gbetaxy ); 
                else gpbetay = ggz[gl+4] * gghyy[gl] /  cos( gbetaxy );


               // Effective basal pressures
               
               if ( sico.LAYERS == 0 ) {

                    if ( gl == 0 ) {

                        gpx = ggamma3 * glocz[gl] * gghxx[gl];
                        gpy = ggamma3 * glocz[gl] * gghyy[gl];

                    } else if ( gl == 1 ) {

                        gpx = -0.5 * ggamma4 * glocz[gl] * ghnorm / galpha[gl] * ggalphax[gl];
                        gpy = -0.5 * ggamma4 * glocz[gl] * ghnorm / galpha[gl] * ggalphay[gl];
                
                    } else {

                        gpx = -0.5 * ggamma5 * glocz[gl] * ghnorm / galpha[gl] * ggalphax[gl];
                        gpy = -0.5 * ggamma5 * glocz[gl] * ghnorm / galpha[gl] * ggalphay[gl];
                    }
                    
                } else {

                    gpx = 0;
                    gpy = 0;
                }
                    
                if ( gcslide == 1 && sflow.DELTA[gl] > 0 ) gslidefractb = fmin(1, fmax(0, pow(fdiv( gtslide - gtsum, gtslide, 0 ), sflow.SLIDEEXP[gl] ))) * sflow.SLIDEDEF[gl];
                else gslidefractb = 0; // factor for global gravity components

                gcompx[gl] = galpha[gl] * ( gangx[gl] * gslidefractb + glocx[gl] * ( 1 - gslidefractb ) - gpbetax - gpx );
                gcompy[gl] = galpha[gl] * ( gangy[gl] * gslidefractb + glocy[gl] * ( 1 - gslidefractb ) - gpbetay - gpy );
                
            } else {
            
                gcompx[gl] = 0;
                gcompy[gl] = 0;
            }
        }
    }


    // Drag terms

    if ( sico.MODEL == 7 ) {

        if ( sico.LAYERS != 0 ) {

            if ( ghflow > sico.HFLOWMIN && ghflow3 > sico.HFLOWMIN ) { gcdrag[0] = 1; gcdrag[1] = 1; gcdrag[2] = 1; }
            else { gcdrag[0] = 0; gcdrag[1] = 0; gcdrag[2] = 0; }

            if ( sflow.DRAG_k > 0 ) {
            
                gcdragx[0] = gcdrag[0] *= sflow.DRAG_k; gcdragx[1] = gcdrag[1] *= sflow.DRAG_k; gcdragx[2] = gcdrag[2] *= sflow.DRAG_k; 
                gcdragy[0] = gcdrag[0] *= sflow.DRAG_k; gcdragy[1] = gcdrag[1] *= sflow.DRAG_k; gcdragy[2] = gcdrag[2] *= sflow.DRAG_k;
                
                
                if ( fsign( gghx1 + gghx2 ) == fsign( gu[2] - gu[0] ) && fabs( gu[0] ) > fabs( gu[2] ) && gh[0] > sico.HFLOWMIN ) 
                    gcdragx[0] *= sin(atan(fabs(gghx1 + gghx2))) * galpha[0]; else gcdragx[0] = 0;
                if ( fsign( gghx1 ) == fsign( gu[1] - gu[0] ) && fabs( gu[0] ) > fabs( gu[1] ) &&  gh[0] > sico.HFLOWMIN ) 
                    gcdragx[1] *= sin(atan(fabs(gghx1))) * galpha[1]; else gcdragx[1] = 0;
                if ( fsign( gghx1 + gghx2 ) == fsign( gu[2] - gu[1] ) && fabs( gu[1] ) > fabs( gu[2] ) &&  gh[0] > sico.HFLOWMIN ) 
                    gcdragx[2] *= sin(atan(fabs(gghx1 + gghx2))) * galpha[0]; else gcdragx[2] = 0;

                if ( fsign( gghy1 + gghy2 ) == fsign( gv[2] - gv[0] ) && fabs( gv[0] ) > fabs( gv[2] ) &&  gh[0] > sico.HFLOWMIN ) 
                    gcdragy[0] *= sin(atan(fabs(gghy1 + gghy2))) * galpha[0]; else gcdragy[0] = 0;
                if ( fsign( gghy1 ) == fsign( gv[1] - gv[0] ) && fabs( gv[0] ) > fabs( gv[1] ) &&  gh[0] > sico.HFLOWMIN ) 
                    gcdragy[1] *= sin(atan(fabs(gghy1))) * galpha[1]; else gcdragy[1] = 0;
                if ( fsign( gghy1 + gghy2 ) == fsign( gv[2] - gv[1] ) && fabs( gv[1] ) > fabs( gv[2] ) &&  gh[0] > sico.HFLOWMIN ) 
                    gcdragy[2] *= sin(atan(fabs(gghy1 + gghy2))) * galpha[0]; else gcdragy[2] = 0;

            } else {            

                gcdragx[0] = sflow.DRAG_k * galpha[0] * -1; gcdragx[1] = sflow.DRAG_k * galpha[1] * -1; gcdragx[2] = sflow.DRAG_k * galpha[2] * -1; 
                gcdragy[0] = sflow.DRAG_k * galpha[0] * -1; gcdragy[1] = sflow.DRAG_k * galpha[1] * -1; gcdragy[2] = sflow.DRAG_k * galpha[2] * -1;
            }  

        } else {

            gcdragx[0] = gcdrag[0]; gcdragx[1] = gcdrag[1]; gcdragx[2] = gcdrag[2]; 
            gcdragy[0] = gcdrag[0]; gcdragy[1] = gcdrag[1]; gcdragy[2] = gcdrag[2];
        }

        gdragx[0] = gcdragx[0] * ( gu[2] - gu[0] ) * pow( fabs( guvw[2] - guvw[0] ), sflow.DRAG_je - 1 );
        gdragx[1] = gcdragx[1] * ( gu[1] - gu[0] ) * pow( fabs( guvw[1] - guvw[0] ), sflow.DRAG_je - 1 );
        gdragx[2] = gcdragx[2] * ( gu[2] - gu[1] ) * pow( fabs( guvw[2] - guvw[1] ), sflow.DRAG_je - 1 );

        gdragy[0] = gcdragy[0] * ( gv[2] - gv[0] ) * pow( fabs( guvw[2] - guvw[0] ), sflow.DRAG_je - 1 );
        gdragy[1] = gcdragy[1] * ( gv[1] - gv[0] ) * pow( fabs( guvw[1] - guvw[0] ), sflow.DRAG_je - 1 );
        gdragy[2] = gcdragy[2] * ( gv[2] - gv[1] ) * pow( fabs( guvw[2] - guvw[1] ), sflow.DRAG_je - 1 );

        if ( sico.LAYERS != 0 && ghflow2 > sico.HFLOWMIN ) { gdragx[0] = 0; gdragy[0] = 0; }

    } else { gdragx[0] = 0; gdragx[1] = 0; gdragx[2] = 0; gdragy[0] = 0; gdragy[1] = 0; gdragy[2] = 0; }


    // Source terms

    if ( sico.SLOMO > 1.0 ) { // equilibrium-of-motion model

        if ( galpha[0] > 0 ) {

            ggs[0] = 0;

            if ( sflow.CDEFORM[0] == 0 ) gslideterm = 0; else gslideterm = 1 / sflow.CDEFORM[0];

            if ( gcslide == 1 && sflow.CDEFORM[0] != 0 ) gslidefractb = fmin(1, fmax(0, pow(fdiv( gtslide - gtsum, gtslide, 0 ), sflow.SLIDEEXP[0] ))) * sflow.SLIDEDEF[0];
            else gslidefractb = 0;
           
            if ( sico.GLACIER == 1 ) {

                ggs[1] = fsign( ggrav[5] ) * fabs( ghnorm * ( 2 * sflow.NY[0] * pow( sflow.RHO1 * ggrav[5] * ghnorm, 3 ) / 5.0 
                    + sflow.CDEFORM[0] * pow( sflow.RHO1 * ggrav[5] * ghnorm, 2 )));
                ggs[2] = fsign( ggrav[6] ) * fabs( ghnorm * ( 2 * sflow.NY[0] * pow( sflow.RHO1 * ggrav[6] * ghnorm, 3 ) / 5.0 
                    + sflow.CDEFORM[0] * pow( sflow.RHO1 * ggrav[6] * ghnorm, 2 )));
                        
            } else if ( sico.LAYERS == 0 && sico.MODEL < 7 ) {

                ggs[1] = sico.SLOMO * ghnorm * galpha[0] * ( ggrav[5] * ghnorm * ghnorm / sflow.NY[0] * 0.5
                    + ( ggrav[2] * ( 1 - gslidefractb ) + gangx[0] * gslidefractb ) * gslideterm ) + ghnorm * ( gdragx[0] + gdragx[1] );
                ggs[2] = sico.SLOMO * ghnorm * galpha[0] * ( ggrav[6] * ghnorm * ghnorm / sflow.NY[0] * 0.5
                    + ( ggrav[3] * ( 1 - gslidefractb ) + gangy[0] * gslidefractb ) * gslideterm ) + ghnorm * ( gdragy[0] + gdragy[1] );         

            } else if ( sico.LAYERS == 0 && sico.MODEL == 7 ) {

                ggs[1] = sico.SLOMO * ghnorm * galpha[0] * ( ggrav[13] * ghnorm * ghnorm / sflow.NY[0] * 0.5
                    + ( ggrav[2] * ( 1 - gslidefractb ) + gangx[0] * gslidefractb ) * gslideterm ) + ghnorm * ( gdragx[0] + gdragx[1] );
                ggs[2] = sico.SLOMO * ghnorm * galpha[0] * ( ggrav[14] * ghnorm * ghnorm/ sflow.NY[0] * 0.5
                    + ( ggrav[3] * ( 1 - gslidefractb ) + gangy[0] * gslidefractb ) * gslideterm ) + ghnorm * ( gdragy[0] + gdragy[1] );  

            } else {
            
                ggs[1] = sico.SLOMO * ghnorm1s * ( ggrav[5] * ghnorm1s * ghnorm1s / sflow.NY[0] * 0.5
                    + ( ggrav[2] * ( 1 - gslidefractb ) + gangx[0] * gslidefractb ) * gslideterm );
                ggs[2] = sico.SLOMO * ghnorm1s * ( ggrav[6] * ghnorm1s * ghnorm1s / sflow.NY[0] * 0.5
                    + ( ggrav[3] * ( 1 - gslidefractb ) + gangy[0] * gslidefractb ) * gslideterm );         
            }
            
        } else { ggs[1] = 0; ggs[2] = 0; }

        if ( sico.MODEL == 7 ) {

            ggs[3] = 0;
            ggs[6] = 0;

            if ( sflow.CDEFORM[1] == 0 ) gslideterm = 0; else gslideterm = 1 / sflow.CDEFORM[1];

            if ( gcslide == 1 && sflow.CDEFORM[1] != 0 ) gslidefractb = fmin(1, fmax(0, pow(fdiv( gtslide - gtsum, gtslide, 0 ), sflow.SLIDEEXP[1] ))) * sflow.SLIDEDEF[1];
            else gslidefractb = 0;

            if ( galpha[1] > 0 && sico.LAYERS == 0 ) {

                ggs[4] = sico.SLOMO * ghnorm * galpha[1] * ( ggrav[13] * ghnorm * ghnorm / sflow.NY[1] * 0.5
                    + ( ggrav[2] * ( 1 - gslidefractb ) + gangx[0] * gslidefractb ) * gslideterm ) - ghnorm * ( 1 / ggamma1 * gdragx[1] - gdragx[2] );
                ggs[5] = sico.SLOMO * ghnorm * galpha[1] * ( ggrav[14] * ghnorm * ghnorm / sflow.NY[1] * 0.5
                    + ( ggrav[3] * ( 1 - gslidefractb ) + gangy[0] * gslidefractb ) * gslideterm ) - ghnorm * ( 1 / ggamma1 * gdragy[1] - gdragy[2] );

            } else if ( galpha[1] > 0 ) {

                ggs[4] = sico.SLOMO * ghnorm2s * ( ggrav[9] * ghnorm2s * ghnorm2s / sflow.NY[1] * 0.5
                    + ( ggrav[7] * ( 1 - gslidefractb ) + gangx[1] * gslidefractb ) * gslideterm );
                ggs[5] = sico.SLOMO * ghnorm2s * ( ggrav[10] * ghnorm2s * ghnorm2s / sflow.NY[1] * 0.5
                    + ( ggrav[8] * ( 1 - gslidefractb ) + gangy[1] * gslidefractb ) * gslideterm );         
             
            } else { ggs[4] = 0; ggs[5] = 0; }

            if ( sflow.CDEFORM[2] == 0 ) gslideterm = 0; else gslideterm = 1 / sflow.CDEFORM[2];

            if ( gcslide == 1 && sflow.CDEFORM[2] != 0 ) gslidefractb = fmin(1, fmax(0, pow(fdiv( gtslide - gtsum, gtslide, 0 ), sflow.SLIDEEXP[2] ))) * sflow.SLIDEDEF[2];
            else gslidefractb = 0;

            if ( galpha[2] > 0 && sico.LAYERS == 0 ) {

                ggs[7] = sico.SLOMO * ghnorm * galpha[2] * ( ggrav[13] * ghnorm * ghnorm / sflow.NY[2] * 0.5
                    + ( ggrav[2] * ( 1 - gslidefractb ) + gangx[0] * gslidefractb ) * gslideterm ) - ghnorm * ( 1 / ggamma0 * gdragx[0] + 1 / ggamma2 * gdragx[2] );
                ggs[8] = sico.SLOMO * ghnorm * galpha[2] * ( ggrav[14] * ghnorm * ghnorm / sflow.NY[2] * 0.5
                    + ( ggrav[3] * ( 1 - gslidefractb ) + gangy[0] * gslidefractb ) * gslideterm ) - ghnorm * ( 1 / ggamma0 * gdragy[0] + 1 / ggamma2 * gdragy[2] );

            } else if ( galpha[2] > 0 ) {
            
                ggs[7] = sico.SLOMO * ghnorm3s * ( ggrav[13] * ghnorm3s * ghnorm3s / sflow.NY[2] * 0.5
                    + ( ggrav[11] * ( 1 - gslidefractb ) + gangx[2] * gslidefractb ) * gslideterm );
                ggs[8] = sico.SLOMO * ghnorm3s * ( ggrav[14] * ghnorm3s * ghnorm3s / sflow.NY[2] * 0.5
                    + ( ggrav[12] * ( 1 - gslidefractb ) + gangy[2] * gslidefractb ) * gslideterm );         

            } else { ggs[7] = 0; ggs[8] = 0; }
        }
        
        ggs[gt] = 0;
        
    } else { // model for extremely rapid flows:

        if ( galpha[0] > 0 ) {
 
            ggs[0] = 0;
            ggs[1] = sflow.LAMBDA_FL * ghnorm * ( gcompx[0] + gdragx[0] + gdragx[1] );
            ggs[2] = sflow.LAMBDA_FL * ghnorm * ( gcompy[0] + gdragy[0] + gdragy[1] );

        } else { ggs[0] = 0; ggs[1] = 0; ggs[2] = 0; }
    }

    if ( sico.MODEL == 7 && sico.SLOMO <= 1 ) {

        if ( galpha[1] > 0 ) {

            ggs[3] = 0; 
            ggs[4] = ghnorm * ( gcompx[1] - 1 / ggamma1 * gdragx[1] + gdragx[2] );
            ggs[5] = ghnorm * ( gcompy[1] - 1 / ggamma1 * gdragy[1] + gdragy[2] );

        } else { ggs[3] = 0; ggs[4] = 0; ggs[5] = 0; }

        if ( galpha[2] > 0 ) {

            ggs[6] = 0; 
            ggs[7] = ghnorm * ( gcompx[2] - 1 / ggamma0 * gdragx[0] - 1 / ggamma2 * gdragx[2] );
            ggs[8] = ghnorm * ( gcompy[2] - 1 / ggamma0 * gdragy[0] - 1 / ggamma2 * gdragy[2] );

        } else { ggs[6] = 0; ggs[7] = 0; ggs[8] = 0; }
    }

    if ( sico.TEMP > 0 ) { // temperature term

        if ( sico.MODEL == 7 ) {

            galphari = galpha[0] + galpha[1];
            gmyri = fdiv( tan( sflow.DELTA[0] ) * galpha[0] + tan( sflow.DELTA[1] ) * galpha[1], galphari, 0 );
            grhori = fdiv( sflow.RHO1 * galpha[0] + sflow.RHO2 * galpha[1], galphari, 0 );
            
            grho = sflow.RHO1 * galpha[0] + sflow.RHO2 * galpha[1] + sflow.RHO3 * galpha[2];
            guvwz = guvw[0] * galpha[0] + guvw[1] * galpha[1] + guvw[2] * galpha[2];
            gub = ( gu[0] * galpha[0] + gu[1] * galpha[1] + gu[2] * galpha[2] ) * sflow.LAMBDA;
            gvb = ( gv[0] * galpha[0] + gv[1] * galpha[1] + gv[2] * galpha[2] ) * sflow.LAMBDA;
            gc = sflow.CR * galpha[0] + sflow.CI * galpha[1] + sflow.CF * galpha[2];
            gkappa = sflow.KAPPAR * galpha[0] + sflow.KAPPAI * galpha[1] + sflow.KAPPAF * galpha[2];
            gtau = grhori * ggrav[4] * ghnorm * gmyri * galphari;

            ggtempx2 = 0; ggtempy2 = 0; ggkappax = 0; ggkappay = 0; ggtempx = 0; ggtempy = 0, gtemp = 0;

            if ( galpha[1] > 0 ) {

                rmelt = sflow.L / sflow.Q * grhori / sflow.RHO2 * ggrav[4] * gmyri * ghnorm * guvwz * galphari;

                ggs[3] -= rmelt;
                ggs[4] -= rmelt * gu[1];
                ggs[5] -= rmelt * gv[1];
                ggs[6] += rmelt / ggamma2;
                ggs[7] += rmelt / ggamma2 * gu[1];
                ggs[8] += rmelt / ggamma2 * gv[1];

                ggs[gt] = ghnorm * ( gkappa / ( grho * gc ) * ( ggtempx2 + ggtempy2 ) + 1 / ( grho * gc ) * ( ggkappax * ggtempx + ggkappay * ggtempy ))
                    - 1 / ( grho * gc ) * sflow.LAMBDAS * sflow.AS * ( gtemp - sflow.TEMPS ) + 1 / ( grho * gc ) * sflow.LAMBDAB * gkappa * ( sflow.TEMPB - gtemp )
                    + 1 / ( grho * gc ) * sflow.L * gtau * ( gub + gvb ) 
                    - 1 / ( grho * gc ) * ( sflow.CL / sflow.Q * grhori / sflow.RHO2 * ggrav[4] * gmyri ) * ghnorm * guvwz * galphari;

            } else {
                
                ggs[gt] = ghnorm * ( gkappa / ( grho * gc ) * ( ggtempx2 + ggtempy2 ) + 1 / ( grho * gc ) * ( ggkappax * ggtempx + ggkappay * ggtempy ))
                    - 1 / ( grho * gc ) * sflow.LAMBDAS * sflow.AS * ( gtemp - sflow.TEMPS ) + 1 / ( grho * gc ) * sflow.LAMBDAB * gkappa * ( sflow.TEMPB - gtemp )
                    + 1 / ( grho * gc ) * sflow.L * gtau * ( gub + gvb );                
            }

        } else {

            galphari = galpha[0];
            gmyri = tan( sflow.DELTA[0] );
            grhori = sflow.RHO1;
            
            grho = sflow.RHO1;
            guvwz = guvw[0];
            gub = gu[0] * sflow.LAMBDA;
            gvb = gv[0] * sflow.LAMBDA;
            gc = sflow.CR * galpha[0];
            gkappa = sflow.KAPPAR * galpha[0];
            gtau = grhori * ggrav[4] * ghnorm * gmyri * galphari;

            ggtempx2 = 0; ggtempy2 = 0; ggkappax = 0; ggkappay = 0; ggtempx = 0; ggtempy = 0, gtemp = 0;

            ggs[gt] = ghnorm * ( gkappa / ( grho * gc ) * ( ggtempx2 + ggtempy2 ) + 1 / ( grho * gc ) * ( ggkappax * ggtempx + ggkappay * ggtempy ))
                - 1 / ( grho * gc ) * sflow.LAMBDAS * sflow.AS * ( gtemp - sflow.TEMPS ) + 1 / ( grho * gc ) * sflow.LAMBDAB * gkappa * ( sflow.TEMPB - gtemp )
                + 1 / ( grho * gc ) * sflow.L * gtau * ( gub + gvb );            
        }

    } else ggs[gt] = 0;   

    return ggs;
}

float *fd ( float *gh, float gghx, float gghy, float *gnuss, float *gnvss, float *gnufs, float *gnvfs,
    float *gnuff, float *gnvff, float *ggalpha10, float *ggalpha20, float *ggalpha30, float gbetax, float gbetay, float gbetaxy, float *gnbetax, float *gnbetay, float *ggz,
    float *gkappau, float gekin, float *ggrav, float *gangx, float *gangy, int gcslide, float gtsum, float gtslide, float gzfrag, struct ico sico, struct flow sflow ) { 
    // function for source terms (decelerating components)

    int gj, gl, gm, gt, gcdelta, gctufri;
    float *ggd, gw[3][6], guvw[3][6], ggalpha[3][9], gfricx[3], gfricy[3], gvisx[3], gvisy[3], gturbx[3], gturby[3], gflufrix[3], gflufriy[3], gambdragx[3], gambdragy[3], 
        gaddfri[3], gdelta, gpbx[3], gpby[3], guratio[3], gvratio[3], gnu[3][9], gnv[3][9], ghnorm, ghnormx[3], gslidefractb = 0, gdx, gdy, ggamma0, ggamma2, gvc, gangz[3], 
        ggux[3][6], gguy[3][6], ggvx[3][6], ggvy[3][6], gnye[3][6], gxterm1, gyterm1, gxterm2, gyterm2, gxterm3, gyterm3, ggzx, ggzy;

    ggd = (float*) calloc( sico.NVECTMIN, sizeof(float));

    gdx = sico.CSZ / cos( gbetax ); // topography-corrected cell sizes and flow height
    gdy = sico.CSZ / cos( gbetay );
    ghnorm = gh[0] * cos( gbetaxy );
    gghx *= cos( gbetaxy ); gghy *= cos( gbetaxy );

    for( gj=0; gj<9; gj++ ) { 

        gnu[0][gj] = gnuss[gj]; gnu[1][gj] = gnufs[gj]; gnu[2][gj] = gnuff[gj]; gnv[0][gj] = gnvss[gj]; gnv[1][gj] = gnvfs[gj]; gnv[2][gj] = gnvff[gj];
        ggalpha[0][gj] = ggalpha10[gj]; ggalpha[1][gj] = ggalpha20[gj]; ggalpha[2][gj] = ggalpha30[gj]; // velocities and phase fractions (renaming for convenience)
        
        for ( gl = 0; gl < sico.PMAX; gl++ ) { gnu[gl][gj] /= cos( gbetaxy ); gnv[gl][gj] /= cos( gbetaxy ); } // topography-corrected velocities
    }

    if ( ggalpha[0][0] > 0 && ggalpha[2][0] > 0 ) ggamma0 = sflow.RHO3 / sflow.RHO1; else ggamma0 = 0;
    if ( ggalpha[1][0] > 0 && ggalpha[2][0] > 0 ) ggamma2 = sflow.RHO3 / sflow.RHO2; else ggamma2 = 1;

    for ( gj=0; gj<6; gj++ ) {

        for ( gl = 0; gl < 3; gl++ ) {

            if ( sico.NONHYDRO > 99 ) gw[0][gj] = gnu[gl][gj] * tan( gnbetax[gj] ) + gnv[gl][gj] * tan( gnbetay[gj] ); // non-hydrostatic effects
            else gw[gl][gj] = 0;

            guvw[gl][gj] = pow( pow( gnu[gl][gj], 2 ) + pow( gnv[gl][gj], 2 ) + pow( gw[gl][gj], 2 ), 0.5 );
            guratio[gl] = fdiv( gnu[gl][0], guvw[gl][0], 0 );
            gvratio[gl] = fdiv( gnv[gl][0], guvw[gl][0], 0 );
        }
    }


    if ( sico.SLOMO > 1.0 ) { // Equilibrium-of-motion model: no deceleration terms

        for ( gl=0; gl<sico.NVECTMIN; gl++ ) ggd[gl] = 0;

    } else { // all other models

        for ( gl = 0; gl < sico.PMAX; gl++ ) {

            if ( ggalpha[gl][0] > 0 ) {
            
                if (sico.LAYERS > 0 ) ghnormx[gl] = ghnorm * ggalpha[gl][0]; else ghnormx[gl] = ghnorm;

                if ( sflow.PHI[gl] < 0 ) {
                
                    gvc = ( sflow.PHI[gl] * -180 / sico.PI ) * pow( -1 * sflow.ADDFRI[gl] * ghnorm, 0.333333333 );
                
                    if ( guvw[gl][0] > gvc ) {
                    
                        gcdelta = 0;
                        gctufri = 1;
                    
                    } else {
                    
                        gcdelta = 1;
                        gctufri = 0;
                    } 
                    
                } else {
                
                    gcdelta = 1;
                    gctufri = 1;
                }
                
                gdelta = gcdelta * sflow.DELTA[gl];

                if ( gcslide == 1 && sflow.DELTA[gl] > 0 ) gslidefractb = fmin(1, fmax(0, pow(fdiv( gtslide - gtsum, gtslide, 0 ), sflow.SLIDEEXP[gl] ))) * sflow.SLIDEDEF[gl];
                else gslidefractb = 0;

                if ( gl == 0 && sico.LAYERS == 0 ) gangz[gl] = ( 1 - ggamma0 ) * gangx[gl+3];
                else if ( gl == 1 && sico.LAYERS == 0 ) gangz[gl] = ggamma2 * gangx[gl+3];
                else gangz[gl] = gangx[gl+3];

                if ( gl == 0 ) gm = 8; else if ( gl == 1 ) gm = 9; else gm = 10;

                if ( sico.CURVCTRL == 1 ) { ggzx = ggz[gm] + gkappau[gl] * pow( gnu[gl][0], 2 ); ggzy = ggz[gm] + gkappau[gl+3] * pow( gnv[gl][0], 2 ); } // curvature effects
                else { ggzx = ggz[gm]; ggzy = ggz[gm]; }

                gpbx[gl] = ggzx * ( 1 - gslidefractb ) + gangz[gl] * gslidefractb; // effective basal pressures
                gpby[gl] = ggzy * ( 1 - gslidefractb ) + gangz[gl] * gslidefractb;

                if ( sflow.ESP != 0 ) {

                    ggux[gl][0] = ( gnu[gl][2] - gnu[gl][5] ) / ( 2 * gdx ); // gradients of flow velocity
                    ggvy[gl][0] = ( gnv[gl][1] - gnv[gl][4] ) / ( 2 * gdy );
                }

                gfricx[gl] = sflow.LAMBDA_FL * ( guratio[gl] * tan( gdelta ) * gpbx[gl] 
                    + sflow.COHESION[gl] / ( ghnormx[gl] ) * fsign( gnu[gl][0] ) - sflow.ESP * tan( sflow.PHI[gl] ) * fabs( sflow.CDEFORM[gl] ) * ggux[gl][0] * gghx
                    + gzfrag * sflow.FRAGML * gnu[gl][0] ); // friction terms
                gfricy[gl] = sflow.LAMBDA_FL * ( gvratio[gl] * tan( gdelta ) * gpby[gl] 
                    + sflow.COHESION[gl] / ( ghnormx[gl] ) * fsign( gnv[gl][0] ) - sflow.ESP * tan( sflow.PHI[gl] ) * fabs( sflow.CDEFORM[gl] ) * ggvy[gl][0] * gghy
                    + gzfrag * sflow.FRAGML * gnv[gl][0] );

                if ( sflow.ADDFRI[gl] > 0 ) { // fluid friction

                    gaddfri[gl] = sflow.ADDFRI[gl];
                
                    gflufrix[gl] = pow( gaddfri[gl], 2 ) / pow( ffmax(sico.HFLOWMIN, ghnormx[gl] ), 1.3333 ) * gnu[gl][0] * guvw[gl][0];
                    gflufriy[gl] = pow( gaddfri[gl], 2 ) / pow( ffmax(sico.HFLOWMIN, ghnormx[gl] ), 1.3333 ) * gnv[gl][0] * guvw[gl][0];
                                
                } else { gflufrix[gl] = 0; gflufriy[gl] = 0; }

                if ( sflow.ADDFRI[gl] < 0 && gctufri == 1 ) { // turbulent friction
      
                    gturbx[gl] = guratio[gl] * ggzx * pow( guvw[gl][0], 2 ) / ( -1 * sflow.ADDFRI[gl] * ghnormx[gl] );
                    gturby[gl] = gvratio[gl] * ggzy * pow( guvw[gl][0], 2 ) / ( -1 * sflow.ADDFRI[gl] * ghnormx[gl] );

                } else { gturbx[gl] = 0; gturby[gl] = 0; }  

                if ( sflow.AMBDRAG > 0 ) {

                    if ( fsign( gnu[0][0] * ggalpha[0][0] + gnu[1][0] * ggalpha[1][0] + gnu[2][0] * ggalpha[2][0] ) == fsign( gghx )) gambdragx[gl] = 0;
                    else gambdragx[gl] = gdx * fabs( gghx ) / ghnormx[gl] * sflow.AMBDRAG;
                    gambdragx[gl] = ffmin( sflow.AMBDRAG, gambdragx[gl] ) * gnu[gl][0] * guvw[gl][0]; // effective ambient drag coefficients

                    if ( fsign( gnv[0][0] * ggalpha[0][0] + gnv[1][0] * ggalpha[1][0] + gnv[2][0] * ggalpha[2][0] ) == fsign( gghy )) gambdragy[gl] = 0;
                    else gambdragy[gl] = gdy * fabs( gghy ) / ghnormx[gl] * sflow.AMBDRAG;
                    gambdragy[gl] = ffmin( sflow.AMBDRAG, gambdragy[gl] ) * gnv[gl][0] * guvw[gl][0];

                } else { gambdragx[gl] = 0; gambdragy[gl] = 0; }

                if ( sflow.NY[gl] > 0.0000001 ) { // viscosity terms

                    ggux[gl][0] = ( gnu[gl][2] - gnu[gl][5] ) / ( 2 * gdx ); // gradients of flow velocity
                    gguy[gl][0] = ( gnu[gl][1] - gnu[gl][4] ) / ( 2 * gdy );
                    ggvx[gl][0] = ( gnv[gl][2] - gnv[gl][5] ) / ( 2 * gdx );
                    ggvy[gl][0] = ( gnv[gl][1] - gnv[gl][4] ) / ( 2 * gdy );

                    ggux[gl][1] = ( gnu[gl][3] - gnu[gl][8] ) / ( 2 * gdx );
                    gguy[gl][1] = ( gnu[gl][1] - gnu[gl][0] ) / ( 1 * gdy );
                    ggvx[gl][1] = ( gnv[gl][3] - gnv[gl][8] ) / ( 2 * gdx );
                    ggvy[gl][1] = ( gnv[gl][1] - gnv[gl][0] ) / ( 1 * gdy );

                    ggux[gl][2] = ( gnu[gl][2] - gnu[gl][0] ) / ( 1 * gdx );
                    gguy[gl][2] = ( gnu[gl][3] - gnu[gl][7] ) / ( 2 * gdy );
                    ggvx[gl][2] = ( gnv[gl][2] - gnv[gl][0] ) / ( 1 * gdx );
                    ggvy[gl][2] = ( gnv[gl][3] - gnv[gl][7] ) / ( 2 * gdy );

                    ggux[gl][4] = ( gnu[gl][7] - gnu[gl][6] ) / ( 2 * gdx );
                    gguy[gl][4] = ( gnu[gl][0] - gnu[gl][4] ) / ( 1 * gdy );
                    ggvx[gl][4] = ( gnv[gl][7] - gnv[gl][6] ) / ( 2 * gdx );
                    ggvy[gl][4] = ( gnv[gl][0] - gnv[gl][4] ) / ( 1 * gdy );

                    ggux[gl][5] = ( gnu[gl][0] - gnu[gl][5] ) / ( 1 * gdx );
                    gguy[gl][5] = ( gnu[gl][8] - gnu[gl][6] ) / ( 2 * gdy );
                    ggvx[gl][5] = ( gnv[gl][0] - gnv[gl][5] ) / ( 1 * gdx );
                    ggvy[gl][5] = ( gnv[gl][8] - gnv[gl][6] ) / ( 2 * gdy );

                    for ( gj=0; gj<6; gj++ ) gnye[gl][gj] = sflow.NY[gl]; // effective kinematic viscosity

                    if ( ggalpha[gl][0] > 0 ) {

                        gxterm1 = 2 / gdx * sflow.LAMBDA_NU * ( gnye[gl][2] * ggux[gl][2] - gnye[gl][5] * ggux[gl][5] ); // components of viscosity terms
                        gxterm2 = 1 / gdy * sflow.LAMBDA_NU * ( gnye[gl][1] * ggvx[gl][1] - gnye[gl][4] * ggvx[gl][4] );
                        gxterm3 = 1 / gdy * sflow.LAMBDA_NU * ( gnye[gl][1] * gguy[gl][1] - gnye[gl][4] * gguy[gl][4] );

                        gyterm1 = 2 / gdy * sflow.LAMBDA_NU * ( gnye[gl][1] * ggvy[gl][1] - gnye[gl][4] * ggvy[gl][4] );
                        gyterm2 = 1 / gdx * sflow.LAMBDA_NU * ( gnye[gl][2] * gguy[gl][2] - gnye[gl][5] * gguy[gl][5] );
                        gyterm3 = 1 / gdx * sflow.LAMBDA_NU * ( gnye[gl][2] * ggvx[gl][2] - gnye[gl][5] * ggvx[gl][5] );

                        gvisx[gl] = -gxterm1 - gxterm2 - gxterm3 + sflow.LAMBDA_NU * gnye[gl][0] * fdiv( sflow.VIS_chi[gl] * gnu[gl][0], sflow.HREF, sico.HFLOWMIN );
                        gvisy[gl] = -gyterm1 - gyterm2 - gyterm3 + sflow.LAMBDA_NU * gnye[gl][0] * fdiv( sflow.VIS_chi[gl] * gnv[gl][0], sflow.HREF, sico.HFLOWMIN );
                            // viscosity terms

                    } else { gvisx[gl] = 0; gvisy[gl] = 0; }
           
                } else { gvisx[gl] = 0; gvisy[gl] = 0; }

                ggd[3*gl] = 0; // deceleration terms
                if ( gnu[gl][0] != 0 ) ggd[3*gl+1] = ghnorm * ggalpha[gl][0] * ( gfricx[gl] + gturbx[gl] + gflufrix[gl] + gambdragx[gl] + gvisx[gl] ); else ggd[3*gl+1] = 0;
                if ( gnv[gl][0] != 0 ) ggd[3*gl+2] = ghnorm * ggalpha[gl][0] * ( gfricy[gl] + gturby[gl] + gflufriy[gl] + gambdragy[gl] + gvisy[gl] ); else ggd[3*gl+2] = 0;
                
            } else { ggd[3*gl+1] = 0; ggd[3*gl+2] = 0; }
        }
    }
 
    if ( sico.MODEL <= 3 ) gt=3; else gt=9;
    ggd[gt] = 0;

    return ggd;
}

float fsigma ( float gv, float gvm, float gvp, float gd, struct ico sico ) { // function for slope of vector

    float gsigma;

    if ( gvp == gv ) gsigma = 0;
    else {

        gsigma = ( gv - gvm ) / ( gvp - gv ); // input for limiter

        if ( sico.LIMITER == 1 ) gsigma = ffmax( 0, ffmin( 1, gsigma )); // Minmod limiter

        else if ( sico.LIMITER == 2 ) {

            gsigma = ffmax( 0, ffmin( 1, 2 * gsigma ));
            gsigma = ffmax( gsigma, ffmin( gsigma, 2 )); // Superbee limiter
            
        } else if ( sico.LIMITER == 3 ) {

            gsigma = ffmin( 2 * gsigma, 0.5 * ( 1 + gsigma ));
            gsigma = ffmax( 0, ffmin ( 2, gsigma )); // Woodward limiter
            
        } else if ( sico.LIMITER == 4 ) gsigma = ( gsigma + fabs( gsigma )) / ( 1 + fabs( gsigma )); // van Leer limiter

        gsigma *= ( gvp - gv ) / gd; // slope of vector
    }
    return gsigma;
}

float fcfl ( float ghflow, float ghflow2, float ghflow3, float gvflowx, float gvflowy, float gvflowx2, float gvflowy2,
    float gvflowx3, float gvflowy3, struct flow sflow, float *gkx, float *gky, float gbetax, float gbetay, float gbetaxy, float gzfrag, float *ggz, 
    float gghx, float gghy, float gghx1, float gghy1, float gghx2, float gghy2, int gi, struct ico sico ) { // function for determining time step length

    int gl, gctrl_deform;
    float gcelx = 0, gcely = 0, gcel, gcelr, gh, gu[3], gv[3], gfact_slidex, gfact_slidey, gghxx[3], gghyy[3], galpha[3];

    gh = ghflow + ghflow2 + ghflow3;
    
    gu[0] = gvflowx / cos( gbetaxy ); gu[1] = gvflowx2 / cos( gbetaxy ); gu[2] = gvflowx3 / cos( gbetaxy );
    gv[0] = gvflowy / cos( gbetaxy ); gv[1] = gvflowy2 / cos( gbetaxy ); gv[2] = gvflowy3 / cos( gbetaxy );

    if ( sico.MODEL <= 3 ) { gghx1 = gghx; gghy1 = gghy; }

    gghx *= cos( gbetaxy ); gghy *= cos( gbetaxy );
    gghx1 *= cos( gbetaxy ); gghy1 *= cos( gbetaxy );
    gghx2 *= cos( gbetaxy ); gghy2 *= cos( gbetaxy );

    if ( sico.LAYERS > 0 ) {

        gghxx[0] = gghx1; gghxx[1] = gghx2 - gghx1; gghxx[2] = gghx - gghx2;
        gghyy[0] = gghy1; gghyy[1] = gghy2 - gghy1; gghyy[2] = gghy - gghy2;

    } else {

        gghxx[0] = gghx; gghxx[1] = gghx; gghxx[2] = gghx;
        gghyy[0] = gghy; gghyy[1] = gghy; gghyy[2] = gghy;
    }

    if ( gh > sico.HFLOWMIN ) {

        if ( ghflow > sico.HFLOWMIN ) galpha[0] = ghflow / gh; else galpha[0] = 0;
        if ( ghflow2 > sico.HFLOWMIN ) galpha[1] = ghflow2 / gh; else galpha[1] = 0;
        if ( ghflow3 > sico.HFLOWMIN ) galpha[2] = ghflow3 / gh; else galpha[2] = 0;

        gcel = 0;

        for ( gl = 0; gl < sico.PMAX; gl++ ) {

            if ( galpha[gl] * gh > sico.HFLOWMIN ) {

                if ( fabs( sflow.CDEFORM[gl] ) < 1 ) gctrl_deform = 1; else gctrl_deform = 0;


                // Hydraulic pressure coefficients

                if ( gctrl_deform == 1 ) {
                
                    if ( sflow.CDEFORM[gl] >= 0 ) {
                
                        gfact_slidex = fmin(100, fmax(0, 1 - (( 1 - sflow.CDEFORM[gl] ) * ( 1 - gzfrag )) * sflow.GZMIN * 0.95 / ( ggz[gl] ) + sflow.FRAGSML * gzfrag ));
                        gfact_slidey = fmin(100, fmax(0, 1 - (( 1 - sflow.CDEFORM[gl] ) * ( 1 - gzfrag )) * sflow.GZMIN * 0.95 / ( ggz[gl+4] ) + sflow.FRAGSML * gzfrag ));
                        
                    } else {

                        gfact_slidex = fmin(100, fmax(0, 1 - ( 1 + sflow.CDEFORM[gl] ) * ( 1 - gzfrag ) + sflow.FRAGSML * gzfrag ));
                        gfact_slidey = fmin(100, fmax(0, 1 - ( 1 + sflow.CDEFORM[gl] ) * ( 1 - gzfrag ) + sflow.FRAGSML * gzfrag ));
                    }

                } else { 
                
                    gfact_slidex = 1;
                    gfact_slidey = 1;
                }
                
                if ( sflow.DELTA[gl] > 0 ) gcelx = 0.5 * gfact_slidex * gkx[gl] * ggz[gl] * gghxx[gl]; 
                else gcelx = 0.5 * ggz[gl] * gghxx[gl];
                

                if ( sflow.DELTA[gl] > 0 ) gcely = 0.5 * gfact_slidey * gky[gl] * ggz[gl+4] * gghyy[gl]; 
                else gcely = 0.5 * ggz[gl+4] * gghyy[gl];

                gcelx = fabs( gu[gl] ) + fabs( gcelx );
                gcely = fabs( gv[gl] ) + fabs( gcely );
        
                if ( gcelx > gcely ) gcelr = gcelx;  else gcelr = gcely;
                if ( gcelr > gcel ) gcel = gcelr;
            }
        }

    } else gcel = 0;

    return gcel;
}

void foutasc ( float **gparam, int *gpx, int *gpy, char *goutmaps, char *gname, int gk, struct ico sico ) { // function for output of ascii raster maps

    FILE *gfascii;
    char gpath[200];
    int gi, gx, gy;
    float gpout;

    sprintf(gpath, "%s%s.asc", goutmaps, gname ); // writing name of ascii file to string
    gfascii=fopen(gpath, "w"); // opening ascii file

    fprintf( gfascii, "ncols %i\n", sico.N ); // writing header to ascii file
    fprintf( gfascii, "nrows %i\n", sico.M );
    fprintf( gfascii, "xllcenter %.2f\n", sico.BDWEST );
    fprintf( gfascii, "yllcenter %.2f\n", sico.BDSOUTH );
    fprintf( gfascii, "cellsize %.2f\n", sico.CSZ );
    fprintf( gfascii, "NODATA_value %.3f\n", sico.UNDEF );

    gi = 0;
  
    for ( gx=0; gx<sico.M; gx++ ) { for ( gy=0; gy<sico.N; gy++ ) {

        if ( gpx[gi] == gx && gpy[gi] == gy ) {

            gpout = gparam[gi][gk];

            fprintf( gfascii, "%.3f ", gpout ); // writing output to ascii file
            gi += 1;
            
        } else fprintf( gfascii, "-9999.000 " ); // writing no data value to file if cell is not in area of interest
        
    } fprintf(gfascii, "\n"); }

    fclose(gfascii);

    return;
}

void foutascind ( int *gparam, int *gpx, int *gpy, char *goutmaps, char *gname, struct ico sico ) { // function for output of ascii raster maps

    FILE *gfascii;
    char gpath[200];
    int gi, gx, gy;
    int gpout;

    sprintf(gpath, "%s%s.asc", goutmaps, gname ); // writing name of ascii file to string
    gfascii=fopen(gpath, "w"); // opening ascii file

    fprintf( gfascii, "ncols %i\n", sico.N ); // writing header to ascii file
    fprintf( gfascii, "nrows %i\n", sico.M );
    fprintf( gfascii, "xllcenter %.2f\n", sico.BDWEST );
    fprintf( gfascii, "yllcenter %.2f\n", sico.BDSOUTH );
    fprintf( gfascii, "cellsize %.2f\n", sico.CSZ );
    fprintf( gfascii, "NODATA_value %.3f\n", sico.UNDEF );

    gi = 0;
    for ( gx=0; gx<sico.M; gx++ ) { for ( gy=0; gy<sico.N; gy++ ) {

        if ( gpx[gi] == gx && gpy[gi] == gy ) {

            gpout = gparam[gi];
            fprintf( gfascii, "%i ", gpout ); // writing output to ascii file
            gi += 1;
            
        } else fprintf( gfascii, "-9999.000 " ); // writing no data value to file if cell is not in area of interest
        
    } fprintf(gfascii, "\n"); }

    fclose(gfascii);

    return;
}

void foutascindf ( float *gparam, int *gpx, int *gpy, char *goutmaps, char *gname, struct ico sico ) { // function for output of ascii raster maps

    FILE *gfascii;
    char gpath[200];
    int gi, gx, gy;
    float gpout;

    sprintf(gpath, "%s%s.asc", goutmaps, gname ); // writing name of ascii file to string
    gfascii=fopen(gpath, "w"); // opening ascii file

    fprintf( gfascii, "ncols %i\n", sico.N ); // writing header to ascii file
    fprintf( gfascii, "nrows %i\n", sico.M );
    fprintf( gfascii, "xllcenter %.2f\n", sico.BDWEST );
    fprintf( gfascii, "yllcenter %.2f\n", sico.BDSOUTH );
    fprintf( gfascii, "cellsize %.2f\n", sico.CSZ );
    fprintf( gfascii, "NODATA_value %.3f\n", sico.UNDEF );

    gi = 0;
    for ( gx=0; gx<sico.M; gx++ ) { for ( gy=0; gy<sico.N; gy++ ) {

        if ( gpx[gi] == gx && gpy[gi] == gy ) {

            gpout = gparam[gi];
            fprintf( gfascii, "%.2f ", gpout ); // writing output to ascii file
            gi += 1;
            
        } else fprintf( gfascii, "-9999.000 " ); // writing no data value to file if cell is not in area of interest
        
    } fprintf(gfascii, "\n"); }

    fclose(gfascii);

    return;
}

void foutdircoord ( FILE *gf_directions, int gintmode, int *gpx, int *gpy, struct ico sico ) { // function for coordinates for output flow direction file

    int gi, git, gy0 = 0, gx, gy, intv, intv2 = 0;
    float gxmetric, gymetric;

    intv = (int)( sico.N / 75 + 1 ); // cell interval for writing output

    if ( gintmode == 1 ) { intv2 = 1; gy0 = 0; } // defining shift depending on type of parameter
    else if ( gintmode == 2 ) { intv2 = 2; gy0 = 0; }
    else if ( gintmode == 3 ) { intv2 = 2; gy0 = intv; }
    else if ( gintmode == 4 ) { intv2 = 3; gy0 = 0; }
    else if ( gintmode == 5 ) { intv2 = 3; gy0 = intv; }
    else if ( gintmode == 6 ) { intv2 = 3; gy0 = 2*intv; }

    gi = 0;
    for ( gx=0; gx<sico.M; gx+=intv ) { 
      for ( gy=gy0; gy<sico.N; gy+=(intv2*intv)) {

        git = gi;
        while (( gpx[git] != gx || gpy[git] != gy ) && git < sico.IMAX ) git+= 1;
        if ( git < sico.IMAX ) {

            gi = git;
            gxmetric = ( float ) gpy[gi] * sico.CSZ + sico.BDWEST;
            fprintf( gf_directions, "%.0f\t", gxmetric ); // writing x coordinate to file
        }
      }

      if ( sico.MODEL <= 3 ) { if ( gy0 == intv ) gy0 = 0; else gy0 = intv; }
      else { if ( gy0 == 2*intv ) gy0 = intv; else if ( gy0 == intv ) gy0 = 0; else gy0 = 2*intv; }
    }

    fprintf( gf_directions, "\n" );

    gi = 0;
    for ( gx=0; gx<sico.M; gx+=intv ) { 
      for ( gy=gy0; gy<sico.N; gy+=(intv2*intv)) {

        git = gi;
        while (( gpx[git] != gx || gpy[git] != gy ) && git < sico.IMAX ) git+= 1;
        if ( git < sico.IMAX ) {

            gi = git;
            gymetric = sico.BDNORTH - ( float ) gpx[gi] * sico.CSZ;
            fprintf( gf_directions, "%.0f\t", gymetric ); // writing y coordinate to file
        }
      }

      if ( sico.MODEL <= 3 ) { if ( gy0 == intv ) gy0 = 0; else gy0 = intv; }
      else { if ( gy0 == 2*intv ) gy0 = intv; else if ( gy0 == intv ) gy0 = 0; else gy0 = 2*intv; }
    }

    fprintf( gf_directions, "\n" );

    return;
}

void foutdir ( FILE *gf_directions, float **gparam, int gl, int gintmode, int *gpx, int *gpy, struct ico sico ) { // function for output flow direction file

    int gi, git, gy0 = 0, gx, gy, intv, intv2 = 0;
    float ghflowi = 0, gparami = 0;

    intv = (int)( sico.N / 75 + 1 ); // cell interval for writing output

    if ( gintmode == 1 ) { intv2 = 1; gy0 = 0; } // defining shift depending on type of parameter
    else if ( gintmode == 2 ) { intv2 = 2; gy0 = 0; }
    else if ( gintmode == 3 ) { intv2 = 2; gy0 = intv; }
    else if ( gintmode == 4 ) { intv2 = 3; gy0 = 0; }
    else if ( gintmode == 5 ) { intv2 = 3; gy0 = intv; }
    else if ( gintmode == 6 ) { intv2 = 3; gy0 = 2*intv; }

    gi = 0;
    for ( gx=0; gx<sico.M; gx+=intv ) {
      for ( gy=gy0; gy<sico.N; gy+=(intv2*intv)) {

        git = gi;
        while (( gpx[git] != gx || gpy[git] != gy ) && git < sico.IMAX ) git+= 1;
        if ( git < sico.IMAX ) {

            gi = git;

            if ( gintmode == 1 ) ghflowi = gparam[gi][0];
            else if ( gintmode == 2 ) ghflowi = gparam[gi][0];
            else if ( gintmode == 3 ) ghflowi = gparam[gi][3];
            else if ( gintmode == 4 ) ghflowi = gparam[gi][0];
            else if ( gintmode == 5 ) ghflowi = gparam[gi][3];
            else if ( gintmode == 6 ) ghflowi = gparam[gi][6];

            if ( ghflowi >= sico.IMPTHR[0] ) {

                if ( gintmode == 1 || gintmode == 2 ) gparami = gparam[gi][gl] / gparam[gi][0];
                else if ( gintmode == 3 ) gparami = gparam[gi][gl] / gparam[gi][3];
                else if ( gintmode == 4 ) gparami = gparam[gi][gl] / gparam[gi][0];
                else if ( gintmode == 5 ) gparami = gparam[gi][gl] / gparam[gi][3];
                else if ( gintmode == 6 ) gparami = gparam[gi][gl] / gparam[gi][6];
                fprintf( gf_directions, "%.2f\t", gparami );
            }
            else fprintf( gf_directions, "0\t" ); // writing parameter to file
        }
      }

      if ( sico.MODEL <= 3 ) { if ( gy0 == intv ) gy0 = 0; else gy0 = intv; }
      else { if ( gy0 == 2*intv ) gy0 = intv; else if ( gy0 == intv ) gy0 = 0; else gy0 = 2*intv; }
    }

    fprintf( gf_directions, "\n" );

    return;
}


#ifdef WITHGRASS


    SEGMENT finrasti ( char *gm, struct ico sico ) { // input of integer GRASS raster maps

        int gf, gx, gy, gin;
        CELL *gc;
        SEGMENT gseg_in;

        Segment_open ( &gseg_in, G_tempfile(), sico.M, sico.N, sico.NSEGRC, sico.NSEGRC, sizeof(int), sico.NSEGS );

        gc = Rast_allocate_c_buf();
        gf = Rast_open_old ( gm, sico.MAINMAPSET );

        for ( gx = 0; gx < sico.M; gx++ ) {

            Rast_get_c_row ( gf, gc, gx );

            for ( gy = 0; gy < sico.N; gy++ ) {

                if ( !Rast_is_c_null_value ( gc + gy ) ) gin = ( int ) gc[gy];
                else gin = (int)sico.UNDEF;
                Segment_put(&gseg_in, &gin, gx, gy);
            }
        }
        Segment_flush(&gseg_in);
        G_free ( gc );
        Rast_close ( gf );

        return gseg_in;
    }

    SEGMENT finrastd ( char *gm, struct ico sico ) { // input of float GRASS raster maps

        int gf, gx, gy;
        float gin;
        DCELL *gc;
        SEGMENT gseg_in;

        Segment_open ( &gseg_in, G_tempfile(), sico.M, sico.N, sico.NSEGRC, sico.NSEGRC, sizeof(float), sico.NSEGS );

        gc = Rast_allocate_d_buf();
        gf = Rast_open_old ( gm, sico.MAINMAPSET );

        for ( gx = 0; gx < sico.M; gx++ ) {

            Rast_get_d_row ( gf, gc, gx );

            for ( gy = 0; gy < sico.N; gy++ ) {

                if ( !Rast_is_d_null_value ( gc + gy ) ) gin = ( float ) gc[gy];
                else gin = sico.UNDEF;
                Segment_put(&gseg_in, &gin, gx, gy);
            }
        }

        Segment_flush(&gseg_in);
        G_free ( gc );
        Rast_close ( gf );

        return gseg_in;
    }

    void foutrast ( char *gm, float ***goutv, struct ico sico, int gk, float gtsum ) { // output of GRASS raster maps

        int gf, gx, gy;
        DCELL *gc;

        if ( gm != NULL ) {
            gc = Rast_allocate_d_buf();
            gf = Rast_open_new( gm, DCELL_TYPE );

            for ( gx = 0; gx < sico.M; gx++ ) {

                if ( gm != NULL ) {

                    for ( gy = 0; gy < sico.N; gy++ ) {

                        gc[gy] = (DCELL) goutv[gx][gy][gk];
                    }
                    Rast_put_d_row( gf, gc );
                }
            }
            Rast_close( gf );
            G_free( gc );

        }

        return;
    }


#else


    float *finaschdr ( char *gindir, char *gname ) { // function for input of ascii raster map header

        FILE *gfascii;
        char gpath[200], *tok;
        char delim[] = " \t\n";
        int gi, gx, gy;
        float gtest3, *gparam, gtestthrsh;
        char *gtest, *gtest2;
        char **strtodhlp=NULL;

        gparam = (float*) calloc( 7, sizeof(float));

        sprintf(gpath, "%s%s", gindir, gname ); // writing name of ascii file to string
        gfascii=fopen(gpath, "r"); // opening ascii file

        // Check if file exists
        if (gfascii == NULL) {
            printf("Could not open file %s", gname);
            return 0;
        } 

        for ( gi=0; gi<6; gi++ ) {
            gtest = fcparam ( gfascii );
            tok=strtok(gtest, delim);
            while ( tok != NULL ) {
                gtest2=tok;
                tok=strtok(NULL, delim);
            }
            gparam[gi] = strtod (gtest2, strtodhlp); // reading header lines and writing to array
            free(gtest);
            free(tok);
        }

        gtestthrsh = gparam[5]; // threshold for validity

        gi = 0;
        for ( gx=0; gx<(int)gparam[1]; gx++ ) { // counting number of cells

            gtest = flparam ( gfascii );

            tok = strtok(gtest, delim);
            for ( gy=0; gy<(int)gparam[0]; gy++ ) {
                gtest2 = tok;
                tok = strtok(NULL, delim);
                gtest3 = strtod (gtest2, strtodhlp);
                if ( gtest3 >= gtestthrsh ) gi += 1;
            }
            free(gtest);
            free(tok);
        }
        gparam[6] = (float)gi;
        fclose( gfascii );
        return gparam;
    }

    int *finascxy ( char *gindir, char *gname, int gcontrol, struct ico sico ) { // function for input of ascii raster x and y coordinates

        FILE *gfascii;
        char *gpath, *tok;
        char delim[] = " \t\n";
        int gi, gx, gy, *gpxy;
        float gtest3, gtestthrsh;
        char **strtodhlp=NULL;

        char *gtest;
        char *gtest2;

        gpath = (char*) calloc ( 1000, sizeof(char));
        gpxy = (int*) calloc ( sico.IMAX, sizeof(int));

        gtestthrsh = sico.UNDEF; // threshold for validity

        sprintf(gpath, "%s%s", gindir, gname ); // writing name of ascii file to string
        gfascii=fopen(gpath, "r"); // opening ascii file

        // Check if file exists
        if (gfascii == NULL) {
            printf("Could not open file %s", gname);
            return 0;
        } 

        for ( gi=0; gi<6; gi++ ) {
            gtest = fcparam ( gfascii ); // reading header lines
            free(gtest);
        }

        gi = 0;
        for ( gx=0; gx<sico.M; gx+=sico.MESH ) {

            gtest = flparam ( gfascii );
            tok = strtok(gtest, delim);
            for ( gy=0; gy<sico.N; gy+=sico.MESH ) {

                gtest2 = tok;
                tok = strtok(NULL, delim);
                gtest3 = strtod (gtest2, strtodhlp);

                if ( gtest3 >= gtestthrsh ) {

                    if ( gcontrol == 1 ) gpxy[gi] = gx; // writing internal coordinates of all rows and columns to array
                    else gpxy[gi] = gy;

                    gi += 1;
                }
            }
            free(gtest);
            free(tok);
        }
        fclose( gfascii );
        free(gpath);
        return gpxy;
    }

    float *finascval ( char *gindir, char *gname, int *gpx, int *gpy, struct ico sico ) { // function for input of ascii raster map values

        FILE *gfascii;
        char *gpath, *tok;
        char delim[] = " \t\n";
        int gi, gx, gy;
        float gtest3, *gparam;
        char **strtodhlp=NULL;

        char *gtest;
        char *gtest2;

        gpath = (char*) calloc ( 1000, sizeof(char));
        gparam = (float*) calloc( sico.IMAX, sizeof(float));

        sprintf(gpath, "%s%s", gindir, gname ); // writing name of ascii file to string
        gfascii=fopen(gpath, "r"); // opening ascii file

        // Check if file exists
        if (gfascii == NULL) {
            printf("Could not open file %s", gname);
           return 0;
        } 
 
        for ( gi=0; gi<6; gi++ ) {
            gtest = fcparam ( gfascii ); // reading header lines
            free(gtest);
        }

        gi = 0;
        for ( gx=0; gx<sico.M; gx+=sico.MESH ) {

            gtest = flparam ( gfascii );
            tok = strtok(gtest, delim);
            for ( gy=0; gy<sico.N; gy+=sico.MESH ) {

                gtest2 = tok;
                tok = strtok(NULL, delim);
                gtest3 = strtod(gtest2, strtodhlp);

                if ( gx == gpx[gi] && gy == gpy[gi] ) {
                    gparam[gi] = gtest3; // writing raster values of all rows and columns to array
                    gi += 1;
                }
            }
            free(gtest);
            free(tok);
        }
        fclose( gfascii );
        free(gpath);
        return gparam;
    }

    int *finascvali ( char *gindir, char *gname, int *gpx, int *gpy, struct ico sico ) { // function for input of ascii raster map values

        FILE *gfascii;
        char *gpath, *tok;
        char delim[] = " \t\n";
        int gi, gx, gy;
        int gtest3, *gparam;
        char **strtodhlp=NULL;

        char *gtest;
        char *gtest2;

        gpath = (char*) calloc ( 1000, sizeof(char));
        gparam = (int*) calloc( sico.IMAX, sizeof(int));

        sprintf(gpath, "%s%s", gindir, gname ); // writing name of ascii file to string
        gfascii=fopen(gpath, "r"); // opening ascii file
        
        // Check if file exists
        if (gfascii == NULL) {
            printf("Could not open file %s", gname);
           return 0;
        }       
      
        for ( gi=0; gi<6; gi++ ) {
            gtest = fcparam ( gfascii ); // reading header lines
            free(gtest);
        }

        gi = 0;
        for ( gx=0; gx<sico.M; gx+=sico.MESH ) {

            gtest = flparam ( gfascii );
            tok = strtok(gtest, delim);
            for ( gy=0; gy<sico.N; gy+=sico.MESH ) {

                gtest2 = tok;
                tok = strtok(NULL, delim);
                gtest3 = strtod(gtest2, strtodhlp);

                if ( gx == gpx[gi] && gy == gpy[gi] ) {
                    gparam[gi] = gtest3; // writing raster values of all rows and columns to array
                    gi += 1;
                }
            }
            free(gtest);
            free(tok);
        }
        fclose( gfascii );
        free(gpath);

        return gparam;
    }


#endif


float *finlistd ( char *gtest ) { // function for input of comma-delimited multiple floating point parameters

    char delim[] = ",";
    int gi;
    float *gparam;
    char *gtest2, *tok;
    char **strtodhlp=NULL;

    gparam = (float*) calloc( 500, sizeof(float));

    tok=strtok(gtest, delim);

    gi=1;
    while ( tok != NULL ) {
        gtest2=tok;
        tok=strtok(NULL, delim);
        gparam[gi] = strtod (gtest2, strtodhlp); // storing component of parameter string
        gi += 1;
    }
    gparam[0] = (float)(gi-1);
        
    return gparam;
}

int fcountlines( char *gindir, char *gfilename ) { // https://www.geeksforgeeks.org/c-program-count-number-lines-file/

FILE *fp;
int count = 0; // Line counter (result)
char *gpath, c;

    gpath = (char*) calloc ( 1000, sizeof(char));

    // Open the file
    sprintf(gpath, "%s%s", gindir, gfilename ); // writing name of text file to string
    fp = fopen(gpath, "r");

    // Check if file exists
    if (fp == NULL) {
        printf("Could not open file %s", gfilename);
        return 0;
    }

    // Extract characters from file and store in character c
    for (c = getc(fp); c != EOF; c = getc(fp))
        if (c == '\n') // Increment count if this character is newline
    count = count + 1;

    // Close the file
    fclose(fp);

    free(gpath);
    return count;
}


// -- STOP --- Functions ----------------------------------------------------------------------------------------


int main ( int argc, char *argv[] ) { // calling main function


// -- START -- Declaration of variables -------------------------------------------------------------------------


    FILE *fparam, *f_paramcomm, *f_mapprof;

    time_t timer = 0;

    char *prefix, *outmaps, *outfiles, *outaimec, *outvr, *outrplots, *parapy, *rscript = 0, *rlibs = 0, *proflist0, *ctrlplist0, *hydrocoordslist0, *hydrographslist0, *hydrographslist[100], 
        *hydreleasename[100], *elevname, *hreleasename = 0, *hreleasename2 = 0, *hreleasename3 = 0, *hentrmaxname = 0, *hentrmaxname2 = 0, *hentrmaxname3 = 0, *vinxname = 0, *vinyname = 0, 
        *vinxname2 = 0, *vinyname2 = 0, *vinxname3 = 0, *vinyname3 = 0, *zonesname = 0, *centrname = 0, *phiname = 0, *phi2name = 0, *phi3name = 0, *deltaname = 0, *delta2name = 0, *delta3name = 0, 
        *addfriname = 0, *addfri2name = 0, *addfri3name = 0, *cohesionname = 0, *cohesion2name = 0, *cohesion3name = 0, *nyname = 0, *ny2name = 0, *ny3name = 0, *cdeformname = 0, 
        *ambdragname = 0, *trans12name = 0, *trans13name = 0, *trans23name = 0, *zfragname = 0, *adaptoname = 0, *frictioname = 0, *transformoname = 0, *treleasename = 0, *trelstopname = 0, 
        *stoptimename = 0, *tslidename = 0, *impactareaname = 0, *hdepositname = 0, *pbg1name = 0, *pbg2name = 0, *pbg3name = 0, *mv, *mv2, **mv0; char *madd = (char*) calloc(1000, sizeof(char));

    char *path = (char*) calloc(1000, sizeof(char));
    char *wkdir = (char*) calloc(800, sizeof(char));
    char *indir = (char*) calloc(800, sizeof(char));

    int x, y, *px, *py, *iin, *inn, *innn, *cin, hydrograph, adaptograph, frictiograph, transformograph, hydnum = 0, hydnin = 0, hydnout, *hydi = 0, hydj = 0, hydk = 0, *hydx = 0, *hydy = 0, hydt = 0, 
        *hydtmax = 0, hydtmaxx = 0, *hydp0 = 0, **hydp = 0, adatmax = 0, adak, adat = 0, fritmax = 0, frik, frit = 0, tratmax = 0, trak, tratx = 0, nvect_all, nvect_red, lmax, xint, ccontinue, csuccess, 
        nsum, nout, time_start, time_stop, ctrlr, ctrlv, ctrlvv, ctrlvvv, ccfl, fj[4][2], i, j, jj, jjj, jmin = 0, k = 0, l, ll, p, z, prec_hflow, prec_vol, prec_ekin, nzones, ctrl_trelease,
        imax = 0, iloop, ctrl_hydout, ix, hydcols, ctrl_release, ctrl_basechange, cflowpre, anid, andist, ansteps, ctr_ansteps, iz, 
        highpass_min, highpass_max, anctr0, cslide, cflow;

    float *pelev, *qelev, *relev, gkx[3], gky[3], *kappau, *vm, *cdrag, *disp, *gze, *gf, *gg, *gs, *flowpar, tout, tmax, vflowx = 0, vflowy = 0, vflowx1, vflowy1, vflowx2, vflowy2, 
        vflowx3, vflowy3, time_elapsed, tlength, tlength0, tlengthx, tlengthpre, tint, tsum, cfl, cflmax, grav[21], dw_dt[10], hflow_maxmax, hflow_max0, hflow_max, hflow_max02, 
        hflow_max2, hflow_max03, hflow_max3, vflow_max, vflow_max2, vflow_max3, **hydhyd0 = 0, ***hydhyd = 0, **adaada = 0, **frifri = 0, **tratra = 0, vol_flow0, vol_flow, vol_flow02, vol_flow2, 
        vol_flow03, vol_flow3, vol_entr, vol_entr2, vol_entr3, vol_edge = 0, vol_edge2 = 0, vol_edge3 = 0, whx, why, whx1 = 0, why1 = 0, whx2 = 0, why2 = 0, vcelr0, vcelr, ctrl_zfrag, *hydelev = 0, 
        *hydalpha = 0, *hydx_metric = 0, *hydy_metric = 0, *hydl = 0, 
        hhyd0 = 0, hhyd = 0, hyde = 0, hyde2 = 0, hyde3 = 0, hydfalpha = 0, hydout_xmin_metric = 0, hydout_xmax_metric = 0, hydout_ymin_metric = 0, hydout_ymax_metric = 0, 
        hydh = 0, hydh2 = 0, hydh3 = 0, hydv = 0, hydv2 = 0, hydv3 = 0, hydq = 0, hydq2 = 0, hydq3 = 0, hydlmax = 0, ekin, ekin_flow, vol_noflux, vol_noflux2, vol_noflux3, carea = 0, 
        hydmalpha, hydnalpha, hydbeta, hydm0, hydmx, hydmy, hydm, hydfcorr, wu[9], wv[9], wu2[9], wv2[9], wu3[9], wv3[9], nbetax[9], nbetay[9], walpha[9], walpha2[9], walpha3[9], 
        vflowxj[9], vflowyj[9], corrfact, qentr, qentr1, qentr2, qentr3, qmelt = 0, qmelt1 = 0, qmelt2 = 0, qmelt3 = 0,
        alpha, alphav, mom, walphax, walphay, walphax2, walphay2, walphax3, walphay3, wbetax, wbetay, wbetaxh, wbetayh, wbetaxy, wbetax2 = 0, wbetay2 = 0, wbetaxy2 = 0, wbetaxh2 = 0, wbetayh2 = 0,
        wbetax3 = 0, wbetay3 = 0, wbetaxy3 = 0, wbetaxh3 = 0, wbetayh3 = 0, wdx, wdy, wgrav[21], vol_hydbef, vol_hydaft, pelevhydtest, vol_hydbef2, vol_hydaft2, vol_hydbef3, vol_hydaft3, 
        vol_hyd, vol_hyd2, vol_hyd3, *gdecel, hflow = 0, hflown = 0, trat, ttrelease, ttrelstop, tfact, trans, wdu[3], wdv[3], xwdu[3], xwdv[3], welev[9], wh[9], wh1[9], wh2[9], 
        whflow = 0, whflow2 = 0, whflow3 = 0, dw_dttest[10], wintbtest, wintctest, awttest[10], hekin, hekin_sum, hekin_max, hydbef[10],
        rhrem, qvol_test, qh_test, qtrelstart = 0, qtrelstop = 0, qtrelspan, ctrl_pressthr, elevtot[9], betav,
        alpha1 = 0, alpha2 = 0, alpha3 = 0, alphab1 = 0, alphab2 = 0, alphab3 = 0, momaddsx, momaddfsx, momaddfx, momaddsy, momaddfsy, momaddfy, momfact = 0, hflowj[9], hentrmaxx,
        phreleaseall = 0, qhreleaseall, vol_flow0all, anslopex[3], anslopey[3], anslopexy[3], angx[6], angy[6], *anpx = 0, *anpy = 0, *anctr = 0, anpx1, anpy1, anpx2, anpy2, anwht[4], ansumh[3], 
        ansumslopex[3], ansumslopey[3], ansumslopexy[3], anavgslopex[3], anavgslopey[3], anavgslopexy[3], kmax = 3, hflownn[3],  
        tflow_maxmax1 = 0, tflow_maxmax2 = 0, tflow_maxmax3 = 0, tflow_maxmax = 0, pflow_maxmax1 = 0, pflow_maxmax2 = 0, 
        pflow_maxmax3 = 0, pflow_maxmax = 0, treachmaxmax = 0, htsunmaxmax = 0, vflow_maxmax1 = 0, vflow_maxmax2 = 0, vflow_maxmax3 = 0, vflow_maxmax = 0, basechange_max = 0, 
        basechange_min = 0, rhrelease1, rhentrmax1, vhrelease, vhentrmax, phexagg = 1.0, hrelease = 0, slidefract = 0, highpass, highpass_halfdiff, highpass_sum, highpass_avg, wtemp;

    float *proflist = 0, *ctrlplist = 0, *profelev = 0, *profdeposit = 0, ***profdata = 0, **profeval = 0, *profnx = 0, *profny = 0, profdiff, profdiffnx, profdiffny, profwhtx1, profwhtx2, 
        profwhty1, profwhty2, *profdiffabs = 0, profhrelease = 0, proflength_hmax = 0, proflength_impactarea = 0, proflength_hdeposit = 0, proflength_simulated = 0, proflength_hmax_ctrl = 0, 
        profdiff_impactarea, profratio_impactarea, profdiff_hdeposit, profratio_hdeposit, evalrp_imp = 0, evalrn_imp = 0, evalrtp_imp = 0, evalrtn_imp = 0, evalrfp_imp = 0, evalrfn_imp = 0, 
        evalrtpp_imp = 0, evalrfpp_imp = 0, evalrtnn_imp = 0, evalrpp_imp = 0, evalrnp_imp = 0, evalrt_imp = 0, evalrf_imp = 0, evalcsi_imp = 0, evalhss_imp = 0, evalauroc_imp = 0, evald2pc_imp = 0, 
        evalfoc_imp = 0, evalspi_imp = 0, evalrp_dep = 0, evalrn_dep = 0, evalrtp_dep = 0, evalrtn_dep = 0, evalrfp_dep = 0, evalrfn_dep = 0, evalrtpp_dep = 0, evalrfpp_dep = 0, evalrtnn_dep = 0, 
        evalrpp_dep = 0, evalrnp_dep = 0, evalrt_dep = 0, evalrf_dep = 0, evalcsi_dep = 0, evalhss_dep = 0, evalauroc_dep = 0, evald2pc_dep = 0, evalfoc_dep = 0, evalspi_dep = 0,
        profxm, profym, *ctrlpxm = 0, *ctrlpym = 0, ***ctrlpdata = 0, *hydrocoordslist = 0, randnum2 = 0, randnum3 = 0, randnum4 = 0, topofact = 1.0;
    int profi, profj, profk = 0, profl, profctrl = 0, *profx = 0, *profy = 0, profn, profdiffx, profdiffy, profmax = 0, evaltp_imp = 0, evaltn_imp = 0, evalfp_imp = 0, evalfn_imp = 0,
        evaltp_dep = 0, evaltn_dep = 0, evalfp_dep = 0, evalfn_dep = 0, evalall = 0, *ctrlpx = 0, *ctrlpy = 0, randmax = 0;

    int paradistort = 0, paracontourshmin = 0, paracontourshmax = 0, paracontourshint = 0, paracontourszmin = 0, paracontourszmax = 0, paracontourszint = 0, paraimp;    
    float paramin = 0, pararef = 0, paratsunref = 0, paraxmetric = 0, paraymetric = 0, paraelev = 0, parahflow = 0, parahrelease = 0, parah = 0, parahs = 0, parahmax = 0, 
        parahflow1 = 0, parahflow2 = 0, parahflow3 = 0, paraalpha = 0, paraalphamax = 0, parafactr = 0, parafactg = 0, parafactb = 0, paracorrtsun = 0, 
        paraaddtsun = 0, paracolfactr = 0, paracolfactb = 0, paracolfactg = 0, parar = 0, parab = 0, parag = 0, parad = 0, parared = 0, paragreen = 0, parablue = 0;

    struct ico sico;
    struct flow sflow;


    #ifdef WITHGRASS


        SEGMENT seg_elev, seg_hrelease, seg_hrelease2, seg_hrelease3, seg_vinx, seg_viny, seg_vinx2, seg_viny2, seg_vinx3, seg_viny3,
        seg_hentrmax, seg_hentrmax2, seg_hentrmax3, seg_zones, seg_centr, seg_phi, seg_phi2, seg_phi3, seg_delta, seg_delta2, seg_delta3, seg_addfri, seg_addfri2, seg_addfri3,
        seg_cohesion, seg_cohesion2, seg_cohesion3, seg_ny, seg_ny2, seg_ny3, seg_cdeform, seg_ambdrag, seg_trans12, seg_trans13, seg_trans23, seg_zfrag, seg_trelease, seg_trelstop, 
        seg_stoptime, seg_tslide, seg_impactarea, seg_hdeposit, seg_pbg1, seg_pbg2, seg_pbg3;
        char *yint = (char*) calloc(1000, sizeof(char));
        int aoi, zones, impactarea, pbg1, pbg2, pbg3;
        float pout, elev, hrelease2, hrelease3, vinx, viny, vinx2, viny2, vinx3, viny3,
        hentrmax, hentrmax2, hentrmax3, zfrag, centr, phi, phi2, phi3, delta, delta2, delta3, addfri, addfri2, addfri3, cohesion, cohesion2, cohesion3, ny, ny2, ny3, cdeform, ambdrag, 
        trans12, trans13, trans23, trelease, trelstop, stoptime, tslide, hdeposit;

        struct GModule *module;
        struct Cell_head cellhd;
        struct Option *input_opt1;


    #else


        float *aschdr;


    #endif


// -- STOP --- Declaration of variables -------------------------------------------------------------------------


    time_start = clock(); // start time


// Preparing environment


    #ifdef WITHGRASS


        xint=1;
        yint=getenv("XINT");
        if ( yint != NULL ) xint=atoi(yint);
        
        G_gisinit( argv[0] );

        module = G_define_module();
        G_add_keyword(_("raster"));
        G_add_keyword(_("landslide"));
        G_add_keyword(_("numerical simulation"));
        module->description = _("The mass flow simulation tool");

        sprintf( wkdir, "%s/%s/.tmp/rtemp", G_location_path(),G_mapset());
        sprintf( path, "%s/param%d.txt", wkdir, xint ); // writing name of parameter file to string

        input_opt1 = G_define_standard_option(G_OPT_F_BIN_INPUT);
        input_opt1->key = "input1";
        input_opt1->description = _("Parameter file 1");
        input_opt1->required = NO;
        input_opt1->answer = path;

        if (G_parser(argc, argv))
            exit(EXIT_FAILURE);

        G_get_set_window( &cellhd );
        
        G_init_tempfile();
        sico.NSEGRC = 16;
        sico.NSEGS = 16;

        fparam=fopen( input_opt1->answer, "r" ); // opening parameter file


    #else


        if ( argc > 1 ) sprintf( wkdir, "%s", argv[1] );
        else {
            printf( "ERROR: Please provide a parameter file.\n" );
            fflush(stdout);
            exit( EXIT_SUCCESS );
        }
        xint=1;

        sprintf( path, "%s/param%d.txt", wkdir, xint ); // writing name of parameter file to string
        fparam=fopen( path, "r" ); // opening parameter file


    #endif


// -- START -- Reading and preprocessing parameter file input ---------------------------------------------------


    sico.PI = 3.1415926536; // pi
    
    sico.ELEV = 0; sico.RELM = 0; sico.RELM2 = 0; sico.RELM3 = 0; sico.RELV = 0; sico.RELV2 = 0; sico.RELV3 = 0; sico.ENTR = 0; sico.ENTR2 = 0; sico.ENTR3 = 0; 
    sico.ZONES = 0; sico.CENTR = 0; sico.PHI = 0; sico.PHI2 = 0; sico.PHI3 = 0; sico.DELTA = 0; sico.DELTA2 = 0; sico.DELTA3 = 0; sico.ADDFRI = 0; sico.ADDFRI2 = 0; sico.ADDFRI3 = 0;  
    sico.COHESION = 0; sico.COHESION2 = 0; sico.COHESION3 = 0;sico.NY = 0; sico.NY2 = 0; sico.NY3 = 0; sico.CDEFORM = 0; sico.AMBDRAG = 0; sico.TRANS12 = 0; sico.TRANS13 = 0; sico.TRANS23 = 0; 
    sico.ZFRAG = 0; sico.TRELEASE = 0; sico.TRELSTOP = 0; sico.STOPTIME = 0; sico.TSLIDE = 0; sico.IMPACTAREA = 0; sico.HDEPOSIT = 0; sico.PBG = 0; // initializing controls for raster maps
    hydrograph = 0; adaptograph = 0; frictiograph = 0; transformograph = 0; sico.PROFILE = 0; sico.CTRLPOINTS = 0;
        // initializing identifiers for use of hydrograph(s), adaptograph, frictiograph, transformograph, profile, and control points (1) or not (0)

    if( !fparam ) { // managing lack of parameter file
        printf( "ERROR: Unable to open parameter file: '%s'\n", path );
        fflush(stdout);
        exit( EXIT_SUCCESS );
    }

    sico.MAINMAPSET = fcparam ( fparam ); // name of main mapset (relevant with GRASS)
    sico.MULT = fiparam ( fparam ); // identifier for single model run (0) or multiple model runs (1)
    
    if ( sico.MULT == 0 ) {

        printf("\nStarting model execution.\n\n");
        fflush(stdout); // forcing immediate display
    }    
    
    sico.MESH = fiparam ( fparam ); // mesh spacing (multiple of ascii raster cell size)


    // Pathes for input and output

    prefix = fcparam ( fparam ); // prefix for output
    indir = fcparam ( fparam ); // name of directory with input files
    outmaps = fcparam ( fparam ); // path and prefix for storing output ascii raster maps
    outfiles = fcparam ( fparam ); // path and prefix for storing output text files
    outaimec = fcparam ( fparam ); // path and prefix for storing output aimec files
    outvr = fcparam ( fparam ); // path and prefix for storing output virtual reality files
    outrplots = fcparam ( fparam ); // path and prefix for storing output R plots
    

    #ifndef WITHGRASS


        sprintf( path, "%s/%sresults/", wkdir, prefix );
        mkdir( path, 0700 );
        sprintf( outmaps, "%s/%sresults/%sascii/", wkdir, prefix, prefix );
        mkdir( outmaps, 0700 );
        sprintf( outfiles, "%s/%sresults/%sfiles/", wkdir, prefix, prefix );
        mkdir( outfiles, 0700 );
        sprintf( outaimec, "%s/%sresults/%saimec/", wkdir, prefix, prefix );
        mkdir( outaimec, 0700 );
        sprintf( outvr, "%s/%sresults/%svr/", wkdir, prefix, prefix );
        mkdir( outvr, 0700 );
        sprintf( path, "%s/%sresults/%svr/data/", wkdir, prefix, prefix );
        mkdir( path, 0700 );
        sprintf( outrplots, "%s/%sresults/%splots/", wkdir, prefix, prefix );
        mkdir( outrplots, 0700 );
        sprintf( path, "%s/%sresults/%splots/%sprofiles_timesteps/", wkdir, prefix, prefix, prefix );
        mkdir( path, 0700 );
        sprintf( path, "%s/%sresults/%splots/%smaps_timesteps/", wkdir, prefix, prefix, prefix );
        mkdir( path, 0700 );
        sprintf( path, "%s/%sresults/%saimec/depth/", wkdir, prefix, prefix );
        mkdir( path, 0700 );
        sprintf( path, "%s/%sresults/%saimec/pressure/", wkdir, prefix, prefix );
        mkdir( path, 0700 );


    #endif


    if ( sico.MULT == 0 ) sprintf(path, "%s%sparamcomm.txt", outfiles, prefix); // commented parameter file
    else sprintf(path, "%s%sparamcomm%d.txt", outfiles, prefix, xint);
    f_paramcomm=fopen(path, "w");

    fprintf( f_paramcomm, "Name of main mapset (relevant with GRASS)\t%s\n", sico.MAINMAPSET );
    fprintf( f_paramcomm, "Identifier for single model run (0) or multiple model runs (1)\t%i\n", sico.MULT );
    fprintf( f_paramcomm, "Mesh spacing (multiple of ascii raster cell size)\t%i\n", sico.MESH );
    fprintf( f_paramcomm, "Prefix for output\t%s\n", prefix );
    fprintf( f_paramcomm, "Name of directory with input files\t%s\n", indir );
    fprintf( f_paramcomm, "Path and prefix for storing output ascii raster maps\t%s\n", outmaps );
    fprintf( f_paramcomm, "Path and prefix for storing output text files\t%s\n", outfiles );
    fprintf( f_paramcomm, "Path and prefix for storing output aimec files\t%s\n", outaimec );
    fprintf( f_paramcomm, "Path and prefix for storing output virtual reality files\t%s\n", outvr );
    fprintf( f_paramcomm, "Path and prefix for storing output R plots\t%s\n", outrplots );


    // General controls
 
    sico.PMAX = fiparam ( fparam ); // number of phases
    fprintf( f_paramcomm, "Number of phases\t%i\n", sico.PMAX );
      
    if ( sico.PMAX == 1 ) sico.MODEL = 1; else sico.MODEL = 7; // model (1=one-phase, 7=multi-phase)

    if ( sico.MODEL <= 3 ) { // size of arrays of state variables

        nvect_all = 20;
        nvect_red = 8;
        sico.NVECTMIN = 4;
        
    } else {

        nvect_all = 49;
        nvect_red = 26;
        sico.NVECTMIN = 10;
    }

    sico.AFLAG = fiparam ( fparam ); // control for additional output raster maps (0=no, 1=yes)
    fprintf( f_paramcomm, "Control for additional output raster maps (0=no, 1=yes)\t%i\n", sico.AFLAG );    

    sico.TSUNAMI = fiparam ( fparam ); // tsunami mode (0=no, 1=yes)
    fprintf( f_paramcomm, "Tsunami mode (0=no, 1=yes)\t%i\n", sico.TSUNAMI );
    
    sico.GRAVITY = fdparam ( fparam ); // gravity (m/s²)
    fprintf( f_paramcomm, "Gravity (m/s²)\t%.2f\n", sico.GRAVITY );    
    
    sico.LIMITER = fiparam ( fparam ); // numerical limiter (1=Minmod, 2=Superbee, 3=Woodward, 4=van Leer)
    fprintf( f_paramcomm, "Numerical limiter (1=Minmod, 2=Superbee, 3=Woodward, 4=van Leer)\t%i\n", sico.LIMITER );     
    
    sico.LAYERS0 = fiparam ( fparam ); // layer mode (0=off, 1=on)
    fprintf( f_paramcomm, "Layer mode (0=no, 1=yes)\t%i\n", sico.LAYERS0 );    
 
    sico.TOPOFACT = fiparam ( fparam ); // input and output of flow heights and velocities (0=vertical, 1=topography-following)
    fprintf( f_paramcomm, "Input and output of flow heights and velocities (0=vertical, 1=topography-following)\t%i\n", sico.TOPOFACT );
  
    sico.SURFACE = fiparam ( fparam ); // surface control (0=off, 1=edge control, 2=force balancing, 3=all)
    fprintf( f_paramcomm, "Surface control (0=no, 1=edge control, 2=force balancing, 3=all)\t%i\n", sico.SURFACE );     
    
    sico.ENTRAINMENT = fiparam ( fparam ); // entrainment control (0=off, 1=on)
    fprintf( f_paramcomm, "Entrainment and deposition control (0=no, 1-4=type of approach)\t%i\n", sico.ENTRAINMENT );      
    
    sico.STOPPING = fiparam ( fparam ); // stopping control (0=off, 1-3=type of approach)
    fprintf( f_paramcomm, "Stopping control (0=no, 1-3=type of approach)\t%i\n", sico.STOPPING );      
    
    sico.NONHYDRO = fiparam ( fparam ); // control for non-hydrostatic effects (0=no, 1=dispersion, 2=enhanced gravity, 2=dispersion and enhanced gravity)
    fprintf( f_paramcomm, "Control for non-hydrostatic effects (0=no, 1=dispersion, 2=enhanced gravity, 2=dispersion and enhanced gravity)\t%i\n", sico.NONHYDRO );     

    sico.TEMP = fiparam ( fparam ); // control for temperature evolution and ice melt (0=off, 1=on)
    fprintf( f_paramcomm, "Control for temperature evolution and ice melt (0=off, 1=on)\t%i\n", sico.TEMP );    

    sico.CURVCTRL = 1; // curvature control (0=off, 1=on)    
    sico.SEPFLUX = 0; // control for phase separation (0=no, 1=yes)
    sico.HYDADD = 2; // control for hydrograph management ( 0=reset flow, 1=impose on flow, 2=impose on centre)
    sico.ORIGINAL = 0; // control for deceleration management (0=suppress friction- or viscosity-induced change of direction, 1=do not suppress)


    // Elevation raster map

    elevname = fcparam ( fparam ); // name of elevation map
    if ( strcmp( elevname, "None" ) != 0 ) sico.ELEV = 1; 
    fprintf( f_paramcomm, "Elevation map\t%s\n", elevname );


    // Setting spatial extent and cell size


    #ifdef WITHGRASS // if GRASS rasters are used


        sico.CSZ = cellhd.ew_res; // cell size
        sico.BDWEST = cellhd.west; // western boundary
        sico.BDNORTH = cellhd.north; // northern boundary
        sico.BDSOUTH = cellhd.south; // southern boundary
        sico.BDEAST = cellhd.east; // eastern boundary
        sico.N = cellhd.cols; // number of cells in y direction (W-E)
        sico.M = cellhd.rows; // number of cells in x direction (S-N)
        sico.UNDEF = -9999; // no data value


    #else // if ascii rasters are used


        aschdr = finaschdr( indir, elevname ); // reading header of ascii raster
        sico.N = (int)aschdr[0] / sico.MESH; // number of cells in y direction (W-E)
        sico.M = (int)aschdr[1] / sico.MESH; // number of cells in x direction (S-N)
        sico.CSZ = aschdr[4] * sico.MESH; // cell size
        sico.BDWEST = aschdr[2] - sico.CSZ / 2; // western boundary
        sico.BDSOUTH = aschdr[3] - sico.CSZ / 2; // southern boundary
        sico.BDNORTH = aschdr[3] + sico.M * sico.CSZ; // northern boundary
        sico.BDEAST = aschdr[2] + sico.N * sico.CSZ; // eastern boundary
        sico.UNDEF = aschdr[5]; // no data value
        sico.IMAX = (int)aschdr[6] / pow( sico.MESH, 2 ); // number of cells

        free(aschdr);


    #endif


    // Raster maps

    hreleasename = fcparam ( fparam ); // name of MIXTURE or PHASE 1 release height map
    if ( strcmp( hreleasename, "None" ) != 0 ) sico.RELM = 1; 
    fprintf( f_paramcomm, "Release mass MIXTURE or PHASE 1 map\t%s\n", hreleasename );
        
    hreleasename2 = fcparam ( fparam ); // name of PHASE 2 release height map
    if ( strcmp( hreleasename2, "None" ) != 0 ) sico.RELM2 = 1; 
    fprintf( f_paramcomm, "Release mass PHASE 2 map\t%s\n", hreleasename2 );
    
    hreleasename3 = fcparam ( fparam ); // name of PHASE 3 release height map
    if ( strcmp( hreleasename3, "None" ) != 0 ) sico.RELM3 = 1; 
    fprintf( f_paramcomm, "Release mass PHASE 3 map\t%s\n", hreleasename3 );    

    rhrelease1 = fdparam ( fparam ); // fraction of PHASE 1 release height
    fprintf( f_paramcomm, "Fraction of PHASE 1 release height\t%.2f\n", rhrelease1 );    
    
    vhrelease = fdparam ( fparam ); // variation of release height 
    fprintf( f_paramcomm, "Variation of release height\t%.2f\n", vhrelease );

    treleasename = fcparam ( fparam ); // name of release time map
    if ( strcmp( treleasename, "None" ) != 0 ) sico.TRELEASE = 1;
    fprintf( f_paramcomm, "Release time map\t%s\n", treleasename );

    trelstopname = fcparam ( fparam ); // name of release stop time map
    if ( strcmp( trelstopname, "None" ) != 0 ) sico.TRELSTOP = 1;
    fprintf( f_paramcomm, "Release stop time map\t%s\n", trelstopname );

    stoptimename = fcparam ( fparam ); // name of stopping time map
    if ( strcmp( stoptimename, "None" ) != 0 ) sico.STOPTIME = 1;
    fprintf( f_paramcomm, "Stopping time map\t%s\n", stoptimename );

    tslidename = fcparam ( fparam ); // name of time until which to apply global gravity map
    if ( strcmp( tslidename, "None" ) != 0 ) sico.TSLIDE = 1;
    fprintf( f_paramcomm, "Time of initial sliding map\t%s\n", tslidename );
   
    zonesname = fcparam ( fparam ); // name of zones map
    if ( strcmp( zonesname, "None" ) != 0 ) sico.ZONES = 1;
    fprintf( f_paramcomm, "Zones map\t%s\n", zonesname );

    impactareaname = fcparam ( fparam ); // name of observed impact area map
    if ( strcmp( impactareaname, "None" ) != 0 ) sico.IMPACTAREA = 1;
    fprintf( f_paramcomm, "Observed impact area map\t%s\n", impactareaname );

    hdepositname = fcparam ( fparam ); // name of observed height of deposit map
    if ( strcmp( hdepositname, "None" ) != 0 ) sico.HDEPOSIT = 1;
    fprintf( f_paramcomm, "Observed height of deposit map\t%s\n", hdepositname );
   
    vinxname = fcparam ( fparam ); // name of MIXTURE or PHASE 1 release velocity in x direction map
    if ( strcmp( vinxname, "None" ) != 0 ) sico.RELV = 1; 
    fprintf( f_paramcomm, "Release velocity x MIXTURE or PHASE 1 map\t%s\n", vinxname );    
    
    vinxname2 = fcparam ( fparam ); // name of PHASE 2 release velocity in x direction map
    if ( strcmp( vinxname2, "None" ) != 0 ) sico.RELV2 = 1; 
    fprintf( f_paramcomm, "Release velocity x PHASE 2 map\t%s\n", vinxname2 );      
    
    vinxname3 = fcparam ( fparam ); // name of PHASE 3 release velocity in x direction map
    if ( strcmp( vinxname3, "None" ) != 0 ) sico.RELV3 = 1; 
    fprintf( f_paramcomm, "Release velocity x PHASE 3 map\t%s\n", vinxname3 );      
    
    vinyname = fcparam ( fparam ); // name of MIXTURE or PHASE 1 release velocity in y direction map
    if ( strcmp( vinyname, "None" ) != 0 ) sico.RELV = 1; 
    fprintf( f_paramcomm, "Release velocity y MIXTURE or PHASE 1 map\t%s\n", vinyname );      
    
    vinyname2 = fcparam ( fparam ); // name of PHASE 2 release velocity in y direction map
    if ( strcmp( vinyname2, "None" ) != 0 ) sico.RELV2 = 1; 
    fprintf( f_paramcomm, "Release velocity y PHASE 2 map\t%s\n", vinyname2 );    
    
    vinyname3 = fcparam ( fparam ); // name of PHASE 3 release velocity in y direction map
    if ( strcmp( vinyname3, "None" ) != 0 ) sico.RELV3 = 1; 
    fprintf( f_paramcomm, "Release velocity y PHASE 3 map\t%s\n", vinyname3 );  
    
    hentrmaxname = fcparam ( fparam ); // name of maximum height of MIXTURE or PHASE 1 entrainment map
    if ( strcmp( hentrmaxname, "None" ) != 0 ) sico.ENTR = 1; 
    fprintf( f_paramcomm, "Maximum height of entrainment MIXTURE or PHASE 1 map\t%s\n", hentrmaxname );    
    
    hentrmaxname2 = fcparam ( fparam ); // name of maximum height of PHASE 2 entrainment map
    if ( strcmp( hentrmaxname2, "None" ) != 0 ) sico.ENTR2 = 1;
    fprintf( f_paramcomm, "Maximum height of entrainment PHASE 2 map\t%s\n", hentrmaxname2 );    
    
    hentrmaxname3 = fcparam ( fparam ); // name of maximum height of PHASE 3 entrainment map
    if ( strcmp( hentrmaxname3, "None" ) != 0 ) sico.ENTR3 = 1;
    fprintf( f_paramcomm, "Maximum height of entrainment PHASE 3 map\t%s\n", hentrmaxname3 );    

    rhentrmax1 = fdparam ( fparam ); // fraction of PHASE 1 maximum height of entrainment
    fprintf( f_paramcomm, "Fraction of PHASE 1 maximum height of entrainment\t%.2f\n", rhentrmax1 );
    
    vhentrmax = fdparam ( fparam ); // variation of maximum height of entrainment
    fprintf( f_paramcomm, "Variation of maximum height of entrainment\t%.2f\n", vhentrmax );

    phiname = fcparam ( fparam ); // name of internal friction angle of PHASE 1 map
    if ( strcmp( phiname, "None" ) != 0 ) sico.PHI = 1;
    fprintf( f_paramcomm, "Internal friction angle of PHASE 1 map\t%s\n", phiname );

    phi2name = fcparam ( fparam ); // name of internal friction angle of PHASE 2 map
    if ( strcmp( phi2name, "None" ) != 0 ) sico.PHI2 = 1;
    fprintf( f_paramcomm, "Internal friction angle of PHASE 2 map\t%s\n", phi2name );

    phi3name = fcparam ( fparam ); // name of internal friction angle of PHASE 3 map
    if ( strcmp( phi3name, "None" ) != 0 ) sico.PHI3 = 1;
    fprintf( f_paramcomm, "Internal friction angle of PHASE 3 map\t%s\n", phi3name );

    deltaname = fcparam ( fparam ); // name of basal friction angle of PHASE 1 map
    if ( strcmp( deltaname, "None" ) != 0 ) sico.DELTA = 1;
    fprintf( f_paramcomm, "Basal friction angle of PHASE 1 map\t%s\n", deltaname );

    delta2name = fcparam ( fparam ); // name of basal friction angle of PHASE 2 map
    if ( strcmp( delta2name, "None" ) != 0 ) sico.DELTA2 = 1;
    fprintf( f_paramcomm, "Basal friction angle of PHASE 2 map\t%s\n", delta2name );

    delta3name = fcparam ( fparam ); // name of basal friction angle of PHASE 3 map
    if ( strcmp( delta3name, "None" ) != 0 ) sico.DELTA3 = 1;
    fprintf( f_paramcomm, "basal friction angle of PHASE 3 map\t%s\n", delta3name );

    addfriname = fcparam ( fparam ); // name of additional friction parameter of PHASE 1 map
    if ( strcmp( addfriname, "None" ) != 0 ) sico.ADDFRI = 1;
    fprintf( f_paramcomm, "Additional friction parameter of PHASE 1 map\t%s\n", addfriname );

    addfri2name = fcparam ( fparam ); // name of additional friction parameter of PHASE 2 map
    if ( strcmp( addfri2name, "None" ) != 0 ) sico.ADDFRI2 = 1;
    fprintf( f_paramcomm, "Additional friction parameter of PHASE 2 map\t%s\n", addfri2name );

    addfri3name = fcparam ( fparam ); // name of additional friction parameter of PHASE 3 map
    if ( strcmp( addfri3name, "None" ) != 0 ) sico.ADDFRI3 = 1;
    fprintf( f_paramcomm, "Additional friction parameter of PHASE 3 map\t%s\n", addfri3name );

    cohesionname = fcparam ( fparam ); // name of cohesion of PHASE 1 map
    if ( strcmp( cohesionname, "None" ) != 0 ) sico.COHESION = 1;
    fprintf( f_paramcomm, "Cohesion of PHASE 1 map\t%s\n", cohesionname );

    cohesion2name = fcparam ( fparam ); // name of cohesion of PHASE 2 map
    if ( strcmp( cohesion2name, "None" ) != 0 ) sico.COHESION2 = 1;
    fprintf( f_paramcomm, "Cohesion of PHASE 2 map\t%s\n", cohesion2name );

    cohesion3name = fcparam ( fparam ); // name of cohesion of PHASE 3 map
    if ( strcmp( cohesion3name, "None" ) != 0 ) sico.COHESION3 = 1;
    fprintf( f_paramcomm, "Cohesion of PHASE 3 map\t%s\n", cohesion3name );

    nyname = fcparam ( fparam ); // name of kinematic viscosity of PHASE 1 map
    if ( strcmp( nyname, "None" ) != 0 ) sico.NY = 1;
    fprintf( f_paramcomm, "Kinematic viscosity of PHASE 1 map\t%s\n", nyname );

    ny2name = fcparam ( fparam ); // name of kinematic viscosity of PHASE 2 map
    if ( strcmp( ny2name, "None" ) != 0 ) sico.NY2 = 1;
    fprintf( f_paramcomm, "Kinematic viscosity of PHASE 2 map\t%s\n", ny2name );

    ny3name = fcparam ( fparam ); // name of kinematic viscosity of PHASE 3 map
    if ( strcmp( ny3name, "None" ) != 0 ) sico.NY3 = 1;
    fprintf( f_paramcomm, "Kinematic viscosity of PHASE 3 map\t%s\n", ny3name );

    cdeformname = fcparam ( fparam ); // name of deformation coefficient coefficient of PHASE 1 map
    if ( strcmp( cdeformname, "None" ) != 0 ) sico.CDEFORM = 1;
    fprintf( f_paramcomm, "Deformation coefficient of PHASE 1 map\t%s\n", cdeformname );

    zfragname = fcparam ( fparam ); // name of fragmentation zones map
    if ( strcmp( zfragname, "None" ) != 0 ) sico.ZFRAG = 1;
    fprintf( f_paramcomm, "Fragmentation zones map\t%s\n", zfragname );

    ambdragname = fcparam ( fparam ); // name of ambient drag map
    if ( strcmp( ambdragname, "None" ) != 0 ) sico.AMBDRAG = 1;
    fprintf( f_paramcomm, "Ambient drag map\t%s\n", ambdragname );

    centrname = fcparam ( fparam ); // name of entrainment coefficient map
    if ( strcmp( centrname, "None" ) != 0 ) sico.CENTR = 1;
    fprintf( f_paramcomm, "Entrainment coefficient map\t%s\n", centrname );

    trans12name = fcparam ( fparam ); // name of PHASE 1 - PHASE 2 transformation coefficient map
    if ( strcmp( trans12name, "None" ) != 0 ) sico.TRANS12 = 1;
    fprintf( f_paramcomm, "PHASE 1 - PHASE 2 transformation coefficient map\t%s\n", trans12name );

    trans13name = fcparam ( fparam ); // name of PHASE 1 - PHASE 3 transformation coefficient map
    if ( strcmp( trans13name, "None" ) != 0 ) sico.TRANS13 = 1;
    fprintf( f_paramcomm, "PHASE 1 - PHASE 3 transformation coefficient map\t%s\n", trans13name );

    trans23name = fcparam ( fparam ); // name of PHASE 2 - PHASE 3 transformation coefficient map
    if ( strcmp( trans23name, "None" ) != 0 ) sico.TRANS23 = 1;
    fprintf( f_paramcomm, "PHASE 2 - PHASE 3 transformation coefficient map\t%s\n", trans23name );

    pbg1name = fcparam ( fparam ); // name of orthophoto channel 1 map
    fprintf( f_paramcomm, "Name of orthophoto channel 1 map\t%s\n", pbg1name );

    pbg2name = fcparam ( fparam ); // name of orthophoto channel 2 map
    fprintf( f_paramcomm, "Name of orthophoto channel 2 map\t%s\n", pbg2name );

    pbg3name = fcparam ( fparam ); // name of orthophoto channel 3 map
    fprintf( f_paramcomm, "Name of orthophoto channel 3 map\t%s\n", pbg3name );

    if ( strcmp( pbg1name, "None" ) != 0 ) sico.PBG = 1;
    if ( sico.PBG == 1 && strcmp( pbg2name, "None" ) == 0 ) strcpy( pbg2name, pbg1name );
    if ( sico.PBG == 1 && strcmp( pbg3name, "None" ) == 0 ) strcpy( pbg3name, pbg1name );


    // Hydrograph(s), adaptograph, frictiograph, and transformograph

    hydrocoordslist0 = fcparam ( fparam ); // parameters of hydrograph profiles
    fprintf( f_paramcomm, "Parameters of hydrograph profiles\t%s\n", hydrocoordslist0 );

    hydrographslist0 = fcparam ( fparam ); // pathes to input hydrograph files
    fprintf( f_paramcomm, "Pathes to input hydrograph files\t%s\n", hydrographslist0 );

    char *htest, *htok;
 
    if ( strcmp( hydrocoordslist0, "None" ) != 0 ) {

        hydrocoordslist = finlistd ( hydrocoordslist0 );
        hydrograph = 1;
        hydnum = (int)( hydrocoordslist[0] / 4 );
        
        hydi = (int*) calloc(( hydnum ), sizeof(int));
            // allocating memory to variable for internal coordinate of hydrograph point
        hydelev = (float*) calloc(( hydnum ), sizeof(float));
            // allocating memory to variable for elevation of hydrograph point
        hydalpha = (float*) calloc(( hydnum ), sizeof(float));
            // allocating memory to variable for aspect of hydrograph point
        hydx_metric = (float*) calloc(( hydnum ), sizeof(float));
            // allocating memory to variable for x coordinate of hydrograph
        hydy_metric = (float*) calloc(( hydnum ), sizeof(float));
            // allocating memory to variable for y coordinate of hydrograph
        hydl = (float*) calloc(( hydnum ), sizeof(float)); // allocating memory to variable for length of hydrograph profile
        hydx = (int*) calloc(( hydnum ), sizeof(int)); // allocating memory to variable for internal x coordinate of hydrograph
        hydy = (int*) calloc(( hydnum ), sizeof(int)); // allocating memory to variable for internal y coordinate of hydrograph

        if ( strcmp( hydrographslist0, "None" ) != 0 ) {

            htok=strtok(hydrographslist0, ",");
            hydrographslist[0] = (char*) calloc( 1000, sizeof(char));
            i=1;
    
            while ( htok != NULL ) {

                hydrographslist[i] = (char*) calloc( 1000, sizeof(char));
                htest=htok;
                htok=strtok(NULL, ",");
                strcpy( hydrographslist[i], htest ); // storing component of parameter string
     
                i += 1;
            }
            sprintf(hydrographslist[0], "%d", i-1 );
            hydnin = atoi( hydrographslist[0] );

        } else hydnin = 0;
        
        hydnout = hydnum - hydnin;
    }

    if ( hydrograph == 1 ) {
    
        hydtmax = (int*) calloc( hydnin, sizeof(int)); // allocating memory to variable for number(s) of hydrograph time steps
        hydtmaxx = 0; // initializing maximum number of hydrograph time steps

        for ( hydj = 0; hydj < hydnin; hydj++ ) { // loop over all input hydrographs

            hydreleasename[hydj] = (char*) calloc( 1000, sizeof(char)); // allocating memory to pointers to input hydrograph files
            strcpy( hydreleasename[hydj], hydrographslist[hydj+1] ); // path to hydrograph file
            
            hydtmax[hydj] = fcountlines( indir, hydreleasename[hydj] ) - 2; // number of hydrograph time steps
            if ( hydtmax[hydj] > hydtmaxx ) hydtmaxx = hydtmax[hydj]; // updating maximum number of hydrograph time steps
        }

        for ( hydj = 0; hydj < hydnum; hydj++ ) { // loop over all hydrographs

            hydx_metric[hydj] = hydrocoordslist[4*hydj+1]; // x coordinate of hydrograph
            hydy_metric[hydj] = hydrocoordslist[4*hydj+2]; // y coordinate of hydrograph
            hydl[hydj] = hydrocoordslist[4*hydj+3]; // length of hydrograph profile
            hydalpha[hydj] = hydrocoordslist[4*hydj+4] * sico.PI / 180; // aspect of hydrograph profile
        }

        hydhyd = alloc_dmatrix3( hydtmaxx+1, 7, hydnin ); // allocating memory to array of input hydrograph data
    
    } else { hydnin = 0; hydnout = 0; }

    adaptoname = fcparam ( fparam ); //path to adaptograph file

    if ( strcmp( adaptoname, "None" ) != 0 ) {
        adaptograph = 1;
        adatmax = fcountlines( indir, adaptoname ) - 2; // number of adaptograph time steps
        adaada = alloc_dmatrix( adatmax+1, 8 ); // allocating memory to array of input adaptograph data
    }

    fprintf( f_paramcomm, "Path to adaptograph file\t%s\n", adaptoname ); 

    frictioname = fcparam ( fparam ); // name of frictiograph file
    if ( strcmp( frictioname, "None" ) != 0 ) {
        frictiograph = 1;
        fritmax = fcountlines( indir, frictioname ) - 2; // number of frictiograph time steps
        frifri = alloc_dmatrix( fritmax+1, 13 ); // allocating memory to array of input frictiograph data
    }

    fprintf( f_paramcomm, "Path to frictiograph file\t%s\n", frictioname ); 

    transformoname = fcparam ( fparam ); // name of transformograph file
    if ( strcmp( transformoname, "None" ) != 0 ) {
        transformograph = 1;
        tratmax = fcountlines( indir, transformoname ) - 2; // number of transformograph time steps
        tratra = alloc_dmatrix( tratmax+1, 7 ); // allocating memory to array of input transformograph data
    }

    fprintf( f_paramcomm, "Path to transformograph file\t%s\n", transformoname ); 


    // Flow parameters

    lmax = fiparam ( fparam ); // number of flow parameters
    fprintf( f_paramcomm, "Number of flow parameters\t%i\n", lmax );
    
    flowpar = (float*) calloc( lmax, sizeof(float)); // allocating memory for storage of flow parameter values

    for ( l=0; l<lmax; l++ ) flowpar[l] = fdparam ( fparam );

    if ( sico.MODEL == 1 ) { // single-phase model

        sflow.RHO1 = flowpar[0]; // density
        fprintf( f_paramcomm, "Density (kg/m³)\t%.2f\n", sflow.RHO1 );

        sflow.PHI0[0] = flowpar[1] * sico.PI / 180; // internal friction angle
        sflow.DELTA0[0] = flowpar[2] * sico.PI / 180; // basal friction angle
        sflow.ADDFRI0[0] = flowpar[3]; // additional friction parameter
        fprintf( f_paramcomm, "Internal friction angle (degrees)\t%.2f\n", sflow.PHI0[0] * 180 / sico.PI );
        fprintf( f_paramcomm, "Basal friction angle (degrees)\t%.2f\n", sflow.DELTA0[0] * 180 / sico.PI );
        fprintf( f_paramcomm, "Additional friction parameter\t%.4f\n", sflow.ADDFRI0[0] );

        sflow.COHESION0[0] = flowpar[4]; // cohesion
        fprintf( f_paramcomm, "Cohesion (N/m²)\t%.2f\n", sflow.COHESION0[0] );

        sflow.NY0[0] = pow(10, flowpar[5]); // viscosity
        fprintf( f_paramcomm, "Log10 of kinematic viscosity (m²/s)\t%.2f\n", log10( sflow.NY0[0] ));

        sflow.CDEFORM0 = flowpar[6]; // deformation coefficient
        fprintf( f_paramcomm, "Deformation coefficient\t%.2f\n", sflow.CDEFORM0 );

        sflow.SLIDEEXP[0] = flowpar[7]; // exponent for evolution of global gravity component
        sflow.SLIDEDEF[0] = flowpar[8]; // fraction of global gravity
        fprintf( f_paramcomm, "Exponent for evolution of sliding component\t%.2f\n", sflow.SLIDEEXP[0] );
        fprintf( f_paramcomm, "Initial fraction of sliding\t%.2f\n", sflow.SLIDEDEF[0] );

        sflow.ESP = flowpar[9]; // strength of frictional energy loss in shearing
        fprintf( f_paramcomm, "Strength of frictional energy loss in shearing\t%.2f\n", sflow.ESP );

        sflow.FRAGSML = flowpar[10]; // fragmentation parameter SMl
        sflow.FRAGML = flowpar[11]; // fragmentation parameter Ml
        fprintf( f_paramcomm, "Fragmentation parameter SMl\t%.4f\n", sflow.FRAGSML ); 
        fprintf( f_paramcomm, "Fragmentation parameter Ml\t%.4f\n", sflow.FRAGML );

        sflow.AMBDRAG0 = flowpar[12]; // ambient drag coefficient
        fprintf( f_paramcomm, "Ambient drag coefficient\t%.4f\n", sflow.AMBDRAG0 ); 

        sflow.CENTR0 = flowpar[13]; // entrainment coefficient
        sflow.CSTOP = flowpar[14]; // stopping criterion
        fprintf( f_paramcomm, "Log 10 of entrainment coefficient\t%.2f\n", sflow.CENTR0 ); 
        fprintf( f_paramcomm, "Stopping criterion\t%.2f\n", sflow.CSTOP ); 
 
    } else if ( sico.MODEL == 7 ) { // multi-phase model

        sflow.RHO1 = flowpar[0]; // PHASE 1 density
        sflow.RHO2 = flowpar[1]; // PHASE 2 density
        sflow.RHO3 = flowpar[2]; // PHASE 3 density
        fprintf( f_paramcomm, "Density of PHASE 1 (kg/m³)\t%.2f\n", sflow.RHO1 );
        fprintf( f_paramcomm, "Density of PHASE 2 (kg/m³)\t%.2f\n", sflow.RHO2 );
        fprintf( f_paramcomm, "Density of PHASE 3 (kg/m³)\t%.2f\n", sflow.RHO3 );

        sflow.PHI0[0] = flowpar[3] * sico.PI / 180; // internal friction angle of PHASE 1
        sflow.PHI0[1] = flowpar[4] * sico.PI / 180; // internal friction angle of PHASE 2
        sflow.PHI0[2] = flowpar[5] * sico.PI / 180; // internal friction angle of PHASE 3
        sflow.DELTA0[0] = flowpar[6] * sico.PI / 180; // basal friction angle of PHASE 1
        sflow.DELTA0[1] = flowpar[7] * sico.PI / 180; // basal friction angle of PHASE 2
        sflow.DELTA0[2] = flowpar[8] * sico.PI / 180; // basal friction angle of PHASE 3
        sflow.ADDFRI0[0] = flowpar[9]; // additional friction parameter of PHASE 1
        sflow.ADDFRI0[1] = flowpar[10]; // additional friction parameter of PHASE 2
        sflow.ADDFRI0[2] = flowpar[11]; // additional friction parameter of PHASE 3
        fprintf( f_paramcomm, "Internal friction angle of PHASE 1 (degrees)\t%.2f\n", sflow.PHI0[0] * 180 / sico.PI );
        fprintf( f_paramcomm, "Internal friction angle of PHASE 2 (degrees)\t%.2f\n", sflow.PHI0[1] * 180 / sico.PI );
        fprintf( f_paramcomm, "Internal friction angle of PHASE 3 (degrees)\t%.2f\n", sflow.PHI0[2] * 180 / sico.PI );
        fprintf( f_paramcomm, "Basal friction angle of PHASE 1 (degrees)\t%.2f\n", sflow.DELTA0[0] * 180 / sico.PI );
        fprintf( f_paramcomm, "Basal friction angle of PHASE 2 (degrees)\t%.2f\n", sflow.DELTA0[1] * 180 / sico.PI );
        fprintf( f_paramcomm, "Basal friction angle of PHASE 3 (degrees)\t%.2f\n", sflow.DELTA0[2] * 180 / sico.PI );
        fprintf( f_paramcomm, "Additional friction parameter of PHASE 1\t%.4f\n", sflow.ADDFRI0[0] );
        fprintf( f_paramcomm, "Additional friction parameter of PHASE 2\t%.4f\n", sflow.ADDFRI0[1] );
        fprintf( f_paramcomm, "Additional friction parameter of PHASE 3\t%.4f\n", sflow.ADDFRI0[2] );

        sflow.COHESION0[0] = flowpar[12]; // cohesion
        sflow.COHESION0[1] = flowpar[13]; // cohesion
        sflow.COHESION0[2] = flowpar[14]; // cohesion
        fprintf( f_paramcomm, "Cohesion of PHASE 1\t%.2f\n", sflow.COHESION0[0] );
        fprintf( f_paramcomm, "Cohesion of PHASE 2\t%.2f\n", sflow.COHESION0[1] );
        fprintf( f_paramcomm, "Cohesion of PHASE 3\t%.2f\n", sflow.COHESION0[2] );

        sflow.NY0[0] = pow(10, flowpar[15]); // viscosity of PHASE 1
        sflow.NY0[1] = pow(10, flowpar[16]); // viscosity of PHASE 2
        sflow.NY0[2] = pow(10, flowpar[17]); // viscosity of PHASE 3
        fprintf( f_paramcomm, "Log10 of kinematic viscosity of PHASE 1 (m²/s)\t%.2f\n", log10( sflow.NY0[0] ));
        fprintf( f_paramcomm, "Log10 of kinematic viscosity of PHASE 2 (m²/s)\t%.2f\n", log10( sflow.NY0[1] ));
        fprintf( f_paramcomm, "Log10 of kinematic viscosity of PHASE 3 (m²/s)\t%.2f\n", log10( sflow.NY0[2] ));

        sflow.CDEFORM0 = flowpar[18]; // deformation coefficient
        sflow.CDEFORM[1] = flowpar[19]; // deformation coefficient
        sflow.CDEFORM[2] = flowpar[20]; // deformation coefficient
        fprintf( f_paramcomm, "Deformation coefficient of PHASE 1\t%.2f\n", sflow.CDEFORM0 );
        fprintf( f_paramcomm, "Deformation coefficient of PHASE 2\t%.2f\n", sflow.CDEFORM[1] );
        fprintf( f_paramcomm, "Deformation coefficient of PHASE 3\t%.2f\n", sflow.CDEFORM[2] );

        sflow.SLIDEEXP[0] = flowpar[21]; // exponent for evolution of global gravity component
        sflow.SLIDEEXP[1] = flowpar[22]; // exponent for evolution of global gravity component
        sflow.SLIDEEXP[2] = flowpar[23]; // exponent for evolution of global gravity component
        sflow.SLIDEDEF[0] = flowpar[24]; // fraction of global gravity
        sflow.SLIDEDEF[1] = flowpar[25]; // fraction of global gravity
        sflow.SLIDEDEF[2] = flowpar[26]; // fraction of global gravity
        
        fprintf( f_paramcomm, "Exponent for evolution of sliding component of PHASE 1\t%.2f\n", sflow.SLIDEEXP[0] );
        fprintf( f_paramcomm, "Exponent for evolution of sliding component of PHASE 2\t%.2f\n", sflow.SLIDEEXP[1] );
        fprintf( f_paramcomm, "Exponent for evolution of sliding component of PHASE 3\t%.2f\n", sflow.SLIDEEXP[2] );        
        fprintf( f_paramcomm, "Initial fraction of sliding of PHASE 1\t%.2f\n", sflow.SLIDEDEF[0] );
        fprintf( f_paramcomm, "Initial fraction of sliding of PHASE 2\t%.2f\n", sflow.SLIDEDEF[1] );
        fprintf( f_paramcomm, "Initial fraction of sliding of PHASE 3\t%.2f\n", sflow.SLIDEDEF[2] );

        sflow.ESP = flowpar[27]; // strength of frictional energy loss in shearing
        fprintf( f_paramcomm, "Strength of frictional energy loss in shearing\t%.2f\n", sflow.ESP );
        
        sflow.FRAGSML = flowpar[28]; // fragmentation parameter SMl
        sflow.FRAGML = flowpar[29]; // fragmentation parameter Ml
        fprintf( f_paramcomm, "Fragmentation parameter SMl\t%.4f\n", sflow.FRAGSML ); 
        fprintf( f_paramcomm, "Fragmentation parameter Ml\t%.4f\n", sflow.FRAGML );

        sflow.AMBDRAG0 = flowpar[30]; // ambient drag coefficient
        fprintf( f_paramcomm, "Ambient drag coefficient\t%.4f\n", sflow.AMBDRAG0 );
       
        sflow.CENTR0 = flowpar[31]; // entrainment coefficient
        sflow.CSTOP = flowpar[32]; // stopping criterion
        fprintf( f_paramcomm, "Log 10 of entrainment coefficient\t%.2f\n", sflow.CENTR0 ); 
        fprintf( f_paramcomm, "Stopping criterion\t%.2f\n", sflow.CSTOP ); 

        sflow.DRAG_k = flowpar[33]; // mass flux parameter for drag (m/s)
        sflow.DRAG_m = flowpar[34]; // exponent for scaling of the fluid-like drag contributions to flow resistance
        sflow.DRAG_n = flowpar[35]; // exponent for scaling of drag parameter with solid fraction
        sflow.DRAG_vterm = flowpar[36]; // terminal velocity (m/s)
        sflow.DRAG_rep = flowpar[37]; // particle Reynolds number
        sflow.DRAG_je = flowpar[38]; // exponent for drag
        fprintf( f_paramcomm, "Mass flux parameter for drag (m/s)\t%2f\n", sflow.DRAG_k );
        fprintf( f_paramcomm, "Exponent for scaling of the fluid-like drag contributions to flow resistance\t%2f\n", sflow.DRAG_m );
        fprintf( f_paramcomm, "Exponent for scaling of drag parameter with solid fraction\t%i\n", sflow.DRAG_n );
        fprintf( f_paramcomm, "Terminal velocity (m/s)\t%2f\n", sflow.DRAG_vterm );
        fprintf( f_paramcomm, "Particle Reynolds number\t%2f\n", sflow.DRAG_rep );
        fprintf( f_paramcomm, "Exponent for drag\t%i\n", sflow.DRAG_je );

        sflow.VM_n0 = flowpar[39]; // virtual mass number
        sflow.VM_l = flowpar[40]; // l parameter related to the virtual mass coefficients
        sflow.VM_n = flowpar[41]; // n parameter related to the virtual mass coefficients
        fprintf( f_paramcomm, "Virtual mass number\t%2f\n", sflow.VM_n0 );
        fprintf( f_paramcomm, "l parameter related to the virtual mass coefficients\t%2f\n", sflow.VM_l );
        fprintf( f_paramcomm, "n parameter related to the virtual mass coefficients\t%2f\n", sflow.VM_n );

        sflow.TRANS120 = flowpar[42]; // transformation coefficient PHASE 1 - PHASE 2
        sflow.TRANS130 = flowpar[43]; // transformation coefficient PHASE 1 - PHASE 3
        sflow.TRANS230 = flowpar[44]; // transformation coefficient PHASE 2 - PHASE 3
        fprintf( f_paramcomm, "Transformation coefficient PHASE 1 - PHASE 2\t%.4f\n", sflow.TRANS120 ); 
        fprintf( f_paramcomm, "Transformation coefficient PHASE 1 - PHASE 3\t%.4f\n", sflow.TRANS130 ); 
        fprintf( f_paramcomm, "Transformation coefficient PHASE 2 - PHASE 3\t%.4f\n", sflow.TRANS230 );
        
        sflow.TEMP0 = flowpar[45]; // temperature of release mass
        sflow.TEMPS = flowpar[46]; // atmospheric temperature
        sflow.TEMPB = flowpar[47]; // ground temperature
        sflow.L = flowpar[48]; // fraction of the shear heat transferred into the avalanche and causing ice melt
        sflow.LAMBDA = flowpar[49]; // fraction of basal sliding
        
        fprintf( f_paramcomm, "Temperature of release mass\t%.1f\n", sflow.TEMP0 ); 
        fprintf( f_paramcomm, "Atmospheric temperature\t%.1f\n", sflow.TEMPS ); 
        fprintf( f_paramcomm, "Ground temperature\t%.1f\n", sflow.TEMPB );
        fprintf( f_paramcomm, "Fraction of the shear heat transferred into the avalanche and causing ice melt\t%.1f\n", sflow.L ); 
        fprintf( f_paramcomm, "Fraction of basal sliding\t%.1f\n", sflow.LAMBDA );
    }


    // Hard-coded flow parameters

    sflow.disp_MULT = 1.0; // multiplication factor for dispersion   
    sflow.CCONST = 5.0; // coefficient for constraining dispersion and phase separation
    sflow.BETAPLAIN = 1.0 * sico.PI / 180; // maximum slope angle considered as plane surface
    sflow.FRICMIN = 0.0 * sico.PI / 180; // minimum value of internal and basal friction
    sflow.EKINCOEF = pow( 10, -6.0 ); // kinetic energy coefficient
    sflow.FRI_exp = 0.0; // phase fraction scaling exponent
    sflow.VIS_chi[0] = 1.0; // vertical PHASE 1 velocity distribution (0=no shearing, 3=parabolic profile)
    sflow.KPMAX[0] = 1.0; // maximum value of PHASE 1 passive earth pressure coefficient
    sflow.HREF = 1.0; // reference flow height for slow flows (with decomposition of viscosity terms)

    if ( sico.MODEL == 7 ) { // multi-phase model

        sflow.THETAS = 0.0; // water content of the entrainable material 
        sflow.VM_fact = 1.0; // reduction factor for virtual mass coefficients        
        sflow.sep_SHEAR = 1.0; // shear factor for phase separation
        sflow.VIS_chi[1] = 1.0; // vertical PHASE 2 velocity distribution (0=no shearing, 3=parabolic profile)
        sflow.VIS_chi[2] = 1.0; // vertical PHASE 3 velocity distribution (0=no shearing, 3=parabolic profile)
        sflow.KPMAX[1] = 1.0; // maximum value of PHASE 2 passive earth pressure coefficient
        sflow.KPMAX[2] = 1.0; // maximum value of PHASE 3 passive earth pressure coefficient
        
        if ( sico.TEMP != 0 ) { // temperature effects
        
            sflow.KAPPAR = 3.36; // thermal conductivities
            sflow.KAPPAI = 2.25;
            sflow.KAPPAF = 0.55;
            sflow.CR = 1000; // heat capacities
            sflow.CI = 2000;
            sflow.CF = 4000;
            sflow.Q = 333.5; // latent heat for ice
            sflow.LAMBDAS = 75; // heat exchange coefficient at the free surface
            sflow.LAMBDAB = 5; // temperature profile dependent parameter at the bottom
            sflow.AS = 0.1; // area exposed to the heat exchange at the free surface
            sflow.CL = 1.0; // ice melt coefficient
        }
    }
  

    // Thresholds and further parameters

    sico.IMPTHR[0] = fdparam ( fparam ); // threshold for display of flow height (m)
    fprintf( f_paramcomm, "Threshold for display of flow height (m)\t%.3f\n", sico.IMPTHR[0] );

    sico.IMPTHR[1] = fdparam ( fparam ); // threshold for display of flow kinetic energy (J)
    fprintf( f_paramcomm, "Threshold for display of flow kinetic energy (J)\t%.3f\n", sico.IMPTHR[1] );

    sico.IMPTHR[2] = fdparam ( fparam ); // threshold for display of flow pressure (Pa)
    fprintf( f_paramcomm, "Threshold for display of flow pressure (Pa)\t%.3f\n", sico.IMPTHR[2] );

    sico.IMPTHR[3] = fdparam ( fparam ); // threshold of flow height for display of velocity (m)
    fprintf( f_paramcomm, "Threshold of flow height for display of velocity (m)\t%.3f\n", sico.IMPTHR[3] );

    sico.HFLOWMIN = fdparam ( fparam ); // threshold flow height for simulation (m)
    fprintf( f_paramcomm, "Threshold flow height for simulation (m)\t%.6f\n", sico.HFLOWMIN );

    tout = fdparam ( fparam ); // time interval for writing output (s)
    fprintf( f_paramcomm, "Time interval for writing output (s)\t%.2f\n", tout );

    tmax = fdparam ( fparam ); // time at which to stop simulation (s)
    fprintf( f_paramcomm, "Time at which to stop simulation (s)\t%.2f\n", tmax );

    trat = round( tmax / tout );
    tmax = trat * tout; // ensuring that time to stop is multiple of output time step length

    sico.SLOMO = fdparam ( fparam ); // factor for slow motion
    fprintf( f_paramcomm, "Factor for slow motion\t%.2f\n", sico.SLOMO );

    if ( sico.SLOMO < 0 ) { // defining glacier mode
    
        sico.GLACIER = 1;
        sico.SLOMO = fabs( sico.SLOMO );
        
    } else { sico.GLACIER = 0; }

    sflow.LAMBDA_NU = fdparam ( fparam ); // viscosity control
    sflow.LAMBDA_FL = fdparam ( fparam ); // flux control
    fprintf( f_paramcomm, "Viscosity control\t%.5f\n", sflow.LAMBDA_NU );
    fprintf( f_paramcomm, "Flux control\t%.5f\n", sflow.LAMBDA_FL );
   
    sico.CFL[0] = fdparam ( fparam ); // CFL criterion
    fprintf( f_paramcomm, "CFL criterion\t%.2f\n", sico.CFL[0] );

    sico.CFL[1] = fdparam ( fparam ); // length of first time step (s)
    fprintf( f_paramcomm, "Length of first time step (s)\t%.4f\n", sico.CFL[1] );

    proflist0 = fcparam ( fparam ); // profile vertices (x and y coordinates)
    fprintf( f_paramcomm, "Profile vertices (x and y coordinates)\t%s\n", proflist0 );
  
    if ( strcmp( proflist0, "None" ) != 0 ) {

        proflist = finlistd ( proflist0 );
        sico.PROFILE = (int)( proflist[0] / 2 );

        if ( sico.MODEL <= 3 ) profmax = 5;
        else if ( sico.MODEL == 7 ) profmax = 15;

        profelev = (float*) calloc( sico.M+sico.N, sizeof(float)); 
        profdeposit = (float*) calloc( sico.M+sico.N, sizeof(float));    
        profdata = alloc_dmatrix3( sico.M+sico.N, (int)(tmax/tout+3), profmax );
        profeval = alloc_dmatrix( sico.M+sico.N, 3 );
        profx = (int*) calloc( sico.PROFILE, sizeof(int));
        profy = (int*) calloc( sico.PROFILE, sizeof(int));
        profnx = (float*) calloc( sico.M+sico.N, sizeof(float));
        profny = (float*) calloc( sico.M+sico.N, sizeof(float));
        profdiffabs = (float*) calloc( sico.M+sico.N, sizeof(float));
         
        sprintf(path, "%s%smapprof.txt", outfiles, prefix); // profile coordinates file
        f_mapprof=fopen(path, "w");
    
        for ( profi=0; profi<sico.PROFILE; profi++ ) {

             profxm = proflist[2*profi+1];
             profy[profi] = (int)(( profxm - sico.BDWEST ) / sico.CSZ + 0.5 );
             if ( profi == 0 ) fprintf( f_mapprof, "%.6f", profxm );
             else fprintf( f_mapprof, "\t%.6f", profxm );
        }

        for ( profi=0; profi<sico.PROFILE; profi++ ) {

             profym = proflist[2*profi+2];
             profx[profi] = (int)(( sico.BDNORTH - profym ) / sico.CSZ + 0.5 );
             if ( profi == 0 ) fprintf( f_mapprof, "\n%.6f", profym );
             else fprintf( f_mapprof, "\t%.6f", profym );
         }
         
         fclose(f_mapprof);
    }

    ctrlplist0 = fcparam ( fparam ); // control points (x and y coordinates)
    fprintf( f_paramcomm, "Control points (x and y coordinates)\t%s\n", ctrlplist0 );
  
    if ( strcmp( ctrlplist0, "None" ) != 0 ) {

        ctrlplist = finlistd ( ctrlplist0 );
        sico.CTRLPOINTS = (int)( ctrlplist[0] / 2 );

        if ( sico.MODEL <= 3 ) k = 2; else if ( sico.MODEL == 7 ) k = 6;

        ctrlpx = (int*) calloc( sico.CTRLPOINTS, sizeof(int)); 
        ctrlpy = (int*) calloc( sico.CTRLPOINTS, sizeof(int));
        
        ctrlpxm = (float*) calloc( sico.CTRLPOINTS, sizeof(float)); 
        ctrlpym = (float*) calloc( sico.CTRLPOINTS, sizeof(float));
        
        ctrlpdata = alloc_dmatrix3( sico.CTRLPOINTS, (int)(tmax/tout+3), k );
 
        for ( i=0; i<sico.CTRLPOINTS; i++ ) {

             ctrlpxm[i] = ctrlplist[2*i+1];
             ctrlpym[i] = ctrlplist[2*i+2];

             ctrlpy[i] = (int)(( ctrlpxm[i] - sico.BDWEST ) / sico.CSZ + 0.5 );
             ctrlpx[i] = (int)(( sico.BDNORTH - ctrlpym[i] ) / sico.CSZ + 0.5 );
        }
    }


    // Parameters for virtual reality and R interfaces

    paradistort = fiparam ( fparam ); // control for distorting background orthophoto (0 = off, 1 = on)
    fprintf( f_paramcomm, "Control for distorting background orthophoto (0 = off, 1 = on)\t%i\n", paradistort ); 

    paramin = fdparam ( fparam ); // minimum flow height for virtual reality (m)
    fprintf( f_paramcomm, "Minimum flow height for virtual reality (m)\t%.2f\n", paramin ); 

    pararef = fdparam ( fparam ); // reference flow height for virtual reality (m)
    fprintf( f_paramcomm, "Reference flow height for virtual reality (m)\t%.2f\n", pararef ); 
        
    paratsunref = fdparam ( fparam ); // reference tsunami height for virtual reality (m)
    fprintf( f_paramcomm, "Reference tsunami height for virtual reality (m)\t%.2f\n", paratsunref );   

    paracontourshmin = fiparam ( fparam ); // minimum level for flow height contours in Paraview (m)
    fprintf( f_paramcomm, "Minimum level for flow height contours in Paraview (m)\t%i\n", paracontourshmin ); 

    paracontourshmax = fiparam ( fparam ); // maximum level for flow height contours in Paraview (m)
    fprintf( f_paramcomm, "Maximum level for flow height contours in Paraview (m)\t%i\n", paracontourshmax );

    paracontourshint = fiparam ( fparam ); // interval for flow height contours in Paraview (m)
    fprintf( f_paramcomm, "Interval for flow height contours in Paraview (m)\t%i\n", paracontourshint );

    paracontourszmin = fiparam ( fparam ); // minimum level for elevation contours in Paraview (m)
    fprintf( f_paramcomm, "Minimum level for elevation contours in Paraview (m)\t%i\n", paracontourszmin ); 

    paracontourszmax = fiparam ( fparam ); // maximum level for elevation contours in Paraview (m)
    fprintf( f_paramcomm, "Maximum level for elevation contours in Paraview (m)\t%i\n", paracontourszmax );

    paracontourszint = fiparam ( fparam ); // interval for elevation contours in Paraview (m)
    fprintf( f_paramcomm, "Interval for elevation contours in Paraview (m)\t%i\n", paracontourszint );

    parar = fdparam ( fparam ); // weight of red colour for virtual reality (neglected for multi-phase model)
    fprintf( f_paramcomm, "Weight of red colour for virtual reality (neglected for multi-phase model)\t%.2f\n", parar );

    parag = fdparam ( fparam ); // weight of green colour for virtual reality (neglected for multi-phase model)
    fprintf( f_paramcomm, "Weight of green colour for virtual reality (neglected for multi-phase model)\t%.2f\n", parag );
        
    parab = fdparam ( fparam ); // weight of blue colour for virtual reality (neglected for multi-phase model)
    fprintf( f_paramcomm, "Weight of blue colour for virtual reality (neglected for multi-phase model)\t%.2f\n", parab );               

    parad = fdparam ( fparam ); // exponent of transparency curve for virtual reality
    fprintf( f_paramcomm, "Exponent of transparency curve for virtual reality\t%.2f\n", parad );

    phexagg = fdparam ( fparam ); // exaggeration factor for flow heights in profile plots
    fprintf( f_paramcomm, "Exaggeration factor for flow heights in profile plots\t%.2f\n", phexagg );
       
    parapy = fcparam ( fparam ); // path to pvpython for visualization in Paraview
    fprintf( f_paramcomm, "Path to pvpython for visualization in Paraview\t%s\n", parapy );        

    rscript = fcparam ( fparam ); // path to Rscript for visualization in R
    fprintf( f_paramcomm, "Path to Rscript for visualization in R\t%s\n", rscript ); 

    rlibs = fcparam ( fparam ); // path to R packages for map plots
    fprintf( f_paramcomm, "Path to R packages for map plots\t%s\n", rlibs ); 

    fclose(fparam); // closing parameter file
    fclose(f_paramcomm); // closing commented parameter file
   

// -- STOP --- Reading and preprocessing parameter file input ---------------------------------------------------


// -- START -- Reading hydrograph, frictiograph, and transformograph data ---------------------------------------


    if ( hydrograph == 1 ) {

        if ( sico.MODEL <= 3 ) hydcols = 3; else hydcols = 7;

        hydlmax = 0;
        for ( hydj = 0; hydj < hydnin + hydnout; hydj++ ) if ( fabs(hydl[hydj]) > hydlmax ) hydlmax = fabs(hydl[hydj]);

        for ( hydj = 0; hydj < hydnin; hydj++ ) {
            hydhyd0 = finhyd( indir, hydreleasename[hydj], hydtmax[hydj], hydtmaxx ); // reading data from all input hydrographs
            for ( i=0; i<hydtmaxx+1; i++ ) { for ( j=0; j<hydcols; j++ ) hydhyd[i][j][hydj] = hydhyd0[i][j]; }
            free( hydhyd0[0] ); free( hydhyd0 );
        }

        for ( hydj = 0; hydj < hydnin + hydnout; hydj++ ) { // internal x and y coordinates of all hydrographs

            hydy[hydj] = (int)(( hydx_metric[hydj] - sico.BDWEST ) / sico.CSZ + 0.5 );
            hydx[hydj] = (int)(( sico.BDNORTH - hydy_metric[hydj] ) / sico.CSZ + 0.5 );
        }
    }

    if ( adaptograph == 1 ) {

        adaada = finhyd( indir, adaptoname, adatmax, adatmax ); // reading data from input adaptograph
    }

    if ( frictiograph == 1 ) {

        frifri = finhyd( indir, frictioname, fritmax, fritmax ); // reading data from input frictiograph
    }

    if ( transformograph == 1 ) {

        tratra = finhyd( indir, transformoname, tratmax, tratmax ); // reading data from input transformograph
    }


// -- STOP --- Reading hydrograph, frictiograph, and transformograph data ---------------------------------------


// -- START -- Reading input raster maps (GRASS) ----------------------------------------------------------------


    #ifdef WITHGRASS


        seg_elev = finrastd( elevname, sico );
        
        if ( sico.RELM == 1 ) seg_hrelease = finrastd( hreleasename, sico );
        if ( sico.MODEL == 7 && sico.RELM2 == 1 ) seg_hrelease2 = finrastd( hreleasename2, sico );
        if ( sico.MODEL == 7 && sico.RELM3 == 1 ) seg_hrelease3 = finrastd( hreleasename3, sico );
        
        if ( sico.RELV == 1 ) seg_vinx = finrastd( vinxname, sico );
        if ( sico.MODEL == 7 && sico.RELV2 == 1 ) seg_vinx2 = finrastd( vinxname2, sico );
        if ( sico.MODEL == 7 && sico.RELV3 == 1 ) seg_vinx3 = finrastd( vinxname3, sico );
        
        if ( sico.RELV == 1 ) seg_viny = finrastd( vinyname, sico );
        if ( sico.MODEL == 7 && sico.RELV2 == 1 ) seg_viny2 = finrastd( vinyname2, sico );
        if ( sico.MODEL == 7 && sico.RELV3 == 1 ) seg_viny3 = finrastd( vinyname3, sico );
        
        if ( sico.ENTR == 1 ) seg_hentrmax = finrastd( hentrmaxname, sico );
        if ( sico.MODEL == 7 && sico.ENTR2 == 1 ) seg_hentrmax2 = finrastd( hentrmaxname2, sico );
        if ( sico.MODEL == 7 && sico.ENTR3 == 1 ) seg_hentrmax3 = finrastd( hentrmaxname3, sico );

        if ( sico.ZONES == 1 ) seg_zones = finrasti( zonesname, sico );
        if ( sico.CENTR == 1 ) seg_centr = finrastd( centrname, sico );
        if ( sico.CDEFORM == 1 ) seg_cdeform = finrastd( cdeformname, sico );
        if ( sico.PHI == 1 ) seg_phi = finrastd( phiname, sico );
        if ( sico.PHI2 == 1 ) seg_phi2 = finrastd( phi2name, sico ); 
        if ( sico.PHI3 == 1 ) seg_phi3 = finrastd( phi3name, sico );
        if ( sico.DELTA == 1 ) seg_delta = finrastd( deltaname, sico );
        if ( sico.DELTA2 == 1 ) seg_delta2 = finrastd( delta2name, sico );
        if ( sico.DELTA3 == 1 ) seg_delta3 = finrastd( delta3name, sico );
        if ( sico.ADDFRI == 1 ) seg_addfri = finrastd( addfriname, sico );
        if ( sico.ADDFRI2 == 1 ) seg_addfri2 = finrastd( addfri2name, sico );
        if ( sico.ADDFRI3 == 1 ) seg_addfri3 = finrastd( addfri3name, sico );
        if ( sico.COHESION == 1 ) seg_cohesion = finrastd( cohesionname, sico );
        if ( sico.COHESION2 == 1 ) seg_cohesion2 = finrastd( cohesion2name, sico );
        if ( sico.COHESION3 == 1 ) seg_cohesion3 = finrastd( cohesion3name, sico );
        if ( sico.NY == 1 ) seg_ny = finrastd( nyname, sico );
        if ( sico.NY2 == 1 ) seg_ny2 = finrastd( ny2name, sico );
        if ( sico.NY3 == 1 ) seg_ny3 = finrastd( ny3name, sico );
        if ( sico.AMBDRAG == 1 ) seg_ambdrag = finrastd( ambdragname, sico );
        if ( sico.TRANS12 == 1 ) seg_trans12 = finrastd( trans12name, sico );
        if ( sico.TRANS13 == 1 ) seg_trans13 = finrastd( trans13name, sico );
        if ( sico.TRANS23 == 1 ) seg_trans23 = finrastd( trans23name, sico );
        if ( sico.ZFRAG == 1 ) seg_zfrag = finrastd( zfragname, sico );
        if ( sico.TRELEASE == 1 ) seg_trelease = finrastd( treleasename, sico );
        if ( sico.TRELSTOP == 1 ) seg_trelstop = finrastd( trelstopname, sico );
        if ( sico.STOPTIME == 1 ) seg_stoptime = finrastd( stoptimename, sico );
        if ( sico.TSLIDE == 1 ) seg_tslide = finrastd( tslidename, sico );
        if ( sico.IMPACTAREA == 1 ) seg_impactarea = finrasti( impactareaname, sico );
        if ( sico.HDEPOSIT == 1 ) seg_hdeposit = finrastd( hdepositname, sico );
        if ( sico.PBG == 1 ) {
            seg_pbg1 = finrasti( pbg1name, sico );
            seg_pbg2 = finrasti( pbg2name, sico );
            seg_pbg3 = finrasti( pbg3name, sico );
        }

        sico.IMAX = 0;
        for ( x = 0; x < sico.M; x++ ) {
            for ( y = 0; y < sico.N; y++ ) sico.IMAX += 1; // number of raster cells
        }


    #endif


// -- STOP --- Reading input raster maps (GRASS) ----------------------------------------------------------------


// -- START -- Preparing arrays for input -----------------------------------------------------------------------


    int *pzones, *pimpactarea = 0, *ppbg1 = 0, *ppbg2 = 0, *ppbg3 = 0, *ppbgt1 = 0, *ppbgt2 = 0, *ppbgt3 = 0, *ppbgx1 = 0, *ppbgx2 = 0, *ppbgx3 = 0;
    float *phrelease = 0, *phrelease2 = 0, *phrelease3 = 0, *qhrelease = 0, *qhrelease2 = 0, *qhrelease3 = 0, *pvinx = 0, *pvinx2 = 0, *pvinx3 = 0, *pviny = 0, *pviny2 = 0, *pviny3 = 0,
    *phentrmax = 0, *phentrmax2 = 0, *phentrmax3 = 0, *pzfrag = 0, *pcentr = 0, *pcdeform = 0, *pphi = 0, *pphi2 = 0, *pphi3 = 0, *pdelta = 0, *pdelta2 = 0, *pdelta3 = 0, 
    *paddfri = 0, *paddfri2 = 0, *paddfri3 = 0, *pcohesion = 0, *pcohesion2 = 0, *pcohesion3 = 0, *pny = 0, *pny2 = 0, *pny3 = 0, *pambdrag = 0, *ptrans12 = 0, *ptrans13 = 0, *ptrans23 = 0, 
    *ptrelease = 0, *ptrelstop = 0, *pstoptime = 0, *ptslide = 0, *phdeposit = 0, qtinit[hydnin+hydnout],  *betax2 = 0, *betay2 = 0, *betaxy2 = 0, *betaxh2 = 0, *betayh2 = 0, 
    *betax3 = 0, *betay3 = 0, *betaxy3 = 0, *betaxh3 = 0, *betayh3 = 0, *htsun = 0, *htsunmax = 0;

    int *cdomain = (int*) calloc( sico.IMAX, sizeof(int));
    int *cdomain2 = (int*) calloc( sico.IMAX, sizeof(int));
    int *cplain = (int*) calloc( sico.IMAX, sizeof(int));
    int *cstopped = (int*) calloc( sico.IMAX, sizeof(int));
    int *pxslide = (int*) calloc( sico.IMAX, sizeof(int));

    qelev = (float*) calloc( sico.IMAX, sizeof(float));
    relev = (float*) calloc( sico.IMAX, sizeof(float));

    float *betax = (float*) calloc( sico.IMAX, sizeof(float));
    float *betay = (float*) calloc( sico.IMAX, sizeof(float));
    float *betaxh = (float*) calloc( sico.IMAX, sizeof(float));
    float *betayh = (float*) calloc( sico.IMAX, sizeof(float));
    float *betaxy = (float*) calloc( sico.IMAX, sizeof(float));
    float *pelev0 = (float*) calloc( sico.IMAX, sizeof(float));
    float *anx = (float*) calloc( sico.IMAX, sizeof(float));
    float *anu = (float*) calloc( sico.IMAX, sizeof(float));

    if ( sico.LAYERS0 > 0 ) {

        betax2 = (float*) calloc( sico.IMAX, sizeof(float));
        betay2 = (float*) calloc( sico.IMAX, sizeof(float));
        betaxy2 = (float*) calloc( sico.IMAX, sizeof(float));
        betaxh2 = (float*) calloc( sico.IMAX, sizeof(float));
        betayh2 = (float*) calloc( sico.IMAX, sizeof(float));
        betax3 = (float*) calloc( sico.IMAX, sizeof(float));
        betay3 = (float*) calloc( sico.IMAX, sizeof(float));
        betaxy3 = (float*) calloc( sico.IMAX, sizeof(float));
    }
    
    if ( sico.LAYERS0 > 0 || sico.SURFACE > 1 ) {   
    
        betaxh3 = (float*) calloc( sico.IMAX, sizeof(float));
        betayh3 = (float*) calloc( sico.IMAX, sizeof(float));
    }

    if ( sico.ZONES != 1 ) pzones = (int*) calloc( sico.IMAX, sizeof(int));

    if ( sico.TRELEASE == 1 ) qhrelease = (float*) calloc( sico.IMAX, sizeof(float));
    
    if ( sico.MODEL == 7 ) {
       
        if ( sico.TRELEASE == 1 ) qhrelease2 = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.TRELEASE == 1 ) qhrelease3 = (float*) calloc( sico.IMAX, sizeof(float));
    }
    

    #ifdef WITHGRASS


        float *v;
        float ***outv = alloc_dmatrix3( sico.M, sico.N, nvect_all );

        v = (float*) calloc( nvect_all, sizeof(float));
        px = (int*) calloc( sico.IMAX, sizeof(int));
        py = (int*) calloc( sico.IMAX, sizeof(int));

        pelev = (float*) calloc( sico.IMAX, sizeof(float));
        
        phrelease = (float*) calloc( sico.IMAX, sizeof(float));
        pvinx = (float*) calloc( sico.IMAX, sizeof(float));
        pviny = (float*) calloc( sico.IMAX, sizeof(float));
        phentrmax = (float*) calloc( sico.IMAX, sizeof(float));

        if ( sico.MODEL == 7 ) {

            phrelease2 = (float*) calloc( sico.IMAX, sizeof(float));
            pvinx2 = (float*) calloc( sico.IMAX, sizeof(float));
            pviny2 = (float*) calloc( sico.IMAX, sizeof(float));
            phentrmax2 = (float*) calloc( sico.IMAX, sizeof(float));

            phrelease3 = (float*) calloc( sico.IMAX, sizeof(float));
            pvinx3 = (float*) calloc( sico.IMAX, sizeof(float));
            pviny3 = (float*) calloc( sico.IMAX, sizeof(float));
            phentrmax3 = (float*) calloc( sico.IMAX, sizeof(float));
        }

        if ( sico.ZONES == 1 ) pzones = (int*) calloc( sico.IMAX, sizeof(int));
        if ( sico.CENTR == 1 ) pcentr = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.CDEFORM == 1 ) pcdeform = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.PHI == 1 ) pphi = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.PHI2 == 1 ) pphi2 = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.PHI3 == 1 ) pphi3 = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.DELTA == 1 ) pdelta = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.DELTA2 == 1 ) pdelta2 = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.DELTA3 == 1 ) pdelta3 = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.ADDFRI == 1 ) paddfri = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.ADDFRI2 == 1 ) paddfri2 = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.ADDFRI3 == 1 ) paddfri3 = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.COHESION == 1 ) pcohesion = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.COHESION2 == 1 ) pcohesion2 = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.COHESION3 == 1 ) pcohesion3 = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.NY == 1 ) pny = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.NY2 == 1 ) pny2 = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.NY3 == 1 ) pny3 = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.AMBDRAG == 1 ) pambdrag = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.TRANS12 == 1 ) ptrans12 = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.TRANS13 == 1 ) ptrans13 = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.TRANS23 == 1 ) ptrans23 = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.ZFRAG == 1 ) pzfrag = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.TRELEASE == 1 ) ptrelease = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.TRELSTOP == 1 ) ptrelstop = (float*) calloc( sico.IMAX, sizeof(float));
        pstoptime = (float*) calloc( sico.IMAX, sizeof(float));
        ptslide = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.IMPACTAREA == 1 ) pimpactarea = (int*) calloc( sico.IMAX, sizeof(int));
        if ( sico.HDEPOSIT == 1 ) phdeposit = (float*) calloc( sico.IMAX, sizeof(float));
        
        ppbg1 = (int*) calloc( sico.IMAX, sizeof(int));
        ppbg2 = (int*) calloc( sico.IMAX, sizeof(int));
        ppbg3 = (int*) calloc( sico.IMAX, sizeof(int));


    #endif
   
    int **icor = alloc_imatrix( sico.M, sico.N );
    int **in = alloc_imatrix( sico.IMAX, 9 );
    int **ibasket = alloc_imatrix( 2, sico.IMAX );
    int **icheck = alloc_imatrix( sico.IMAX, 2 );
    int *ib = (int*) calloc( 2, sizeof(int));

    float *asigma_xelev = calloc( sico.IMAX, sizeof(float));
    float *asigma_yelev = calloc( sico.IMAX, sizeof(float));

    kappau = (float*) calloc( 6, sizeof(float));

    float **aw = alloc_dmatrix( sico.IMAX, sico.NVECTMIN );
    float **awt = alloc_dmatrix( sico.IMAX, nvect_all );
    float **af = alloc_dmatrix( sico.IMAX, sico.NVECTMIN );
    float **ag = alloc_dmatrix( sico.IMAX, sico.NVECTMIN );
    float **as = alloc_dmatrix( sico.IMAX, sico.NVECTMIN );
    float **ad = alloc_dmatrix( sico.IMAX, sico.NVECTMIN );

    float **asigma_x = alloc_dmatrix( sico.IMAX, sico.NVECTMIN );
    float **asigma_y = alloc_dmatrix( sico.IMAX, sico.NVECTMIN );
    float **asigma_f = alloc_dmatrix( sico.IMAX, sico.NVECTMIN );
    float **asigma_g = alloc_dmatrix( sico.IMAX, sico.NVECTMIN );

    float **wintelev = alloc_dmatrix( sico.IMAX, 4 );
    float *wintelevd = (float*) calloc( sico.IMAX, sizeof(float));
    float ***winta = alloc_dmatrix3( sico.IMAX, 4, sico.NVECTMIN );
    float ***wintb = alloc_dmatrix3( sico.IMAX, 6, sico.NVECTMIN );
    float ***wintc = alloc_dmatrix3( sico.IMAX, 4, sico.NVECTMIN );
    float **wintd = alloc_dmatrix( sico.IMAX, sico.NVECTMIN );
    float **wintdtest = alloc_dmatrix( sico.IMAX, sico.NVECTMIN );

    float **f = alloc_dmatrix( 4, sico.NVECTMIN );
    float **g = alloc_dmatrix( 4, sico.NVECTMIN );
    float **s = alloc_dmatrix( 4, sico.NVECTMIN );
    float ***d = alloc_dmatrix3( sico.IMAX, 6, sico.NVECTMIN );

    float *dx = (float*) calloc( sico.IMAX, sizeof(float));
    float *dy = (float*) calloc( sico.IMAX, sizeof(float));
    
    if ( sico.PBG == 1 ) {
            
        ppbgt1 = (int*) calloc( sico.IMAX, sizeof(int));
        ppbgt2 = (int*) calloc( sico.IMAX, sizeof(int));
        ppbgt3 = (int*) calloc( sico.IMAX, sizeof(int));
        ppbgx1 = (int*) calloc( sico.IMAX, sizeof(int));
        ppbgx2 = (int*) calloc( sico.IMAX, sizeof(int));
        ppbgx3 = (int*) calloc( sico.IMAX, sizeof(int));
        anpx = (float*) calloc( sico.IMAX, sizeof(float));
        anpy = (float*) calloc( sico.IMAX, sizeof(float));
        anctr = (float*) calloc( sico.IMAX, sizeof(float));
        
                     
        timer=time( NULL ); // timer for randomization
        srand( timer ); // seed for random number
        randmax = 1000; // maximum for random number
    }
    

// -- STOP --- Preparing arrays for input -----------------------------------------------------------------------


    if ( sico.TSUNAMI != 0 ) {
               
            htsun = (float*) calloc( sico.IMAX, sizeof(float));
            htsunmax = (float*) calloc( sico.IMAX, sizeof(float));
    } 

    mv0 = (char**) calloc( 1000, sizeof(char*));
    mv = (char*) calloc( 1000, sizeof(char));
    mv2 = (char*) calloc( 1000, sizeof(char));

    if ( sico.MODEL <= 3 ) {

        mv0[0] = "hflow"; mv0[1] = "vflowx"; mv0[2] = "vflowy"; mv0[3] = "temperature"; mv0[4] = "basechange"; mv0[5] = "vflow";
        mv0[6] = "tflow"; mv0[7] = "pflow"; mv0[8] = "hflow_max"; mv0[9] = "vflow_max";
        mv0[10] = "tflow_max"; mv0[11] = "pflow_max"; mv0[12] = "vfront"; mv0[13] = "r1front"; mv0[14] = "r3front"; mv0[15] = "r1max"; mv0[16] = "r3max"; mv0[17] = "vhmax"; mv0[18] = "treach";
        
    } else if ( sico.MODEL == 7 ) {

        mv0[0]  = "hflow1";     mv0[1]  = "vflowx1";    mv0[2]  = "vflowy1";    mv0[3]  = "hflow2";      mv0[4]  = "vflowx2";     mv0[5] = "vflowy2";
        mv0[6]  = "hflow3";     mv0[7]  = "vflowx3";    mv0[8]  = "vflowy3";    mv0[9]  = "temperature"; mv0[10] = "basechange1"; mv0[11] = "basechange2"; mv0[12] = "basechange3";
        mv0[13] = "vflow1";     mv0[14] = "vflow2";     mv0[15] = "vflow3";     mv0[16] = "hflow";       mv0[17] = "tflow1";      mv0[18] = "tflow2";
        mv0[19] = "tflow3";     mv0[20] = "tflow";      mv0[21] = "pflow1";     mv0[22] = "pflow2";      mv0[23] = "pflow3";      mv0[24] = "pflow";
        mv0[25] = "basechange"; mv0[26] = "hflow1_max"; mv0[27] = "vflow1_max"; mv0[28] = "hflow2_max";  mv0[29] = "vflow2_max";  mv0[30] = "hflow3_max";
        mv0[31] = "vflow3_max"; mv0[32] = "hflow_max";  mv0[33] = "tflow1_max"; mv0[34] = "pflow1_max";  mv0[35] = "tflow2_max";  mv0[36] = "pflow2_max";
        mv0[37] = "tflow3_max"; mv0[38] = "pflow3_max"; mv0[39] = "tflow_max";  mv0[40] = "pflow_max";   mv0[41] = "vfront";      mv0[42] = "r1front"; 
        mv0[43] = "r3front";    mv0[44] = "r1max";      mv0[45] = "r3max";      mv0[46] = "vhmax";       mv0[47] = "treach"; // variables for names of output raster maps
    }


// -- START -- Defining, opening, and pre-processing output files -----------------------------------------------


    FILE *f_summary, *f_profile = 0, *f_profile_aimec = 0, *f_ctrlpoints = 0, *f_evaluation = 0, *f_evaluationh = 0, *f_aimec = 0, *f_aimech = 0, *f_volumes, 
        *f_directions = 0, *f_directions2 = 0, *f_directions3 = 0, *f_nout, *f_hydout = 0, *f_hydinfo[hydnin+hydnout], *f_hydtrans[hydnin+hydnout], 
        *f_vr[(int)(tmax/tout+3)], *f_paraviewi, *f_paraviewc, *f_blenderi, *f_unreali, *f_rmultval, *f_rroc, *f_rhydrograph, *f_chydrograph, *f_rprofile, *f_cprofile, *f_rmap, *f_cmap;

    if ( sico.MULT == 0 ) sprintf(path, "%s%ssummary.txt", outfiles, prefix); // summary file
    else sprintf(path, "%s%ssummary%d.txt", outfiles, prefix, xint);
    f_summary=fopen(path, "w");

    if ( sico.CTRLPOINTS > 0 ) {

        if ( sico.MULT == 0 ) sprintf(path, "%s%sctrlpoints.txt", outfiles, prefix); // control point data file
        else sprintf(path, "%s%sctrlpoints%d.txt", outfiles, prefix, xint);
        f_ctrlpoints=fopen(path, "w");
    }

    if ( sico.MULT == 0 ) sprintf(path, "%s%svolumes.txt", outfiles, prefix); // volumes file
    else sprintf(path, "%s%svolumes%d.txt", outfiles, prefix, xint);      
    f_volumes=fopen(path, "w");

    if ( hydrograph == 1 ) {

        sprintf(path, "%s%shydprofiles.txt", outfiles, prefix); // hydrograph profiles file
        f_hydout=fopen(path, "w");

        fprintf(f_hydout, "ID\tx1\txC\tx2\ty1\tyC\ty2\n"); // printing header of hydrograph profiles file

        if ( sico.MULT == 0 ) {

            for ( i = 0; i < hydnin + hydnout; i++ ) { // loop over all hydrographs

                sprintf(path, "%s%shydinfo%i.txt", outfiles, prefix, i+1); // hydrograph info file
                f_hydinfo[i]=fopen(path, "w");

                if ( i >= hydnin ) fprintf(f_hydinfo[i],
                    "T\tH1\tV1\tE1\tQ1\tH2\tV2\tE2\tQ2\tH3\tV3\tE3\tQ3\n0.0\t0.00\t0.00\t0.00\t0.00\t0.00\t0.00\t0.00\t0.00\t0.00\t0.00\t0.00\t0.00\n");
                else fprintf(f_hydinfo[i], "T\tH1\tV1\tE1\tQ1\tH2\tV2\tE2\tQ2\tH3\tV3\tE3\tQ3\n"); // printing header of hydrograph info file
                
                sprintf(path, "%s%shydtrans%i.txt", outfiles, prefix, i+1); // hydrograph transfer file
                f_hydtrans[i]=fopen(path, "w");

                fprintf(f_hydtrans[i], "T\tQ1\tV1\tQ2\tV2\tQ3\tV3\n"); // printing header of hydrograph transfer file
            }

        }
    }

    if ( sico.MULT == 0 ) {
    
        sprintf(path, "%s%sdirections1.txt", outfiles, prefix); // PHASE 1 flow direction file (for display as arrows)
        f_directions=fopen(path, "w");

        if ( sico.MODEL == 7 ) {
        
            sprintf(path, "%s%sdirections2.txt", outfiles, prefix); // PHASE 2 flow direction file (for display as arrows)
            f_directions2=fopen(path, "w");

            sprintf(path, "%s%sdirections3.txt", outfiles, prefix); // PHASE 3 flow direction file (for display as arrows)
            f_directions3=fopen(path, "w");
        }
        
    } else if ( xint == 1 ) {
        
        sprintf(path, "%s%sevaluationh.txt", outfiles, prefix); // evaluation file (header for multiple model runs)
        f_evaluationh=fopen(path, "w");
        
        if ( sico.MODEL == 1 ) fprintf( f_evaluationh, "nrun\tvhr\tvhem\trho1\tphi1\tdelta1\taddfri1\tcohesion1\tny1\tcdeform1\tslideexp1\tslidefract1\tesp\tfragsml\tfragml\tcambdrag\tcentr\tcstop\tVE1\tVE2\tVE3\tLio\tLdo\tLdm\tdLi\trLi\tLdm\tdLd\trLd\trTPi\trTNi\trFPi\trFNi\tCSIi\tHSSi\tAUROCi\tD2PCi\tFoCi\tSPIi\trTPd\trTNd\trFPd\trFNd\tCSId\tHSSd\tAUROCd\tD2PCd\tFoCd\tSPId\n" );
        else if ( sico.MODEL == 7 ) fprintf( f_evaluationh, "nrun\tvhr\trhrs\tvhem\trhes\trho1\trho2\trho3\tphi1\tphi2\tphi3\tdelta1\tdelta2\tdelta3\taddfri1\taddfri2\taddfri3\tcohesion1\tcohesion2\tcohesion3\tny1\tny2\tny3\tcdeform1\tcdeform2\tcdeform3\tslideexp1\tslideexp2\tslideexp3\tslidefract1\tslidefract2\tslidefract3\tesp\tfragsml\tfragml\tcambdrag\tcentr\tcstop\tdrk\tdrm\tdrn\tdrvterm\tdrrep\tdrje\tvmn0\tvml\tvmn\tctrans12\tctrans13\tctrans23\tmtempm\tmtempa\tmtempb\tmfmelt\tmfslide\tVE1\tVE2\tVE3\tLio\tLdo\tLdm\tdLi\trLi\tLdm\tdLd\trLd\trTPi\trTNi\trFPi\trFNi\tCSIi\tHSSi\tAUROCi\tD2PCi\tFoCi\tSPIi\trTPd\trTNd\trFPd\trFNd\tCSId\tHSSd\tAUROCd\tD2PCd\tFoCd\tSPId\n" );
        
        fclose(f_evaluationh);

        sprintf(path, "%s%saimech.txt", outaimec, prefix); // aimec file (header for multiple model runs)
        f_aimech=fopen(path, "w");
        
        if ( sico.MODEL == 1 ) fprintf( f_aimech, "nrun, vhr, vhem, rho1, phi1, delta1, addfri1, cohesion1, ny1, cdeform1, slideexp1, slidefract1, esp, fragsml, fragml, cambdrag, centr, cstop\n");
        else if ( sico.MODEL == 7 ) fprintf( f_aimech, "nrun, vhr, rhrs, vhem, rhes, rho1, rho2, rho3, phi1, phi2, phi3, delta1, delta2, delta3, addfri1, addfri2, addfri3, cohesion1, cohesion2, cohesion3, ny1, ny2, ny3, cdeform1, cdeform2, cdeform3, slideexp1, slideexp2, slideexp3, slidefract1, slidefract2, slidefract3, esp, fragsml, fragml, cambdrag, centr, cstop, drk, drm, drn, drvterm, drrep, drje, vmn0, vml, vmn, ctrans12, ctrans13, ctrans23. mtempm, mtempa, mtempb, mfmelt, mfslide\n");
        
        fclose(f_aimech);
    }


// -- STOP --- Defining, opening, and pre-processing output files -----------------------------------------------


// -- START -- Preparing arrays (if GRASS is used) --------------------------------------------------------------


    #ifdef WITHGRASS


        i = 0;

        for ( x = 0; x < sico.M; x++ ) { 
          for ( y = 0; y < sico.N; y++ ) {

            aoi = 1;
            Segment_get(&seg_elev, &elev, x, y); // reading data from segmentation files
            
            if ( sico.RELM == 1 ) Segment_get(&seg_hrelease, &hrelease, x, y);
            if ( sico.MODEL == 7 && sico.RELM2 == 1 ) Segment_get(&seg_hrelease2, &hrelease2, x, y);
            if ( sico.MODEL == 7 && sico.RELM3 == 1 ) Segment_get(&seg_hrelease3, &hrelease3, x, y);
            
            if ( sico.RELV == 1 ) Segment_get(&seg_vinx, &vinx, x, y);
            if ( sico.MODEL == 7 && sico.RELV2 == 1 ) Segment_get(&seg_vinx2, &vinx2, x, y);
            if ( sico.MODEL == 7 && sico.RELV3 == 1 ) Segment_get(&seg_vinx3, &vinx3, x, y);
            
            if ( sico.RELV == 1 ) Segment_get(&seg_viny, &viny, x, y);
            if ( sico.MODEL == 7 && sico.RELV2 == 1 ) Segment_get(&seg_viny2, &viny2, x, y);
            if ( sico.MODEL == 7 && sico.RELV3 == 1 ) Segment_get(&seg_viny3, &viny3, x, y);
            
            if ( sico.ENTR == 1 ) Segment_get(&seg_hentrmax, &hentrmax, x, y);
            if ( sico.MODEL == 7 && sico.ENTR2 == 1 ) Segment_get(&seg_hentrmax2, &hentrmax2, x, y);
            if ( sico.MODEL == 7 && sico.ENTR3 == 1 ) Segment_get(&seg_hentrmax3, &hentrmax3, x, y);

            if ( sico.ZONES == 1 ) Segment_get(&seg_zones, &zones, x, y);     
            if ( sico.CENTR == 1 ) Segment_get(&seg_centr, &centr, x, y);
            if ( sico.CDEFORM == 1 ) Segment_get(&seg_cdeform, &cdeform, x, y);
            if ( sico.PHI == 1 ) Segment_get(&seg_phi, &phi, x, y);
            if ( sico.PHI2 == 1 ) Segment_get(&seg_phi2, &phi2, x, y);
            if ( sico.PHI3 == 1 ) Segment_get(&seg_phi3, &phi3, x, y);
            if ( sico.DELTA == 1 ) Segment_get(&seg_delta, &delta, x, y);
            if ( sico.DELTA2 == 1 ) Segment_get(&seg_delta2, &delta2, x, y);
            if ( sico.DELTA3 == 1 ) Segment_get(&seg_delta3, &delta3, x, y);
            if ( sico.ADDFRI == 1 ) Segment_get(&seg_addfri, &addfri, x, y);
            if ( sico.ADDFRI2 == 1 ) Segment_get(&seg_addfri2, &addfri2, x, y);
            if ( sico.ADDFRI3 == 1 ) Segment_get(&seg_addfri3, &addfri3, x, y);
            if ( sico.COHESION == 1 ) Segment_get(&seg_cohesion, &cohesion, x, y);
            if ( sico.COHESION2 == 1 ) Segment_get(&seg_cohesion2, &cohesion2, x, y);
            if ( sico.COHESION3 == 1 ) Segment_get(&seg_cohesion3, &cohesion3, x, y);
            if ( sico.NY == 1 ) Segment_get(&seg_ny, &ny, x, y);
            if ( sico.NY2 == 1 ) Segment_get(&seg_ny2, &ny2, x, y);
            if ( sico.NY3 == 1 ) Segment_get(&seg_ny3, &ny3, x, y);
            if ( sico.AMBDRAG == 1 ) Segment_get(&seg_ambdrag, &ambdrag, x, y);
            if ( sico.TRANS12 == 1 ) Segment_get(&seg_trans12, &trans12, x, y);
            if ( sico.TRANS13 == 1 ) Segment_get(&seg_trans13, &trans13, x, y);
            if ( sico.TRANS23 == 1 ) Segment_get(&seg_trans23, &trans23, x, y);
            if ( sico.ZFRAG == 1 ) Segment_get(&seg_zfrag, &zfrag, x, y);
            if ( sico.TRELEASE == 1 ) Segment_get(&seg_trelease, &trelease, x, y);
            if ( sico.TRELSTOP == 1 ) Segment_get(&seg_trelstop, &trelstop, x, y);
            if ( sico.STOPTIME == 1 ) Segment_get(&seg_stoptime, &stoptime, x, y);
            if ( sico.TSLIDE == 1 ) Segment_get(&seg_tslide, &tslide, x, y);
            if ( sico.IMPACTAREA == 1 ) Segment_get(&seg_impactarea, &impactarea, x, y);
            if ( sico.HDEPOSIT == 1 ) Segment_get(&seg_hdeposit, &hdeposit, x, y);
            if ( sico.PBG == 1 ) {
                Segment_get(&seg_pbg1, &pbg1, x, y);
                Segment_get(&seg_pbg2, &pbg2, x, y);
                Segment_get(&seg_pbg3, &pbg3, x, y);                                          
            }

            if ( aoi >= 1 ) {

                px[i] = x; // internal x and y coordinates of cell
                py[i] = y;
                icor[x][y] = i;

                pelev[i] = elev; pelev0[i] = elev;

                if ( sico.RELM == 1 && hrelease != sico.UNDEF ) phrelease[i] = hrelease; // writing data to arrays
                else phrelease[i] = 0;
                if ( sico.MODEL == 7 && sico.RELM2 == 1  && hrelease2 != sico.UNDEF ) phrelease2[i] = hrelease2;
                else if ( sico.MODEL == 7 ) phrelease2[i] = 0;
                if ( sico.MODEL == 7 && sico.RELM3 == 1  && hrelease3 != sico.UNDEF ) phrelease3[i] = hrelease3;
                else if ( sico.MODEL == 7 ) phrelease3[i] = 0;
                if ( sico.RELV == 1 && vinx != sico.UNDEF ) pvinx[i] = vinx;
                else pvinx[i] = 0;
                if ( sico.MODEL == 7 && sico.RELV2 == 1  && vinx2 != sico.UNDEF ) pvinx2[i] = vinx2;
                else if ( sico.MODEL == 7 ) pvinx2[i] = 0;
                if ( sico.MODEL == 7 && sico.RELV3 == 1  && vinx3 != sico.UNDEF ) pvinx3[i] = vinx3;
                else if ( sico.MODEL == 7 ) pvinx3[i] = 0;
                if ( sico.RELV == 1 && viny != sico.UNDEF ) pviny[i] = viny;
                else pviny[i] = 0;
                if ( sico.MODEL == 7 && sico.RELV2 == 1  && viny2 != sico.UNDEF ) pviny2[i] = viny2;
                else if ( sico.MODEL == 7 ) pviny2[i] = 0;
                if ( sico.MODEL == 7 && sico.RELV3 == 1  && viny3 != sico.UNDEF ) pviny3[i] = viny3;
                else if ( sico.MODEL == 7 ) pviny3[i] = 0;
                if ( sico.ENTR == 1 ) phentrmax[i] = hentrmax;
                else phentrmax[i] = 9999;
                if ( sico.MODEL == 7 && sico.ENTR2 == 1 && hentrmax2 != sico.UNDEF ) phentrmax2[i] = hentrmax2;
                else if ( sico.MODEL == 7 ) phentrmax2[i] = 0;
                if ( sico.MODEL == 7 && sico.ENTR3 == 1 && hentrmax3 != sico.UNDEF ) phentrmax3[i] = hentrmax3;
                else if ( sico.MODEL == 7 ) phentrmax3[i] = 0;
                if ( sico.ZONES == 1 && zones != sico.UNDEF ) pzones[i] = zones;
                else if ( sico.ZONES== 1 ) pzones[i] = 0;        
                if ( sico.CENTR == 1 && centr != sico.UNDEF ) pcentr[i] = centr;
                else if ( sico.CENTR == 1 ) pcentr[i] = sico.UNDEF;
                if ( sico.PHI == 1 && phi != sico.UNDEF ) pphi[i] = phi * sico.PI / 180;
                else if ( sico.PHI == 1 ) pphi[i] = sico.UNDEF;
                if ( sico.PHI2 == 1 && phi2 != sico.UNDEF ) pphi2[i] = phi2 * sico.PI / 180;
                else if ( sico.PHI2 == 1 ) pphi2[i] = sico.UNDEF;
                if ( sico.PHI3 == 1 && phi3 != sico.UNDEF ) pphi3[i] = phi3 * sico.PI / 180;
                else if ( sico.PHI3 == 1 ) pphi3[i] = sico.UNDEF;
                if ( sico.DELTA == 1 && delta != sico.UNDEF ) pdelta[i] = delta * sico.PI / 180;
                else if ( sico.DELTA == 1 ) pdelta[i] = sico.UNDEF;
                if ( sico.DELTA2 == 1 && delta2 != sico.UNDEF ) pdelta2[i] = delta2 * sico.PI / 180;
                else if ( sico.DELTA2 == 1 ) pdelta2[i] = sico.UNDEF;
                if ( sico.DELTA3 == 1 && delta3 != sico.UNDEF ) pdelta3[i] = delta3 * sico.PI / 180;
                else if ( sico.DELTA3 == 1 ) pdelta3[i] = sico.UNDEF;
                if ( sico.ADDFRI == 1 && addfri != sico.UNDEF ) paddfri[i] = addfri;
                else if ( sico.ADDFRI == 1 ) paddfri[i] = sico.UNDEF;
                if ( sico.ADDFRI2 == 1 && addfri2 != sico.UNDEF ) paddfri2[i] = addfri2;
                else if ( sico.ADDFRI2 == 1 ) paddfri2[i] = sico.UNDEF;
                if ( sico.ADDFRI3 == 1 && addfri3 != sico.UNDEF ) paddfri3[i] = addfri3;
                else if ( sico.ADDFRI3 == 1 ) paddfri3[i] = sico.UNDEF;
                if ( sico.COHESION == 1 && cohesion != sico.UNDEF ) pcohesion[i] = cohesion;
                else if ( sico.COHESION == 1 ) pcohesion[i] = sico.UNDEF;
                if ( sico.COHESION2 == 1 && cohesion2 != sico.UNDEF ) pcohesion2[i] = cohesion2;
                else if ( sico.COHESION2 == 1 ) pcohesion2[i] = sico.UNDEF;
                if ( sico.COHESION3 == 1 && cohesion3 != sico.UNDEF ) pcohesion3[i] = cohesion3;
                else if ( sico.COHESION3 == 1 ) pcohesion3[i] = sico.UNDEF;
                if ( sico.NY == 1 && ny != sico.UNDEF ) pny[i] = pow( 10, ny );
                else if ( sico.NY == 1 ) pny[i] = sico.UNDEF;
                if ( sico.NY2 == 1 && ny2 != sico.UNDEF ) pny2[i] = pow( 10, ny2 );
                else if ( sico.NY2 == 1 ) pny2[i] = sico.UNDEF;
                if ( sico.NY3 == 1 && ny3 != sico.UNDEF ) pny3[i] = pow( 10, ny3 );
                else if ( sico.NY3 == 1 ) pny3[i] = sico.UNDEF;
                if ( sico.CDEFORM == 1 && cdeform != sico.UNDEF ) pcdeform[i] = cdeform;
                else if ( sico.CDEFORM == 1 ) pcdeform[i] = sico.UNDEF;
                if ( sico.AMBDRAG == 1 && ambdrag != sico.UNDEF ) pambdrag[i] = ambdrag;
                else if ( sico.AMBDRAG == 1 ) pambdrag[i] = sico.UNDEF;
                if ( sico.TRANS12 == 1 && trans12 != sico.UNDEF ) ptrans12[i] = trans12;
                else if ( sico.TRANS12 == 1 ) ptrans12[i] = sico.UNDEF;
                if ( sico.TRANS13 == 1 && trans13 != sico.UNDEF ) ptrans13[i] = trans13;
                else if ( sico.TRANS13 == 1 ) ptrans13[i] = sico.UNDEF;
                if ( sico.TRANS23 == 1 && trans23 != sico.UNDEF ) ptrans23[i] = trans23;
                else if ( sico.TRANS23 == 1 ) ptrans23[i] = sico.UNDEF;
                if ( sico.ZFRAG == 1 && zfrag != sico.UNDEF ) pzfrag[i] = zfrag;
                else if ( sico.ZFRAG == 1 ) pzfrag[i] = sico.UNDEF;
                if ( sico.TRELEASE == 1 && trelease != sico.UNDEF ) ptrelease[i] = trelease;
                else if ( sico.TRELEASE == 1 ) ptrelease[i] = sico.UNDEF;
                if ( sico.TRELSTOP == 1 && trelstop != sico.UNDEF ) ptrelstop[i] = trelstop;
                else if ( sico.TRELSTOP == 1 ) ptrelstop[i] = sico.UNDEF;
                if ( sico.STOPTIME == 1 && stoptime != sico.UNDEF ) pstoptime[i] = stoptime;
                else pstoptime[i] = 5;
                if ( sico.TSLIDE == 1 && tslide != sico.UNDEF ) ptslide[i] = tslide;
                else ptslide[i] = sico.UNDEF;
                if ( sico.IMPACTAREA == 1 && impactarea != sico.UNDEF ) pimpactarea[i] = impactarea;
                else if ( sico.IMPACTAREA == 1 ) pimpactarea[i] = sico.UNDEF;
                if ( sico.HDEPOSIT == 1 && hdeposit != sico.UNDEF ) phdeposit[i] = hdeposit;
                else if ( sico.HDEPOSIT == 1 ) phdeposit[i] = sico.UNDEF;
                
                if ( sico.PBG == 1 && pbg1 != sico.UNDEF ) ppbg1[i] = pbg1;
                else if ( sico.PBG == 1 ) ppbg1[i] = sico.UNDEF;
                if ( sico.PBG == 1 && pbg2 != sico.UNDEF ) ppbg2[i] = pbg2;
                else if ( sico.PBG == 1 ) ppbg2[i] = sico.UNDEF;
                if ( sico.PBG == 1 && pbg3 != sico.UNDEF ) ppbg3[i] = pbg3;
                else if ( sico.PBG == 1 ) ppbg3[i] = sico.UNDEF;
                               
                i += 1;
            }
          }
        }


// -- STOP --- Preparing arrays (if GRASS is used) --------------------------------------------------------------


    #else


// -- START -- Reading ascii rasters and pre-processing arrays (if GRASS is not used) ---------------------------


        px = finascxy( indir, elevname, 1, sico ); // internal x and y coordinates of cell
        py = finascxy( indir, elevname, 2, sico );

        pelev = finascval( indir, elevname, px, py, sico ); // reading data from ascii rasters
        
        if ( sico.RELM == 1 ) phrelease = finascval( indir, hreleasename, px, py, sico );
        else phrelease = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.MODEL == 7 && sico.RELM2 == 1 ) phrelease2 = finascval( indir, hreleasename2, px, py, sico );
        else if ( sico.MODEL == 7 ) phrelease2 = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.MODEL == 7 && sico.RELM3 == 1 ) phrelease3 = finascval( indir, hreleasename3, px, py, sico );
        else if ( sico.MODEL == 7 ) phrelease3 = (float*) calloc( sico.IMAX, sizeof(float));

        if ( sico.RELV == 1 ) pvinx = finascval( indir, vinxname, px, py, sico );
        else pvinx = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.MODEL == 7 && sico.RELV2 == 1 ) pvinx2 = finascval( indir, vinxname2, px, py, sico );
        else if ( sico.MODEL == 7 ) pvinx2 = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.MODEL == 7 && sico.RELV3 == 1 ) pvinx3 = finascval( indir, vinxname3, px, py, sico );
        else if ( sico.MODEL == 7 ) pvinx3 = (float*) calloc( sico.IMAX, sizeof(float));

        if ( sico.RELV == 1 ) pviny = finascval( indir, vinyname, px, py, sico );
        else pviny = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.MODEL == 7 && sico.RELV2 == 1 ) pviny2 = finascval( indir, vinyname2, px, py, sico ); 
        else if ( sico.MODEL == 7 ) pviny2 = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.MODEL == 7 && sico.RELV3 == 1 ) pviny3 = finascval( indir, vinyname3, px, py, sico );
        else if ( sico.MODEL == 7 ) pviny3 = (float*) calloc( sico.IMAX, sizeof(float));

        if ( sico.ENTR == 1 ) phentrmax = finascval( indir, hentrmaxname, px, py, sico );
        else phentrmax = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.MODEL == 7 && sico.ENTR2 == 1 ) phentrmax2 = finascval( indir, hentrmaxname2, px, py, sico );
        else if ( sico.MODEL == 7 ) phentrmax2 = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.MODEL == 7 && sico.ENTR3 == 1 ) phentrmax3 = finascval( indir, hentrmaxname3, px, py, sico );
        else if ( sico.MODEL == 7 ) phentrmax3 = (float*) calloc( sico.IMAX, sizeof(float));

        if ( sico.ZONES == 1 ) pzones = finascvali( indir, zonesname, px, py, sico );
        if ( sico.CENTR == 1 ) pcentr = finascval( indir, centrname, px, py, sico );
        if ( sico.PHI == 1 ) pphi = finascval( indir, phiname, px, py, sico );
        if ( sico.PHI2 == 1 ) pphi2 = finascval( indir, phi2name, px, py, sico );
        if ( sico.PHI3 == 1 ) pphi3 = finascval( indir, phi3name, px, py, sico );
        if ( sico.DELTA == 1 ) pdelta = finascval( indir, deltaname, px, py, sico );
        if ( sico.DELTA2 == 1 ) pdelta2 = finascval( indir, delta2name, px, py, sico );
        if ( sico.DELTA3 == 1 ) pdelta3 = finascval( indir, delta3name, px, py, sico );
        if ( sico.ADDFRI == 1 ) paddfri = finascval( indir, addfriname, px, py, sico );
        if ( sico.ADDFRI2 == 1 ) paddfri2 = finascval( indir, addfri2name, px, py, sico );
        if ( sico.ADDFRI3 == 1 ) paddfri3 = finascval( indir, addfri3name, px, py, sico );
        if ( sico.COHESION == 1 ) pcohesion = finascval( indir, cohesionname, px, py, sico );
        if ( sico.COHESION2 == 1 ) pcohesion2 = finascval( indir, cohesion2name, px, py, sico );
        if ( sico.COHESION3 == 1 ) pcohesion3 = finascval( indir, cohesion3name, px, py, sico );
        if ( sico.NY == 1 ) pny = finascval( indir, nyname, px, py, sico );
        if ( sico.NY2 == 1 ) pny2 = finascval( indir, ny2name, px, py, sico );
        if ( sico.NY3 == 1 ) pny3 = finascval( indir, ny3name, px, py, sico );
        if ( sico.CDEFORM == 1 ) pcdeform = finascval( indir, cdeformname, px, py, sico );
        if ( sico.AMBDRAG == 1 ) pambdrag = finascval( indir, ambdragname, px, py, sico );
        if ( sico.TRANS12 == 1 ) ptrans12 = finascval( indir, trans12name, px, py, sico );
        if ( sico.TRANS13 == 1 ) ptrans13 = finascval( indir, trans13name, px, py, sico );
        if ( sico.TRANS23 == 1 ) ptrans23 = finascval( indir, trans23name, px, py, sico );
        if ( sico.ZFRAG == 1 ) pzfrag = finascval( indir, zfragname, px, py, sico );
        if ( sico.TRELEASE == 1 ) ptrelease = finascval( indir, treleasename, px, py, sico );
        if ( sico.TRELSTOP == 1 ) ptrelstop = finascval( indir, trelstopname, px, py, sico );
        if ( sico.STOPTIME == 1 ) pstoptime = finascval( indir, stoptimename, px, py, sico );
        else pstoptime = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.TSLIDE == 1 ) ptslide = finascval( indir, tslidename, px, py, sico );
        else ptslide = (float*) calloc( sico.IMAX, sizeof(float));
        if ( sico.IMPACTAREA == 1 ) pimpactarea = finascvali( indir, impactareaname, px, py, sico );
        if ( sico.HDEPOSIT == 1 ) phdeposit = finascval( indir, hdepositname, px, py, sico );
        if ( sico.PBG == 1 ) ppbg1 = finascvali( indir, pbg1name, px, py, sico );
        if ( sico.PBG == 1 ) ppbg2 = finascvali( indir, pbg2name, px, py, sico );
        if ( sico.PBG == 1 ) ppbg3 = finascvali( indir, pbg3name, px, py, sico );

        for ( i=0; i<sico.IMAX; i++ ) { // finalizing input arrays

            icor[px[i]][py[i]] = i;
            pelev0[i] = pelev[i];
            
            if ( sico.RELM != 1 || phrelease[i] == sico.UNDEF ) phrelease[i] = 0;
            if ( sico.MODEL == 7 && ( sico.RELM2 != 1 || phrelease2[i] == sico.UNDEF )) phrelease2[i] = 0;
            if ( sico.MODEL == 7 && ( sico.RELM3 != 1 || phrelease3[i] == sico.UNDEF )) phrelease3[i] = 0;

            if ( sico.RELV != 1 || pvinx[i] == sico.UNDEF ) pvinx[i] = 0;
            if ( sico.MODEL == 7 && ( sico.RELV2 != 1 || pvinx2[i] == sico.UNDEF )) pvinx2[i] = 0;
            if ( sico.MODEL == 7 && ( sico.RELV3 != 1 || pvinx3[i] == sico.UNDEF )) pvinx3[i] = 0;

            if ( sico.RELV != 1 || pviny[i] == sico.UNDEF ) pviny[i] = 0;
            if ( sico.MODEL == 7 && ( sico.RELV2 != 1 || pviny2[i] == sico.UNDEF )) pviny2[i] = 0;
            if ( sico.MODEL == 7 && ( sico.RELV3 != 1 || pviny3[i] == sico.UNDEF )) pviny3[i] = 0;

            if ( sico.ENTR != 1 || phentrmax[i] == sico.UNDEF ) phentrmax[i] = 9999;
            if ( sico.MODEL == 7 && ( sico.ENTR2 != 1 || phentrmax2[i] == sico.UNDEF )) phentrmax2[i] = 0;
            if ( sico.MODEL == 7 && ( sico.ENTR3 != 1 || phentrmax3[i] == sico.UNDEF )) phentrmax3[i] = 0;

            if ( sico.ZONES == 1 && pzones[i] != sico.UNDEF ) pzones[i] = pzones[i];
            else if ( sico.ZONES == 1 ) pzones[i] = 0;
            if ( sico.CENTR == 1 && pcentr[i] != sico.UNDEF ) pcentr[i] = pcentr[i];
            else if ( sico.CENTR == 1 ) pcentr[i] = sico.UNDEF;
            if ( sico.PHI == 1 && pphi[i] != sico.UNDEF ) pphi[i] = pphi[i] * sico.PI / 180;
            else if ( sico.PHI == 1 ) pphi[i] = sico.UNDEF;
            if ( sico.PHI2 == 1 && pphi2[i] != sico.UNDEF ) pphi2[i] = pphi2[i] * sico.PI / 180;
            else if ( sico.PHI2 == 1 ) pphi2[i] = sico.UNDEF;
            if ( sico.PHI3 == 1 && pphi3[i] != sico.UNDEF ) pphi3[i] = pphi3[i] * sico.PI / 180;
            else if ( sico.PHI3 == 1 ) pphi3[i] = sico.UNDEF;
            if ( sico.DELTA == 1 && pdelta[i] != sico.UNDEF ) pdelta[i] = pdelta[i] * sico.PI / 180;
            else if ( sico.DELTA == 1 ) pdelta[i] = sico.UNDEF;
            if ( sico.DELTA2 == 1 && pdelta2[i] != sico.UNDEF ) pdelta2[i] = pdelta2[i] * sico.PI / 180;
            else if ( sico.DELTA2 == 1 ) pdelta2[i] = sico.UNDEF;
            if ( sico.DELTA3 == 1 && pdelta3[i] != sico.UNDEF ) pdelta3[i] = pdelta3[i] * sico.PI / 180;
            else if ( sico.DELTA3 == 1 ) pdelta3[i] = sico.UNDEF;
            if ( sico.ADDFRI == 1 && paddfri[i] != sico.UNDEF ) paddfri[i] = paddfri[i];
            else if ( sico.ADDFRI == 1 ) paddfri[i] = sico.UNDEF;
            if ( sico.ADDFRI2 == 1 && paddfri2[i] != sico.UNDEF ) paddfri2[i] = paddfri2[i];
            else if ( sico.ADDFRI2 == 1 ) paddfri2[i] = sico.UNDEF;
            if ( sico.ADDFRI3 == 1 && paddfri3[i] != sico.UNDEF ) paddfri3[i] = paddfri3[i];
            else if ( sico.ADDFRI3 == 1 ) paddfri3[i] = sico.UNDEF;
            if ( sico.COHESION == 1 && pcohesion[i] != sico.UNDEF ) pcohesion[i] = pcohesion[i];
            else if ( sico.COHESION == 1 ) pcohesion[i] = sico.UNDEF;
            if ( sico.COHESION2 == 1 && pcohesion2[i] != sico.UNDEF ) pcohesion2[i] = pcohesion2[i];
            else if ( sico.COHESION2 == 1 ) pcohesion2[i] = sico.UNDEF;
            if ( sico.COHESION3 == 1 && pcohesion3[i] != sico.UNDEF ) pcohesion3[i] = pcohesion3[i];
            else if ( sico.COHESION3 == 1 ) pcohesion3[i] = sico.UNDEF;
            if ( sico.NY == 1 && pny[i] != sico.UNDEF ) pny[i] = pow( 10, pny[i] );
            else if ( sico.NY == 1 ) pny[i] = sico.UNDEF;
            if ( sico.NY2 == 1 && pny2[i] != sico.UNDEF ) pny2[i] = pow( 10, pny2[i] );
            else if ( sico.NY2 == 1 ) pny2[i] = sico.UNDEF;
            if ( sico.NY3 == 1 && pny3[i] != sico.UNDEF ) pny3[i] = pow( 10, pny3[i] );
            else if ( sico.NY3 == 1 ) pny3[i] = sico.UNDEF;
            if ( sico.CDEFORM == 1 && pcdeform[i] != sico.UNDEF ) pcdeform[i] = pcdeform[i];
            else if ( sico.CDEFORM == 1 ) pcdeform[i] = sico.UNDEF;
            if ( sico.AMBDRAG == 1 && pambdrag[i] != sico.UNDEF ) pambdrag[i] = pambdrag[i];
            else if ( sico.AMBDRAG == 1 ) pambdrag[i] = sico.UNDEF;
            if ( sico.TRANS12 == 1 && ptrans12[i] != sico.UNDEF ) ptrans12[i] = ptrans12[i];
            else if ( sico.TRANS12 == 1 ) ptrans12[i] = sico.UNDEF;
            if ( sico.TRANS13 == 1 && ptrans13[i] != sico.UNDEF ) ptrans13[i] = ptrans13[i];
            else if ( sico.TRANS13 == 1 ) ptrans13[i] = sico.UNDEF;
            if ( sico.TRANS23 == 1 && ptrans23[i] != sico.UNDEF ) ptrans23[i] = ptrans23[i];
            else if ( sico.TRANS23 == 1 ) ptrans23[i] = sico.UNDEF;
            if ( sico.ZFRAG == 1 && pzfrag[i] != sico.UNDEF ) pzfrag[i] = pzfrag[i];
            else if ( sico.ZFRAG == 1 ) pzfrag[i] = sico.UNDEF;
            if ( sico.TRELEASE == 1 && ptrelease[i] != sico.UNDEF ) ptrelease[i] = ptrelease[i];
            else if ( sico.TRELEASE == 1 ) ptrelease[i] = sico.UNDEF;
            if ( sico.TRELSTOP == 1 && ptrelstop[i] != sico.UNDEF ) ptrelstop[i] = ptrelstop[i];
            else if ( sico.TRELSTOP == 1 ) ptrelstop[i] = sico.UNDEF;
            if ( sico.STOPTIME == 1 && pstoptime[i] != sico.UNDEF ) pstoptime[i] = pstoptime[i];
            else pstoptime[i] = 5;
            if ( sico.TSLIDE == 1 && ptslide[i] != sico.UNDEF ) ptslide[i] = ptslide[i];
            else ptslide[i] = sico.UNDEF;
            if ( sico.IMPACTAREA == 1 && pimpactarea[i] != sico.UNDEF ) pimpactarea[i] = pimpactarea[i];
            else if ( sico.IMPACTAREA == 1 ) pimpactarea[i] = sico.UNDEF;
            if ( sico.HDEPOSIT == 1 && phdeposit[i] != sico.UNDEF ) phdeposit[i] = phdeposit[i];
            else if ( sico.HDEPOSIT == 1 ) phdeposit[i] = sico.UNDEF;
            
            if ( sico.PBG == 1 && ppbg1[i] != sico.UNDEF ) ppbg1[i] = ppbg1[i];
            else if ( sico.PBG == 1 ) ppbg1[i] = sico.UNDEF;
            if ( sico.PBG == 1 && ppbg2[i] != sico.UNDEF ) ppbg2[i] = ppbg2[i];
            else if ( sico.PBG == 1 ) ppbg2[i] = sico.UNDEF;
            if ( sico.PBG == 1 && ppbg3[i] != sico.UNDEF ) ppbg3[i] = ppbg3[i];
            else if ( sico.PBG == 1 ) ppbg3[i] = sico.UNDEF;
            
        }
        
        if ( sico.PBG == 0 ) {
               
            ppbg1 = (int*) calloc( sico.IMAX, sizeof(int));
            ppbg2 = (int*) calloc( sico.IMAX, sizeof(int));
            ppbg3 = (int*) calloc( sico.IMAX, sizeof(int));
        }


    #endif


// -- STOP --- Reading ascii rasters and pre-processing arrays (if GRASS is not used) ---------------------------


// -- START -- Computing maximum release heights and release volumes --------------------------------------------


    hflow_max0 = 0; // initializing maximum PHASE 1 release height
    vol_flow0 = 0; // initializing PHASE 1 release volume
    
    if ( sico.MODEL == 7 ) { 
    
        hflow_max02 = 0; hflow_max03 = 0; // initializing maximum PHASE 2 and PHASE 3 release heights
        vol_flow02 = 0; vol_flow03 = 0; // initializing PHASE 2 and PHASE 3 release volumes
    }

    for ( i=0; i<sico.IMAX; i++ ) { // updating release height maxima and volumes:

        if ( sico.MODEL == 7 && sico.RELM3 == 0 && rhrelease1 < 1 ) {
        
            phrelease[i] = phrelease[i] * rhrelease1;
            phrelease2[i] = 0;
            phrelease3[i] = phrelease[i] * ( 1 - rhrelease1 ) / rhrelease1;
        }

        if ( sico.MODEL == 7 && sico.ENTR3 == 0 && rhentrmax1 < 1 ) {
        
            phentrmax[i] = phentrmax[i] * rhentrmax1;
            phentrmax2[i] = 0;
            phentrmax3[i] = phentrmax[i] * ( 1 - rhentrmax1 ) / rhentrmax1;        
        }

        if ( sico.RELM == 1 ) {
            phrelease[i] *= vhrelease;
            if ( phrelease[i] > hflow_max0 ) hflow_max0 = phrelease[i]; // updating maximum PHASE 1 release height
            vol_flow0 += fmax(0, phrelease[i]) * pow ( sico.CSZ, 2 ); // updating PHASE 1 release volume
        }
        if ( sico.MODEL == 7 && sico.RELM2 == 1 ) {
            phrelease2[i] *= vhrelease;
            if ( phrelease2[i] > hflow_max02 ) hflow_max02 = phrelease2[i]; // updating maximum PHASE 2 release height
            vol_flow02 += fmax(0, phrelease2[i]) * pow ( sico.CSZ, 2 ); // updating PHASE 2 release volume
        }
        if ( sico.MODEL == 7 && ( sico.RELM3 == 1 || rhrelease1 < 1 )) {
            phrelease3[i] *= vhrelease;
            if ( phrelease3[i] > hflow_max03 ) hflow_max03 = phrelease3[i]; // updating maximum PHASE 3 release height
            vol_flow03 += fmax(0, phrelease3[i]) * pow ( sico.CSZ, 2 ); // updating PHASE 3 release volume
        }
        
        if ( sico.ENTR == 1 ) phentrmax[i] *= vhentrmax;
        if ( sico.MODEL == 7 && sico.ENTR2 == 1 ) phentrmax2[i] *= vhentrmax;
        if ( sico.MODEL == 7 && ( sico.ENTR3 == 1 || rhentrmax1 < 1 )) phentrmax3[i] *= vhentrmax;
    }
    
    
// -- STOP --- Computing maximum release heights and release volumes --------------------------------------------


// -- START -- Definition of computational domains --------------------------------------------------------------


    // *** cdomain[i] = 0: 1st row edge cells
    // *** cdomain[i] = 1: all other cells
    // *** cdomain2[i] = 1: 1st, 2nd or 3rd row edge cells
    // *** cdomain2[i] = 0: all other cells
    // *** Further refinement and dynamic adaptation of cdomain[i] at the beginning of each time step

    ib[0] = 0;

    for ( i=0; i<sico.IMAX; i++ ) {

        iin = fin( i, pelev, sico ); // cell neighbourhood
        for ( j=0; j<9; j++ ) in[i][j] = iin[j];
        ctrlv = 1; ctrlvv = 1; ctrlvvv = 1; // resetting controls
        for ( j=1; j<9; j++ ) { // loop over all neighbour cells

            if ( in[i][j] < 0 || in[i][j] >= sico.IMAX ) ctrlv = 0; // 1st row edge cells

            inn = fin( in[i][j], pelev, sico ); // neighbourhood of neighbour cell
            for ( jj=1; jj<9; jj++ ) { // loop over neighbourhood
                if ( inn[jj] < 0 || inn[jj] >= sico.IMAX ) ctrlvv = 0; // 2nd row edge cells

                innn = fin( inn[jj], pelev, sico ); // neighbourhood of neighbour cell
                for ( jjj=1; jjj<9; jjj++ ) // loop over neighbourhood
                    if ( innn[jjj] < 0 || innn[jjj] >= sico.IMAX ) ctrlvvv = 0; // 3rd row edge cells

                free( innn );
            }
            free( inn );
        }
        free( iin );

        icheck[i][0] = 0;
        icheck[i][1] = 0;

        if ( ctrlv == 0 ) { cdomain[i] = 0; cdomain2[i] = 1; } // setting domain controls
        else if ( ctrlvv == 0 || ctrlvvv == 0 ) { cdomain[i] = 1; cdomain2[i] = 1; } //!!!CHECK 3rd row edge cells
        else { cdomain[i] = 1; cdomain2[i] = 0; }

        if ( cdomain2[i] == 0 ) {

            ibasket[0][ib[0]] = i;
            ib[0] += 1;
        }
    }
    
    nzones = 0; // zones for material budgets
    
    for ( i=0; i<sico.IMAX; i++ ) {
    
        if ( sico.ZONES != 1 ) pzones[i] = 0;
        if ( pzones[i] > nzones ) nzones = pzones[i];
    }
    
    nzones += 1;
    
    float vol_zone1[nzones], vol_zone2[nzones], vol_zone3[nzones], vol_czone1[nzones], vol_czone2[nzones], vol_czone3[nzones];


// -- STOP --- Definition of computational domains --------------------------------------------------------------


// -- START -- Computing slopes and topography-following cell sizes ---------------------------------------------


    // *** Slopes are computed from the 3x3 environment of each raster cell
    // *** Cell sizes applied in the simulation depend on the definition of the height-to-depth conversion

    for ( i=0; i<sico.IMAX; i++ ) { // loop over all cells:

        qelev[i] = pelev[i];
        relev[i] = 0;

        if ( cdomain[i] != 0 ) { // if cell is no 1st row edge cell:

            betax[i] = fbeta( pelev[in[i][5]], pelev[in[i][2]], 2.0, sico ); // slopes
            betay[i] = fbeta( pelev[in[i][4]], pelev[in[i][1]], 2.0, sico );
            betaxy[i] = fbetaxy( betax[i], betay[i] );

        } else { betax[i] = 0; betay[i] = 0; betaxy[i] = 0; }

        dx[i] = sico.CSZ;
        dy[i] = sico.CSZ;
    }


// -- STOP --- Computing slopes and topography-following cell sizes ---------------------------------------------


// -- START -- Converting heights into depths and deriving release momenta --------------------------------------


    for ( i=0; i<sico.IMAX; i++ ) {

        sico.COLLAPSE = 0;
        if ( sico.TRELSTOP == 1 ) { if ( ptrelstop[i] < 0 && ptrelstop[i] != sico.UNDEF ) { sico.COLLAPSE = 1; ptrelstop[i] *= -1; }} // checking whether to use collapse mode

        ctrl_trelease = 0;
        if ( sico.TRELEASE == 1 ) { if ( ptrelease[i] > 0 ) ctrl_trelease = 1; }

        if ( sico.TOPOFACT == 0 ) topofact = 1.0;
        else topofact = cos( betaxy[i] );

        if ( sico.MODEL <= 3 ) {

            if ( ctrl_trelease == 0 && sico.COLLAPSE == 0 ) {

                if ( phrelease[i] > 0 ) aw[i][0] = fmax(0, phrelease[i]) / topofact; // release depth !!!CHECK
                aw[i][1] = aw[i][0] * pviny[i] * -1 * topofact; // release momenta in x and y directions
                aw[i][2] = aw[i][0] * pvinx[i] * topofact;
                
            } else { aw[i][0] = 0; aw[i][1] = 0; aw[i][2] = 0; }

            aw[i][3] = sflow.TEMP0; // temperature
            
        } else if ( sico.MODEL == 7 ) {

            if ( ctrl_trelease == 0 && sico.COLLAPSE == 0 ) {

                aw[i][0] = fmax(0, phrelease[i]) / topofact; // PHASE 1 release depth
                aw[i][3] = fmax(0, phrelease2[i]) / topofact; // PHASE 2 release depth
                aw[i][6] = fmax(0, phrelease3[i]) / topofact; // PHASE 3 release depth

                aw[i][1] = aw[i][0] * pviny[i] * -1 * topofact; // PHASE 1 release momenta in x and y directions
                aw[i][2] = aw[i][0] * pvinx[i] * topofact;
                aw[i][4] = aw[i][3] * pviny2[i] * -1 * topofact; // PHASE 2 release momenta in x and y directions
                aw[i][5] = aw[i][3] * pvinx2[i] * topofact;
                aw[i][7] = aw[i][6] * pviny3[i] * -1 * topofact; // PHASE 3 release momenta in x and y directions
                aw[i][8] = aw[i][6] * pvinx3[i] * topofact;

            } else { aw[i][0] = 0; aw[i][1] = 0; aw[i][2] = 0; aw[i][3] = 0; aw[i][4] = 0; aw[i][5] = 0; aw[i][6] = 0; aw[i][7] = 0; aw[i][8] = 0; }
            
            aw[i][9] = sflow.TEMP0; // temperature
        }
    }


// -- STOP --- Converting heights into depths and deriving release momenta --------------------------------------


// -- START -- Computing slopes and topography-following cell sizes ---------------------------------------------


    // *** Slopes are computed from the 3x3 environment of each raster cell
    // *** Cell sizes applied in the simulation depend on the definition of the height-to-depth conversion

    for ( i=0; i<sico.IMAX; i++ ) { // loop over all cells:

        if ( cdomain[i] != 0 ) { // if cell is no 1st row edge cell:

            cin = fnoosc( aw, in, i, sico );

            betax[i] = fbeta( pelev[cin[0]], pelev[cin[1]], (float)(cin[4]), sico ); // slopes
            betay[i] = fbeta( pelev[cin[2]], pelev[cin[3]], (float)(cin[5]), sico );
            betaxh[i] = fbeta( pelev[cin[0]]+aw[cin[0]][0], pelev[cin[1]]+aw[cin[1]][0], (float)(cin[4]), sico ); // slopes including flow heights
            betayh[i] = fbeta( pelev[cin[2]]+aw[cin[2]][0], pelev[cin[3]]+aw[cin[3]][0], (float)(cin[5]), sico );
            betaxy[i] = fbetaxy( betax[i], betay[i] );

            if ( sico.LAYERS0 > 0 ) {
                            
                betax2[i] = fbeta( pelev[cin[0]]+aw[cin[0]][0], pelev[cin[1]]+aw[cin[1]][0], (float)(cin[4]), sico );
                betay2[i] = fbeta( pelev[cin[2]]+aw[cin[2]][0], pelev[cin[3]]+aw[cin[3]][0], (float)(cin[5]), sico );

                betaxy2[i] = fbetaxy( betax2[i], betay2[i] );
                                                               
                betaxh2[i] = fbeta( pelev[cin[0]]+aw[cin[0]][0]+aw[cin[0]][3], pelev[cin[1]]+aw[cin[1]][0]+aw[cin[1]][3], (float)(cin[4]), sico );
                betayh2[i] = fbeta( pelev[cin[2]]+aw[cin[2]][0]+aw[cin[2]][3], pelev[cin[3]]+aw[cin[3]][0]+aw[cin[3]][3], (float)(cin[5]), sico );
                                    
                betax3[i] = fbeta( pelev[cin[0]]+aw[cin[0]][0]+aw[cin[0]][3], pelev[cin[1]]+aw[cin[1]][0]+aw[cin[1]][3], (float)(cin[4]), sico );
                betay3[i] = fbeta( pelev[cin[2]]+aw[cin[2]][0]+aw[cin[2]][3], pelev[cin[3]]+aw[cin[3]][0]+aw[cin[3]][3], (float)(cin[5]), sico );
                
                betaxy3[i] = fbetaxy( betax3[i], betay3[i] );
            }

            if ( sico.LAYERS0 > 0 || sico.SURFACE > 1 ) {
                                     
                betaxh3[i] = fbeta( pelev[cin[0]]+aw[cin[0]][0]+aw[cin[0]][3]+aw[cin[0]][6], pelev[cin[1]]+aw[cin[1]][0]+aw[cin[1]][3]+aw[cin[1]][6], (float)(cin[4]), sico );
                betayh3[i] = fbeta( pelev[cin[2]]+aw[cin[2]][0]+aw[cin[2]][3]+aw[cin[2]][6], pelev[cin[3]]+aw[cin[3]][0]+aw[cin[3]][3]+aw[cin[3]][6], (float)(cin[5]), sico );
            }
     
            free( cin );

        } else { betax[i] = 0; betay[i] = 0; betaxy[i] = 0; }
    }

    if ( sico.SURFACE > 1 ) {
     
        for ( i=0; i<sico.IMAX; i++ ) {
              
            cplain[i] = 1; // control for plain area
                
            for ( l=1; l<9; l++ ) {

                if ( sico.MODEL <= 3 ) hflow = aw[in[i][0]][0];
                else if ( sico.MODEL == 7 ) hflow = aw[in[i][0]][0] + aw[in[i][0]][3] + aw[in[i][0]][6];
            
                if ( fabs(betaxh3[in[i][0]]) > sflow.BETAPLAIN || fabs(betayh3[in[i][0]]) > sflow.BETAPLAIN || hflow <= sico.HFLOWMIN ) cplain[i] = 0;
            }
        }
    }
    

// -- STOP --- Computing slopes and topography-following cell sizes ---------------------------------------------


// -- START -- Initializing vector elements, volumes, and maximum depths ----------------------------------------


    hflow_max = 0; // initializing maximum mixture or PHASE 1 flow depth
    hflow_maxmax = 0; // initializing absolute maximum flow depth
    vol_flow = 0; // initializing mixture or PHASE 1 flow volume
    vol_entr = 0; // initializing mixture or PHASE 1 entrained or deposited volume
    
    for ( z=0; z<nzones; z++ ) vol_zone1[z] = 0; // zone-specific volumes

    if ( sico.MODEL == 7 ) {
    
        hflow_max2 = 0; hflow_max3 = 0; // initializing maximum PHASE 2 and PHASE 3 flow depths
        vol_flow2 = 0; vol_flow3 = 0;  // initializing PHASE 2 and PHASE 3 flow volumes
        vol_entr2 = 0; vol_entr3 = 0; // initializing PHASE 2 and PHASE 3 entrained or deposited volumes
        
        for ( z=0; z<nzones; z++ ) { vol_zone2[z] = 0; vol_zone3[z] = 0; } // zone-specific volumes
    }

    for ( i=0; i<sico.IMAX; i++ ) {

        if ( sico.PBG == 1 && paradistort >= 1 ) { 
        
            ppbgt1[i] = 0; ppbgt2[i] = 0; ppbgt3[i] = 0; ppbgx1[i] = ppbg1[i]; ppbgx2[i] = ppbg2[i]; ppbgx3[i] = ppbg3[i]; anctr[i] = 0; anpx[i] = px[i]; anpy[i] = py[i];
        }

        if ( sico.TOPOFACT == 0 ) topofact = 1.0;
        else topofact = cos( betaxy[i] );

        carea = pow ( sico.CSZ, 2 ); // area of cell

        awt[i][0] = aw[i][0] * topofact;

        if ( aw[i][0] * topofact > hflow_max ) hflow_max = aw[i][0] * topofact; // updating maximum mixture or PHASE 1 release depth
        vol_flow += aw[i][0] * carea; // updating mixture or PHASE 1 release volume
        vol_zone1[pzones[i]] += aw[i][0] * carea;

        if ( sico.MODEL <= 3 ) { // one-phase model

            for ( k=4; k<8; k++ ) awt[i][k] = 0;
            if ( sico.RELM == 1 ) awt[i][8] = aw[i][0] * topofact; else awt[i][8] = 0;
            for ( k=9; k<12; k++ ) awt[i][k] = 0;
            
            if ( hflow_max > hflow_maxmax ) hflow_maxmax = hflow_max; // updating absolute maximum flow depth
            
        }  else if ( sico.MODEL == 7 ) { // multi-phase model

            awt[i][3] = aw[i][3] * topofact;
            awt[i][6] = aw[i][6] * topofact;

            for ( k=10; k<15; k++ ) awt[i][k] = 0;
            awt[i][16] = ( aw[i][0] + aw[i][3] + aw[i][6] ) * topofact;
            for ( k=17; k<25; k++ ) awt[i][k] = 0;
            if ( sico.RELM == 1 ) awt[i][26] = aw[i][0] * topofact; else awt[i][26] = 0;
            awt[i][27] = 0;
            if ( sico.RELM2 == 1 ) awt[i][28] = aw[i][3] * topofact; else awt[i][28] = 0;
            awt[i][29] = 0;
            if ( sico.RELM3 == 1 || rhrelease1 < 1 ) awt[i][30] = aw[i][6] * topofact; else awt[i][30] = 0;
            awt[i][31] = 0;
            awt[i][32] = ( aw[i][0] + aw[i][3] + aw [i][6] ) * topofact;
            for ( k=33; k<41; k++ ) awt[i][k] = 0;

            if ( aw[i][3] * topofact > hflow_max2 ) hflow_max2 = aw[i][3] * topofact; // updating maximum PHASE 2 release depth
            if ( aw[i][6] * topofact > hflow_max3 ) hflow_max3 = aw[i][6] * topofact; // updating maximum PHASE 3 release depth
            if ( hflow_max + hflow_max2 + hflow_max3 > hflow_maxmax ) hflow_maxmax = hflow_max + hflow_max2 + hflow_max3; // updating absolute maximum flow depth
            
            vol_flow2 += aw[i][3] * carea; // updating PHASE 2 release volume
            vol_flow3 += aw[i][6] * carea; // updating PHASE 3 release volume
            vol_zone2[pzones[i]] += aw[i][3] * carea;
            vol_zone3[pzones[i]] += aw[i][6] * carea;
        }

        for ( l=-8; l<=-3; l++ ) awt[i][nvect_all+l] = sico.UNDEF; // velocities, phase fractions, time of reach
        awt[i][nvect_all-2] = sico.UNDEF * -100;
        awt[i][nvect_all-1] = sico.UNDEF;

        cstopped[i] = 0; // control for stopping
    }

    for ( i=0; i<sico.IMAX; i++ ) {
    
        if ( cdomain2[i] == 0 ) {

            asigma_xelev[i] = fsigma( pelev[i], pelev[in[i][5]], pelev[in[i][2]], sico.CSZ, sico ); // gradients of elevation
            asigma_yelev[i] = fsigma( pelev[i], pelev[in[i][4]], pelev[in[i][1]], sico.CSZ, sico );
        } else {
        
            asigma_xelev[i] = 0; asigma_yelev[i] = 0;
        }
    }

    fj[0][0]=1; fj[0][1]=1; fj[1][0]=1; fj[1][1]=-1; fj[2][0]=-1; fj[2][1]=1; fj[3][0]=-1; fj[3][1]=-1; // factors for positive or negative gradients

    for ( i=0; i<sico.IMAX; i++ ) {

        for ( j=0; j<4; j++ ) {
        
            if ( cdomain2[i] == 0 ) wintelev[i][j]=pelev[in[i][j]]+fj[j][0]*0.25*dx[i]*asigma_xelev[in[i][j]]+fj[j][1]*0.25*dy[i]*asigma_yelev[in[i][j]];
            else wintelev[i][j] = pelev[i];
        }
    }

    if ( sico.PROFILE > 0 ) {

        for ( profctrl = 0; profctrl < sico.M+sico.N; profctrl++ ) {
            for ( profi = 0; profi < (int)(tmax/tout+3); profi++ ) {     
                for ( profj = 0; profj < profmax; profj++ ) profdata[profctrl][profi][profj] = 0;
            }
        }

        profctrl = 0;
        profnx[profctrl] = profx[0];
        profny[profctrl] = profy[0];
        profdiffabs[0] = 0;

        for ( profi=1; profi<sico.PROFILE; profi++ ) {
                           
            profdiffx = profx[profi] - profx[profi-1];
            profdiffy = profy[profi] - profy[profi-1];
            profdiff = pow( pow( profdiffx, 2 ) + pow( profdiffy, 2 ), 0.5 );
            profn = (int)( profdiff + 0.5 );
                               
            profdiffnx = (float)profdiffx / (float)profn;
            profdiffny = (float)profdiffy / (float)profn;
                               
            for ( profj = 0; profj < profn; profj++ ) {
                               
                profctrl += 1;
                profnx[profctrl] = profnx[profctrl-1] + profdiffnx;
                profny[profctrl] = profny[profctrl-1] + profdiffny;
                profdiffabs[profctrl] = profdiffabs[profctrl-1] + sico.CSZ * pow( pow( profdiffnx, 2 ) + pow( profdiffny, 2 ), 0.5 );
            }
        }
            
        for ( profi = 0; profi <= profctrl; profi++ ) {

            profwhtx1 = 1 - ( profnx[profi] - (int)profnx[profi] );
            profwhtx2 = 1 - ((int)profnx[profi] + 1 - profnx[profi] );
            profwhty1 = 1 - ( profny[profi] - (int)profny[profi] );
            profwhty2 = 1 - ((int)profny[profi] + 1 - profny[profi] );
              
            profelev[profi] = pelev[icor[(int)profnx[profi]][(int)profny[profi]]] * profwhtx1 * profwhty1
                + pelev[icor[(int)profnx[profi]][(int)profny[profi]+1]] * profwhtx1 * profwhty2
                + pelev[icor[(int)profnx[profi]+1][(int)profny[profi]]] * profwhtx2 * profwhty1
                + pelev[icor[(int)profnx[profi]+1][(int)profny[profi]+1]] * profwhtx2 * profwhty2;

            if ( sico.HDEPOSIT == 1 ) profdeposit[profi] = phdeposit[icor[(int)profnx[profi]][(int)profny[profi]]] * profwhtx1 * profwhty1
                + phdeposit[icor[(int)profnx[profi]][(int)profny[profi]+1]] * profwhtx1 * profwhty2
                + phdeposit[icor[(int)profnx[profi]+1][(int)profny[profi]]] * profwhtx2 * profwhty1
                + phdeposit[icor[(int)profnx[profi]+1][(int)profny[profi]+1]] * profwhtx2 * profwhty2;
                
            for ( profl = 0; profl < profmax; profl++ ) {

                if ( sico.MODEL <= 3 ) {
                     
                    if ( profl == 0 ) profk = 0; else if ( profl == 1 ) profk = 5; else if ( profl == 2 ) profk = 6; else if ( profl == 3 ) profk = 7; else profk = 4;
                         
                } else if ( sico.MODEL == 7 ) { 

                    if ( profl == 0 ) profk = 0; else if ( profl == 1 ) profk = 3; else if ( profl == 2 ) profk = 6; else if ( profl == 3 ) profk = 13; else if ( profl == 4 ) profk = 14;
                    else if ( profl == 5 ) profk = 15; else if ( profl == 6 ) profk = 17; else if ( profl == 7 ) profk = 18; else if ( profl == 8 ) profk = 19; 
                    else if ( profl == 9 ) profk = 21; else if ( profl == 10 ) profk = 22; else if ( profl == 11 ) profk = 23; else if ( profl == 12 ) profk = 10; 
                    else if ( profl == 13 ) profk = 11; else profk = 12;
                }
                
                profdata[profi][0][profl] = awt[icor[(int)profnx[profi]][(int)profny[profi]]][profk] * profwhtx1 * profwhty1
                    + awt[icor[(int)profnx[profi]][(int)profny[profi]+1]][profk] * profwhtx1 * profwhty2
                    + awt[icor[(int)profnx[profi]+1][(int)profny[profi]]][profk] * profwhtx2 * profwhty1
                    + awt[icor[(int)profnx[profi]+1][(int)profny[profi]+1]][profk] * profwhtx2 * profwhty2;
            }
        }
    }

    if ( sico.CTRLPOINTS > 0 ) {

        for ( j=0; j<sico.CTRLPOINTS; j++ ) {

            if ( sico.TOPOFACT == 0 ) topofact = 1.0;
            else topofact = cos( betaxy[icor[ctrlpx[j]][ctrlpy[j]]] );

            if ( sico.MODEL <= 3 ) {

                ctrlpdata[j][0][0] = awt[icor[ctrlpx[j]][ctrlpy[j]]][0] * topofact;
                ctrlpdata[j][0][1] = awt[icor[ctrlpx[j]][ctrlpy[j]]][4];
                    
            } else if ( sico.MODEL == 7 ) {

                ctrlpdata[j][0][0] = awt[icor[ctrlpx[j]][ctrlpy[j]]][0] * topofact;
                ctrlpdata[j][0][1] = awt[icor[ctrlpx[j]][ctrlpy[j]]][3] * topofact;
                ctrlpdata[j][0][2] = awt[icor[ctrlpx[j]][ctrlpy[j]]][6] * topofact; 
                ctrlpdata[j][0][3] = awt[icor[ctrlpx[j]][ctrlpy[j]]][10];
                ctrlpdata[j][0][4] = awt[icor[ctrlpx[j]][ctrlpy[j]]][11];
                ctrlpdata[j][0][5] = awt[icor[ctrlpx[j]][ctrlpy[j]]][12];
            }
        }
    }


// -- STOP --- Initializing vector elements, volumes, and maximum depths ----------------------------------------


// -- START -- Definition of hydrograph profiles ----------------------------------------------------------------


    // *** Profile parameters needed for adding material through input hydrographs and for quantifying material at output hydrographs at each time step (if applicable)
    // *** Centre and end points of each hydrograph profile are written to text files

    if ( hydrograph == 1 ) { // if at least one hydrograph is provided
    
        for ( i=0; i<sico.IMAX; i++ ) { 

            for ( hydj = 0; hydj < ( hydnin + hydnout ); hydj++ ) {
                if ( px[i] == hydx[hydj] && py[i] == hydy[hydj] ) {

                    hydi[hydj] = i; // coordinates and elevation of centres of hydrographs
                    hydelev[hydj] = pelev[i];
                }
            }
            qtinit[hydj] = sico.UNDEF; // initializing start time of hydrograph discharge
        }

        hydp = alloc_imatrix( (int)(2*hydlmax/sico.CSZ+2), hydnin+hydnout ); // allocating memory to hydrograph profiles

        for ( hydj = 0; hydj < ( hydnin + hydnout ); hydj++ ) { // loop over all hydrographs

            if ( hydalpha[hydj] < 0 ) hydalpha[hydj] = falpha ( betax[hydi[hydj]], betay[hydi[hydj]], sico ) + sico.PI * 0.5;
                // cross-slope direction
            if ( hydalpha[hydj] > 2 * sico.PI ) hydalpha[hydj] -= 2 * sico.PI;
            hydp0 = fhydp ( hydx[hydj], hydy[hydj], fabs(hydl[hydj]), px, py, hydalpha[hydj], hydlmax, sico );
            for ( i=0; i<(int)(2*hydlmax/sico.CSZ+2); i++ ) hydp[i][hydj] = hydp0[i]; // hydrograph profile
            free( hydp0 );

            if ( hydl[hydj] > 0 ) { // Adjusting central point according to lowest elevation along hydrograph profile

                pelevhydtest = sico.UNDEF * -1;

                for ( hydk=1; hydk<= hydp[0][hydj]; hydk++ ) { // loop over all cells of hydrograph profile

                    if ( pelev[hydp[hydk][hydj]] < pelevhydtest ) {

                        pelevhydtest = pelev[hydp[hydk][hydj]];
                        hydi[hydj] = hydp[hydk][hydj];
                        hydx[hydj] = px[hydp[hydk][hydj]];
                        hydy[hydj] = py[hydp[hydk][hydj]];
                    }
                }

                hydelev[hydj] = pelevhydtest;
            }

            hydx_metric[hydj] = ( float ) hydy[hydj] * sico.CSZ + sico.BDWEST;
            hydy_metric[hydj] = sico.BDNORTH - ( float ) hydx[hydj] * sico.CSZ;

            hydout_xmin_metric = ( float ) py[hydp[1][hydj]] * sico.CSZ + sico.BDWEST; // external coordinates of profile terminal points
            hydout_xmax_metric = ( float ) py[hydp[hydp[0][hydj]][hydj]] * sico.CSZ + sico.BDWEST;
            hydout_ymin_metric = sico.BDNORTH - ( float ) px[hydp[1][hydj]] * sico.CSZ;
            hydout_ymax_metric = sico.BDNORTH - ( float ) px[hydp[hydp[0][hydj]][hydj]] * sico.CSZ;

            if ( hydj < hydnin ) fprintf(f_hydout, "I%i\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", hydj+1,
                hydout_xmin_metric, hydx_metric[hydj], hydout_xmax_metric, hydout_ymin_metric, hydy_metric[hydj], hydout_ymax_metric);
            else fprintf(f_hydout, "O%i\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", hydj-hydnin+1,
                hydout_xmin_metric, hydx_metric[hydj], hydout_xmax_metric, hydout_ymin_metric, hydy_metric[hydj], hydout_ymax_metric);
                // writing hydrograph profile terminal points and centre to text file

        }
    }


// -- STOP --- Definition of hydrograph profiles ----------------------------------------------------------------


// -- START -- Definition / initialization of control variables -------------------------------------------------


    tlength0 = sico.CFL[1]; // length of first time step
    tlengthpre = sico.CFL[1]; // length of previous time step
    tsum = 0; // total time
    tint = 0; // time since last output
    cflmax = 0; // maximum CFL value
    nsum = 1; // number of time steps
    nout = 1; // number of output time steps
    ccontinue = 1; // control whether simulation should be continued (1 or 0)
    csuccess = 1; // control for success of simulation (1 or 0)
    //coverflow = 0; // control for overflow i.e. part of the flow leaving the study area (1 or 0), deactivated
    if ( sico.MODEL <= 3 ) qh_test = hflow_max0; else qh_test = hflow_max0 + hflow_max02 + hflow_max03;  // control release height for collapse
    cflowpre = 0; // control for activation of flow model at previous time step
    vol_hyd = 0; vol_hyd2 = 0; vol_hyd3 = 0; ctrl_hydout = 1; // hydrograph volumes and control for hydrograph output
    hekin_max = 0; // initialization of maximum value of flow kinetic energy
    ctr_ansteps = 1; // counter for output time steps for refesthing dynamic evolution of orthophoto
    ansteps = (int)( fmax( 1, 100 / tout )); // interval of output time steps for refreshing dynamic evolution of orthophoto
    highpass = 0.5; // colour range factor for high-pass filter


// -- STOP --- Definition / initialization of control variables -------------------------------------------------


// -- START -- Summary of release and initial state of the flow -------------------------------------------------


    // *** Information is written to the display (immediately) and to the summary file

    if ( sico.MODEL <= 3 ) { // one-phase models

        if ( hflow_max0 >= 1000 ) prec_hflow = 0; // precision for summary of maximum flow heights and depths
        else if ( hflow_max0 >= 100 ) prec_hflow = 1;
        else if ( hflow_max0 >= 10 ) prec_hflow = 2;
        else if ( hflow_max0 >= 1 ) prec_hflow = 3;
        else if ( hflow_max0 >= 0.1 ) prec_hflow = 4;
        else if ( hflow_max0 == 0 ) prec_hflow = 2;
        else prec_hflow = 5;

        if ( vol_flow0 >= 100000000 ) prec_vol = 0; // precision for summary of flow volumes
        else if ( vol_flow0 >= 1000000 ) prec_vol = 1;
        else if ( vol_flow0 >= 100000 ) prec_vol = 2;
        else if ( vol_flow0 >= 10000 ) prec_vol = 3;
        else if ( vol_flow0 >= 1000 ) prec_vol = 4;
        else if ( vol_flow0 == 0 ) prec_vol = 3;
        else prec_vol = 5;

        prec_ekin = prec_vol-3; if ( prec_ekin < 0 ) prec_ekin = 0; // precision for summary of flow kinetic energies

        if ( sico.MULT == 0 ) {

            printf("   nout\tnsum\tcfl\ttlength\ttsum\tdmax\tvmax\tvolume\tekin\n");
            printf("   R\t-----\t-----\t-----\t-----\t*%.*f\t-----\t*%.*f\t-----\n",
                prec_hflow, hflow_max0, prec_vol, vol_flow0/1000 ); // displaying release parameters

            printf("   0\t0\t-----\t-----\t0.0\t%.*f\t0.00\t%.*f\t-----\n",
                prec_hflow, hflow_max, prec_vol, vol_flow/1000 ); // displaying initial parameters

            fflush(stdout); // forcing immediate display
        }

        fprintf(f_summary, "nout\tnsum\tcfl\ttlength\ttsum\tdmax\tvmax\tvolume\tekin\n");
        fprintf(f_summary, "R\t-----\t-----\t-----\t-----\t*%.*f\t-----\t*%.*f\t-----\n", prec_hflow, hflow_max0, prec_vol, vol_flow0 ); // writing release parameters to summary file
        fprintf(f_summary, "0\t0\t-----\t-----\t0.0\t%.*f\t0.00\t%.*f\t\n", prec_hflow, hflow_max, prec_vol, vol_flow ); // writing initial parameters to summary file
        
    } else if ( sico.MODEL == 7 ) { // multi-phase model

        if ( hflow_max0 >= 1000 || hflow_max02 >= 1000 || hflow_max03 >= 1000 ) prec_hflow = 0;  // precision for summary of maximum flow heights and depths
        else if ( hflow_max0 >= 100 || hflow_max02 >= 100 || hflow_max03 >= 100 ) prec_hflow = 1;
        else if ( hflow_max0 >= 10 || hflow_max02 >= 10 || hflow_max03 >= 10 ) prec_hflow = 2;
        else if ( hflow_max0 >= 1 || hflow_max02 >= 1 || hflow_max03 >= 1 ) prec_hflow = 3;
        else if ( hflow_max0 >= 0.1 || hflow_max02 >= 0.1 || hflow_max03 >= 0.1 ) prec_hflow = 4;
        else if ( hflow_max0 == 0 && hflow_max02 == 0 && hflow_max03 == 0 ) prec_hflow = 2;
        else prec_hflow = 5;

        if ( vol_flow0 >= 100000000 || vol_flow02 >= 10000000 || vol_flow03 >= 10000000 ) prec_vol = 0; // precision for summary of flow volumes
        else if ( vol_flow0 >= 1000000 || vol_flow02 >= 1000000 || vol_flow03 >= 1000000 ) prec_vol = 1;
        else if ( vol_flow0 >= 100000 || vol_flow02 >= 100000 || vol_flow03 >= 100000 ) prec_vol = 2;
        else if ( vol_flow0 >= 10000 || vol_flow02 >= 10000 || vol_flow03 >= 10000 ) prec_vol = 3;
        else if ( vol_flow0 >= 1000 || vol_flow02 >= 1000  || vol_flow03 >= 1000 ) prec_vol = 4;
        else if ( vol_flow0 == 0 && vol_flow02 == 0 && vol_flow03 == 0 ) prec_vol = 3;
        else prec_vol = 5;

        prec_ekin = prec_vol-3; if ( prec_ekin < 0 ) prec_ekin = 0; // precision for summary of flow kinetic energies

        if ( sico.MULT == 0 ) {

            printf("   nout\tnsum\tcfl\ttlength\ttsum\tdmax1\tvmax1\tdmax2\tvmax2\tdmax3\tvmax3\tvolume1\tvolume2\tvolume3\tekin\n");
            printf("   R\t-----\t-----\t-----\t-----\t*%.*f\t-----\t*%.*f\t-----\t*%.*f\t-----\t*%.*f\t*%.*f\t*%.*f\t-----\n",
                prec_hflow, hflow_max0, prec_hflow, hflow_max02, prec_hflow, hflow_max03, prec_vol, vol_flow0/1000, prec_vol, vol_flow02/1000, prec_vol, vol_flow03/1000);
                    // displaying release parameters
            printf("   0\t0\t-----\t-----\t0.0\t%.*f\t0.00\t%.*f\t0.00\t%.*f\t0.00\t%.*f\t%.*f\t%.*f\t-----\n", prec_hflow, hflow_max, prec_hflow, hflow_max2,
                     prec_hflow, hflow_max3, prec_vol, vol_flow/1000, prec_vol, vol_flow2/1000, prec_vol, vol_flow3/1000); // displaying initial parameters

            fflush(stdout); // forcing immediate display
        }

        fprintf(f_summary, "nout\tnsum\tcfl\ttlength\ttsum\tdmax1\tvmax1\tdmax2\tvmax2\tdmax3\tvmax3\tvolume1\tvolume2\tvolume3\tekin\n");
        fprintf(f_summary, "R\t-----\t-----\t-----\t-----\t*%.*f\t-----\t*%.*f\t-----\t*%.*f\t-----\t*%.*f\t*%.*f\t*%.*f\t-----\n",
            prec_hflow, hflow_max0, prec_hflow, hflow_max02, prec_hflow, hflow_max03, prec_vol, vol_flow0, prec_vol, vol_flow02, prec_vol, vol_flow03 ); // writing release parameters to summary file
        fprintf(f_summary, "0\t0\t-----\t-----\t0.0\t%.*f\t0.00\t%.*f\t0.00\t%.*f\t0.00\t%.*f\t%.*f\t%.*f\t-----\n",
            prec_hflow, hflow_max, prec_hflow, hflow_max2, prec_hflow, hflow_max3, prec_vol, vol_flow, prec_vol, vol_flow2, prec_vol, vol_flow3 ); // writing initial parameters to summary file
    }
    
    fprintf(f_volumes, "nout\tnsum"); // volumes file
    
    for ( z=0; z<nzones; z++ ) {
    
        fprintf(f_volumes, "\tz%i_vflow1\tz%ivchange1", z, z );
        if ( sico.MODEL == 7 ) fprintf(f_volumes, "\tz%i_vflow2\tz%ivchange2\tz%i_vflow3\tz%ivchange3", z, z, z, z );
    }

    fprintf(f_volumes, "\n0\t0"); // volumes file
    
    for ( z=0; z<nzones; z++ ) {
    
        fprintf(f_volumes, "\t%.3f\t0", vol_zone1[z] );
        if ( sico.MODEL == 7 ) fprintf(f_volumes, "\t%.3f\t0\t%.3f\t0", vol_zone2[z], vol_zone3[z] );
    }
    
    fprintf(f_volumes, "\n");
   
    
// -- STOP --- Summary of release and initial state of the flow -------------------------------------------------


// -- START -- Preparing and writing output raster maps of initial values ---------------------------------------


    #ifdef WITHGRASS


        for ( i=0; i<sico.IMAX; i++ ) {

            if ( sico.TOPOFACT == 0 ) topofact = 1.0;
            else topofact = cos( betaxy[i] );

            v[0] = awt[i][0]; // mixture or PHASE 1 flow height
            if ( aw[i][0] > sico.IMPTHR[3] ) v[1] = fdiv( aw[i][2], aw[i][0] * topofact, sico.IMPTHR[3] ); else v[1] = 0;
            if ( aw[i][0] > sico.IMPTHR[3] ) v[2] = fdiv( -aw[i][1], aw[i][0] * topofact, sico.IMPTHR[3] ); else v[2] = 0; // PHASE 1 x and y velocities

            if ( sico.MODEL <= 3 ) { // one-phase models

                if ( aw[i][0] > sico.IMPTHR[3] ) v[3] = aw[i][3] / aw[i][0]; else v[3] = 0; // temperature
                v[4] = awt[i][4]; // change of basal surface
                v[5] = awt[i][5]; // velocity
                v[6] = awt[i][6]; // flow kinetic energy
                v[7] = awt[i][7]; // flow pressure
                        
            } else if ( sico.MODEL == 7 ) { // multi-phase model

                v[3] = aw[i][3]; // PHASE 2 flow height
                if ( aw[i][3] > sico.IMPTHR[3] ) v[4] = fdiv( aw[i][5], aw[i][3] * topofact, sico.IMPTHR[3] ); else v[4] = 0;
                if ( aw[i][3] > sico.IMPTHR[3] ) v[5] = fdiv( -aw[i][4], aw[i][3] * topofact, sico.IMPTHR[3] ); else v[5] = 0; // PHASE 2 x and y velocities
                v[6] = aw[i][6]; // PHASE 3 flow height
                if ( aw[i][6] > sico.IMPTHR[3] ) v[7] = fdiv( aw[i][5], aw[i][6] * topofact, sico.IMPTHR[3] ); else v[7] = 0;
                if ( aw[i][6] > sico.IMPTHR[3] ) v[8] = fdiv( -aw[i][4], aw[i][6] * topofact, sico.IMPTHR[3] ); else v[8] = 0; // PHASE 3 x and y velocities
                if ( aw[i][0] + aw[i][3] + aw[i][6] > sico.IMPTHR[3] ) v[9] = aw[i][9] / ( aw[i][0] + aw[i][3] + aw[i][6] ); else v[9] = 0; // temperature
                v[10] = awt[i][10];
                v[11] = awt[i][11];
                v[12] = awt[i][12]; // change of basal surface
                v[13] = awt[i][13];
                v[14] = awt[i][14];
                v[15] = awt[i][15]; // flow velocities
                v[16] = awt[i][16]; // total flow height
                v[17] = awt[i][17];
                v[18] = awt[i][18];
                v[19] = awt[i][19];
                v[20] = awt[i][20]; // flow kinetic energies
                v[21] = awt[i][21];
                v[22] = awt[i][22];
                v[23] = awt[i][23];
                v[24] = awt[i][24]; // flow pressures
                v[25] = awt[i][25]; //hentri + hentri2 + hentri3; // total change of basal surface !!!CHECK
            }

            for ( k=0; k<nvect_red; k++ ) outv[px[i]][py[i]][k] = v[k];
        }


    #endif


    if ( sico.PBG == 0 ) {

         for ( i=0; i<sico.IMAX; i++ ) {

             alpha = falpha( betax[i], betay[i], sico );
             ppbg1[i] = (int)( 255.0 * ( cos( 45 * sico.PI / 180.0 ) * cos( betaxy[i] ) + sin( 45 * sico.PI / 180.0 ) * sin( betaxy[i] ) * cos( 135 * sico.PI / 180.0 - alpha )));
         }

         sprintf( mv, "%shillshade0000", prefix );
         foutascind( ppbg1, px, py, outmaps, mv, sico ); // writing ascii raster map
    }

    if ( sico.MULT == 0 ) { // for single model run

        if ( sico.TSUNAMI != 0 ) {

             for ( i=0; i<sico.IMAX; i++ ) {

                 htsun[i] = 0;
             }

             sprintf( mv, "%shtsun0000", prefix );
             foutascindf( htsun, px, py, outmaps, mv, sico ); // writing ascii raster map
        }

        if ( nout < 10 ) sprintf( madd, "000"); // fill string (for maintaining correct order in list of maps)
        else if ( nout < 100 ) sprintf( madd, "00");
        else if ( nout < 1000 ) sprintf( madd, "0");

        for ( k=0; k<nvect_red; k++ ) {
                
            sprintf( mv, "%s%s%s%i", prefix, mv0[k], madd, nout-1 ); // names of output raster maps

            if ( sico.AFLAG == 1 || ( sico.MODEL <= 3 && ( k==0 || k==3 || k==4 )) || ( sico.MODEL == 7 && ( k==0 || k==3 || k==6 || k==9 || k==10 || k==11 ||  k == 12 || k == 16 || k == 25 ))) {


                #ifdef WITHGRASS


                    foutrast ( mv, outv, sico, k, tsum ); // if GRASS is used, writing GRASS raster maps


                #endif


                if ( sico.MODEL <= 3 ) foutasc( awt, px, py, outmaps, mv, k, sico ); // writing ascii raster maps (one-phase models)
                else if ( sico.MODEL == 7 ) foutasc( awt, px, py, outmaps, mv, k, sico ); // writing ascii raster maps (multi-phase model)
            }
        }

        if ( sico.MODEL <= 3 ) { // ascii raster maps of maximum flow height at time step
            sprintf( mv, "%s%s%s%i", prefix, mv0[8], madd, nout-1 );
            foutasc( awt, px, py, outmaps, mv, 8, sico );
                    
        } else if ( sico.MODEL == 7 ) {
            sprintf( mv, "%s%s%s%i", prefix, mv0[32], madd, nout-1 );
            foutasc( awt, px, py, outmaps, mv, 32, sico );
        }
    }

    sprintf( mv, "%selev", prefix ); // name of output elevation raster map   


    #ifdef WITHGRASS


        for ( i=0; i<sico.IMAX; i++ ) outv[px[i]][py[i]][0] = pelev[i];
        foutrast ( mv, outv, sico, 0, tsum ); // writing GRASS raster map


    #endif


    for ( i=0; i<sico.IMAX; i++ ) awt[i][0] = pelev[i];
    foutasc( awt, px, py, outmaps, mv, 0, sico ); // writing ascii raster map

    if ( sico.IMPACTAREA != 0 ) {

        sprintf( mv, "%simpactarea", prefix ); // name of output impact area raster map  
        sprintf( mv2, "%sobsbini", prefix ); 

        #ifdef WITHGRASS


            for ( i=0; i<sico.IMAX; i++ ) { outv[px[i]][py[i]][0] = pimpactarea[i]; if( pimpactarea[i] > 0 ) outv[px[i]][py[i]][1] = 1; else outv[px[i]][py[i]][1] = 0; }
            foutrast ( mv, outv, sico, 0, tsum ); // writing GRASS raster maps
            foutrast ( mv2, outv, sico, 1, tsum );

        #endif


        for ( i=0; i<sico.IMAX; i++ ) { awt[i][0] = pimpactarea[i]; if( pimpactarea[i] > 0 ) awt[i][1] = 1; else awt[i][1] = 0; }
        foutasc( awt, px, py, outmaps, mv, 0, sico ); // writing ascii raster maps
        foutasc( awt, px, py, outmaps, mv2, 1, sico );
    }

    if ( sico.HDEPOSIT != 0 ) {

        sprintf( mv, "%shdeposit", prefix ); // name of output height of deposition raster map 
        sprintf( mv2, "%sobsbind", prefix );

        #ifdef WITHGRASS


            for ( i=0; i<sico.IMAX; i++ ) { outv[px[i]][py[i]][0] = phdeposit[i]; if( phdeposit[i] > 0 ) outv[px[i]][py[i]][1] = 1; else outv[px[i]][py[i]][1] = 0; }
            foutrast ( mv, outv, sico, 0, tsum ); // writing GRASS raster maps
            foutrast ( mv2, outv, sico, 1, tsum );

        #endif


        for ( i=0; i<sico.IMAX; i++ ) { awt[i][0] = phdeposit[i]; if( phdeposit[i] > 0 ) awt[i][1] = 1; else awt[i][1] = 0; }
        foutasc( awt, px, py, outmaps, mv, 0, sico ); // writing ascii raster maps
        foutasc( awt, px, py, outmaps, mv2, 1, sico );
    }

    
// -- STOP --- Preparing and writing output raster maps of initial values ---------------------------------------


    
// *** Start of loop over time steps until break criterion is met: ----------------------------------------------


    while ( tsum <= tmax && ccontinue == 1 ) {
    
        if ( tint >= tout ) tint -= tout; // resetting time interval for output, if required


// *** Start of loop over two steps, each moving the system half of a cell (NOC scheme) -------------------------


        for ( p=0; p<2; p++ ) {
        
            cslide = 0; cflow = 0;
            for ( i=0; i<sico.IMAX; i++ ) {
                for ( k=0; k<sico.NVECTMIN; k++ ) {

                    awt[i][k] = 0; // resetting temporary state variables
                    wintd[i][k] = 0;
                    wintdtest[i][k] = 0;
                }
              
                pxslide[i] = 0; anctr0 = 0;

                if ( ptslide[i] > tsum || ptslide[i] == -777 ) cslide = 1; // resetting controls for relevance of sliding and flowing
                if ( ptslide[i] <= tsum || ptslide[i] == -777 || cslide < 2 ) cflow = 1;

                if ( ptslide[i] == -777 && sico.MODEL == 7 && cdomain[i] != 0 ) {
                
                    anctr0 = 1;
                    for ( j=0; j<9; j++ ) { if ( aw[in[i][j]][6] > sico.HFLOWMIN ) anctr0 = 0; }
                }
                
                if ( ptslide[i] > tsum || anctr0 == 1 ) pxslide[i] = 1;
            }

            if ( cflow == cflowpre ) tlength = tlength0; // time step length
            else tlength = 0.001;
            cflowpre = cflow;
            ib[1] = 0;


// -- START -- Preparing data for progressive collapse (constant volume) ----------------------------------------


            // *** Options trelease and trelstop, negative value of trelstop activates this release mode
            // *** Release of mass starts from the top of the release area
            // *** Similar volumes are released at each time step, until all release material has been released
            // *** The basal surface is updated at each time step (in contrast to all other release types, the initial terrain includes the release mass)
            // *** Applicable with the one-phase models and the multi-phase model

            if ( sico.COLLAPSE == 1 && sico.TRELEASE == 1 && sico.TRELSTOP == 1 ) {

                if ( adaptograph == 1 ) {

                    for ( adak=0; adak<=adatmax; adak++ ) { // identifying relevant line of input adaptograph
                        if ( adaada[adak][0] <= tsum ) adat=adak;
                        else break;
                    }
                }

                qtrelspan = 0;

                for ( i=0; i<sico.IMAX; i++ ) {

                    if ( tsum >= ptrelease[i] && tsum < ptrelstop[i] ) {

                        if ( ptrelstop[i] - ptrelease[i] > qtrelspan ) {

                            qtrelspan = ptrelstop[i] - ptrelease[i]; // start, stop, and duration of progressive collapse
                            qtrelstart = ptrelease[i];
                            qtrelstop = ptrelstop[i];
                        }
                    }
                }

                rhrem = tlengthpre / qtrelspan; // fraction of input hydrograph to be released at current time step

                if ( adaptograph == 1 ) {

                    if ( rhrem * adaada[adat][1] + adaada[adat][2] > 0 ) rhrem = rhrem * adaada[adat][1] + adaada[adat][2]; else rhrem = 0;
                }

                if ( tsum >= qtrelstart && tsum < qtrelstop && qh_test > 0 ) {

                    qvol_test = 0; // initializing target release volume

                    if ( sico.MODEL <= 3 ) vol_flow0all = vol_flow0;
                    else vol_flow0all = vol_flow0 + vol_flow02 + vol_flow03;

                    while ( qvol_test < vol_flow0all * rhrem ) { // loop over test release heights until target release volume is met

                        qvol_test = 0;
                        qh_test -= vol_flow0all * pow(10, -10 ); // updating test release height
                        
                        for ( i=0; i<sico.IMAX; i++ )  {

                            if ( tsum >= ptrelease[i] && tsum < ptrelstop[i] ) {
                            
                                if ( sico.MODEL <= 3 ) phreleaseall = phrelease[i]; // total release height
                                else phreleaseall = phrelease[i] + phrelease2[i] + phrelease3[i];
                            
                                qvol_test += sico.CSZ * sico.CSZ * ffmax( 0, phreleaseall - qh_test ); // updating target volume
                            }
                        }
                    }

                    for ( i=0; i<sico.IMAX; i++ ) {

                        if ( tsum >= ptrelease[i] && tsum < ptrelstop[i] && qh_test > 0 ) {

                            if ( sico.MODEL <= 3 ) { phreleaseall = phrelease[i]; alpha1 = 1; } // total release height and phase fractions
                            else { phreleaseall = phrelease[i] + phrelease2[i] + phrelease3[i];
                                if ( phreleaseall > 0 ) { alpha1 = phrelease[i] / phreleaseall; alpha2 = phrelease2[i] / phreleaseall; }
                            }

                            if ( phreleaseall > 0 ) {

                                qhreleaseall = ffmax( 0, phreleaseall - qh_test ); // release hydrograph height
                                qhrelease[i] = qhreleaseall * alpha1; phrelease[i] -= qhrelease[i]; // phase-specific release hydrograph heights and remaining release heights
                                
                                if ( sico.MODEL == 7 ) { 

                                    qhrelease2[i] = qhreleaseall * alpha2; phrelease2[i] -= qhrelease2[i];
                                    qhrelease3[i] = qhreleaseall * ( 1 - alpha1 - alpha2 ); phrelease3[i] -= qhrelease3[i];
                                }
                                
                                pelev[i] -= qhreleaseall; // updating terrain elevation

                                if ( sico.MODEL <= 3 ) awt[i][4] = awt[i][4] - qhrelease[i]; // updating changes of basal topography
                                else { awt[i][10] = awt[i][10] - qhrelease[i]; awt[i][11] = awt[i][11] - qhrelease2[i]; awt[i][12] = awt[i][12] - qhrelease3[i]; }
                            }
                        }
                    }

                } else {

                    for ( i=0; i<sico.IMAX; i++ ) {

                        qhrelease[i] = 0; // setting release hydrograph height to zero outside of release time span
                        if ( sico.MODEL == 7 ) { qhrelease2[i] = 0; qhrelease3[i] = 0; }
                    }
                }
                
                for ( ix=0; ix<ib[0]; ix++ ) { // updating slopes and topography-following cell sizes

                    i = ibasket[0][ix];

                    cin = fnoosc( aw, in, i, sico );

                    betax[i] = fbeta( pelev[cin[0]], pelev[cin[1]], (float)(cin[4]), sico ); // slopes
                    betay[i] = fbeta( pelev[cin[2]], pelev[cin[3]], (float)(cin[5]), sico );
                    betaxh[i] = fbeta( pelev[cin[0]]+aw[cin[0]][0], pelev[cin[1]]+aw[cin[1]][0], (float)(cin[4]), sico ); // slopes including flow heights
                    betayh[i] = fbeta( pelev[cin[2]]+aw[cin[2]][0], pelev[cin[3]]+aw[cin[3]][0], (float)(cin[5]), sico );
                    betaxy[i] = fbetaxy( betax[i], betay[i] );

                    if ( sico.LAYERS0 > 0 ) {
                            
                        betax2[i] = fbeta( pelev[cin[0]]+aw[cin[0]][0], pelev[cin[1]]+aw[cin[1]][0], (float)(cin[4]), sico );
                        betay2[i] = fbeta( pelev[cin[2]]+aw[cin[2]][0], pelev[cin[3]]+aw[cin[3]][0], (float)(cin[5]), sico );

                        betaxy2[i] = fbetaxy( betax2[i], betay2[i] );
                                                               
                        betaxh2[i] = fbeta( pelev[cin[0]]+aw[cin[0]][0]+aw[cin[0]][3], pelev[cin[1]]+aw[cin[1]][0]+aw[cin[1]][3], (float)(cin[4]), sico );
                        betayh2[i] = fbeta( pelev[cin[2]]+aw[cin[2]][0]+aw[cin[2]][3], pelev[cin[3]]+aw[cin[3]][0]+aw[cin[3]][3], (float)(cin[5]), sico );
                                    
                        betax3[i] = fbeta( pelev[cin[0]]+aw[cin[0]][0]+aw[cin[0]][3], pelev[cin[1]]+aw[cin[1]][0]+aw[cin[1]][3], (float)(cin[4]), sico );
                        betay3[i] = fbeta( pelev[cin[2]]+aw[cin[2]][0]+aw[cin[2]][3], pelev[cin[3]]+aw[cin[3]][0]+aw[cin[3]][3], (float)(cin[5]), sico );
                        
                        betaxy3[i] = fbetaxy( betax3[i], betay3[i] );
                    }

                    if ( sico.LAYERS0 > 0 || sico.SURFACE > 1 ) {
                                     
                        betaxh3[i] = fbeta( pelev[cin[0]]+aw[cin[0]][0]+aw[cin[0]][3]+aw[cin[0]][6], pelev[cin[1]]+aw[cin[1]][0]+aw[cin[1]][3]+aw[cin[1]][6], (float)(cin[4]), sico );
                        betayh3[i] = fbeta( pelev[cin[2]]+aw[cin[2]][0]+aw[cin[2]][3]+aw[cin[2]][6], pelev[cin[3]]+aw[cin[3]][0]+aw[cin[3]][3]+aw[cin[3]][6], (float)(cin[5]), sico );
                    }

                    free( cin );
                }

            } else if ( sico.TRELEASE == 1 ) { // with release time but no progressive collapse, setting release hydrograph height to release height for further use

                if ( adaptograph == 1 ) {

                    for ( adak=0; adak<=adatmax; adak++ ) { // identifying relevant line of input adaptograph
                        if ( adaada[adak][0] <= tsum ) adat=adak;
                        else break;
                    }
                }

                for ( i=0; i<sico.IMAX; i++ ) { // setting release hydrograph height, avoiding negative values

                    if ( adaptograph == 1 ) {

                        if ( phrelease[i] * adaada[adat][1] + adaada[adat][2] > 0 ) qhrelease[i] = phrelease[i] * adaada[adat][1] + adaada[adat][2]; else qhrelease[i] = 0;
                        if ( sico.MODEL == 7 ) { 
                            if ( phrelease2[i] * adaada[adat][3] + adaada[adat][4] > 0 ) qhrelease2[i] = phrelease2[i] * adaada[adat][3] + adaada[adat][4]; else qhrelease2[i] = 0;
                            if ( phrelease3[i] * adaada[adat][5] + adaada[adat][6] > 0 ) qhrelease3[i] = phrelease3[i] * adaada[adat][5] + adaada[adat][6]; else qhrelease3[i] = 0;
                        }
                    } else {
                    
                        if ( phrelease[i] > 0 ) qhrelease[i] = phrelease[i]; else qhrelease[i] = 0;
                        if ( sico.MODEL == 7 ) { 
                            if ( phrelease2[i] > 0 ) qhrelease2[i] = phrelease2[i]; else qhrelease2[i] = 0;
                            if ( phrelease3[i] > 0 ) qhrelease3[i] = phrelease3[i]; else qhrelease3[i] = 0;
                        }      
                    }
                }
            }
            
            
// -- STOP --- Preparing data for progressive collapse (constant volume) ----------------------------------------


// -- START -- Preparing data for extrusive release -------------------------------------------------------------


            // *** Options trelease and trelstop, positive value of trelstop activates this release mode
            // *** Release height has to be provided as rate of extrusion (m/s)
            // *** Similar volumes are released at each time step, until all release material has been released
            // *** Applicable with the one-phase models and the multi-phase model

            if ( sico.TRELSTOP == 1 && sico.COLLAPSE == 0 ) tfact = tlengthpre; else tfact = 1;

            for ( i=0; i<sico.IMAX; i++ ) {

                if ( sico.TRELEASE == 1 ) ttrelease = ptrelease[i]; // release time
                else ttrelease = 0;

                if ( sico.TRELSTOP == 1 && ptrelstop[i] > 0 ) ttrelstop = ptrelstop[i]; // release hydrograph stop time
                else ttrelstop = ttrelease;

                if ( sico.TRELEASE == 1 && tsum >= ttrelease && ( ptrelease[i] != 1000 * sico.UNDEF || tsum < ttrelstop )) {

                    if ( sico.MODEL <= 3 ) {

                        aw[i][0] += qhrelease[i] * tfact; // release depth

                        aw[i][1] += qhrelease[i] * tfact * pviny[i] * -1;
                            // release momentum in x direction
                        aw[i][2] += qhrelease[i] * tfact * pvinx[i];
                            // release momentum in y direction
                            
                    } else if ( sico.MODEL == 7 ) {

                        aw[i][0] += qhrelease[i] * tfact; // PHASE 1 release depth
                        aw[i][3] += qhrelease2[i] * tfact; // PHASE 2 release depth
                        aw[i][6] += qhrelease3[i] * tfact; // PHASE 3 release depth

                        aw[i][1] += qhrelease[i] * tfact * pviny[i] * -1; // PHASE 1 release momenta in x and y directions
                        aw[i][2] += qhrelease[i] * tfact * pvinx[i];

                        aw[i][4] += qhrelease2[i] * tfact * pviny2[i] * -1; // PHASE 2 release momenta in x and y directions
                        aw[i][5] += qhrelease2[i] * tfact * pvinx2[i];

                        aw[i][7] += qhrelease3[i] * tfact * pviny3[i] * -1; // PHASE 3 release momenta in x and y directions
                        aw[i][8] += qhrelease3[i] * tfact * pvinx3[i];
                    }

                    if ( sico.COLLAPSE == 0 ) ptrelease[i] = 1000 * sico.UNDEF;
                }
            }


// -- STOP --- Preparing data for extrusive release -------------------------------------------------------------


// -- START -- Updating flow depths and velocities according to input hydrographs -------------------------------


            vol_hydbef = 0; vol_hydaft = 0;
            if ( sico.MODEL == 7 ) { vol_hydbef2 = 0; vol_hydaft2 = 0; }
            if ( sico.MODEL == 7 ) { vol_hydbef3 = 0; vol_hydaft3 = 0; }

            if ( hydrograph == 1 ) {

                if ( tsum == 0 || sico.HYDADD == 0 ) tlengthx = 1;
                else tlengthx = tlengthpre;

                for ( hydj = 0; hydj < hydnin; hydj++ ) { // loop over all input hydrographs

                    for ( hydk=0; hydk<=hydtmax[hydj]; hydk++ ) { // identifying relevant line of input hydrograph
                        if ( hydhyd[hydk][0][hydj] <= tsum ) hydt=hydk;
                        else break;
                    }

                    if ( sico.MODEL <= 3 ) hhyd0 = hydhyd[hydt][1][hydj]; // total flow height at centre of hydrograph
                    else if ( sico.MODEL == 7 ) hhyd0 = hydhyd[hydt][1][hydj] + hydhyd[hydt][3][hydj] + hydhyd[hydt][5][hydj];

                    if ( hhyd0 > sico.HFLOWMIN && sico.HYDADD < 2 ) {

                        for ( hydk=1; hydk<= hydp[0][hydj]; hydk++ ) { // loop over all cells of hydrograph profile

                            carea = pow ( sico.CSZ, 2 ); // area of cell

                            hhyd = hhyd0 - ( pelev[hydp[hydk][hydj]] - hydelev[hydj] ); // total flow height at cell

                            for ( k=0; k<sico.NVECTMIN; k++ ) {
                                if ( sico.HYDADD == 1 ) hydbef[k] = aw[hydp[hydk][hydj]][k]; else hydbef[k] = 0;
                            }

                            vol_hydbef += aw[hydp[hydk][hydj]][0] * carea;
                            if ( sico.MODEL == 7 ) vol_hydbef2 += aw[hydp[hydk][hydj]][3] * carea;
                            if ( sico.MODEL == 7 ) vol_hydbef3 += aw[hydp[hydk][hydj]][6] * carea;

                            if ( hhyd > sico.HFLOWMIN ) {

                                aw[hydp[hydk][hydj]][0] = hydbef[0] + tlengthx * hydhyd[hydt][1][hydj] * hhyd / hhyd0 * cos( betaxy[hydp[hydk][hydj]] );
                                    // mixture or PHASE 1 flow depth
                                aw[hydp[hydk][hydj]][1] = hydbef[1] + tlengthx * hydhyd[hydt][1][hydj] * hhyd / hhyd0 * hydhyd[hydt][2][hydj]
                                    * cos( hydalpha[hydj] - sico.PI * 0.5 )  * cos( betaxy[hydp[hydk][hydj]] ); // mixture or PHASE 1 x momentum
                                aw[hydp[hydk][hydj]][2] = hydbef[2] + tlengthx * hydhyd[hydt][1][hydj] * hhyd / hhyd0 * hydhyd[hydt][2][hydj]
                                    * sin( hydalpha[hydj] - sico.PI * 0.5 )  * cos( betaxy[hydp[hydk][hydj]] ); // mixture or PHASE 1 y momentum
                            }

                            if ( hhyd > sico.HFLOWMIN && sico.MODEL == 7 ) {

                                aw[hydp[hydk][hydj]][3] = hydbef[3] + tlengthx * hydhyd[hydt][3][hydj] * hhyd / hhyd0  * cos( betaxy[hydp[hydk][hydj]] ); // PHASE 2 flow depth
                                aw[hydp[hydk][hydj]][4] = hydbef[4] + tlengthx * hydhyd[hydt][3][hydj] * hhyd / hhyd0 * hydhyd[hydt][4][hydj]
                                    * cos( hydalpha[hydj] - sico.PI * 0.5 ) * cos( betaxy[hydp[hydk][hydj]] ); // PHASE 2 x momentum
                                aw[hydp[hydk][hydj]][5] = hydbef[5] + tlengthx * hydhyd[hydt][3][hydj] * hhyd / hhyd0 * hydhyd[hydt][4][hydj]
                                    * sin( hydalpha[hydj] - sico.PI * 0.5 ) * cos( betaxy[hydp[hydk][hydj]] ); // PHASE 2 y momentum
                            }

                            if ( hhyd > sico.HFLOWMIN && sico.MODEL == 7 ) {

                                aw[hydp[hydk][hydj]][6] = hydbef[6] + tlengthx * hydhyd[hydt][5][hydj] * hhyd / hhyd0  * cos( betaxy[hydp[hydk][hydj]] ); // PHASE 3 flow depth
                                aw[hydp[hydk][hydj]][7] = hydbef[7] + tlengthx * hydhyd[hydt][5][hydj] * hhyd / hhyd0 * hydhyd[hydt][6][hydj]
                                    * cos( hydalpha[hydj] - sico.PI * 0.5 ) * cos( betaxy[hydp[hydk][hydj]] ); // PHASE 3 x momentum
                                aw[hydp[hydk][hydj]][8] = hydbef[8] + tlengthx * hydhyd[hydt][5][hydj] * hhyd / hhyd0 * hydhyd[hydt][6][hydj]
                                    * sin( hydalpha[hydj] - sico.PI * 0.5 ) * cos( betaxy[hydp[hydk][hydj]] ); // PHASE 3 y momentum
                            }

                            vol_hydaft += aw[hydp[hydk][hydj]][0] * carea;
                            if ( sico.MODEL == 7 ) vol_hydaft2 += aw[hydp[hydk][hydj]][3] * carea;
                            if ( sico.MODEL == 7 ) vol_hydaft3 += aw[hydp[hydk][hydj]][6] * carea;

                            if ( fvalid( aw[hydp[hydk][hydj]][0], aw[hydp[hydk][hydj]][3], aw[hydp[hydk][hydj]][6], sico ) == 1
                                && cstopped[hydp[hydk][hydj]] == 0 && cdomain[hydp[hydk][hydj]] > 0 ) { // flow cells

                                for ( j=0; j<9; j++ ) { // 1st row flow boundary cells

                                    if ( icheck[in[hydp[hydk][hydj]][j]][1] == 0 && cstopped[in[hydp[hydk][hydj]][j]] == 0
                                        && cdomain[in[hydp[hydk][hydj]][j]] > 0 ) {

                                        ibasket[1][ib[1]] = in[hydp[hydk][hydj]][j];
                                        ib[1] += 1;
                                        icheck[in[hydp[hydk][hydj]][j]][1] = 1;
                                    }
                                }
                            }
                        }

                    } else if ( hhyd0 / ( sico.CSZ * sico.CSZ ) > sico.HFLOWMIN && sico.HYDADD == 2 ) {

                        carea = pow ( sico.CSZ, 2 ); // area of cell

                        vol_hydbef += aw[hydi[hydj]][0] * carea;
                        if ( sico.MODEL == 7 ) vol_hydbef2 += aw[hydi[hydj]][3] * carea;
                        if ( sico.MODEL == 7 ) vol_hydbef3 += aw[hydi[hydj]][6] * carea;

                        aw[hydi[hydj]][0] += tlengthx * hydhyd[hydt][1][hydj] / carea;
                            // mixture or PHASE 1 flow depth
                        aw[hydi[hydj]][1] += tlengthx * hydhyd[hydt][1][hydj] / carea * hydhyd[hydt][2][hydj]
                            * cos( hydalpha[hydj] - sico.PI * 0.5 ); // mixture or PHASE 1 x momentum
                        aw[hydi[hydj]][2] += tlengthx * hydhyd[hydt][1][hydj] / carea * hydhyd[hydt][2][hydj]
                            * sin( hydalpha[hydj] - sico.PI * 0.5 ); // mixture or PHASE 1 y momentum

                        if ( sico.MODEL == 7 ) {

                            aw[hydi[hydj]][3] += tlengthx * hydhyd[hydt][3][hydj] / carea;
                                // PHASE 2 flow depth
                            aw[hydi[hydj]][4] += tlengthx * hydhyd[hydt][3][hydj] / carea * hydhyd[hydt][4][hydj]
                                * cos( hydalpha[hydj] - sico.PI * 0.5 ); // PHASE 2 x momentum
                            aw[hydi[hydj]][5] += tlengthx * hydhyd[hydt][3][hydj] / carea * hydhyd[hydt][4][hydj]
                                * sin( hydalpha[hydj] - sico.PI * 0.5 ); // PHASE 2 y momentum
                        }

                        if ( sico.MODEL == 7 ) {

                            aw[hydi[hydj]][6] += tlengthx * hydhyd[hydt][5][hydj] / carea;
                                // PHASE 3 flow depth
                            aw[hydi[hydj]][7] += tlengthx * hydhyd[hydt][5][hydj] / carea * hydhyd[hydt][6][hydj]
                                * cos( hydalpha[hydj] - sico.PI * 0.5 ); // PHASE 3 x momentum
                            aw[hydi[hydj]][8] += tlengthx * hydhyd[hydt][5][hydj] / carea * hydhyd[hydt][6][hydj]
                                * sin( hydalpha[hydj] - sico.PI * 0.5 ); // PHASE 3 y momentum
                        }

                        vol_hydaft += aw[hydi[hydj]][0] * carea;
                        if ( sico.MODEL == 7 ) vol_hydaft2 += aw[hydi[hydj]][3] * carea;
                        if ( sico.MODEL == 7 ) vol_hydaft3 += aw[hydi[hydj]][6] * carea;
                        
                        if ( fvalid( aw[hydi[hydj]][0], aw[hydi[hydj]][3], aw[hydi[hydj]][6], sico ) == 1
                            && cstopped[hydi[hydj]] == 0 && cdomain[hydi[hydj]] > 0 ) { // flow cells

                            for ( j=0; j<9; j++ ) { // 1st row flow boundary cells

                                if ( icheck[in[hydi[hydj]][j]][1] == 0 && cstopped[in[hydi[hydj]][j]] == 0
                                    && cdomain[in[hydi[hydj]][j]] > 0 ) {

                                    ibasket[1][ib[1]] = in[hydi[hydj]][j];
                                    ib[1] += 1;
                                    icheck[in[hydi[hydj]][j]][1] = 1;
                                }
                            }
                        }                        
                    }
                }

                vol_hyd += ( vol_hydaft - vol_hydbef ); // volume added through hydrograph
                if ( sico.MODEL == 7 ) vol_hyd2 += ( vol_hydaft2 - vol_hydbef2 );
                if ( sico.MODEL == 7 ) vol_hyd3 += ( vol_hydaft3 - vol_hydbef3 );
            }


// -- STOP --- Updating flow depths and velocities according to input hydrographs -------------------------------


// -- START -- Writing hydrograph infos to files ----------------------------------------------------------------


            if ( p == 1 && hydrograph == 1 ) {

                for ( hydj = 0; hydj < hydnin; hydj++ ) { // loop over all input hydrographs

                    if ( sico.MULT == 0 && ctrl_hydout == 1 ) { // for first time step after output

                        hydh = aw[hydi[hydj]][0] / cos( betaxy[hydi[hydj]] ); // mixture or PHASE 1 flow height
                        if ( hydh > sico.HFLOWMIN ) hydv = pow( pow( aw[hydi[hydj]][1], 2 ) + pow( aw[hydi[hydj]][2], 2 ), 0.5 ) / aw[hydi[hydj]][0];
                        else hydv= 0; // mixture or PHASE 1 flow velocity

                        if ( sico.MODEL <= 3 ) {

                            hyde = aw[hydi[hydj]][4]; // entrained height
                            hydh2 = 0; hydv2 = 0; hyde2 = 0; hydh3 = 0; hydv3 = 0; hyde3 = 0;
                            
                        } else if ( sico.MODEL == 7 ) {

                            hydh2 = aw[hydi[hydj]][3] / cos( betaxy[hydi[hydj]] ); // PHASE 2 flow height
                            if ( hydh2 > sico.HFLOWMIN ) hydv2 = pow( pow( aw[hydi[hydj]][4], 2 ) + pow( aw[hydi[hydj]][5], 2 ), 0.5 ) / aw[hydi[hydj]][3];
                            else hydv2= 0; // PHASE 2 flow velocity
                            hydh3 = aw[hydi[hydj]][6] / cos( betaxy[hydi[hydj]] ); // PHASE 3 flow height
                            if ( hydh3 > sico.HFLOWMIN ) hydv3 = pow( pow( aw[hydi[hydj]][7], 2 ) + pow( aw[hydi[hydj]][8], 2 ), 0.5 ) / aw[hydi[hydj]][6];
                            else hydv3= 0; // PHASE 3 flow velocity
                            hyde = aw[hydi[hydj]][10];
                            hyde2 = aw[hydi[hydj]][11];
                            hyde3 = aw[hydi[hydj]][12]; // changes of basal topography
                        }

                        hydq = 0; hydq2 = 0; hydq3 = 0; // resetting discharges

                        if ( hydalpha[hydj] == 0 || hydalpha[hydj] == sico.PI * 0.5 || hydalpha[hydj] == sico.PI || hydalpha[hydj] == 3 * sico.PI * 0.5 ) hydfalpha = 1;
                        else if ( fabs( 1 / sin( hydalpha[hydj] )) < fabs( 1 / cos( hydalpha[hydj] ))) hydfalpha = fabs( 1 / sin( hydalpha[hydj] ));
                        else hydfalpha = fabs( 1 / cos( hydalpha[hydj] )); // correction factor for profile direction

                        for ( hydk = 1; hydk <= hydp[0][hydj]; hydk++ ) { // loop over all cells of profile

                            alpha = falpha( betax[hydp[hydk][hydj]], betay[hydp[hydk][hydj]], sico ); // aspect
                            hydbeta = atan ( tan( betaxy[hydp[hydk][hydj]] ) * cos ( alpha - hydalpha[hydj] + sico.PI * 0.5 )); // corrected slope

                            hydfcorr = sico.CSZ * hydfalpha / cos( hydbeta ); // reference length for discharge
                            hydm0 = pow( pow(aw[hydp[hydk][hydj]][1], 2 ) + pow( aw[hydp[hydk][hydj]][2], 2), 0.5 );

                            if ( aw[hydp[hydk][hydj]][0] > sico.HFLOWMIN && hydm0 > 0 ) {

                                hydmx = aw[hydp[hydk][hydj]][1];
                                hydmy = aw[hydp[hydk][hydj]][2];

                                if ( hydmy > 0 ) hydmalpha = acos( hydmx / hydm0 );
                                else hydmalpha = 2 * sico.PI - acos( hydmx / hydm0 );
                                hydnalpha = hydalpha[hydj] - sico.PI * 0.5; if ( hydnalpha < 0 ) hydnalpha += 2 * sico.PI;
                                hydm = hydm0 * cos ( hydnalpha - hydmalpha ); // momentum in direction perpendicular to profile
                                hydq += hydm * hydfcorr; // updating mixture or PHASE 1 discharge
                            }

                            if ( sico.MODEL == 7 ) {

                                hydm0 = pow( pow(aw[hydp[hydk][hydj]][4], 2 ) + pow( aw[hydp[hydk][hydj]][5], 2), 0.5 );

                                if ( aw[hydp[hydk][hydj]][3] > sico.HFLOWMIN && hydm0 > 0 ) {

                                    hydmx = aw[hydp[hydk][hydj]][4];
                                    hydmy = aw[hydp[hydk][hydj]][5];

                                    if ( hydmy > 0 ) hydmalpha = acos( hydmx / hydm0 );
                                    else hydmalpha = 2 * sico.PI - acos( hydmx / hydm0 );
                                    hydnalpha = hydalpha[hydj] - sico.PI * 0.5; if ( hydnalpha < 0 ) hydnalpha += 2 * sico.PI;
                                    hydm = hydm0 * cos ( hydnalpha - hydmalpha ); // momentum in direction perpendicular to profile
                                    hydq2 += hydm * hydfcorr; // updating PHASE 2 discharge
                                }

                                hydm0 = pow( pow(aw[hydp[hydk][hydj]][7], 2 ) + pow( aw[hydp[hydk][hydj]][8], 2), 0.5 );

                                if ( aw[hydp[hydk][hydj]][6] > sico.HFLOWMIN && hydm0 > 0 ) {

                                    hydmx = aw[hydp[hydk][hydj]][7];
                                    hydmy = aw[hydp[hydk][hydj]][8];

                                    if ( hydmy > 0 ) hydmalpha = acos( hydmx / hydm0 );
                                    else hydmalpha = 2 * sico.PI - acos( hydmx / hydm0 );
                                    hydnalpha = hydalpha[hydj] - sico.PI * 0.5; if ( hydnalpha < 0 ) hydnalpha += 2 * sico.PI;
                                    hydm = hydm0 * cos ( hydnalpha - hydmalpha ); // momentum in direction perpendicular to profile
                                    hydq3 += hydm * hydfcorr; // updating PHASE 3 discharge
                                }
                            }
                        }

                        fprintf(f_hydinfo[hydj], "%.1f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
                            tsum, hydh, hydv, hyde, hydq, hydh2, hydv2, hyde2, hydq2, hydh3, hydv3, hyde3, hydq3); // writing hydrograph info to file
                    }
                }
                
                ctrl_hydout = 0; // resetting control for hydrograph output
            }


// -- STOP --- Writing hydrograph infos to files ----------------------------------------------------------------


// -- START -- Updating computational domains -------------------------------------------------------------------


            // *** cdomain[i] = 0: 1st row edge cell
            // *** cdomain[i] = 1: no-flux cell (suppression of oscillations)
            // *** cdomain[i] = 3: included cells, flow depth < minimum at & around cell
            // *** cdomain[i] = 4: included cells, flow depth >= minimum around cell (flow boundary)
            // *** cdomain[i] = 5: included cells, flow depth >= minimum (flow)

            if ( cslide < 2 ) iloop = ib[0]; else iloop = sico.IMAX;
            for ( ix=0; ix<iloop; ix++ ) {

                if ( cslide < 2 ) i = ibasket[0][ix]; else i = ix;

                if (( fvalid( aw[i][0], aw[i][3], aw[i][6], sico ) == 1 || ( sico.TRELEASE == 1 && ( ptrelease[i] > 0 || sico.COLLAPSE != 0 )))
                    && cstopped[i] == 0 && cdomain[i] > 0 ) { // flow cells

                    for ( j=0; j<9; j++ ) { // 1st row flow boundary cells

                        if ( icheck[in[i][j]][1] == 0 && cstopped[in[i][j]] == 0 && cdomain[in[i][j]] > 0 ) {

                            ibasket[1][ib[1]] = in[i][j];
                            ib[1] += 1;
                            icheck[in[i][j]][1] = 1;
                        }
                    }
                }
            }

            ib[0] = 0;

            for ( ix=0; ix<ib[1]; ix++ ) { // 2nd row flow boundary cells

                i = ibasket[1][ix];

                for ( j=0; j<9; j++ ) {

                    if ( icheck[in[i][j]][0] == 0 && cstopped[in[i][j]] == 0 && cdomain[in[i][j]] > 0 ) {

                        ibasket[0][ib[0]] = in[i][j];
                        ib[0] += 1;
                        icheck[in[i][j]][0] = 1;
                    }
                }
            }

            for ( ix=0; ix<ib[0]; ix++ ) { // updating domains and resetting checks

                i = ibasket[0][ix];

                if ( sico.MODEL <= 3 ) hflow = aw[i][0];
                else if ( sico.MODEL == 7 ) hflow = aw[i][0] + aw[i][3] + aw[i][6];

                if ( cdomain[i] != 0 ) { // for all non-edge cells

                    ctrlv = 0; ctrlr = 1; // resetting controls
                    for ( j=1; j<9; j++ ) { // loop over neighbourhood

                        if ( sico.MODEL <= 3 ) hflown = aw[in[i][j]][0];
                        else if ( sico.MODEL == 7 ) hflown = aw[in[i][j]][0] + aw[in[i][j]][3] + aw[in[i][j]][6];

                        if ( hflown > sico.HFLOWMIN || hflow > sico.HFLOWMIN ) ctrlv = 1;
                            // setting control for cell to positive if not at or surrounded by no-flux area

                        if ( cstopped[i] == 1 || cstopped[in[i][j]] == 1 ) ctrlr = 0;
                            // setting control for cell to negative if at or adjacent to stopped area
                    }
                    if ( ctrlv == 1 && ctrlr == 1 ) cdomain[i] = 3; // updating computational domain
                }

                if ( cdomain[i] == 3 ) { // if cell is part of computational domain

                    ctrlr = 0; // resetting control

                    for ( j=1; j<9; j++ ) { // loop over neighbourhood

                        if ( sico.MODEL <= 3 ) hflown = aw[in[i][j]][0];
                        else if ( sico.MODEL == 7 ) hflown = aw[in[i][j]][0] + aw[in[i][j]][3] + aw[in[i][j]][6];
                        if ( hflown > sico.HFLOWMIN ) ctrlr = 1; // if neighbour cell is flow cell, setting control to positive
                    }
                    if ( ctrlr == 1 ) cdomain[i] = 4; // flow depth > minimum around cell
                    if ( hflow > sico.HFLOWMIN ) cdomain[i] = 5; // flow depth > minimum
                }

                icheck[i][0] = 0;
                icheck[i][1] = 0;
            }


// -- STOP --- Updating computational domains -------------------------------------------------------------------


// *** Start of loop over various time step lengths until CFL criterion is met ----------------------------------


            ccfl = 0; // resetting control for fulfilment of CFL criterion
            while ( ccfl == 0 ) {


// -- START -- Applying initial block sliding -------------------------------------------------------------------


                // *** Global weighted averages of gravity components are computed and used instead of or additionally to local values
                // *** Duration before slide changes to flow locally defined through parameter tslide (raster map) or, with tslide=-777, through presence of fluid

                if ( cslide > 0 ) {
                    
                    for ( k=0; k<kmax; k++ ) { ansumh[k] = 0; ansumslopex[k] = 0; ansumslopey[k] = 0; ansumslopexy[k] = 0; } // initializing cumulative flow height and slopes
                   
                    for ( ix=0; ix<ib[0]; ix++ ) {

                        i = ibasket[0][ix];
                        if ( pxslide[i] > 0 && cdomain[i] == 5 ) {
                            
                            if ( sico.MODEL <= 3 ) hflownn[0] = aw[i][0];
                            else if ( sico.LAYERS == 0 ) hflownn[0] = aw[i][0] + aw[i][3] + aw[i][6];
                            else {
                                
                                hflownn[0] = aw[i][0];
                                hflownn[1] = aw[i][3];
                                hflownn[2] = aw[i][6];                               
                            }

                            cin = fnoosc( aw, in, i, sico );

                            if ( sico.MODEL <= 3 || sico.LAYERS == 0 ) {
                                                    
                                anslopex[0] = fbeta( pelev[cin[0]], pelev[cin[1]], (float)(cin[4]), sico ); // slopes
                                anslopey[0] = fbeta( pelev[cin[2]], pelev[cin[3]], (float)(cin[5]), sico );
                                anslopexy[0] = fbetaxy( anslopex[0], anslopey[0] );

                                anslopex[1] = anslopex[0]; anslopex[2] = anslopex[0];
                                anslopey[1] = anslopey[0]; anslopey[2] = anslopey[0];
                                anslopexy[1] = anslopexy[0]; anslopexy[2] = anslopexy[0];
                                        
                            } else {

                                anslopex[0] = fbeta( pelev[cin[0]], pelev[cin[1]], (float)(cin[4]), sico );
                                anslopey[0] = fbeta( pelev[cin[2]], pelev[cin[3]], (float)(cin[5]), sico );
                                anslopexy[0] = fbetaxy( anslopex[0], anslopey[0] );                  
                                            
                                anslopex[1] = fbeta( pelev[cin[0]] + aw[cin[0]][0], pelev[cin[1]] + aw[cin[1]][0], (float)(cin[4]), sico );
                                anslopey[1] = fbeta( pelev[cin[2]] + aw[cin[2]][0], pelev[cin[3]] + aw[cin[3]][0], (float)(cin[5]), sico );
                                anslopexy[1] = fbetaxy( anslopex[1], anslopey[1] ); 

                                anslopex[2] = fbeta( pelev[cin[0]] + aw[cin[0]][0] + aw[cin[0]][3], pelev[cin[1]] + aw[cin[1]][0] + aw[cin[1]][3], (float)(cin[4]), sico );
                                anslopey[2] = fbeta( pelev[cin[2]] + aw[cin[2]][0] + aw[cin[2]][3], pelev[cin[3]] + aw[cin[3]][0] + aw[cin[3]][3], (float)(cin[5]), sico );
                                anslopexy[2] = fbetaxy( anslopex[2], anslopey[2] );
                            }

                            free( cin );

                            for ( k=0; k<kmax; k++ ) {

                                ansumh[k] += hflownn[k]; // cumulative flow height
                                ansumslopex[k] += hflownn[k] * tan(anslopex[k]); // weighted cumulative x and y slopes
                                ansumslopey[k] += hflownn[k] * tan(anslopey[k]);
                                ansumslopexy[k] += hflownn[k] * tan(anslopexy[k]);
                            }
                        }
                    }

                    for ( ix=0; ix<ib[0]; ix++ ) {

                        i = ibasket[0][ix];
                        if ( pxslide[i] > 0 && cdomain[i] == 5 ) {

                            for ( k=0; k<kmax; k++ ) { 

                                anavgslopex[k] = atan( ansumslopex[k] / ansumh[k] ); // average x and y slopes
                                anavgslopey[k] = atan( ansumslopey[k] / ansumh[k] );
                                anavgslopexy[k] = atan( ansumslopexy[k] / ansumh[k] );
                                angx[k] = sico.GRAVITY * sin( anavgslopex[k] );
                                angy[k] = sico.GRAVITY * sin( anavgslopey[k] );
                                angx[k+3] = sico.GRAVITY * cos( anavgslopexy[k] );
                                angy[k+3] = sico.GRAVITY * cos( anavgslopexy[k] );
                            }
                        }
                    }
                    
                } else {
                
                    for ( k=0; k<kmax; k++ ) {                 
                        angx[k] = 0;
                        angy[k] = 0;
                        angx[k+3] = 0;
                        angy[k+3] = 0;
                    }
                }
                
                sflow.GZMIN = sico.GRAVITY; // minimum gravity component in slope-normal direction

                for ( ix=0; ix<ib[0]; ix++ ) {

                    i = ibasket[0][ix];

                    if ( cdomain[i] <= 3 ) {

                        if ( sico.MODEL <= 3 ) hflownn[0] = aw[i][0];
                        else hflownn[0] = aw[i][0] + aw[i][3] + aw[i][6];

                        if ( sico.GRAVITY * cos( betaxy[i] ) < sflow.GZMIN && hflownn[0] > sico.HFLOWMIN ) sflow.GZMIN = sico.GRAVITY * cos( betaxy[i] );
                    }
                }


// -- STOP --- Applying initial block sliding -------------------------------------------------------------------


// -- START -- Flow propagation with NOC scheme -----------------------------------------------------------------


                // *** High-resolution scheme, system is moved half a cell size each time step
                // *** Technically, the system is shifted one cell each second time step
                // *** At the end, the new state variables are written to a temporary vector (awt[i][k])
                // *** The values from the temporary vector are later written to the permanent vector (aw[i][k])
                // *** as soon as it has been verified that the time step length is sufficiently short to fulfil the CFL criterion
                // *** Otherwise, the procedure is repeated with shorter time step length


                // Fluxes, source terms, and gradients (original coordinate system)

                if ( cflow == 1 ) {

                  for ( ix=0; ix<ib[0]; ix++ ) {

                    i = ibasket[0][ix];

                    if ( cdomain[i] <= 3 ) {

                        for ( k=0; k<sico.NVECTMIN; k++ ) { af[i][k]=0; ag[i][k]=0; as[i][k]=0; asigma_x[i][k]=0; asigma_y[i][k]=0; }

                    } else {

                        if ( sico.PHI == 1 ) {
                            if ( pphi[i] != sico.UNDEF ) sflow.PHI[0] = pphi[i]; // internal friction angle of mixture or PHASE 1
                            else sflow.PHI[0] = sflow.PHI0[0];
                        } else sflow.PHI[0] = sflow.PHI0[0];

                        if ( sico.PHI2 == 1 ) {
                            if ( pphi2[i] != sico.UNDEF ) sflow.PHI[1] = pphi2[i]; // internal friction angle of PHASE 2
                            else sflow.PHI[1] = sflow.PHI0[1];
                        } else sflow.PHI[1] = sflow.PHI0[1];

                        if ( sico.PHI3 == 1 ) {
                            if ( pphi3[i] != sico.UNDEF ) sflow.PHI[2] = pphi3[i]; // internal friction angle of PHASE 3
                            else sflow.PHI[2] = sflow.PHI0[2];
                        } else sflow.PHI[2] = sflow.PHI0[2];

                        if ( sico.DELTA == 1 ) {
                            if ( pdelta[i] != sico.UNDEF ) sflow.DELTA[0] = pdelta[i]; // basal friction angle of mixture or PHASE 1
                            else sflow.DELTA[0] = sflow.DELTA0[0];
                        } else sflow.DELTA[0] = sflow.DELTA0[0];

                        if ( sico.DELTA2 == 1 ) {
                            if ( pdelta2[i] != sico.UNDEF ) sflow.DELTA[1] = pdelta2[i]; // basal friction angle of PHASE 2
                            else sflow.DELTA[1] = sflow.DELTA0[1];
                        } else sflow.DELTA[1] = sflow.DELTA0[1];

                        if ( sico.DELTA3 == 1 ) {
                            if ( pdelta3[i] != sico.UNDEF ) sflow.DELTA[2] = pdelta3[i]; // basal friction angle of PHASE 3
                            else sflow.DELTA[2] = sflow.DELTA0[2];
                        } else sflow.DELTA[2] = sflow.DELTA0[2];

                        if ( sico.ADDFRI == 1 ) {
                            if ( paddfri[i] != sico.UNDEF ) sflow.ADDFRI[0] = paddfri[i]; // additional friction parameter of PHASE 1
                            else sflow.ADDFRI[0] = sflow.ADDFRI0[0];
                        } else sflow.ADDFRI[0] = sflow.ADDFRI0[0];
                        
                        if ( sico.ADDFRI2 == 1 ) {
                            if ( paddfri2[i] != sico.UNDEF ) sflow.ADDFRI[1] = paddfri2[i]; // additional friction parameter of PHASE 2
                            else sflow.ADDFRI[1] = sflow.ADDFRI0[1];
                        } else sflow.ADDFRI[1] = sflow.ADDFRI0[1];

                        if ( sico.ADDFRI3 == 1 ) {
                            if ( paddfri3[i] != sico.UNDEF ) sflow.ADDFRI[2] = paddfri3[i]; // additional friction parameter of PHASE 3
                            else sflow.ADDFRI[2] = sflow.ADDFRI0[2];
                        } else sflow.ADDFRI[2] = sflow.ADDFRI0[2];

                        if ( sico.COHESION == 1 ) {
                            if ( pcohesion[i] != sico.UNDEF ) sflow.COHESION[0] = pcohesion[i]; // cohesion of PHASE 1
                            else sflow.COHESION[0] = sflow.COHESION0[0];
                        } else sflow.COHESION[0] = sflow.COHESION0[0];
                        
                        if ( sico.COHESION2 == 1 ) {
                            if ( pcohesion2[i] != sico.UNDEF ) sflow.COHESION[1] = pcohesion2[i]; // cohesion of PHASE 2
                            else sflow.COHESION[1] = sflow.COHESION0[1];
                        } else sflow.COHESION[1] = sflow.COHESION0[1];

                        if ( sico.COHESION3 == 1 ) {
                            if ( pcohesion3[i] != sico.UNDEF ) sflow.COHESION[2] = pcohesion3[i]; // cohesion of PHASE 3
                            else sflow.COHESION[2] = sflow.COHESION0[2];
                        } else sflow.COHESION[2] = sflow.COHESION0[2];

                        if ( sico.NY == 1 ) {
                            if ( pny[i] != sico.UNDEF ) sflow.NY[0] = pny[i]; // viscosity of mixture or PHASE 1
                            else sflow.NY[0] = sflow.NY0[0];
                        } else sflow.NY[0] = sflow.NY0[0];

                        if ( sico.NY2 == 1 ) {
                            if ( pny2[i] != sico.UNDEF ) sflow.NY[1] = pny2[i]; // viscosity of PHASE 2
                            else sflow.NY[1] = sflow.NY0[1];
                        } else sflow.NY[1] = sflow.NY0[1];

                        if ( sico.NY3 == 1 ) {
                            if ( pny3[i] != sico.UNDEF ) sflow.NY[2] = pny3[i]; // viscosity of PHASE 3
                            else sflow.NY[2] = sflow.NY0[2];
                        } else sflow.NY[2] = sflow.NY0[2];

                        if ( sico.CDEFORM == 1 ) {
                            if ( pcdeform[i] != sico.UNDEF ) sflow.CDEFORM[0] = pcdeform[i]; // deformation coefficient of PHASE 1
                            else sflow.CDEFORM[0] = sflow.CDEFORM0;
                        } else sflow.CDEFORM[0] = sflow.CDEFORM0;

                        if ( sico.AMBDRAG == 1 ) {
                            if ( pambdrag[i] != sico.UNDEF ) sflow.AMBDRAG = pambdrag[i]; // ambient drag coefficient
                            else sflow.AMBDRAG = sflow.AMBDRAG0;
                        } else sflow.AMBDRAG = sflow.AMBDRAG0;

                        if ( sico.ZONES == 1 ) {
                            if ( pzones[i] != sico.UNDEF ) sico.LAYERS = pzones[i]; // layer mode
                            else sico.LAYERS = sico.LAYERS0;
                        } else sico.LAYERS = sico.LAYERS0;

                        if ( frictiograph == 1 ) {

                            for ( frik=0; frik<=fritmax; frik++ ) { // identifying relevant line of input frictiograph
                                if ( frifri[frik][0] <= tsum ) frit=frik;
                                else break;
                            }

                            sflow.PHI[0] += frifri[frit][1] * sico.PI / 180; 
                            sflow.DELTA[0] += frifri[frit][2] * sico.PI / 180;
                            sflow.NY[0] += pow(10, frifri[frit][3]); 
                            sflow.CDEFORM[0] = frifri[frit][4];

                            if ( sico.MODEL == 7 ) { 

                                sflow.PHI[1] += frifri[frit][5] * sico.PI / 180;
                                sflow.DELTA[1] += frifri[frit][6] * sico.PI / 180;
                                sflow.NY[1] += pow(10, frifri[frit][7]); 
                                sflow.CDEFORM[1] = frifri[frit][8];                                
                                sflow.PHI[2] += frifri[frit][9] * sico.PI / 180; 
                                sflow.DELTA[2] += frifri[frit][10] * sico.PI / 180;
                                sflow.NY[2] += pow(10, frifri[frit][11]); 
                                sflow.CDEFORM[2] = frifri[frit][12];
                            }
                        }


                        // Components of gravity

                        grav[0] = sico.GRAVITY * cos( betax[i] );
                        grav[1] = sico.GRAVITY * cos( betay[i] );
                        grav[2] = sico.GRAVITY * sin( betax[i] );
                        grav[3] = sico.GRAVITY * sin( betay[i] );
                        grav[4] = sico.GRAVITY * cos( betaxy[i] );
                        grav[5] = sico.GRAVITY * sin( betaxh[i] );
                        grav[6] = sico.GRAVITY * sin( betayh[i] );

                        if ( sico.LAYERS > 0 ) {
                            
                            grav[7] = sico.GRAVITY * sin( betax2[i] );
                            grav[8] = sico.GRAVITY * sin( betay2[i] );                            
                            grav[9] = sico.GRAVITY * sin( betaxh2[i] );
                            grav[10] = sico.GRAVITY * sin( betayh2[i] );
                            grav[11] = sico.GRAVITY * sin( betax3[i] );
                            grav[12] = sico.GRAVITY * sin( betay3[i] );
                            grav[13] = sico.GRAVITY * sin( betaxh3[i] );
                            grav[14] = sico.GRAVITY * sin( betayh3[i] );
                            grav[15] = sico.GRAVITY * cos( betax2[i] );
                            grav[16] = sico.GRAVITY * cos( betay2[i] );                            
                            grav[17] = sico.GRAVITY * cos( betax3[i] );
                            grav[18] = sico.GRAVITY * cos( betay3[i] );
                            grav[19] = sico.GRAVITY * cos( betaxy2[i] );
                            grav[20] = sico.GRAVITY * cos( betaxy3[i] );
                        }


                        // Flow velocities in x and y directions

                        vflowx = fdiv( aw[i][1], aw[i][0], sico.HFLOWMIN );
                        vflowy = fdiv( aw[i][2], aw[i][0], sico.HFLOWMIN );

                        if ( sico.MODEL == 7 ) vflowx2 = fdiv( aw[i][4], aw[i][3], sico.HFLOWMIN ); else vflowx2 = 0;
                        if ( sico.MODEL == 7 ) vflowy2 = fdiv( aw[i][5], aw[i][3], sico.HFLOWMIN ); else vflowy2 = 0;

                        if ( sico.MODEL == 7 ) vflowx3 = fdiv( aw[i][7], aw[i][6], sico.HFLOWMIN ); else vflowx3 = 0;
                        if ( sico.MODEL == 7 ) vflowy3 = fdiv( aw[i][8], aw[i][6], sico.HFLOWMIN ); else vflowy3 = 0;
                        
                        
                        // Total flow depths, velocities, fractions, and slopes of neighbouring cells

                        welev[0] = pelev[i];

                        for ( l=0; l<9; l++ ) {

                            welev[l] = pelev[in[i][l]];

                            nbetax[l] = betax[in[i][l]];
                            nbetay[l] = betay[in[i][l]];

                            if ( sico.MODEL <= 3 ) wh[l] = aw[in[i][l]][0];
                            else if ( sico.MODEL == 7 ) {
                            
                                wh[l] = aw[in[i][l]][0] + aw[in[i][l]][3] + aw[in[i][l]][6];
                                
                                if ( sico.LAYERS > 0 ) {
                                
                                    wh1[l] = aw[in[i][l]][0];
                                    wh2[l] = aw[in[i][l]][0] + aw[in[i][l]][3];
                                }                            
                            }

                            wu[l] = fdiv( aw[in[i][l]][1], aw[in[i][l]][0], sico.HFLOWMIN );
                            wv[l] = fdiv( aw[in[i][l]][2], aw[in[i][l]][0], sico.HFLOWMIN );

                            if ( sico.MODEL == 7 ) wu2[l] = fdiv( aw[in[i][l]][4], aw[in[i][l]][3], sico.HFLOWMIN ); else wu2[l] = 0;
                            if ( sico.MODEL == 7 ) wv2[l] = fdiv( aw[in[i][l]][5], aw[in[i][l]][3], sico.HFLOWMIN ); else wv2[l] = 0;

                            if ( sico.MODEL == 7 ) wu3[l] = fdiv( aw[in[i][l]][7], aw[in[i][l]][6], sico.HFLOWMIN ); else wu3[l] = 0;
                            if ( sico.MODEL == 7 ) wv3[l] = fdiv( aw[in[i][l]][8], aw[in[i][l]][6], sico.HFLOWMIN ); else wv3[l] = 0;

                            walpha[l] = fdiv( aw[in[i][l]][0], wh[l], sico.HFLOWMIN );
                            if ( sico.MODEL == 7 ) walpha2[l] = fdiv( aw[in[i][l]][3], wh[l], sico.HFLOWMIN ); else walpha2[l] = 0;
                            if ( sico.MODEL == 7 ) walpha3[l] = fdiv( aw[in[i][l]][6], wh[l], sico.HFLOWMIN ); else walpha3[l] = 0;
                        }


                        // Gradients of total flow depth, flow velocities, and fractions

                        wdx = dx[i];
                        wdy = dy[i];

                        cin = f0noosc( wh, sico );

                        whx = ( wh[cin[0]] - wh[cin[1]] ) / ((float)(cin[4]) * dx[i] );
                        why = ( wh[cin[2]] - wh[cin[3]] ) / ((float)(cin[5]) * dy[i] );

                        if ( sico.LAYERS > 0 ) {
                            
                            whx1 = ( wh1[cin[0]] - wh1[cin[1]] ) / ((float)(cin[4]) * dx[i] );
                            whx2 = ( wh2[cin[0]] - wh2[cin[1]] ) / ((float)(cin[4]) * dx[i] );
                            why1 = ( wh1[cin[2]] - wh1[cin[3]] ) / ((float)(cin[5]) * dy[i] );
                            why2 = ( wh2[cin[2]] - wh2[cin[3]] ) / ((float)(cin[5]) * dy[i] );
                        }
                        
                        free( cin );

                        wdu[0] = ( wu[2] - wu[5] ) / ( 2 * wdx );
                        wdv[0] = ( wv[1] - wv[4] ) / ( 2 * wdy );

                        if ( sico.MODEL == 7 ) wdu[1] = ( wu2[2] - wu2[5] ) / ( 2 * wdx ); else wdu[1] = 0;
                        if ( sico.MODEL == 7 ) wdv[1] = ( wv2[1] - wv2[4] ) / ( 2 * wdy ); else wdv[1] = 0;

                        if ( sico.MODEL == 7 ) wdu[2] = ( wu3[2] - wu3[5] ) / ( 2 * wdx ); else wdu[2] = 0;
                        if ( sico.MODEL == 7 ) wdv[2] = ( wv3[1] - wv3[4] ) / ( 2 * wdy ); else wdv[2] = 0;

                        if ( walpha[2] ) xwdv[0] = ( wv[2] - wv[5] ) / ( 2 * wdx );
                        if ( walpha[1] ) xwdu[0] = ( wu[1] - wu[4] ) / ( 2 * wdy );

                        if ( sico.MODEL == 7 ) xwdv[1] = ( wv2[2] - wv2[5] ) / ( 2 * wdx ); else xwdv[1] = 0;
                        if ( sico.MODEL == 7 ) xwdu[1] = ( wu2[1] - wu2[4] ) / ( 2 * wdy ); else xwdu[1] = 0;

                        if ( sico.MODEL == 7 ) xwdv[2] = ( wv3[2] - wv3[5] ) / ( 2 * wdx ); else xwdv[2] = 0;
                        if ( sico.MODEL == 7 ) xwdu[2] = ( wu3[1] - wu3[4] ) / ( 2 * wdy ); else xwdu[2] = 0;

                        walphax = ( walpha[2] - walpha[5] ) / ( 2 * wdx );
                        walphay = ( walpha[1] - walpha[4] ) / ( 2 * wdy );

                        if ( sico.MODEL == 7 ) walphax2 = ( walpha2[2] - walpha2[5] ) / ( 2 * wdx ); else walphax2 = 0;
                        if ( sico.MODEL == 7 ) walphay2 = ( walpha2[1] - walpha2[4] ) / ( 2 * wdy ); else walphay2 = 0;

                        if ( sico.MODEL == 7 ) walphax3 = ( walpha3[2] - walpha3[5] ) / ( 2 * wdx ); else walphax3 = 0;
                        if ( sico.MODEL == 7 ) walphay3 = ( walpha3[1] - walpha3[4] ) / ( 2 * wdy ); else walphay3 = 0;


                        // Earth pressure coefficients

                        if ( sico.MODEL <= 3 ) { whflow = aw[i][0], whflow2 = 0, whflow3 = 0; }
                        else if ( sico.MODEL == 7 ) { whflow = aw[i][0], whflow2 = aw[i][3], whflow3 = aw[i][6]; }

                        hekin = 0.5 * whflow * sflow.RHO1 * ( pow( vflowx, 2 ) + pow( vflowy, 2 ));
                        if ( sico.MODEL == 7 ) hekin += 0.5 * whflow2 * sflow.RHO2 * ( pow( vflowx2, 2 ) + pow( vflowy2, 2 ));
                        if ( sico.MODEL == 7 ) hekin += 0.5 * whflow3 * sflow.RHO3 * ( pow( vflowx3, 2 ) + pow( vflowy3, 2 ));

                        for ( l=0; l<sico.PMAX; l++ ) {

                            if ( sflow.DELTA[l] > 0 ) {

                                gkx[l] = fk( whflow+whflow2+whflow3, aw[i][3*l], vflowx, vflowy, wdu[l], hekin, l, sico, sflow );
                                gky[l] = fk( whflow+whflow2+whflow3, aw[i][3*l], vflowx, vflowy, wdv[l], hekin, l, sico, sflow );
                            }
                        }


                        // Curvature, flux, source, and deceleration terms

                        if ( sico.ZFRAG == 1 ) ctrl_zfrag = pzfrag[i]; else ctrl_zfrag = 0; 

                        for ( l=0; l<sico.PMAX; l++ ) {

                            kappau[l] = ( betax[in[i][2]] - betax[in[i][4]] ) * 0.5 * cos ( betax[i] ) / sico.CSZ;  
                            kappau[l+3] = ( betay[in[i][1]] - betay[in[i][5]] ) * 0.5 * cos ( betay[i] ) / sico.CSZ; 
                        }

                        if ( sico.MODEL <= 3 ) wtemp = aw[i][3]; else wtemp = aw[i][9];

                        vm = fvm( whflow, whflow2, whflow3, vflowx, vflowx2, vflowx3, vflowy, vflowy2, vflowy3, sico, sflow );
                        cdrag = fdrag( whflow, whflow2, whflow3, sico, sflow );
                        gze = fgze( whflow, whflow2, whflow3, vflowx, vflowx2, vflowx3, vflowy, vflowy2, vflowy3, wh, whx, why,
                            wdu, wdv, xwdu, xwdv, wu, wv, wu2, wv2, wu3, wv3, walpha, walpha2, walpha3, wdx, wdy,
                            grav, betax[i], betay[i], cdrag, tsum, sico, sflow );
                        disp = fdisp( whflow, whflow2, whflow3, vflowx, vflowx2, vflowx3, vflowy, vflowy2, vflowy3, wdu, wdv, xwdu, xwdv, 
                            betax[i], betay[i], betaxy[i], gkx, gky, vm, cdrag, tsum, sico, sflow );

                        gf = ff( wh, whflow, whflow2, whflow3, vflowx, vflowx2, vflowx3, vflowy, vflowy2, vflowy3,
                            gkx, gze, vm, whx, walphax, cdrag, grav, disp, betax[i], betay[i], betaxy[i], ctrl_zfrag, wtemp, sflow, sico );
                        gg = fg( wh, whflow, whflow2, whflow3, vflowx, vflowx2, vflowx3, vflowy, vflowy2, vflowy3,
                            gky, gze, vm, why, walphay, cdrag, grav, disp, betax[i], betay[i], betaxy[i], ctrl_zfrag, wtemp, sflow, sico );
                        gs = fs( wh, whflow, whflow2, whflow3, vflowx, vflowx2, vflowx3, vflowy, vflowy2, vflowy3,
                            whx, why, whx1, why1, whx2, why2, walphax, walphay, walphax2, walphay2, walphax3, walphay3, grav, gze, betax[i], betay[i], betaxy[i], 
                            cdrag, gkx, gky, pxslide[i], angx, angy, tsum, ptslide[i], ctrl_zfrag, wtemp, sico, sflow );

                        if ( vflowx == 0 ) wu[0] = fdiv( 0.5 * gs[1], whflow, sico.HFLOWMIN );
                        if ( vflowy == 0 ) wv[0] = fdiv( 0.5 * gs[2], whflow, sico.HFLOWMIN );
                        
                        if ( sico.MODEL == 7 ) {
                            
                            if ( vflowx2 == 0 ) wu2[0] = fdiv( 0.5 * gs[4], whflow2, sico.HFLOWMIN );
                            if ( vflowy2 == 0 ) wv2[0] = fdiv( 0.5 * gs[5], whflow2, sico.HFLOWMIN );
                            if ( vflowx3 == 0 ) wu3[0] = fdiv( 0.5 * gs[7], whflow3, sico.HFLOWMIN );
                            if ( vflowy3 == 0 ) wv3[0] = fdiv( 0.5 * gs[8], whflow3, sico.HFLOWMIN );
                        }
                        
                        sico.XDIST = fabs(sico.XREL - px[i]) * sico.CSZ;
                        sico.YDIST = fabs(sico.YREL - py[i]) * sico.CSZ;
                        
                        gdecel = fd( wh, whx, why, wu, wv, wu2, wv2, wu3, wv3, walpha, walpha2, walpha3, betax[i], betay[i], betaxy[i], nbetax, nbetay, gze, 
                            kappau, hekin, grav, angx, angy, pxslide[i], tsum, ptslide[i], ctrl_zfrag, sico, sflow );

                        for ( k=0; k<sico.NVECTMIN; k++ ) { 
                            af[i][k] = gf[k]; ag[i][k] = gg[k]; as[i][k] = gs[k];
                        }
                        for ( k=0; k<sico.NVECTMIN; k++ ) ad[i][k] = gdecel[k];

                        free( vm ); free( cdrag ); free( gze ); free( disp ); free( gf ); free( gg ); free( gs ); free( gdecel );


                        // Slopes of the vector components                        

                        asigma_xelev[i] = fsigma( pelev[i], pelev[in[i][5]], pelev[in[i][2]], sico.CSZ, sico );
                        asigma_yelev[i] = fsigma( pelev[i], pelev[in[i][4]], pelev[in[i][1]], sico.CSZ, sico );

                        for (k=0; k<sico.NVECTMIN; k++) {

                            if ( k == 0 || k == 3 || k == 6 ) {

                                asigma_x[i][k] = fsigma( aw[i][k], aw[in[i][5]][k], aw[in[i][2]][k], sico.CSZ, sico );
                                asigma_y[i][k] = fsigma( aw[i][k], aw[in[i][4]][k], aw[in[i][1]][k], sico.CSZ, sico );
                        
                            } else {
                        
                                asigma_x[i][k] = fsigma( aw[i][k], aw[in[i][5]][k], aw[in[i][2]][k], sico.CSZ, sico );
                                asigma_y[i][k] = fsigma( aw[i][k], aw[in[i][4]][k], aw[in[i][1]][k], sico.CSZ, sico );
                            }                        
                        }
                    }
                }


                // Slopes of the fluxes

                for ( ix=0; ix<ib[0]; ix++ ) {

                    i = ibasket[0][ix];

                    if ( cdomain[i] <= 3 ) {

                        for ( k=0; k<sico.NVECTMIN; k++) { asigma_f[i][k]=0; asigma_g[i][k]=0; }
                         
                    } else {

                        for ( k=0; k<sico.NVECTMIN; k++ ) {

                            asigma_f[i][k] = fsigma( af[i][k], af[in[i][5]][k], af[in[i][2]][k], sico.CSZ, sico );
                            asigma_g[i][k] = fsigma( ag[i][k], ag[in[i][4]][k], ag[in[i][1]][k], sico.CSZ, sico );
                        }
                    }
                }


                // Values of vector at quarter of cell after half time step

                for ( ix=0; ix<ib[0]; ix++ ) {

                    i = ibasket[0][ix];

                    for ( j=0; j<4; j++ ) { 
                    
                        if ( cdomain[i] > 3 ) wintelev[i][j]=pelev[in[i][j]]+fj[j][0]*0.25*dx[i]*asigma_xelev[in[i][j]]+fj[j][1]*0.25*dy[i]*asigma_yelev[in[i][j]];
                        else wintelev[i][j] = pelev[i];
                    }
                }

                for ( ix=0; ix<ib[0]; ix++ ) {

                    i = ibasket[0][ix];

                    if ( cdomain[i] > 3 ) {

                        for ( j=0; j<6; j++ ) { // loop over adjacent cells

                            for ( k=0; k<sico.NVECTMIN; k++ ) { // loop over all components of the vector

                                dw_dt[k]=-asigma_f[in[i][j]][k]-asigma_g[in[i][j]][k]+as[in[i][j]][k]; // gradient of vector
                                dw_dttest[k]=-asigma_f[in[i][j]][k]-asigma_g[in[i][j]][k]+(as[in[i][j]][k]-ad[in[i][j]][k]); // gradient of vector
                            }

                            for ( k=0; k<sico.NVECTMIN; k++ ) { // loop over all components of the vector

                                if ( j < 4 ) winta[i][j][k]=aw[in[i][j]][k]+fj[j][0]*0.25*dx[i]*asigma_x[in[i][j]][k]+fj[j][1]
                                    *0.25*dy[i]*asigma_y[in[i][j]][k]; // value at quarter of cell

                                wintb[i][j][k]=aw[in[i][j]][k]+0.5*tlength*dw_dt[k];
                                wintbtest=aw[in[i][j]][k]+0.5*tlength*dw_dttest[k];

                                if (( fsign(wintb[i][j][k]) == fsign(wintbtest) && fabs(wintb[i][j][k]) > fabs(wintbtest)) || sico.ORIGINAL == 1 ) wintb[i][j][k] = wintbtest;
                                else if ( fsign(wintb[i][j][k]) == fsign(wintbtest)) wintb[i][j][k] = wintb[i][j][k];
                                else if ( k==0 || k==3 || k==6 ) wintb[i][j][k] = aw[in[i][j]][k];
                                else wintb[i][j][k] = 0;

                                if ( j < 4 ) {

                                    wintc[i][j][k]=aw[in[i][j]][k]+0.5*tlength*dw_dt[k]+fj[j][0]*0.25*dx[i]*asigma_x[in[i][j]][k]
                                    +fj[j][1]*0.25*dy[i]*asigma_y[in[i][j]][k];

                                    wintctest=aw[in[i][j]][k]+0.5*tlength*dw_dttest[k]+fj[j][0]*0.25*dx[i]*asigma_x[in[i][j]][k]
                                    +fj[j][1]*0.25*dy[i]*asigma_y[in[i][j]][k];

                                    if (( fsign(wintc[i][j][k]) == fsign(wintctest) && fabs(wintc[i][j][k]) > fabs(wintctest)) || sico.ORIGINAL==1 ) wintc[i][j][k] = wintctest;
                                    else if ( fsign(wintc[i][j][k]) == fsign(wintctest)) wintc[i][j][k] = wintc[i][j][k];
                                    else if ( k==0 || k==3 || k==6 ) wintc[i][j][k] = aw[in[i][j]][k];
                                    else wintc[i][j][k] = 0;

                                }
                            }
                        }

                    } else {

                        for ( j=0; j<4; j++ ) {
                            for ( k=0; k<sico.NVECTMIN; k++ ) {

                                if ( j < 4 ) winta[i][j][k]=0; //!!!CHECK j
                                wintb[i][j][k]=0; 
                                if ( j < 4 ) wintc[i][j][k]=0;
                            }
                        }
                    }
                }


                // Fluxes and source terms (shifted coordinate system)

                vcelr = 0;

                for ( ix=0; ix<ib[0]; ix++ ) {

                    i = ibasket[0][ix];

                    if ( cdomain[i] > 3 ) {

                        if ( sico.PHI == 1 ) {
                            if ( pphi[i] != sico.UNDEF ) sflow.PHI[0] = pphi[i]; // internal friction angle of mixture or PHASE 1
                            else sflow.PHI[0] = sflow.PHI0[0];
                        } else sflow.PHI[0] = sflow.PHI0[0];

                        if ( sico.PHI2 == 1 ) {
                            if ( pphi2[i] != sico.UNDEF ) sflow.PHI[1] = pphi2[i]; // internal friction angle of PHASE 2
                            else sflow.PHI[1] = sflow.PHI0[1];
                        } else sflow.PHI[1] = sflow.PHI0[1];

                        if ( sico.PHI3 == 1 ) {
                            if ( pphi3[i] != sico.UNDEF ) sflow.PHI[2] = pphi3[i]; // internal friction angle of PHASE 3
                            else sflow.PHI[2] = sflow.PHI0[2];
                        } else sflow.PHI[2] = sflow.PHI0[2];

                        if ( sico.DELTA == 1 ) {
                            if ( pdelta[i] != sico.UNDEF ) sflow.DELTA[0] = pdelta[i]; // basal friction angle of mixture or PHASE 1
                            else sflow.DELTA[0] = sflow.DELTA0[0];
                        } else sflow.DELTA[0] = sflow.DELTA0[0];

                        if ( sico.DELTA2 == 1 ) {
                            if ( pdelta2[i] != sico.UNDEF ) sflow.DELTA[1] = pdelta2[i]; // basal friction angle of PHASE 2
                            else sflow.DELTA[1] = sflow.DELTA0[1];
                        } else sflow.DELTA[1] = sflow.DELTA0[1];

                        if ( sico.DELTA3 == 1 ) {
                            if ( pdelta3[i] != sico.UNDEF ) sflow.DELTA[2] = pdelta3[i]; // basal friction angle of PHASE 3
                            else sflow.DELTA[2] = sflow.DELTA0[2];
                        } else sflow.DELTA[2] = sflow.DELTA0[2];

                        if ( sico.ADDFRI == 1 ) {
                            if ( paddfri[i] != sico.UNDEF ) sflow.ADDFRI[0] = paddfri[i]; // additional friction parameter of PHASE 1
                            else sflow.ADDFRI[0] = sflow.ADDFRI0[0];
                        } else sflow.ADDFRI[0] = sflow.ADDFRI0[0];
                        
                        if ( sico.ADDFRI2 == 1 ) {
                            if ( paddfri2[i] != sico.UNDEF ) sflow.ADDFRI[1] = paddfri2[i]; // additional friction parameter of PHASE 2
                            else sflow.ADDFRI[1] = sflow.ADDFRI0[1];
                        } else sflow.ADDFRI[1] = sflow.ADDFRI0[1];

                        if ( sico.ADDFRI3 == 1 ) {
                            if ( paddfri3[i] != sico.UNDEF ) sflow.ADDFRI[2] = paddfri3[i]; // additional friction parameter of PHASE 3
                            else sflow.ADDFRI[2] = sflow.ADDFRI0[2];
                        } else sflow.ADDFRI[2] = sflow.ADDFRI0[2];

                        if ( sico.COHESION == 1 ) {
                            if ( pcohesion[i] != sico.UNDEF ) sflow.COHESION[0] = pcohesion[i]; // cohesion of PHASE 1
                            else sflow.COHESION[0] = sflow.COHESION0[0];
                        } else sflow.COHESION[0] = sflow.COHESION0[0];
                        
                        if ( sico.COHESION2 == 1 ) {
                            if ( pcohesion2[i] != sico.UNDEF ) sflow.COHESION[1] = pcohesion2[i]; // cohesion of PHASE 2
                            else sflow.COHESION[1] = sflow.COHESION0[1];
                        } else sflow.COHESION[1] = sflow.COHESION0[1];

                        if ( sico.COHESION3 == 1 ) {
                            if ( pcohesion3[i] != sico.UNDEF ) sflow.COHESION[2] = pcohesion3[i]; // cohesion of PHASE 3
                            else sflow.COHESION[2] = sflow.COHESION0[2];
                        } else sflow.COHESION[2] = sflow.COHESION0[2];

                        if ( sico.NY == 1 ) {
                            if ( pny[i] != sico.UNDEF ) sflow.NY[0] = pny[i]; // viscosity of PHASE 1
                            else sflow.NY[0] = sflow.NY0[0];
                        } else sflow.NY[0] = sflow.NY0[0];

                        if ( sico.NY2 == 1 ) {
                            if ( pny2[i] != sico.UNDEF ) sflow.NY[1] = pny2[i]; // viscosity of PHASE 2
                            else sflow.NY[1] = sflow.NY0[1];
                        } else sflow.NY[1] = sflow.NY0[1];

                        if ( sico.NY3 == 1 ) {
                            if ( pny3[i] != sico.UNDEF ) sflow.NY[2] = pny3[i]; // viscosity of PHASE 3
                            else sflow.NY[2] = sflow.NY0[2];
                        } else sflow.NY[2] = sflow.NY0[2];

                        if ( sico.CDEFORM == 1 ) {
                            if ( pcdeform[i] != sico.UNDEF ) sflow.CDEFORM[0] = pcdeform[i]; // deformation coefficient of PHASE 1
                            else sflow.CDEFORM[0] = sflow.CDEFORM0;
                        } else sflow.CDEFORM[0] = sflow.CDEFORM0;

                        if ( sico.AMBDRAG == 1 ) {
                            if ( pambdrag[i] != sico.UNDEF ) sflow.AMBDRAG = pambdrag[i]; // ambient drag coefficient
                            else sflow.AMBDRAG = sflow.AMBDRAG0;
                        } else sflow.AMBDRAG = sflow.AMBDRAG0;

                        if ( sico.ZONES == 1 ) {
                            if ( pzones[i] != sico.UNDEF ) sico.LAYERS = pzones[i]; // layer mode
                            else sico.LAYERS = sico.LAYERS0;
                        } else sico.LAYERS = sico.LAYERS0;

                        if ( frictiograph == 1 ) {

                            for ( frik=0; frik<=fritmax; frik++ ) { // identifying relevant line of input frictiograph
                                if ( frifri[frik][0] <= tsum ) frit=frik;
                                else break;
                            }

                            sflow.PHI[0] += frifri[frit][1] * sico.PI / 180; 
                            sflow.DELTA[0] += frifri[frit][2] * sico.PI / 180;
                            sflow.NY[0] += pow(10, frifri[frit][3]); 
                            sflow.CDEFORM[0] = frifri[frit][4];

                            if ( sico.MODEL == 7 ) { 
                            
                                sflow.PHI[1] += frifri[frit][5] * sico.PI / 180;
                                sflow.DELTA[1] += frifri[frit][6] * sico.PI / 180;
                                sflow.NY[1] += pow(10, frifri[frit][7]); 
                                sflow.CDEFORM[1] = frifri[frit][8];                                
                                sflow.PHI[2] += frifri[frit][9] * sico.PI / 180; 
                                sflow.DELTA[2] += frifri[frit][10] * sico.PI / 180;
                                sflow.NY[2] += pow(10, frifri[frit][11]); 
                                sflow.CDEFORM[2] = frifri[frit][12];
                            }
                        }

                        for ( j=0; j<4; j++ ) {


                            // Slopes, topography-following cell sizes, and gravity components at quarter of cell

                            if ( cdomain2[i] == 1 ) {

                                wbetax = betax[i];
                                wbetay = betay[i];
                                wbetaxh = betaxh[i];
                                wbetayh = betayh[i];
                                wbetaxy = betaxy[i];
                                
                                if ( sico.LAYERS > 0 ) {
                            
                                    wbetax2 = betax2[i];
                                    wbetay2 = betay2[i];
                                    wbetaxy2 = betaxy2[i];                          
                                    wbetaxh2 = betaxh2[i];
                                    wbetayh2 = betayh2[i];
                                    wbetax3 = betax3[i];
                                    wbetay3 = betay3[i];
                                    wbetaxy3 = betaxy3[i];
                                    wbetaxh3 = betaxh3[i];
                                    wbetayh3 = betayh3[i];                         
                                }

                            } else {

                                cin = fwnoosc( wintc, in, i, j, sico );

                                wbetax = fbeta( wintelev[cin[0]][j], wintelev[cin[1]][j], (float)(cin[4]), sico );
                                wbetay = fbeta( wintelev[cin[2]][j], wintelev[cin[3]][j], (float)(cin[5]), sico );
                                wbetaxh = fbeta( wintelev[cin[0]][j]+wintc[cin[0]][j][0], wintelev[cin[1]][j]+wintc[cin[1]][j][0], (float)(cin[4]), sico );
                                wbetayh = fbeta( wintelev[cin[2]][j]+wintc[cin[2]][j][0], wintelev[cin[3]][j]+wintc[cin[3]][j][0], (float)(cin[5]), sico );
                                wbetaxy = fbetaxy( wbetax, wbetay );

                                if ( sico.LAYERS > 0 ) {
                            
                                    wbetax2 = fbeta( wintelev[cin[0]][j]+wintc[cin[0]][j][0], wintelev[cin[1]][j]+wintc[cin[1]][j][0], (float)(cin[4]), sico );
                                    wbetay2 = fbeta( wintelev[cin[2]][j]+wintc[cin[2]][j][0], wintelev[cin[3]][j]+wintc[cin[3]][j][0], (float)(cin[5]), sico );

                                    wbetaxy2 = fbetaxy( wbetax2, wbetay2 );
                           
                                    wbetaxh2 = fbeta( wintelev[cin[0]][j]+wintc[cin[0]][j][0]+wintc[cin[0]][j][3], 
                                        wintelev[cin[1]][j]+wintc[cin[1]][j][0]+wintc[cin[1]][j][3], (float)(cin[4]), sico );
                                    wbetayh2 = fbeta( wintelev[cin[2]][j]+wintc[cin[2]][j][0]+wintc[cin[2]][j][3], 
                                        wintelev[cin[3]][j]+wintc[cin[3]][j][0]+wintc[cin[3]][j][3], (float)(cin[5]), sico );
                                    
                                    wbetax3 = fbeta( wintelev[cin[0]][j]+wintc[cin[0]][j][0]+wintc[cin[0]][j][3], 
                                        wintelev[cin[1]][j]+wintc[cin[1]][j][0]+wintc[cin[1]][j][3], (float)(cin[4]), sico );
                                    wbetay3 = fbeta( wintelev[cin[2]][j]+wintc[cin[2]][j][0]+wintc[cin[2]][j][3], 
                                        wintelev[cin[3]][j]+wintc[cin[3]][j][0]+wintc[cin[3]][j][3], (float)(cin[5]), sico );

                                    wbetaxy3 = fbetaxy( wbetax3, wbetay3 );
                                    
                                    wbetaxh3 = fbeta( wintelev[cin[0]][j]+wintc[cin[0]][j][0]+wintc[cin[0]][j][3]+wintc[cin[0]][j][6], 
                                        wintelev[cin[1]][j]+wintc[cin[1]][j][0]+wintc[cin[1]][j][3]+wintc[cin[1]][j][6], (float)(cin[4]), sico );
                                    wbetayh3 = fbeta( wintelev[cin[2]][j]+wintc[cin[2]][j][0]+wintc[cin[2]][j][3]+wintc[cin[2]][j][6], 
                                        wintelev[cin[3]][j]+wintc[cin[3]][j][0]+wintc[cin[3]][j][3]+wintc[cin[3]][j][6], (float)(cin[5]), sico );                       
                                }
                                
                                free( cin );
                            }

                            wdx = dx[i];
                            wdy = dy[i];

                            wgrav[0] = sico.GRAVITY * cos( wbetax );
                            wgrav[1] = sico.GRAVITY * cos( wbetay );
                            wgrav[2] = sico.GRAVITY * sin( wbetax );
                            wgrav[3] = sico.GRAVITY * sin( wbetay );
                            wgrav[4] = sico.GRAVITY * cos( wbetaxy );
                            wgrav[5] = sico.GRAVITY * sin( wbetaxh );
                            wgrav[6] = sico.GRAVITY * sin( wbetayh );
                            
                            if ( sico.LAYERS > 0 ) {
                            
                                wgrav[7] = sico.GRAVITY * sin( wbetax2 );
                                wgrav[8] = sico.GRAVITY * sin( wbetay2 );                            
                                wgrav[9] = sico.GRAVITY * sin( wbetaxh2 );
                                wgrav[10] = sico.GRAVITY * sin( wbetayh2 );
                                wgrav[11] = sico.GRAVITY * sin( wbetax3 );
                                wgrav[12] = sico.GRAVITY * sin( wbetay3 );
                                wgrav[13] = sico.GRAVITY * sin( wbetaxh3 );
                                wgrav[14] = sico.GRAVITY * sin( wbetayh3 );
                                wgrav[15] = sico.GRAVITY * cos( wbetax2 );
                                wgrav[16] = sico.GRAVITY * cos( wbetay2 );                            
                                wgrav[17] = sico.GRAVITY * cos( wbetax3 );
                                wgrav[18] = sico.GRAVITY * cos( wbetay3 );
                                wgrav[19] = sico.GRAVITY * cos( wbetaxy2 );
                                wgrav[20] = sico.GRAVITY * cos( wbetaxy3 );                        
                            }

                            // Total flow depths

                            for ( l=1; l<9; l++ ) {

                                if ( sico.MODEL <= 3 ) wh[l] = wintb[in[i][l]][j][0];
                                else if ( sico.MODEL == 7 ) {
                                
                                    wh[l] = wintb[in[i][l]][j][0] + wintb[in[i][l]][j][3] + wintb[in[i][l]][j][6];
                                  
                                    if ( sico.LAYERS > 0 ) {
                                
                                        wh1[l] = wintb[in[i][l]][j][0];
                                        wh2[l] = wintb[in[i][l]][j][0] + wintb[in[i][l]][j][3];
                                    }
                                }
                            }


                            // Gradients of total flow depth

                            cin = f0noosc( wh, sico );

                            whx = ( wh[cin[0]] - wh[cin[1]] ) / ((float)(cin[4]) * dx[i] );
                            why = ( wh[cin[2]] - wh[cin[3]] ) / ((float)(cin[5]) * dy[i] );

                            if ( sico.LAYERS > 0 ) {
                            
                                whx1 = ( wh1[cin[0]] - wh1[cin[1]] ) / ((float)(cin[4]) * dx[i] );
                                whx2 = ( wh2[cin[0]] - wh2[cin[1]] ) / ((float)(cin[4]) * dx[i] );
                                why1 = ( wh1[cin[2]] - wh1[cin[3]] ) / ((float)(cin[5]) * dy[i] );
                                why2 = ( wh2[cin[2]] - wh2[cin[3]] ) / ((float)(cin[5]) * dy[i] );

                            }
                        
                            free( cin );


                            // Flow velocities in x and y directions

                            vflowx = fdiv( wintb[i][j][1], wintb[i][j][0], sico.HFLOWMIN );
                            vflowy = fdiv( wintb[i][j][2], wintb[i][j][0], sico.HFLOWMIN );

                            if ( sico.MODEL == 7 ) vflowx2 = fdiv( wintb[i][j][4], wintb[i][j][3], sico.HFLOWMIN ); else vflowx2 = 0;
                            if ( sico.MODEL == 7 ) vflowy2 = fdiv( wintb[i][j][5], wintb[i][j][3], sico.HFLOWMIN ); else vflowy2 = 0;

                            if ( sico.MODEL == 7 ) vflowx3 = fdiv( wintb[i][j][7], wintb[i][j][6], sico.HFLOWMIN ); else vflowx3 = 0;
                            if ( sico.MODEL == 7 ) vflowy3 = fdiv( wintb[i][j][8], wintb[i][j][6], sico.HFLOWMIN ); else vflowy3 = 0;


                            // Gradients of flow velocities

                            for ( l=1; l<9; l++ ) {

                                wu[l] = fdiv( wintb[in[i][l]][j][1], wintb[in[i][l]][j][0], sico.HFLOWMIN );
                                wv[l] = fdiv( wintb[in[i][l]][j][2], wintb[in[i][l]][j][0], sico.HFLOWMIN );

                                if ( sico.MODEL == 7 ) wu2[l] = fdiv( wintb[in[i][l]][j][4], wintb[in[i][l]][j][3], sico.HFLOWMIN ); else wu2[l] = 0;
                                if ( sico.MODEL == 7 ) wv2[l] = fdiv( wintb[in[i][l]][j][5], wintb[in[i][l]][j][3], sico.HFLOWMIN ); else wv2[l] = 0;

                                if ( sico.MODEL == 7 ) wu3[l] = fdiv( wintb[in[i][l]][j][7], wintb[in[i][l]][j][6], sico.HFLOWMIN ); else wu3[l] = 0;
                                if ( sico.MODEL == 7 ) wv3[l] = fdiv( wintb[in[i][l]][j][8], wintb[in[i][l]][j][6], sico.HFLOWMIN ); else wv3[l] = 0;

                                walpha[l] = fdiv( wintb[in[i][l]][j][0], wh[l], sico.HFLOWMIN );
                                if ( sico.MODEL == 7 ) walpha2[l] = fdiv( wintb[in[i][l]][j][3], wh[l], sico.HFLOWMIN ); else walpha2[l] = 0;
                                if ( sico.MODEL == 7 ) walpha3[l] = fdiv( wintb[in[i][l]][j][6], wh[l], sico.HFLOWMIN ); else walpha3[l] = 0;
                            }

                            wdu[0] = ( wu[2] - wu[5] ) / ( 2 * wdx );
                            wdv[0] = ( wv[1] - wv[4] ) / ( 2 * wdy );

                            if ( sico.MODEL == 7  ) wdu[1] = ( wu2[2] - wu2[5] ) / ( 2 * wdx ); else wdu[1] = 0;
                            if ( sico.MODEL == 7  ) wdv[1] = ( wv2[1] - wv2[4] ) / ( 2 * wdy ); else wdv[1] = 0;

                            if ( sico.MODEL == 7 ) wdu[2] = ( wu3[2] - wu3[5] ) / ( 2 * wdx ); else wdu[2] = 0;
                            if ( sico.MODEL == 7 ) wdv[2] = ( wv3[1] - wv3[4] ) / ( 2 * wdy ); else wdv[2] = 0;

                            xwdv[0] = ( wv[2] - wv[5] ) / ( 2 * wdx );
                            xwdu[0] = ( wu[1] - wu[4] ) / ( 2 * wdy );

                            if ( sico.MODEL == 7 ) xwdv[1] = ( wv2[2] - wv2[5] ) / ( 2 * wdx ); else xwdv[1] = 0;
                            if ( sico.MODEL == 7 ) xwdu[1] = ( wu2[1] - wu2[4] ) / ( 2 * wdy ); else xwdu[1] = 0;

                            if ( sico.MODEL == 7 ) xwdv[2] = ( wv3[2] - wv3[5] ) / ( 2 * wdx ); else xwdv[2] = 0;
                            if ( sico.MODEL == 7 ) xwdu[2] = ( wu3[1] - wu3[4] ) / ( 2 * wdy ); else xwdu[2] = 0;

                            walphax = ( walpha[2] - walpha[5] ) / ( 2 * wdx );
                            walphay = ( walpha[1] - walpha[4] ) / ( 2 * wdy );

                            if ( sico.MODEL == 7 ) walphax2 = ( walpha2[2] - walpha2[5] ) / ( 2 * wdx ); else walphax2 = 0;
                            if ( sico.MODEL == 7 ) walphay2 = ( walpha2[1] - walpha2[4] ) / ( 2 * wdy ); else walphay2 = 0;

                            if ( sico.MODEL == 7 ) walphax3 = ( walpha3[2] - walpha3[5] ) / ( 2 * wdx ); else walphax3 = 0;
                            if ( sico.MODEL == 7 ) walphay3 = ( walpha3[1] - walpha3[4] ) / ( 2 * wdy ); else walphay3 = 0;

 
                            // Earth pressure coefficients

                            if ( sico.MODEL <= 3 ) { whflow = wintb[i][j][0], whflow2 = 0, whflow3 = 0; }
                            else if ( sico.MODEL == 7 ) { whflow = wintb[i][j][0], whflow2 = wintb[i][j][3], whflow3 = wintb[i][j][6]; }

                            hekin = 0.5 * whflow * sflow.RHO1 * ( pow( vflowx, 2 ) + pow( vflowy, 2 ));
                            if ( sico.MODEL == 7 ) hekin += 0.5 * whflow2 * sflow.RHO2 * ( pow( vflowx2, 2 ) + pow( vflowy2, 2 ));
                            if ( sico.MODEL == 7 ) hekin += 0.5 * whflow3 * sflow.RHO3 * ( pow( vflowx3, 2 ) + pow( vflowy3, 2 ));

                            for ( ll=0; ll<sico.PMAX; ll++ ) {

                                if ( sflow.DELTA[ll] > 0 ) {

                                    gkx[ll] = fk( whflow+whflow2+whflow3, wintb[i][j][3*ll], vflowx, vflowy, wdu[ll], hekin, ll, sico, sflow );
                                    gky[ll] = fk( whflow+whflow2+whflow3, wintb[i][j][3*ll], vflowx, vflowy, wdv[ll], hekin, ll, sico, sflow );
                                }
                            }
                            
                            if ( sico.MODEL <= 3 ) wtemp = wintb[i][j][3]; else wtemp = wintb[i][j][9];


                            // Flux terms

                            if ( sico.ZFRAG == 1 ) ctrl_zfrag = pzfrag[i]; else ctrl_zfrag = 0; 

                            vm = fvm( whflow, whflow2, whflow3, vflowx, vflowx2, vflowx3, vflowy, vflowy2, vflowy3, sico, sflow );
                            cdrag = fdrag( whflow, whflow2, whflow3, sico, sflow );
                            gze = fgze( whflow, whflow2, whflow3, vflowx, vflowx2, vflowx3, vflowy, vflowy2, vflowy3, wh, whx, why, wdu, wdv, xwdu, xwdv, 
                                wu, wv, wu2, wv2, wu3, wv3, walpha, walpha2, walpha3, wdx, wdy, wgrav, wbetax, wbetay, cdrag, tsum, sico, sflow );
                            disp = fdisp( whflow, whflow2, whflow3, vflowx, vflowx2, vflowx3, vflowy, vflowy2, vflowy3, wdu, wdv, xwdu, xwdv, 
                                wbetax, wbetay, wbetaxy, gkx, gky, vm, cdrag, tsum, sico, sflow );

                            gf = ff( wh, whflow, whflow2, whflow3, vflowx, vflowx2, vflowx3, vflowy, vflowy2, vflowy3,
                                gkx, gze, vm, whx, walphax, cdrag, wgrav, disp, wbetax, wbetay, wbetaxy, ctrl_zfrag, wtemp, sflow, sico );
                            gg = fg( wh, whflow, whflow2, whflow3, vflowx, vflowx2, vflowx3, vflowy, vflowy2, vflowy3,
                                gky, gze, vm, why, walphay, cdrag, wgrav, disp, wbetax, wbetay, wbetaxy, ctrl_zfrag, wtemp, sflow, sico );

                            for ( k=0; k<sico.NVECTMIN; k++ ) { f[j][k] = gf[k]; g[j][k] = gg[k]; }
                            free( vm ); free( cdrag ); free( gze ); free( disp ); free( gf ); free( gg );


                            // Flow velocities in x and y directions

                            vflowx = fdiv( wintc[i][j][1], wintc[i][j][0], sico.HFLOWMIN );
                            vflowy = fdiv( wintc[i][j][2], wintc[i][j][0], sico.HFLOWMIN );

                            if ( sico.MODEL == 7 ) vflowx2 = fdiv( wintc[i][j][4], wintc[i][j][3], sico.HFLOWMIN ); else vflowx2 = 0;
                            if ( sico.MODEL == 7 ) vflowy2 = fdiv( wintc[i][j][5], wintc[i][j][3], sico.HFLOWMIN ); else vflowy2 = 0;

                            if ( sico.MODEL == 7 ) vflowx3 = fdiv( wintc[i][j][7], wintc[i][j][6], sico.HFLOWMIN ); else vflowx3 = 0;
                            if ( sico.MODEL == 7 ) vflowy3 = fdiv( wintc[i][j][8], wintc[i][j][6], sico.HFLOWMIN ); else vflowy3 = 0;


                            // Total flow depths, velocities, fractions, and slopes of neighbouring cells

                            for ( l=0; l<9; l++ ) {

                                welev[l] = wintelev[in[i][l]][j];
                                
                                if ( cdomain2[in[i][l]] == 1 ) {

                                    nbetax[l] = betax[in[i][l]];
                                    nbetay[l] = betay[in[i][l]];

                                } else {

                                    nbetax[l] = fbeta( wintelev[in[in[i][l]][5]][j], wintelev[in[in[i][l]][2]][j], 2.0, sico );
                                    nbetay[l] = fbeta( wintelev[in[in[i][l]][4]][j], wintelev[in[in[i][l]][1]][j], 2.0, sico );

                                }

                                if ( sico.MODEL <= 3 ) wh[l] = wintc[in[i][l]][j][0];
                                else if ( sico.MODEL == 7 ) {
                                
                                    wh[l] = wintc[in[i][l]][j][0] + wintc[in[i][l]][j][3] + wintc[in[i][l]][j][6];

                                    if ( sico.LAYERS > 0 ) {
                                
                                        wh1[l] = wintc[in[i][l]][j][0];
                                        wh2[l] = wintc[in[i][l]][j][0] + wintc[in[i][l]][j][3];
                                    }
                                }

                                wu[l] = fdiv( wintc[in[i][l]][j][1], wintc[in[i][l]][j][0], sico.HFLOWMIN );
                                wv[l] = fdiv( wintc[in[i][l]][j][2], wintc[in[i][l]][j][0], sico.HFLOWMIN );

                                if ( sico.MODEL == 7 ) wu2[l] = fdiv( wintc[in[i][l]][j][4], wintc[in[i][l]][j][3], sico.HFLOWMIN ); else wu2[l] = 0;
                                if ( sico.MODEL == 7 ) wv2[l] = fdiv( wintc[in[i][l]][j][5], wintc[in[i][l]][j][3], sico.HFLOWMIN ); else wv2[l] = 0;

                                if ( sico.MODEL == 7 ) wu3[l] = fdiv( wintc[in[i][l]][j][7], wintc[in[i][l]][j][6], sico.HFLOWMIN ); else wu3[l] = 0;
                                if ( sico.MODEL == 7 ) wv3[l] = fdiv( wintc[in[i][l]][j][8], wintc[in[i][l]][j][6], sico.HFLOWMIN ); else wv3[l] = 0;

                                walpha[l] = fdiv( wintc[in[i][l]][j][0], wh[l], sico.HFLOWMIN );
                                if ( sico.MODEL == 7 ) walpha2[l] = fdiv( wintc[in[i][l]][j][3], wh[l], sico.HFLOWMIN ); else walpha2[l] = 0;
                                if ( sico.MODEL == 7 ) walpha3[l] = fdiv( wintc[in[i][l]][j][6], wh[l], sico.HFLOWMIN ); else walpha3[l] = 0;
                            }


                            // Gradients of total flow depth, flow velocities, and fractions

                            cin = fw0noosc( wh, sico );

                            wbetax = fbeta( welev[cin[0]], welev[cin[1]], (float)(cin[4]), sico );
                            wbetay = fbeta( welev[cin[2]], welev[cin[3]], (float)(cin[5]), sico );
                            
                            free( cin );

                            wdx = dx[i];
                            wdy = dy[i];

                            wgrav[0] = sico.GRAVITY * cos( wbetax );
                            wgrav[1] = sico.GRAVITY * cos( wbetay );
                            wgrav[2] = sico.GRAVITY * sin( wbetax );
                            wgrav[3] = sico.GRAVITY * sin( wbetay );
                            wgrav[4] = sico.GRAVITY * cos( wbetaxy );

                            cin = f0noosc( wh, sico );

                            whx = ( wh[cin[0]] - wh[cin[1]] ) / ((float)(cin[4]) * dx[i] );
                            why = ( wh[cin[2]] - wh[cin[3]] ) / ((float)(cin[5]) * dy[i] );

                            if ( sico.LAYERS > 0 ) {
                            
                                whx1 = ( wh1[cin[0]] - wh1[cin[1]] ) / ((float)(cin[4]) * dx[i] );
                                whx2 = ( wh2[cin[0]] - wh2[cin[1]] ) / ((float)(cin[4]) * dx[i] );
                                why1 = ( wh1[cin[2]] - wh1[cin[3]] ) / ((float)(cin[5]) * dy[i] );
                                why2 = ( wh2[cin[2]] - wh2[cin[3]] ) / ((float)(cin[5]) * dy[i] );
                            }
                        
                            free( cin );

                            wdu[0] = ( wu[2] - wu[5] ) / ( 2 * wdx );
                            wdv[0] = ( wv[1] - wv[4] ) / ( 2 * wdy );

                            if ( sico.MODEL == 7 ) wdu[1] = ( wu2[2] - wu2[5] ) / ( 2 * wdx ); else wdu[1] = 0;
                            if ( sico.MODEL == 7 ) wdv[1] = ( wv2[1] - wv2[4] ) / ( 2 * wdy ); else wdv[1] = 0;

                            if ( sico.MODEL == 7 ) wdu[2] = ( wu3[2] - wu3[5] ) / ( 2 * wdx ); else wdu[2] = 0;
                            if ( sico.MODEL == 7 ) wdv[2] = ( wv3[1] - wv3[4] ) / ( 2 * wdy ); else wdv[2] = 0;

                            xwdv[0] = ( wv[2] - wv[5] ) / ( 2 * wdx );
                            xwdu[0] = ( wu[1] - wu[4] ) / ( 2 * wdy );

                            if ( sico.MODEL == 7 ) xwdv[1] = ( wv2[2] - wv2[5] ) / ( 2 * wdx ); else xwdv[1] = 0;
                            if ( sico.MODEL == 7 ) xwdu[1] = ( wu2[1] - wu2[4] ) / ( 2 * wdy ); else xwdu[1] = 0;

                            if ( sico.MODEL == 7 ) xwdv[2] = ( wv3[2] - wv3[5] ) / ( 2 * wdx ); else xwdv[2] = 0;
                            if ( sico.MODEL == 7 ) xwdu[2] = ( wu3[1] - wu3[4] ) / ( 2 * wdy ); else xwdu[2] = 0;

                            walphax = ( walpha[2] - walpha[5] ) / ( 2 * wdx );
                            walphay = ( walpha[1] - walpha[4] ) / ( 2 * wdy );

                            if ( sico.MODEL == 7 ) walphax2 = ( walpha2[2] - walpha2[5] ) / ( 2 * wdx ); else walphax2 = 0;
                            if ( sico.MODEL == 7 ) walphay2 = ( walpha2[1] - walpha2[4] ) / ( 2 * wdy ); else walphay2 = 0;

                            if ( sico.MODEL == 7 ) walphax3 = ( walpha3[2] - walpha3[5] ) / ( 2 * wdx ); else walphax3 = 0;
                            if ( sico.MODEL == 7 ) walphay3 = ( walpha3[1] - walpha3[4] ) / ( 2 * wdy ); else walphay3 = 0;


                            // Source terms

                            if ( sico.MODEL <= 3 ) { whflow = wintc[i][j][0], whflow2 = 0, whflow3 = 0; }
                            else if ( sico.MODEL == 7 ) { whflow = wintc[i][j][0], whflow2 = wintc[i][j][3], whflow3 = wintc[i][j][6]; }

                            hekin = 0.5 * whflow * sflow.RHO1 * ( pow( vflowx, 2 ) + pow( vflowy, 2 ));
                            if ( sico.MODEL == 7 ) hekin += 0.5 * whflow2 * sflow.RHO2 * ( pow( vflowx2, 2 ) + pow( vflowy2, 2 ));
                            if ( sico.MODEL == 7 ) hekin += 0.5 * whflow3 * sflow.RHO3 * ( pow( vflowx3, 2 ) + pow( vflowy3, 2 ));

                            for ( ll=0; ll<sico.PMAX; ll++ ) {

                                if ( sflow.DELTA[ll] > 0 ) {

                                    gkx[ll] = fk( whflow+whflow2+whflow3, wintc[i][j][3*ll], vflowx, vflowy, wdu[ll], hekin, ll, sico, sflow );
                                    gky[ll] = fk( whflow+whflow2+whflow3, wintc[i][j][3*ll], vflowx, vflowy, wdv[ll], hekin, ll, sico, sflow );
                                }
                            }

                            for ( l=0; l<sico.PMAX; l++ ) {

                                kappau[l] = ( betax[in[i][2]] - betax[in[i][4]] ) * 0.5 * cos ( wbetax ) / sico.CSZ;  
                                kappau[l+3] = ( betay[in[i][1]] - betay[in[i][5]] ) * 0.5 * cos ( wbetay ) / sico.CSZ; 
                            }

                            if ( sico.MODEL <= 3 ) wtemp = wintc[i][j][3]; else wtemp = wintc[i][j][9];

                            vm = fvm( whflow, whflow2, whflow3, vflowx, vflowx2, vflowx3, vflowy, vflowy2, vflowy3, sico, sflow );
                            cdrag = fdrag( whflow, whflow2, whflow3, sico, sflow );

                            gze = fgze( whflow, whflow2, whflow3, vflowx, vflowx2, vflowx3, vflowy, vflowy2, vflowy3, wh, whx, why, wdu, wdv, xwdu, xwdv, 
                                wu, wv, wu2, wv2, wu3, wv3, walpha, walpha2, walpha3, wdx, wdy, wgrav, wbetax, wbetay, cdrag, tsum, sico, sflow );

                            gs = fs( wh, whflow, whflow2, whflow3, vflowx, vflowx2, vflowx3, vflowy, vflowy2, vflowy3,
                                whx, why, whx1, why1, whx2, why2, walphax, walphay, walphax2, walphay2, walphax3, walphay3, wgrav, gze, wbetax, wbetay, wbetaxy, 
                                cdrag, gkx, gky, pxslide[i], angx, angy, tsum, ptslide[i], ctrl_zfrag, wtemp, sico, sflow );

                            if ( vflowx == 0 ) wu[0] = fdiv( 0.5 * gs[1], whflow, sico.HFLOWMIN );
                            if ( vflowy == 0 ) wv[0] = fdiv( 0.5 * gs[2], whflow, sico.HFLOWMIN );
                        
                            if ( sico.MODEL == 7 ) {
                            
                                if ( vflowx2 == 0 ) wu2[0] = fdiv( 0.5 * gs[4], whflow2, sico.HFLOWMIN );
                                if ( vflowy2 == 0 ) wv2[0] = fdiv( 0.5 * gs[5], whflow2, sico.HFLOWMIN );
                                if ( vflowx3 == 0 ) wu3[0] = fdiv( 0.5 * gs[7], whflow3, sico.HFLOWMIN );
                                if ( vflowy3 == 0 ) wv3[0] = fdiv( 0.5 * gs[8], whflow3, sico.HFLOWMIN );
                            }

                            for ( k=0; k<sico.NVECTMIN; k++ ) s[j][k] = gs[k];
                            free( gs );


                            // Deceleration terms

                            sico.XDIST = fabs(sico.XREL - px[i]) * sico.CSZ;
                            sico.YDIST = fabs(sico.YREL - py[i]) * sico.CSZ;

                            gdecel = fd( wh, whx, why, wu, wv, wu2, wv2, wu3, wv3, walpha, walpha2, walpha3, wbetax, wbetay, wbetaxy, nbetax, nbetay, gze, 
                                kappau, hekin, wgrav, angx, angy, pxslide[i], tsum, ptslide[i], ctrl_zfrag, sico, sflow );

                            for ( k=0; k<sico.NVECTMIN; k++ ) d[i][j][k] = gdecel[k];
                            free( cdrag ); free( vm ); free( gdecel );

                            vcelr0 = fcfl( whflow, whflow2, whflow3, vflowx, vflowy, vflowx2, vflowy2, vflowx3, vflowy3, sflow, gkx, gky, betax[i], betay[i], betaxy[i], 
                                ctrl_zfrag, gze, whx, why, whx1, why1, whx2, why2, i, sico ); // flow velocity plus wave speed
                        
                            free( gze ); 
                        
                            if ( vcelr0 > vcelr ) vcelr = vcelr0;
                        }


                        // Value of vector at top right corner of the cell

                        wintelevd[i] = 0.25 * ( wintelev[i][0] + wintelev[i][1] + wintelev[i][2] + wintelev[i][3] );

                        for ( k=0; k<sico.NVECTMIN; k++ ) {

                            if ( sico.SLOMO > 1 && ( k == 1 || k == 2 || k == 4 || k == 5 || k == 7 || k == 8 ))

                                wintd[i][k]=0.25*(s[0][k]+s[1][k]+s[2][k]+s[3][k]);
                                
                            else wintd[i][k]=0.25*(winta[i][0][k]+winta[i][1][k]+winta[i][2][k]+winta[i][3][k])
                                -tlength/dx[i]*(0.5*(f[2][k]+f[3][k])-0.5*(f[0][k]+f[1][k]))
                                -tlength/dy[i]*(0.5*(g[1][k]+g[3][k])-0.5*(g[0][k]+g[2][k]))
                                +0.25*tlength*(s[0][k]+s[1][k]+s[2][k]+s[3][k]);

                            if ( sico.SLOMO > 1 && ( k == 1 || k == 2 || k == 4 || k == 5 || k == 7 || k == 8 ))

                                wintdtest[i][k]=0.25*(s[0][k]+s[1][k]+s[2][k]+s[3][k]);
                                
                            else wintdtest[i][k]=0.25*(winta[i][0][k]+winta[i][1][k]+winta[i][2][k]+winta[i][3][k])
                                -tlength/dx[i]*(0.5*(f[2][k]+f[3][k])-0.5*(f[0][k]+f[1][k]))
                                -tlength/dy[i]*(0.5*(g[1][k]+g[3][k])-0.5*(g[0][k]+g[2][k]))
                                +0.25*tlength*(s[0][k]+s[1][k]+s[2][k]+s[3][k]-d[i][0][k]-d[i][1][k]-d[i][2][k]-d[i][3][k]);
                        }
                        
                    } else {

                        wintelevd[i] = pelev[i];

                        for ( k=0; k<sico.NVECTMIN; k++ ) { 
                    
                            wintd[i][k] = aw[i][k];
                            wintdtest[i][k] = aw[i][k];
                        }
                    }
                }


                // Moving vector if second sub-timestep and writing values to temporary vector

                for ( ix=0; ix<ib[0]; ix++ ) {

                    i = ibasket[0][ix];

                    cin = fnoosc( aw, in, i, sico );

                    if ( pxslide[i] < 2 ) {

                        ctrlv = 1; if ( in[i][6] == 0 ) ctrlv = 0;
                        if ( ctrlv == 1 ) {
                            ctrlr = 0;
                            if ( sico.MODEL <= 3 ) {

                                if ( wintd[i][0] > 0  ) ctrlr = 1;
                                if ( wintd[in[i][6]][0] > 0 ) ctrlr = 1;

                            } else if ( sico.MODEL == 7 ) {

                                if ( wintd[i][0] > 0 || wintd[i][3] > 0 || wintd[i][6] > 0 ) ctrlr = 1;
                                if ( wintd[in[i][6]][0] > 0 || wintd[in[i][6]][3] > 0 || wintd[in[i][6]][6] > 0 ) ctrlr = 1;
                            }
                        }

                        if ( ctrlv == 1 && ctrlr == 1 && cdomain[i] > 3 && cdomain2[i] == 0 ) {

                            if ( sico.MODEL <= 3 ) hflow = wintd[i][0];
                            else if ( sico.MODEL == 7 ) hflow = wintd[i][0] + wintd[i][3] + wintd[i][6];

                            if ( sico.MODEL <= 3 ) hflown = wintd[in[i][6]][0];
                            else if ( sico.MODEL == 7 ) hflown = wintd[in[i][6]][0] + wintd[in[i][6]][3] + wintd[in[i][6]][6];

                            if ( p == 1 && cdomain[in[i][6]] > 3 && cplain[i] == 1 && sico.SURFACE > 1 ) qelev[i] = wintelevd[in[i][6]];
                            else if ( p != 1 && cdomain[i] > 3 && cplain[i] == 1 && sico.SURFACE > 1 ) qelev[i] = wintelevd[i];
                            else qelev[i] = pelev[i];

                            for ( k=0; k<sico.NVECTMIN; k++ ) {
                            
                                if ( p == 1 && cdomain[in[i][6]] > 3 ) {
                            
                                    awt[i][k] = wintd[in[i][6]][k];
                                    awttest[k] = wintdtest[in[i][6]][k];

                                } else if ( p != 1 && cdomain[i] > 3 ) {

                                    awt[i][k] = wintd[i][k];
                                    awttest[k] = wintdtest[i][k];
                                
                                } else if ( k<sico.NVECTMIN ) { awt[i][k] = aw[i][k]; awttest[k] = aw[i][k];
                                } else { awt[i][k] = aw[i][k]; awttest[k] = aw[i][k]; }

                                if (( fsign(awt[i][k]) == fsign(awttest[k]) && fabs(awt[i][k]) > fabs(awttest[k])) || sico.ORIGINAL == 1 ) awt[i][k] = awttest[k];
                                else if ( fsign(awt[i][k]) == fsign(awttest[k])) awt[i][k] = awt[i][k];
                                else if ( k==0 || k==3 || k==6 ) awt[i][k] = aw[i][k];
                                else awt[i][k] = 0;
                            }

                        } else {
                        
                            qelev[i] = pelev[i];
                            if ( cstopped[i] != 1 ) { for ( k=0; k<sico.NVECTMIN; k++ ) awt[i][k] = 0;
                            } else { for ( k=0; k<sico.NVECTMIN; k++ ) awt[i][k] = aw[i][k]; }
                        }

                        if ( awt[i][0] < 0 ) { awt[i][0] = 0; awt[i][1] = 0; awt[i][2] = 0; }
                        if ( sico.MODEL == 7 && awt[i][3] < 0 ) { awt[i][3] = 0; awt[i][4] = 0; awt[i][5] = 0; }
                        if ( sico.MODEL == 7 && awt[i][6] < 0 ) { awt[i][6] = 0; awt[i][7] = 0; awt[i][8] = 0; }
                        for ( k=0; k<sico.NVECTMIN; k++ ) { if ( isnan( awt[i][k] ) != 0 ) awt[i][k] = 0; }

                        if ( cdomain2[i] == 1 ) { for ( k=0; k<sico.NVECTMIN; k++ ) awt[i][k] = aw[i][k]; }
                    }
                    
                    free ( cin );
                  }
                }
                

// -- STOP --- Flow propagation with NOC scheme -----------------------------------------------------------------


// -- START -- Evaluating time step length and validity of time step --------------------------------------------


                cfl = vcelr * tlength / sico.CSZ; // updating CFL value

                if ( cfl <= sico.CFL[0] || vcelr == 0 ) { // if time step length meets CFL criterion

                    ccfl = 1; // setting control for fulfilment of CFL criterion to positive
                    if ( cfl > cflmax ) cflmax = cfl; // updating maximum CFL value, if necessary
                    tlength = tlength / sflow.LAMBDA_FL; // adjusting length to flux control

                    tlengthpre = tlength; // updating length of previous time step

                    tsum += tlength; // updating time since release
                    tint += tlength; // updating time since last output

                    if ( tlength <= pow( 10, -7 )) { ccontinue = 0; csuccess = 0; } // defining numerical failure

                    if ( vcelr > 0 ) tlength0 = sico.CFL[0] * 0.90 * sico.CSZ / vcelr; // updating length of time step according to CFL criterion
                    else  tlength0 = sico.CFL[1]; 

                    vcelr = 0; // resetting flow velocity plus wave speed
                    
                } else { // if time step is too long too meet CFL criterion, repeating time step with reduced length

                    tlength = 0.75 * tlength * sico.CFL[0] / cfl; // defining new time step length
                    vcelr = 0; // resetting flow velocity plus wave speed
                }


// -- STOP --- Evaluating time step length and validity of time step --------------------------------------------


                if ( sico.MULT == 0 ) {

                    printf("   %i\t%i\t%.3f\t%.1f  \t%.1f\t...\r", nout, nsum, cflmax, 100*tlength, tsum ); // progress
                    fflush(stdout); // forcing immediate display
                }


// *** End of loop over time step lengths -----------------------------------------------------------------------


            }

            if ( cslide < 2 ) iloop = ib[0]; else iloop = sico.IMAX;
            for ( ix=0; ix<iloop; ix++ ) {

                if ( cslide < 2 ) i = ibasket[0][ix]; else i = ix;
                
                relev[i] = qelev[i] - pelev[i];
                pelev[i] = qelev[i];
                for ( k=0; k<sico.NVECTMIN; k++ ) aw[i][k] = awt[i][k]; // applying values of temporary vectors


// -- START -- Phase transformations ----------------------------------------------------------------------------
 

                if ( sico.MODEL == 7 ) {

                    if ( sico.TRANS12 == 1 ) { // PHASE 1 - PHASE 2 transformation coefficient
                        if ( ptrans12[i] != sico.UNDEF ) sflow.TRANS12 = ptrans12[i];
                        else sflow.TRANS12 = sflow.TRANS120;
                    } else sflow.TRANS12 = sflow.TRANS120;

                    if ( sico.TRANS13 == 1 ) { // PHASE 1 - PHASE 3 transformation coefficient
                        if ( ptrans13[i] != sico.UNDEF ) sflow.TRANS13 = ptrans13[i];
                        else sflow.TRANS13 = sflow.TRANS130;
                    } else sflow.TRANS13 = sflow.TRANS130;

                    if ( sico.TRANS23 == 1 ) { // PHASE 2 - PHASE 3 transformation coefficient
                        if ( ptrans23[i] != sico.UNDEF ) sflow.TRANS23 = ptrans23[i];
                        else sflow.TRANS23 = sflow.TRANS230;
                    } else sflow.TRANS23 = sflow.TRANS230;

                    if ( transformograph == 1 ) { // reading transformograph, if available

                        for ( trak=0; trak<=tratmax; trak++ ) {
                            if ( tratra[trak][0] <= tsum ) tratx=trak;
                            else break;
                        }

                        sflow.TRANS12 = tratra[tratx][1];
                        sflow.TRANS13 = tratra[tratx][2];
                        sflow.TRANS23 = tratra[tratx][3];
                    }

                    if ( aw[i][0] > 0 ) { vflowx = aw[i][1] / aw[i][0]; vflowy = aw[i][2] / aw[i][0]; } else { vflowx = 0; vflowy = 0; } // flow velocities
                    if ( aw[i][3] > 0 ) { vflowx2 = aw[i][4] / aw[i][3]; vflowy2 = aw[i][5] / aw[i][3]; } else { vflowx2 = 0; vflowy2 = 0; }
                    if ( aw[i][6] > 0 ) { vflowx3 = aw[i][7] / aw[i][6]; vflowy3 = aw[i][8] / aw[i][6]; } else { vflowx3 = 0; vflowy3 = 0; }

                    if ( transformograph == 1 ) ekin = 1; // transformograph does not use flow kinetic energy
                    else ekin = ( aw[i][0] * sflow.RHO1 * ( pow( vflowx, 2 ) + pow( vflowy, 2 ))
                        + aw[i][3] * sflow.RHO2 * ( pow( vflowx2, 2 ) + pow( vflowy2, 2 ))
                        + aw[i][6] * sflow.RHO3 * ( pow( vflowx3, 2 ) + pow( vflowy3, 2 ))) * 0.5; // flow kinetic energy

                    if ( sflow.TRANS13 > 0 ) { // PHASE 1 to PHASE 3 transformation

                        trans = ffmin( aw[i][0], tlength * pow( 10, -fabs( sflow.TRANS13 )) * ekin );

                        aw[i][0] -= trans;
                        aw[i][1] -= trans * vflowx;
                        aw[i][2] -= trans * vflowy;

                        aw[i][6] += trans * sflow.RHO1 / sflow.RHO3;
                        aw[i][7] += trans * sflow.RHO1 / sflow.RHO3 * vflowx;
                        aw[i][8] += trans * sflow.RHO1 / sflow.RHO3 * vflowy;
                    }

                    if ( sflow.TRANS13 < 0 ) { // PHASE 3 to PHASE 1 transformation

                        trans = ffmax( -aw[i][6], -tlength * pow( 10, -fabs( sflow.TRANS13 )) * ekin );

                        aw[i][0] -= trans * sflow.RHO3 / sflow.RHO1;
                        aw[i][1] -= trans * vflowx3 * sflow.RHO3 / sflow.RHO1;
                        aw[i][2] -= trans * vflowy3 * sflow.RHO3 / sflow.RHO1;

                        aw[i][6] += trans;
                        aw[i][7] += trans * vflowx3;
                        aw[i][8] += trans * vflowy3;
                    }

                    if ( sflow.TRANS12 > 0 ) { // PHASE 1 to PHASE 2 transformation

                        trans = ffmin( aw[i][0], tlength * pow( 10, -fabs( sflow.TRANS12 )) * ekin );

                        aw[i][0] -= trans;
                        aw[i][1] -= trans * vflowx;
                        aw[i][2] -= trans * vflowy;

                        aw[i][3] += trans * sflow.RHO1 / sflow.RHO2;
                        aw[i][4] += trans * sflow.RHO1 / sflow.RHO2 * vflowx;
                        aw[i][5] += trans * sflow.RHO1 / sflow.RHO2 * vflowy;
                    }

                    if ( sflow.TRANS12 < 0 ) { // PHASE 2 to PHASE 1 transformation

                        trans = ffmax( -aw[i][3], -tlength * pow( 10, -fabs( sflow.TRANS12 )) * ekin );

                        aw[i][0] -= trans * sflow.RHO2 / sflow.RHO1;
                        aw[i][1] -= trans * vflowx2 * sflow.RHO2 / sflow.RHO1;
                        aw[i][2] -= trans * vflowy2 * sflow.RHO2 / sflow.RHO1;

                        aw[i][3] += trans;
                        aw[i][4] += trans * vflowx2;
                        aw[i][5] += trans * vflowy2;
                    }         

                    if ( sflow.TRANS23 > 0 ) { // PHASE 2 to PHASE 3 transformation

                        trans = ffmin( aw[i][3], tlength * pow( 10, -fabs( sflow.TRANS23 )) * ekin );

                        aw[i][3] -= trans;
                        aw[i][4] -= trans * vflowx2;
                        aw[i][5] -= trans * vflowy2;

                        aw[i][6] += trans * sflow.RHO2 / sflow.RHO3;
                        aw[i][7] += trans * sflow.RHO2 / sflow.RHO3 * vflowx2;
                        aw[i][8] += trans * sflow.RHO2 / sflow.RHO3 * vflowy2;
                    }

                    if ( sflow.TRANS23 < 0 ) { // PHASE 3 to PHASE 2 transformation

                        trans = ffmax( -aw[i][6], -tlength * pow( 10, -fabs( sflow.TRANS23 )) * ekin );

                        aw[i][3] -= trans * sflow.RHO3 / sflow.RHO2;
                        aw[i][4] -= trans * vflowx3 * sflow.RHO3 / sflow.RHO2;
                        aw[i][5] -= trans * vflowy3 * sflow.RHO3 / sflow.RHO2;

                        aw[i][6] += trans;
                        aw[i][7] += trans * vflowx3;
                        aw[i][8] += trans * vflowy3;
                    }
                }
            }


// -- STOP --- Phase transformations ----------------------------------------------------------------------------


            if ( cslide < 2 ) iloop = ib[0]; else iloop = sico.IMAX;
            for ( ix=0; ix<iloop; ix++ ) {

                if ( cslide < 2 ) i = ibasket[0][ix]; else i = ix;
                
                if ( cdomain2[i] == 0 ) {
                
                    cin = fnoosc( aw, in, i, sico );

                    betax[i] = fbeta( pelev[cin[0]], pelev[cin[1]], (float)(cin[4]), sico ); // slopes
                    betay[i] = fbeta( pelev[cin[2]], pelev[cin[3]], (float)(cin[5]), sico );
                    betaxh[i] = fbeta( pelev[cin[0]]+aw[cin[0]][0], pelev[cin[1]]+aw[cin[1]][0], (float)(cin[4]), sico ); // slopes including flow heights
                    betayh[i] = fbeta( pelev[cin[2]]+aw[cin[2]][0], pelev[cin[3]]+aw[cin[3]][0], (float)(cin[5]), sico );
                    betaxy[i] = fbetaxy( betax[i], betay[i] );

                    if ( sico.LAYERS0 > 0 ) {
                            
                        betax2[i] = fbeta( pelev[cin[0]]+aw[cin[0]][0], pelev[cin[1]]+aw[cin[1]][0], (float)(cin[4]), sico );
                        betay2[i] = fbeta( pelev[cin[2]]+aw[cin[2]][0], pelev[cin[3]]+aw[cin[3]][0], (float)(cin[5]), sico );

                        betaxy2[i] = fbetaxy( betax2[i], betay2[i] );
                                                               
                        betaxh2[i] = fbeta( pelev[cin[0]]+aw[cin[0]][0]+aw[cin[0]][3], pelev[cin[1]]+aw[cin[1]][0]+aw[cin[1]][3], (float)(cin[4]), sico );
                        betayh2[i] = fbeta( pelev[cin[2]]+aw[cin[2]][0]+aw[cin[2]][3], pelev[cin[3]]+aw[cin[3]][0]+aw[cin[3]][3], (float)(cin[5]), sico );
                                    
                        betax3[i] = fbeta( pelev[cin[0]]+aw[cin[0]][0]+aw[cin[0]][3], pelev[cin[1]]+aw[cin[1]][0]+aw[cin[1]][3], (float)(cin[4]), sico );
                        betay3[i] = fbeta( pelev[cin[2]]+aw[cin[2]][0]+aw[cin[2]][3], pelev[cin[3]]+aw[cin[3]][0]+aw[cin[3]][3], (float)(cin[5]), sico );
                        
                        betaxy3[i] = fbetaxy( betax3[i], betay3[i] );
                    }
                
                    if ( sico.LAYERS0 > 0 || sico.SURFACE > 1 ) {                
                                   
                        betaxh3[i] = fbeta( pelev[cin[0]]+aw[cin[0]][0]+aw[cin[0]][3]+aw[cin[0]][6], pelev[cin[1]]+aw[cin[1]][0]+aw[cin[1]][3]+aw[cin[1]][6], 
                            (float)(cin[4]), sico );
                        betayh3[i] = fbeta( pelev[cin[2]]+aw[cin[2]][0]+aw[cin[2]][3]+aw[cin[2]][6], pelev[cin[3]]+aw[cin[3]][0]+aw[cin[3]][3]+aw[cin[3]][6], 
                            (float)(cin[5]), sico );
                    }
     
                    free( cin );

                    if ( sico.MODEL <= 3 ) hflow = aw[i][0];
                    else if ( sico.MODEL == 7 ) hflow = aw[i][0] + aw[i][3] + aw[i][6];

                    if ( sico.SURFACE > 1 ) {
                       
                        for ( l=1; l<9; l++ ) {

                            if ( sico.MODEL <= 3 ) hflow = aw[in[i][0]][0];
                            else if ( sico.MODEL == 7 ) hflow = aw[in[i][0]][0] + aw[in[i][0]][3] + aw[in[i][0]][6];
            
                            if ( hflow <= sico.HFLOWMIN ) cplain[i] = 0;
                        }
                    }
                }
            }


// -- START -- Entrainment --------------------------------------------------------------------------------------


            // *** 1 = Entrainment coefficient multiplied with flow momentum


            if ( sico.ENTRAINMENT != 0 ) {

                if ( cslide < 2 ) iloop = ib[0]; else iloop = sico.IMAX;
                for ( ix=0; ix<iloop; ix++ ) {

                    if ( cslide < 2 ) i = ibasket[0][ix]; else i = ix;
                    qentr = 0;

                    if ( sico.CENTR == 1 ) {
                        if ( pcentr[i] != sico.UNDEF ) sflow.CENTR = pcentr[i]; // entrainment coefficient
                        else sflow.CENTR = sflow.CENTR0;
                    } else sflow.CENTR = sflow.CENTR0;

                    if ( sico.MODEL <= 3 ) { hflow = aw[i][0]; alpha1 = 1; } // flow heights and phase fractions
                    else { hflow = aw[i][0] + aw[i][3] + aw[i][6]; alpha1 = aw[i][0] / hflow; alpha2 = aw[i][3] / hflow; alpha3 = aw[i][6] / hflow; }
                    
                    if ( sico.MODEL <= 3 ) { hentrmaxx = phentrmax[i]; alphab1 = 1; } // entrainable heights and phase fractions
                    else {
                     
                        hentrmaxx = phentrmax[i] + phentrmax2[i] + phentrmax3[i]; 
                        if ( hentrmaxx > 0 ) { alphab1 = phentrmax[i] / hentrmaxx; alphab2 = phentrmax2[i] / hentrmaxx; alphab3 = phentrmax3[i] / hentrmaxx; }
                        else { alphab1 = 1; alphab2 = 0; alphab3 = 0; }
                    }

                    if ( aw[i][0] > 0 ) { vflowx1 = aw[i][1]/aw[i][0]; vflowy1 = aw[i][2]/aw[i][0]; } else { vflowx1 = 0; vflowy1 = 0; } // flow velocities of individual phases
                    if ( sico.MODEL == 7 && aw[i][3] > 0 ) { vflowx2 = aw[i][4]/aw[i][3]; vflowy2 = aw[i][5]/aw[i][3]; } else { vflowx2 = 0; vflowy2 = 0; }
                    if ( sico.MODEL == 7 && aw[i][6] > 0 ) { vflowx3 = aw[i][7]/aw[i][6]; vflowy3 = aw[i][8]/aw[i][6]; } else { vflowx3 = 0; vflowy3 = 0; }

                    for ( j=0; j<9; j++ ) {

                        if ( sico.MODEL <= 3 && hflow > sico.HFLOWMIN ) hflowj[j] = aw[in[i][j]][0];
                        else if ( hflow > sico.HFLOWMIN ) hflowj[j] = aw[in[i][j]][0] + aw[in[i][j]][3] + aw[in[i][j]][6];
                        elevtot[j] = pelev[in[i][j]]; // elevation of neighbour cell
                        
                        if ( sico.MODEL <= 3 && hflowj[j] > sico.HFLOWMIN && hflow > sico.HFLOWMIN ) {
                            
                            vflowxj[j] = aw[in[i][j]][1] / hflowj[j]; // flow velocities of neighbour cell
                            vflowyj[j] = aw[in[i][j]][2] / hflowj[j];
                            
                        } else if ( sico.MODEL == 7 && hflowj[j] > sico.HFLOWMIN && hflow > sico.HFLOWMIN ) {
                            
                            vflowxj[j] = ( aw[in[i][j]][1] + aw[in[i][j]][4] + aw[in[i][j]][5] ) / hflowj[j];
                            vflowyj[j] = ( aw[in[i][j]][2] + aw[in[i][j]][5] + aw[in[i][j]][8] ) / hflowj[j]; 

                        } else { vflowxj[j] = 0; vflowyj[j] = 0; }
                    }

                    wbetax = fbeta( elevtot[5], elevtot[2], 2.0, sico );
                    wbetay = fbeta( elevtot[4], elevtot[1], 2.0, sico );
                    wbetaxy = fbetaxy( wbetax, wbetay ); // slope, including flow height

                    alphav = falphav( vflowxj[0], vflowyj[0], sico ); // direction of movement
                    alpha = falpha( wbetax, wbetay, sico ); // aspect

                    if ( sico.ENTRAINMENT == 1 && hflow > sico.HFLOWMIN && hentrmaxx > 0 ) { // first model

                        mom = aw[i][0] * sflow.RHO1 * pow( pow( vflowx1, 2 ) + pow( vflowy1, 2 ), 0.5 ); // flow momenta
                        if ( sico.MODEL == 7 ) mom += aw[i][3] * sflow.RHO2 * pow( pow( vflowx2, 2 ) + pow( vflowy2, 2 ), 0.5 );
                        if ( sico.MODEL == 7 ) mom += aw[i][6] * sflow.RHO3 * pow( pow( vflowx3, 2 ) + pow( vflowy3, 2 ), 0.5 );

                        if ( sflow.CENTR > 0 ) betav = ffmax(0, tan( wbetaxy * cos( alpha - alphav ))); // movement-following slope
                        else betav = 1;
              
                        qentr = tlength * pow( 10, -fabs( sflow.CENTR )) * mom * betav; // entrainment rate
                        momfact = 0; // mobility generator
                  
                    } else if ( sico.ENTRAINMENT == 9 && hflow > 0 ) { // melt model !!!CHECK momentum for extremely rapid flows !!!CHECK combination with entrainment
                                                
                        if ( adaptograph == 1 ) {
                        
                            if ( adaada[adat][7] * phrelease[i] * adaada[adat][1] + adaada[adat][2] < 0 ) 
                                qmelt1 = ffmax( -aw[i][0], tlength * ( adaada[adat][7] * phrelease[i] * adaada[adat][1] + adaada[adat][2] )); else qmelt1 = 0;
                            if ( sico.MODEL == 7 && adaada[adat][7] * phrelease2[i] * adaada[adat][3] + adaada[adat][4] < 0 ) 
                                qmelt2 = ffmax( -aw[i][3], tlength * ( adaada[adat][7] * phrelease2[i] * adaada[adat][3] + adaada[adat][4] )); else qmelt2 = 0;
                            if ( sico.MODEL == 7 && adaada[adat][7] * phrelease3[i] * adaada[adat][5] + adaada[adat][6] < 0 ) 
                                qmelt3 = ffmax( -aw[i][6], tlength * ( adaada[adat][7] * phrelease3[i] * adaada[adat][5] + adaada[adat][6] )); else qmelt3 = 0;
                                
                        } else {

                            if ( phrelease[i] < 0 ) 
                                qmelt1 = ffmax( -aw[i][0], tlength * phrelease[i] ); else qmelt1 = 0;
                            if ( sico.MODEL == 7 && phrelease2[i] < 0 ) 
                                qmelt2 = ffmax( -aw[i][3], tlength * phrelease2[i] ); else qmelt2 = 0;
                            if ( sico.MODEL == 7 && phrelease3[i] < 0 ) 
                                qmelt3 = ffmax( -aw[i][6], tlength * phrelease3[i] ); else qmelt3 = 0;
                        }
                        
                        qmelt = qmelt1 + qmelt2 + qmelt3;
                    }

                    if ( qentr != 0 ) {
 
                        if ( qentr > 0 ) qentr1 = ffmin( phentrmax[i], alphab1 * qentr ); // constraining entrainment rate
                        else qentr1 = ffmax( -aw[i][0], alpha1 * qentr );
                        phentrmax[i] -= ffmax( 0, qentr1 ); // updating entrainable height
                        
                        momaddsx = momfact * vflowx1 * qentr1; // contribution of entrainment to momentum
                        momaddsy = momfact * vflowy1 * qentr1;

                        aw[i][0] = aw[i][0] + qentr1; // updating flow depths and momenta (mixture or PHASE 1)
                        aw[i][1] = aw[i][1] + momaddsx;
                        aw[i][2] = aw[i][2] + momaddsy;
                        
                        if ( sico.MODEL == 7 ) {
                        
                            if ( qentr > 0 ) qentr2 = ffmin( phentrmax2[i], alphab2 * qentr ); // constraining entrainment rates
                            else qentr2 = ffmax( -aw[i][3], alpha2 * qentr );
                            if ( qentr > 0 ) qentr3 = ffmin( phentrmax3[i], alphab3 * qentr );
                            else {
                            
                                qentr3 = ffmax( -aw[i][6], alpha3 * qentr );
                                qentr3 = ffmax( qentr3, sflow.THETAS * ( qentr1 + qentr2 ));
                            }
                            phentrmax2[i] -= ffmax( 0, qentr2 ); // updating entrainable heights
                            phentrmax3[i] -= ffmax( 0, qentr3 );
                            
                            momaddfsx = momfact * vflowx2 * qentr2;  // contribution of entrainment to momenta
                            momaddfx = momfact * vflowx3 * qentr3;
                            momaddfsy = momfact * vflowy2 * qentr2;
                            momaddfy = momfact * vflowy3 * qentr3;

                            awt[i][10] = awt[i][10] - qentr1; // updating change of basal topography (PHASE 1)

                            aw[i][3] = aw[i][3] + qentr2; // updating flow depths, momenta, and change of basal topography (PHASE 2)
                            aw[i][4] = aw[i][4] + momaddfsx;
                            aw[i][5] = aw[i][5] + momaddfsy;
                            awt[i][11] = awt[i][11] - qentr2; 

                            aw[i][6] = aw[i][6] + qentr3; // updating flow depths, momenta, and change of basal topography (PHASE 3)
                            aw[i][7] = aw[i][7] + momaddfx;
                            aw[i][8] = aw[i][8] + momaddfy;
                            awt[i][12] = awt[i][12] - qentr3;

                        } else awt[i][4] = awt[i][4] - qentr1; // updating basal topography
                    }
                    
                    if ( sico.ENTRAINMENT == 9 && hflow > 0 && qmelt < 0 ) { //melt model !!!CHECK momentum for extremely rapid flows
                        
                        aw[i][0] = aw[i][0] + qmelt1;

                        if ( sico.MODEL == 7 ) {
                        
                            aw[i][3] = aw[i][3] + qmelt2;
                            aw[i][6] = aw[i][6] + qmelt3;
                        }
                    }
                }
            }


// -- STOP --- Entrainment --------------------------------------------------------------------------------------


// -- START -- Stopping of flow ---------------------------------------------------------------------------------


            // *** 1 = Stopping based on fraction of maximum kinetic energy, material is deposited and simulation is terminated
            // ***     when stopping occurs
            // *** 2 = Stopping based on fraction of maximum momentum, material is deposited and simulation is terminated
            // ***     when stopping occurs
            // *** 3 = Stopping based on flow pressure threshold, material is deposited 
            // ***     and simulation is terminated when stopping occurs
            // *** Negative numbers allow for the release of the fluid material along with an equal amount of solid material (multi-phase model only) 


            if ( tsum > 0 && sico.STOPPING != 0 && p == 1 ) { // stopping through fraction of maximum kinetic energy or momentum, or minimum pressure criterion

                hekin_sum = 0; // resetting sum of momentum or kinetic energy or momentum
                ctrl_pressthr = 0; // resetting minimum pressure criterion

                if ( cslide < 2 ) iloop = ib[0]; else iloop = sico.IMAX;
                for ( ix=0; ix<iloop; ix++ ) {

                    if ( cslide < 2 ) i = ibasket[0][ix]; else i = ix;

                    if ( cdomain2[i] == 0 && tsum >= fabs( pstoptime[i] )) {

                        if ( aw[i][0] > sico.HFLOWMIN ) { vflowx = aw[i][1]/aw[i][0]; vflowy = aw[i][2]/aw[i][0]; } else { vflowx = 0; vflowy = 0; } // flow velocities
                        if ( sico.MODEL == 7 && aw[i][3] > sico.HFLOWMIN ) { vflowx2 = aw[i][4]/aw[i][3]; vflowy2 = aw[i][5]/aw[i][3]; } else { vflowx2 = 0; vflowy2 = 0; }
                        if ( sico.MODEL == 7 && aw[i][6] > sico.HFLOWMIN ) { vflowx3 = aw[i][7]/aw[i][6]; vflowy3 = aw[i][8]/aw[i][6]; } else { vflowx3 = 0; vflowy3 = 0; }

                        if ( fabs( sico.STOPPING ) == 1 ) { // fraction of kinetic energy

                            hekin = 0.5 * aw[i][0] * sflow.RHO1 * ( pow( vflowx, 2 ) + pow( vflowy, 2 )); // flow kinetic energies
                            if ( sico.MODEL == 7 ) hekin += 0.5 * aw[i][3] * sflow.RHO2 * ( pow( vflowx2, 2 ) + pow( vflowy2, 2 ));
                            if ( sico.MODEL == 7 ) hekin += 0.5 * aw[i][6] * sflow.RHO3 * ( pow( vflowx3, 2 ) + pow( vflowy3, 2 ));
                            hekin_sum += hekin; // updating sum of kinetic energy

                        } else if ( fabs( sico.STOPPING ) == 2 ) { // fraction of momentum

                            hekin = aw[i][0] * sflow.RHO1 * pow( pow( vflowx, 2 ) + pow( vflowy, 2 ), 0.5 ); // flow momenta
                            if ( sico.MODEL == 7 ) hekin += aw[i][3] * sflow.RHO2 * pow( pow( vflowx2, 2 ) + pow( vflowy2, 2 ), 0.5 );
                            if ( sico.MODEL == 7 ) hekin += aw[i][6] * sflow.RHO3 * pow( pow( vflowx3, 2 ) + pow( vflowy3, 2 ), 0.5 );
                            hekin_sum += hekin; // updating sum of momentum
                            
                        } else if ( fabs( sico.STOPPING ) == 3 ) { // minimum pressure criterion                        
                        
                            if ( sico.MODEL <= 3 ) hflow = aw[i][0]; // flow heights
                            else if ( sico.MODEL == 7 ) hflow = aw[i][0] + aw[i][3] + aw[i][6];

                            hekin = sflow.RHO1 * ( pow( vflowx, 2 ) + pow( vflowy, 2 )); // flow pressures
                            if ( sico.MODEL == 7 ) hekin += sflow.RHO2 * ( pow( vflowx2, 2 ) + pow( vflowy2, 2 ));
                            if ( sico.MODEL == 7 ) hekin += sflow.RHO3 * ( pow( vflowx3, 2 ) + pow( vflowy3, 2 ));
                            if ( hflow > sico.IMPTHR[0] ) { if ( hekin > sflow.CSTOP ) ctrl_pressthr = 1; } // updating minimum pressure criterion
                        }
                    }
                }

                if ( hekin_sum > hekin_max ) hekin_max = hekin_sum; // updating maximum momentum or kinetic energy

                if ( cslide < 2 ) iloop = ib[0]; else iloop = sico.IMAX;
                for ( ix=0; ix<iloop; ix++ ) {

                    if ( cslide < 2 ) i = ibasket[0][ix]; else i = ix;

                    if ( tsum >= fabs( pstoptime[i] )) {

                        if ( pstoptime[i] < 0 || ( fabs( sico.STOPPING ) == 3 && ctrl_pressthr == 0 ) || ( fabs( sico.STOPPING ) != 3 && hekin_sum < hekin_max * sflow.CSTOP )) { 
                            // applying stopping criterion

                            cstopped[i] = 1; // updating control for stopping

                            if ( sico.STOPPING < 0 && sico.MODEL == 7 ) { // stopping with material release
                            
                                for ( j=0; j<imax; j++ ) {

                                    aw[in[i][j]][0] = aw[in[i][j]][6]; aw[in[i][j]][1] = aw[in[i][j]][7]; aw[in[i][j]][2] = aw[in[i][j]][8]; // constraining PHASE 1 !!!CHECK one-phase
                                }
                            }
                        }
                    }
                }
            }


// -- STOP --- Stopping of flow ---------------------------------------------------------------------------------


// -- START -- Updating vectors of composite, maximum and cumulative values, and basal topography ---------------


            if ( cslide < 2 ) iloop = ib[0]; else iloop = sico.IMAX;
            for ( ix=0; ix<iloop; ix++ ) {

                if ( cslide < 2 ) i = ibasket[0][ix]; else i = ix;

                corrfact = 1; // correction factor

                if ( sico.TOPOFACT == 0 ) topofact = 1.0;
                else topofact = cos( betaxy[i] );

                if ( sico.MODEL <= 3 ) { // one-phase model

                    awt[i][0] = aw[i][0] * topofact;
                    awt[i][4] += relev[i] / corrfact;

                    pelev[i] = pelev0[i] + awt[i][4] * corrfact; // correcting elevation for entrainment and deposition

                    awt[i][5] = pow( pow( fdiv( aw[i][1] / cos( betaxy[i] ), aw[i][0] * topofact, sico.IMPTHR[3] ), 2 ) 
                                   + pow( fdiv( aw[i][2] / cos( betaxy[i] ), aw[i][0] * topofact, sico.IMPTHR[3] ), 2 ) , 0.5 ); // velocity
                    awt[i][6] = 0.5 * sflow.RHO1 * aw[i][0] * topofact * pow( awt[i][5], 2 ); // flow kinetic energy
                    awt[i][7] = sflow.RHO1 * pow( awt[i][5], 2 ); // flow pressure

                    if ( awt[i][0] > awt[i][8] ) awt[i][8] = awt[i][0]; // maximum flow height
                    if ( awt[i][5] > awt[i][9] && awt[i][0] >= sico.IMPTHR[3] ) awt[i][9] = awt[i][5]; // maximum flow velocity
                    if ( awt[i][6] > awt[i][10] ) awt[i][10] = awt[i][6]; //maximum flow kinetic energy
                    if ( awt[i][7] > awt[i][11] ) awt[i][11] = awt[i][7]; // maximum flow pressure

                    if ( awt[i][9] > vflow_maxmax ) vflow_maxmax = awt[i][9]; // overall maximum flow velocity
                    if ( awt[i][6] > tflow_maxmax ) tflow_maxmax = awt[i][6]; // overall maximum flow kinetic energy
                    if ( awt[i][7] > pflow_maxmax ) pflow_maxmax = awt[i][7]; // overall maximum flow pressure
                    if ( awt[i][4] > basechange_max ) basechange_max = awt[i][4]; // overall maximum change of basal topography
                    if ( awt[i][4] < basechange_min ) basechange_min = awt[i][4]; // overall minimum change of basal topography
                
                } else if ( sico.MODEL == 7 ) { // multi-phase model

                    awt[i][0] = aw[i][0] * topofact;
                    awt[i][3] = aw[i][3] * topofact;
                    awt[i][6] = aw[i][6] * topofact;

                    awt[i][10] += relev[i] / corrfact;

                    pelev[i] = pelev0[i] + ( awt[i][10] + awt[i][11] + awt[i][12] ) * corrfact; // correcting elevation for entrainment and deposition

                    awt[i][13] = pow( pow( fdiv( aw[i][1] / cos( betaxy[i] ), aw[i][0] * topofact, sico.IMPTHR[3] ), 2 ) 
                                    + pow( fdiv( aw[i][2] / cos( betaxy[i] ), aw[i][0] * topofact, sico.IMPTHR[3] ), 2 ) , 0.5 );
                    awt[i][14] = pow( pow( fdiv( aw[i][4] / cos( betaxy[i] ), aw[i][3] * topofact, sico.IMPTHR[3] ), 2 ) 
                                    + pow( fdiv( aw[i][5] / cos( betaxy[i] ), aw[i][3] * topofact, sico.IMPTHR[3] ), 2 ) , 0.5 );
                    awt[i][15] = pow( pow( fdiv( aw[i][7] / cos( betaxy[i] ), aw[i][6] * topofact, sico.IMPTHR[3] ), 2 ) 
                                    + pow( fdiv( aw[i][8] / cos( betaxy[i] ), aw[i][6] * topofact, sico.IMPTHR[3] ), 2 ) , 0.5 ); // velocities
                    
                    awt[i][16] = awt[i][0] + awt[i][3] + awt[i][6]; // total flow depth
                    
                    awt[i][17] = 0.5 * sflow.RHO1 * aw[i][0] * topofact * pow( awt[i][13], 2 );
                    awt[i][18] = 0.5 * sflow.RHO2 * aw[i][3] * topofact * pow( awt[i][14], 2 );
                    awt[i][19] = 0.5 * sflow.RHO3 * aw[i][6] * topofact * pow( awt[i][15], 2 );
                    awt[i][20] = awt[i][17] + awt[i][18] + awt[i][19]; // flow kinetic energies
                    
                    if ( awt[i][0] > sico.HFLOWMIN ) awt[i][21] = 2.0 * awt[i][17] / ( aw[i][0] * topofact ); else awt[i][21] = 0;
                    if ( awt[i][3] > sico.HFLOWMIN ) awt[i][22] = 2.0 * awt[i][18] / ( aw[i][3] * topofact ); else awt[i][22] = 0;
                    if ( awt[i][6] > sico.HFLOWMIN ) awt[i][23] = 2.0 * awt[i][19] / ( aw[i][6] * topofact ); else awt[i][23] = 0;
                    if ( awt[i][16] > sico.HFLOWMIN ) awt[i][24] = awt[i][20] / awt[i][16]; else awt[i][24] = 0; // flow pressures
                    
                    awt[i][25] = awt[i][10] + awt[i][11] + awt[i][12]; // total entrained or deposited depth

                    if ( awt[i][0] > awt[i][26] ) awt[i][26] = awt[i][0];
                    if ( awt[i][13] > awt[i][27] && awt[i][0] >= sico.IMPTHR[3] ) awt[i][27] = awt[i][13];
                    if ( awt[i][3] > awt[i][28] ) awt[i][28] = awt[i][3];
                    if ( awt[i][14] > awt[i][29] && awt[i][3] >= sico.IMPTHR[3] ) awt[i][29] = awt[i][14];
                    if ( awt[i][6] > awt[i][30] ) awt[i][30] = awt[i][6];
                    if ( awt[i][15] > awt[i][31] && awt[i][6] >= sico.IMPTHR[3] ) awt[i][31] = awt[i][15]; // maximum flow depths and velocities

                    if ( pow( pow( awt[i][1] + awt[i][4] + awt[i][7], 2 ) + pow( awt[i][2] + awt[i][5] + awt[i][8], 2 ), 0.5 ) > awt[i][48] ) {
                    
                        awt[i][48] = pow( pow( awt[i][1] + awt[i][4] + awt[i][7], 2 ) + pow( awt[i][2] + awt[i][5] + awt[i][8], 2 ), 0.5 ); // maximum total momentum
                        awt[i][46] = pow( pow(( awt[i][1] + awt[i][4] + awt[i][7] ) / ( awt[i][0] + awt[i][3] + awt[i][6] ), 2 )
                            + pow(( awt[i][2] + awt[i][5] + awt[i][8] ) / ( awt[i][0] + awt[i][3] + awt[i][6] ), 2 ), 0.5 ); // flow velocity at maximum total momentum
                        awt[i][44] = awt[i][0] / awt[i][16];
                        awt[i][45] = awt[i][6] / awt[i][16]; // phase fractions at maximum total momentum
                    }

                    if ( awt[i][16] > awt[i][32] ) awt[i][32] = awt[i][16]; // maximum total flow depth

                    if ( awt[i][17] > awt[i][33] ) awt[i][33] = awt[i][17];
                    if ( awt[i][21] > awt[i][34] ) awt[i][34] = awt[i][21];
                    if ( awt[i][18] > awt[i][35] ) awt[i][35] = awt[i][18];
                    if ( awt[i][22] > awt[i][36] ) awt[i][36] = awt[i][22];
                    if ( awt[i][19] > awt[i][37] ) awt[i][37] = awt[i][19];
                    if ( awt[i][23] > awt[i][38] ) awt[i][38] = awt[i][23];
                    if ( awt[i][20] > awt[i][39] ) awt[i][39] = awt[i][20];
                    if ( awt[i][24] > awt[i][40] ) awt[i][40] = awt[i][24]; // maximum kinetic energies and flow pressures

                    if ( awt[i][27] > vflow_maxmax1 ) vflow_maxmax1 = awt[i][27];
                    if ( awt[i][29] > vflow_maxmax2 ) vflow_maxmax2 = awt[i][29];
                    if ( awt[i][31] > vflow_maxmax3 ) vflow_maxmax3 = awt[i][31]; // overall maximum velocities

                    if ( awt[i][17] > tflow_maxmax1 ) tflow_maxmax1 = awt[i][17];
                    if ( awt[i][18] > tflow_maxmax2 ) tflow_maxmax2 = awt[i][18];             
                    if ( awt[i][19] > tflow_maxmax3 ) tflow_maxmax3 = awt[i][19]; // overall maximum flow kinetic energies

                    if ( awt[i][21] > pflow_maxmax1 ) pflow_maxmax1 = awt[i][21];
                    if ( awt[i][22] > pflow_maxmax2 ) pflow_maxmax2 = awt[i][22];
                    if ( awt[i][23] > pflow_maxmax3 ) pflow_maxmax3 = awt[i][23]; // overall maximum flow pressures

                    if ( awt[i][25] > basechange_max ) basechange_max = awt[i][25]; // overall maximum change of basal topography
                    if ( awt[i][25] < basechange_min ) basechange_min = awt[i][25]; // overall minimum change of basal topography
                    
                    if ( awt[i][20] > tflow_maxmax ) tflow_maxmax = awt[i][20];
                    if ( awt[i][24] > pflow_maxmax ) pflow_maxmax = awt[i][24];
                    vflow_maxmax = fmax( fmax( awt[i][27], awt[i][29] ), fmax( awt[i][31], vflow_maxmax ));
                }
            }

            /*for ( i=0; i<sico.IMAX; i++ ) { //!!!CHECK edge control

                if ( cdomain2[i] == 1 ) {
                
                    for ( k=0; k<sico.NVECTMIN; k++ ) {
                
                        awsum = 0;
                        elevsum = 0;
                        awnum = 0;
                
                        for ( l=1; l<9; l++ ) {
                    
                            if ( px[i] > 0 && px[i] < sico.M-1 && py[i] > 0 && py[i] < sico.N-1 ) {
                            
                                if ( cdomain2[in[i][l]] == 0 ) { 
                                
                                    awsum += aw[in[i][l]][k];
                                    elevsum += pelev[in[i][l]]; 
                                    awnum += 1;
                                    
                                }
                            }
                        }
                        
                        if ( awnum > 0 ) {
                        
                            aw[i][k] = (float)awsum / (float)awnum;
                            pelev[i] = (float)elevsum / (float)awnum;
                        }
                    }
                }
            }*/


// -- STOP --- Updating vectors of composite, maximum and cumulative values, and basal topography ---------------


// -- START -- Preparing dynamic evolution of orthophoto --------------------------------------------------------


            if ( sico.PBG == 1 ) {
               
                if ( cslide < 2 ) iloop = ib[0]; else iloop = sico.IMAX;
                for ( ix=0; ix<iloop; ix++ ) {

                    i = ibasket[0][ix];
                        
                    if ( cdomain[i] == 5 ) {

                        if ( sico.MODEL <= 3 ) {

                            vflowx = fdiv( aw[i][1], aw[i][0], sico.HFLOWMIN );
                            vflowy = fdiv( aw[i][2], aw[i][0], sico.HFLOWMIN );
                        
                        } else if ( sico.LAYERS == 0 ) {
                        
                            vflowx = fdiv( aw[i][1] + aw[i][4] + aw[i][7], aw[i][0] + aw[i][3] + aw[i][6], sico.HFLOWMIN );
                            vflowy = fdiv( aw[i][2] + aw[i][5] + aw[i][8], aw[i][0] + aw[i][3] + aw[i][6], sico.HFLOWMIN );
                        
                        } else if ( aw[i][6] > sico.HFLOWMIN ) {
                        
                            vflowx = fdiv( aw[i][7], aw[i][6], sico.HFLOWMIN );
                            vflowy = fdiv( aw[i][8], aw[i][6], sico.HFLOWMIN );
                        
                        } else if ( aw[i][3] > sico.HFLOWMIN ) {
                        
                            vflowx = fdiv( aw[i][4], aw[i][3], sico.HFLOWMIN );
                            vflowy = fdiv( aw[i][5], aw[i][3], sico.HFLOWMIN );
                        
                        } else {
                        
                            vflowx = fdiv( aw[i][1], aw[i][0], sico.HFLOWMIN );
                            vflowy = fdiv( aw[i][2], aw[i][0], sico.HFLOWMIN );
                        }                       

                        anpx[i] += vflowx * tlength / sico.CSZ; // updating sliding distances of raster cells in x and y direction
                        anpy[i] += vflowy * tlength / sico.CSZ;
                    }
                }
            }                        


// -- STOP --- Preparing dynamic evolution of orthophoto --------------------------------------------------------


// *** End of loop over two steps, each moving the system half of a cell (NOC scheme) ---------------------------


        }

        if ( cslide < 2 ) iloop = ib[0]; else iloop = sico.IMAX;
        for ( ix=0; ix<iloop; ix++ ) {

            if ( cslide < 2 ) i = ibasket[0][ix]; else i = ix;

            betax[i] = fbeta( pelev[in[i][5]], pelev[in[i][2]], 2.0, sico ); // slopes
            betay[i] = fbeta( pelev[in[i][4]], pelev[in[i][1]], 2.0, sico );
            betaxy[i] = fbetaxy( betax[i], betay[i] );
        }


// -- START -- Updating maximum values of flow parameters and volumes, time of reach ----------------------------


        hflow_max = 0; vflow_max = 0; vol_flow = 0; vol_edge = 0; vol_entr = 0; vol_noflux = 0; ekin_flow = 0;
        hflow_max2 = 0; vflow_max2 = 0; vol_flow2 = 0; vol_edge2 = 0; vol_entr2 = 0; vol_noflux2 = 0;
        hflow_max3 = 0; vflow_max3 = 0; vol_flow3 = 0; vol_edge3 = 0; vol_entr3 = 0; vol_noflux3 = 0;

        for ( z=0; z<nzones; z++ ) { vol_zone1[z] = 0; vol_zone2[z] = 0; vol_zone3[z] = 0; vol_czone1[z] = 0; vol_czone2[z] = 0; vol_czone3[z] = 0; } // zone-specific volumes

        for ( i=0; i<sico.IMAX; i++ ) {

            if ( cdomain2[i] == 0 ) {

                carea = pow ( sico.CSZ, 2 ); // area of cell

                if ( awt[i][0] > hflow_max ) hflow_max = awt[i][0]; // maximum mixture or PHASE 1 flow depth
                if ( cstopped[i] == 0 && aw[i][0] > sico.HFLOWMIN ) { vol_flow += aw[i][0] * carea; vol_zone1[pzones[i]] += aw[i][0] * carea; } 
                    // mixture or PHASE 1 flow volume
                if ( cstopped[i] == 1 && aw[i][0] > sico.HFLOWMIN ) vol_noflux += aw[i][0] * carea; // mixture or PHASE 1 volume without fluxes

                if ( sico.MODEL <= 3 ) { // one-phase models

                    if ( hflow_max > hflow_maxmax ) hflow_maxmax = hflow_max; // absolute maximum flow depth
                    if ( awt[i][0] >= sico.IMPTHR[3] && awt[i][5] > vflow_max ) vflow_max = awt[i][5]; // maximum flow velocity
                    if ( cdomain[i] != 0 && ( sico.ENTRAINMENT != 0 || sico.STOPPING != 0 )) { vol_entr += awt[i][4] * carea; vol_czone1[pzones[i]] += awt[i][4] * carea; }
                        // entrained or deposited volume

                    ekin_flow += awt[i][6] * carea; // kinetic energy of flow
                    
                } else if ( sico.MODEL == 7 ) { // multi-phase model

                    if ( awt[i][3] > hflow_max2 ) hflow_max2 = awt[i][3];
                    if ( awt[i][6] > hflow_max3 ) hflow_max3 = awt[i][6];
                    if ( hflow_max + hflow_max2 + hflow_max3 > hflow_maxmax ) hflow_maxmax = hflow_max + hflow_max2 + hflow_max3; // absolute maximum flow depths

                    if ( cdomain[i] != 0 && awt[i][0] >= sico.IMPTHR[3] && awt[i][13] > vflow_max ) vflow_max = awt[i][13];
                    if ( cdomain[i] != 0 && awt[i][3] >= sico.IMPTHR[3] && awt[i][14] > vflow_max2 ) vflow_max2 = awt[i][14];
                    if ( cdomain[i] != 0 && awt[i][6] >= sico.IMPTHR[3] && awt[i][15] > vflow_max3 ) vflow_max3 = awt[i][15]; // maximum velocities

                    if ( cstopped[i] == 0 && aw[i][3] > sico.HFLOWMIN ) { vol_flow2 += aw[i][3] * carea; vol_zone2[pzones[i]] += aw[i][3] * carea; }
                    if ( cstopped[i] == 0 && aw[i][6] > sico.HFLOWMIN ) { vol_flow3 += aw[i][6] * carea; vol_zone3[pzones[i]] += aw[i][6] * carea; } // flow volumes

                    if ( cstopped[i] == 1 && aw[i][3] > sico.HFLOWMIN ) vol_noflux2 += aw[i][3] * carea;
                    if ( cstopped[i] == 1 && aw[i][6] > sico.HFLOWMIN ) vol_noflux3 += aw[i][6] * carea; // volumes without fluxes

                    if ( cdomain[i] != 0 && ( sico.ENTRAINMENT != 0 || sico.STOPPING != 0 )) { vol_entr += awt[i][10] * carea; vol_czone1[pzones[i]] += awt[i][10] * carea; }
                    if ( cdomain[i] != 0 && ( sico.ENTRAINMENT != 0 || sico.STOPPING != 0 )) { vol_entr2 += awt[i][11] * carea; vol_czone2[pzones[i]] += awt[i][11] * carea; }
                    if ( cdomain[i] != 0 && ( sico.ENTRAINMENT != 0 || sico.STOPPING != 0 )) { vol_entr3 += awt[i][12] * carea; vol_czone3[pzones[i]] += awt[i][12] * carea; }
                        // entrained or deposited volumes

                    ekin_flow += awt[i][20] * carea; // kinetic energy of flow
                }
            }

            if ( cdomain2[i] == 1 && aw[i][0] > sico.HFLOWMIN ) vol_edge += aw[i][0] * carea;
            if ( sico.MODEL == 7 && cdomain2[i] == 1 && aw[i][3] > sico.HFLOWMIN ) vol_edge2 += aw[i][3] * carea;
            if ( sico.MODEL == 7 && cdomain2[i] == 1 && aw[i][6] > sico.HFLOWMIN ) vol_edge3 += aw[i][6] * carea; // flow volumes leaving area of interest
        }

        for ( i=0; i<sico.IMAX; i++ ) { // updating time of reach !!!CHECK

            if ( cdomain2[i] == 0 ) {

                if ( sico.MODEL <= 3 ) hflow = awt[i][0];
                else if ( sico.MODEL == 7 ) hflow = awt[i][0] + awt[i][3] + awt[i][6];
                if ( awt[i][nvect_all-2] == sico.UNDEF * -100 && hflow >= sico.IMPTHR[0] ) {

                    awt[i][nvect_all-8] = pow( pow(( awt[i][1] + awt[i][4] + awt[i][7] ) / ( awt[i][0] + awt[i][3] + awt[i][6] ), 2 )
                        + pow(( awt[i][2] + awt[i][5] + awt[i][8] ) / ( awt[i][0] + awt[i][3] + awt[i][6] ), 2 ), 0.5 ); 
                    awt[i][nvect_all-7] = awt[i][0] / hflow;
                    awt[i][nvect_all-6] = awt[i][3] / hflow; 
                    awt[i][nvect_all-2] = tsum;
                    
                    if ( tsum > treachmaxmax ) treachmaxmax = tsum;
                }
            }
        }


// -- STOP --- Updating maximum values of flow parameters and volumes, time of reach ----------------------------


// -- START -- Managing stopped flows and numerical failures ----------------------------------------------------


        if (( sico.MODEL <= 3 && vol_flow == 0 ) || ( sico.MODEL == 7 && vol_flow == 0 && vol_flow2 == 0 && vol_flow3 == 0 )) {
            // stopping simulation if the entire flow has stopped moving or left the area of interest (successful simulation)

            if ( sico.MULT == 0 ) {

                printf("\n                       \n                       ------------------ FLOW STOPPED ------------------\n\n");
                fflush(stdout); // forcing immediate display
            }
            
            fprintf(f_summary, "\n                       \n                       ------------------ FLOW STOPPED ------------------\n\n");
            ccontinue = 0; // setting control for continuation to negative
        }

        if ( ccontinue == 0 && csuccess == 0 ) {
            // stopping simulation if time step length indicates numerical failure (unsuccessful simulation)

            if ( sico.MULT == 0 ) {

                printf("\n                       \n                       --------------- NUMERICAL FAILURE ----------------\n\n");
                fflush(stdout); // forcing immediate display
            }
            
            fprintf(f_summary, "\n                       \n                       --------------- NUMERICAL FAILURE ----------------\n\n");
        }

        if ( sico.STOPPING != 0 && ccontinue == 0 && csuccess == 1 ) { // correcting deposited depth and basal topography for stopping

            if ( cslide < 2 ) iloop = ib[0]; else iloop = sico.IMAX;
            for ( ix=0; ix<iloop; ix++ ) {

                if ( cslide < 2 ) i = ibasket[0][ix]; else i = ix;

                corrfact = 1; // correction factor

                if ( sico.MODEL <= 3 ) { // one-phase models

                    awt[i][4] += cstopped[i] * aw[i][0]; // change of basal surface
                    if ( awt[i][4] > basechange_max ) basechange_max = awt[i][4]; // overall maximum change of basal topography
                    
                    pelev[i] += cstopped[i] * aw[i][0] * corrfact; // elevation
                    aw[i][0] -= cstopped[i] * aw[i][0]; // flow depth
                    awt[i][0] -= cstopped[i] * awt[i][0]; // flow thickness
                    
                } else if ( sico.MODEL == 7 ) { // multi-phase model

                    awt[i][10] += cstopped[i] * aw[i][0];
                    awt[i][11] += cstopped[i] * aw[i][3];
                    awt[i][12] += cstopped[i] * aw[i][6];
                    awt[i][25] = awt[i][10] + awt[i][11] + awt[i][12]; // change of basal surface
                    if ( awt[i][25] > basechange_max ) basechange_max = awt[i][25]; // overall maximum change of basal topography
                    
                    pelev[i] += ( cstopped[i] * aw[i][0] + cstopped[i] * aw[i][3] + cstopped[i] * aw[i][6] ) * corrfact; // elevation
                    aw[i][0] -= cstopped[i] * aw[i][0];
                    aw[i][3] -= cstopped[i] * aw[i][3];
                    aw[i][6] -= cstopped[i] * aw[i][6];

                    awt[i][0] -= cstopped[i] * awt[i][0];
                    awt[i][3] -= cstopped[i] * awt[i][3];
                    awt[i][6] -= cstopped[i] * awt[i][6]; // flow thicknesses
                }
            }
        }


// -- STOP --- Managing stopped flows and numerical failures ----------------------------------------------------


        if ( (int)(round( 1000 * tint )) >= (int)( round( 1000 * tout )) || ccontinue == 0 ) { // if defined interval or last time step is reached


// -- START -- Writing hydrograph infos to files ----------------------------------------------------------------


            if ( sico.MULT == 0 && hydrograph == 1 ) {

                ctrl_hydout = 1; // control for hydrograph output

                for ( hydj = hydnin; hydj < hydnin + hydnout; hydj++ ) {  // loop over all output hydrographs

                    hydh = awt[hydi[hydj]][0]; // PHASE 1 flow height
                    if ( hydh > sico.HFLOWMIN ) hydv = pow( pow( aw[hydi[hydj]][1] / aw[hydi[hydj]][0] / cos( betaxy[hydi[hydj]] ), 2 ) 
                            + pow( aw[hydi[hydj]][2] / aw[hydi[hydj]][0] / cos( betaxy[hydi[hydj]] ), 2 ), 0.5 );
                    else hydv= 0; // PHASE 1 flow velocity

                    if ( sico.MODEL <= 3 ) {

                        hyde = awt[hydi[hydj]][4]; // entrained height
                        hydh2 = 0; hydv2 = 0; hyde2 = 0; hydh3 = 0; hydv3 = 0; hyde3 = 0;
                        
                    } else if ( sico.MODEL == 7 ) {

                        hydh2 = awt[hydi[hydj]][3]; // PHASE 2 flow height
                        hydh3 = awt[hydi[hydj]][6]; // PHASE 3 flow height
                        if ( hydh2 > sico.HFLOWMIN ) hydv2 = pow( pow( aw[hydi[hydj]][4] / aw[hydi[hydj]][3] / cos( betaxy[hydi[hydj]] ), 2 ) 
                            + pow( aw[hydi[hydj]][5] / aw[hydi[hydj]][3] / cos( betaxy[hydi[hydj]] ), 2 ), 0.5 );
                        else hydv2= 0; // PHASE 2 flow velocity
                        if ( hydh3 > sico.HFLOWMIN ) hydv3 = pow( pow( aw[hydi[hydj]][7] / aw[hydi[hydj]][6] / cos( betaxy[hydi[hydj]] ), 2 ) 
                            + pow( aw[hydi[hydj]][8] / aw[hydi[hydj]][6] / cos( betaxy[hydi[hydj]] ), 2 ), 0.5 );
                        else hydv3= 0; // PHASE 3 flow velocity
                        hyde = awt[hydi[hydj]][10]; // change of basal surface mixture or PHASE 1
                        hyde2 = awt[hydi[hydj]][11]; // change of basal surface PHASE 2
                        hyde3 = awt[hydi[hydj]][12]; // change of basal surface PHASE 3
                    }

                    hydq = 0; hydq2 = 0; hydq3 = 0; // resetting discharges

                    if ( hydalpha[hydj] == 0 || hydalpha[hydj] == sico.PI * 0.5 || hydalpha[hydj] == sico.PI || hydalpha[hydj] == 3 * sico.PI * 0.5 ) hydfalpha = 1;
                    else if ( fabs( 1 / sin( hydalpha[hydj] )) < fabs( 1 / cos( hydalpha[hydj] ))) hydfalpha = fabs( 1 / sin( hydalpha[hydj] ));
                    else hydfalpha = fabs( 1 / cos( hydalpha[hydj] )); // correction factor for profile direction

                    for ( hydk = 1; hydk <= hydp[0][hydj]; hydk++ ) { // loop over all cells of profile

                        alpha = falpha( betax[hydp[hydk][hydj]], betay[hydp[hydk][hydj]], sico ); // aspect
                        hydbeta = atan ( tan( betaxy[hydp[hydk][hydj]] ) * cos ( alpha - hydalpha[hydj] + sico.PI * 0.5 )); // corrected slope

                        hydfcorr = sico.CSZ * hydfalpha / cos( hydbeta ); // reference length for discharge

                        hydm0 = pow( pow(aw[hydp[hydk][hydj]][1], 2 ) + pow( aw[hydp[hydk][hydj]][2], 2), 0.5 );

                        if ( aw[hydp[hydk][hydj]][0] > sico.HFLOWMIN && hydm0 > 0 ) {

                            hydmx = aw[hydp[hydk][hydj]][1];
                            hydmy = aw[hydp[hydk][hydj]][2];

                            if ( hydmy > 0 ) hydmalpha = acos( hydmx / hydm0 ); else hydmalpha = 2 * sico.PI - acos( hydmx / hydm0 );
                            hydnalpha = hydalpha[hydj] - sico.PI * 0.5; if ( hydnalpha < 0 ) hydnalpha += 2 * sico.PI;
                            hydm = hydm0 * cos( hydnalpha - hydmalpha ); // momentum in direction perpendicular to profile
                            hydq += hydm * hydfcorr; // updating mixture or PHASE 1 discharge
                        }

                        if ( sico.MODEL == 7 ) {

                            hydm0 = pow( pow(aw[hydp[hydk][hydj]][4], 2 ) + pow( aw[hydp[hydk][hydj]][5], 2), 0.5 );

                            if ( aw[hydp[hydk][hydj]][3] > sico.HFLOWMIN && hydm0 > 0 ) {

                                hydmx = aw[hydp[hydk][hydj]][4];
                                hydmy = aw[hydp[hydk][hydj]][5];

                                if ( hydmy > 0 ) hydmalpha = acos( hydmx / hydm0 ); else hydmalpha = 2 * sico.PI - acos( hydmx / hydm0 );
                                hydnalpha = hydalpha[hydj] - sico.PI * 0.5; if ( hydnalpha < 0 ) hydnalpha += 2 * sico.PI;
                                hydm = hydm0 * cos ( hydnalpha - hydmalpha ); // momentum in direction perpendicular to profile
                                hydq2 += hydm * hydfcorr; // updating PHASE 2 discharge
                            }

                            hydm0 = pow( pow(aw[hydp[hydk][hydj]][7], 2 ) + pow( aw[hydp[hydk][hydj]][8], 2), 0.5 );

                            if ( aw[hydp[hydk][hydj]][6] > sico.HFLOWMIN && hydm0 > 0 ) {

                                hydmx = aw[hydp[hydk][hydj]][7];
                                hydmy = aw[hydp[hydk][hydj]][8];

                                if ( hydmy > 0 ) hydmalpha = acos( hydmx / hydm0 ); else hydmalpha = 2 * sico.PI - acos( hydmx / hydm0 );
                                hydnalpha = hydalpha[hydj] - sico.PI * 0.5; if ( hydnalpha < 0 ) hydnalpha += 2 * sico.PI;
                                hydm = hydm0 * cos ( hydnalpha - hydmalpha ); // momentum in direction perpendicular to profile
                                hydq3 += hydm * hydfcorr; // updating PHASE 3 discharge
                            }
                        }
                    }

                    fprintf(f_hydinfo[hydj], "%.1f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
                        tsum, hydh, hydv, hyde, hydq, hydh2, hydv2, hyde2, hydq2, hydh3, hydv3, hyde3, hydq3); // writing hydrograph info to file

                    if ( hydq != 0 || hydq2 != 0 || hydq3 != 0 ) {

                        if ( qtinit[hydj] == sico.UNDEF ) qtinit[hydj] = tsum;

                        fprintf(f_hydtrans[hydj], "%.1f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
                        tsum - qtinit[hydj], hydq, hydv, hydq2, hydv2, hydq3, hydv3 ); // writing hydrograph info to file in input hydrograph format
                    }
                }
            }


// -- STOP --- Writing hydrograph infos to files ----------------------------------------------------------------


// -- START -- Writing profile ----------------------------------------------------------------------------------


            if ( sico.PROFILE > 0 ) {
            
                 for ( profi = 0; profi <= profctrl; profi++ ) {

                     profwhtx1 = 1 - ( profnx[profi] - (int)profnx[profi] );
                     profwhtx2 = 1 - ((int)profnx[profi] + 1 - profnx[profi] );
                     profwhty1 = 1 - ( profny[profi] - (int)profny[profi] );
                     profwhty2 = 1 - ((int)profny[profi] + 1 - profny[profi] );
              
                     for ( profl = 0; profl < profmax; profl++ ) {

                         if ( sico.MODEL <= 3 ) {
                     
                             if ( profl == 0 ) profk = 0; else if ( profl == 1 ) profk = 5; else if ( profl == 2 ) profk = 6; else if ( profl == 3 ) profk = 7; else profk = 4;
                         
                         } else if ( sico.MODEL == 7 ) { 

                             if ( profl == 0 ) profk = 0; else if ( profl == 1 ) profk = 3; else if ( profl == 2 ) profk = 6; else if ( profl == 3 ) profk = 13; 
                             else if ( profl == 4 ) profk = 14;
                             else if ( profl == 5 ) profk = 15; else if ( profl == 6 ) profk = 17; else if ( profl == 7 ) profk = 18; else if ( profl == 8 ) profk = 19; 
                             else if ( profl == 9 ) profk = 21; else if ( profl == 10 ) profk = 22; else if ( profl == 11 ) profk = 23; else if ( profl == 12 ) profk = 10; 
                             else if ( profl == 13 ) profk = 11; else profk = 12;
                         }
                         
                         profdata[profi][nout][profl] = awt[icor[(int)profnx[profi]][(int)profny[profi]]][profk] * profwhtx1 * profwhty1
                             + awt[icor[(int)profnx[profi]][(int)profny[profi]+1]][profk] * profwhtx1 * profwhty2
                             + awt[icor[(int)profnx[profi]+1][(int)profny[profi]]][profk] * profwhtx2 * profwhty1
                             + awt[icor[(int)profnx[profi]+1][(int)profny[profi]+1]][profk] * profwhtx2 * profwhty2;
                     }
                 }
             }


// -- STOP --- Writing profile ----------------------------------------------------------------------------------


// -- START -- Writing control point data -----------------------------------------------------------------------


            if ( sico.CTRLPOINTS > 0 ) {

                for ( j=0; j<sico.CTRLPOINTS; j++ ) {

                    if ( sico.MODEL <= 3 ) {

                        ctrlpdata[j][nout][0] = awt[icor[ctrlpx[j]][ctrlpy[j]]][0];
                        ctrlpdata[j][nout][1] = awt[icor[ctrlpx[j]][ctrlpy[j]]][4];
                    
                    } else if ( sico.MODEL == 7 ) {

                        ctrlpdata[j][nout][0] = awt[icor[ctrlpx[j]][ctrlpy[j]]][0];
                        ctrlpdata[j][nout][1] = awt[icor[ctrlpx[j]][ctrlpy[j]]][3];
                        ctrlpdata[j][nout][2] = awt[icor[ctrlpx[j]][ctrlpy[j]]][6]; 
                        ctrlpdata[j][nout][3] = awt[icor[ctrlpx[j]][ctrlpy[j]]][10];
                        ctrlpdata[j][nout][4] = awt[icor[ctrlpx[j]][ctrlpy[j]]][11];
                        ctrlpdata[j][nout][5] = awt[icor[ctrlpx[j]][ctrlpy[j]]][12];
                    }
                }
            }


// -- START -- Writing csv files for virtual reality ------------------------------------------------------------


            if ( sico.MULT == 0 ) {

                if ( nout == 1 ) jmin = 0; else jmin = 1;

                for ( j=jmin; j<2; j++ ) {

                    if ( sico.PBG == 1 && paradistort >= 1 ) {
                
                        if ( cslide < 2 ) iloop = ib[0]; else iloop = sico.IMAX;

                        for ( ix=0; ix<iloop; ix++ ) {

                            i = ibasket[0][ix];
                        
                            if ( sico.MODEL <= 3 ) hflownn[0] = aw[i][0];
                            else hflownn[0] = aw[i][0] + aw[i][3] + aw[i][6];
                        
                            if ( hflownn[0] > paramin ) {

                                if ( j > 0 ) {
                                
                                    for ( l=1; l<6; l++ ) {
                                    
                                        if ( sico.MODEL <= 3 || ( sico.LAYERS > 0 && aw[in[i][l]][3] < sico.IMPTHR[0] && aw[in[i][l]][6] < sico.IMPTHR[0] ) 
                                            || ( sico.LAYERS > 0 && sflow.PHI0[1] == 0 && sflow.PHI0[2] == 0 )) {
                        
                                            wu[l] = fdiv( aw[in[i][l]][1], aw[in[i][l]][0], sico.HFLOWMIN );
                                            wv[l] = fdiv( aw[in[i][l]][2], aw[in[i][l]][0], sico.HFLOWMIN );
                                            
                                        } else if ( sico.LAYERS == 0 ) {

                                            wu[l] = fdiv( aw[in[i][l]][1] + aw[in[i][l]][4] + aw[in[i][l]][7], aw[in[i][l]][0] + aw[in[i][l]][3] + aw[in[i][l]][6], 
                                                sico.HFLOWMIN );
                                            wv[l] = fdiv( aw[in[i][l]][2] + aw[in[i][l]][5] + aw[in[i][l]][8], aw[in[i][l]][0] + aw[in[i][l]][3] + aw[in[i][l]][6], 
                                                sico.HFLOWMIN );                                        
                                        
                                        } else if ( aw[in[i][l]][6] >= sico.IMPTHR[0] && sflow.PHI0[2] != 0 ) {

                                            wu[l] = fdiv( aw[in[i][l]][7], aw[in[i][l]][6], sico.HFLOWMIN );
                                            wv[l] = fdiv( aw[in[i][l]][8], aw[in[i][l]][6], sico.HFLOWMIN );                                        
                                        
                                        } else {
                                        
                                            wu[l] = fdiv( aw[in[i][l]][4], aw[in[i][l]][3], sico.HFLOWMIN );
                                            wv[l] = fdiv( aw[in[i][l]][5], aw[in[i][l]][3], sico.HFLOWMIN );
                                        }
                                    }
                        
                                    slidefract = pow( pow(( wu[5] - wu[2] ) * betax[i], 2 ) + pow(( wv[4] - wv[1] ) * betay[i], 2 ), 1.0 ) / sico.CSZ;

                                    if ( paradistort == 2 ) {

                                        randnum2 = 0;
                                        randnum3 = 0;                                    
                                        randnum4 = 0;
                               
                                    } else {
                               
                                        randnum2 = 0.0002 * ((rand()%randmax) - 500 );
                                        randnum3 = 0.0001 * ((rand()%randmax) - 500 );
                                        randnum4 = fmax(0, fmin( 1, ( slidefract - 0.02 ) * 2.0 * tout + randnum2 ));
                                    }

                                    ppbg1[i] = ( 1 - randnum4 ) * ppbg1[i] + ( randnum4 ) * (int)( 255 * ( parar + randnum3 ));
                                    ppbg2[i] = ( 1 - randnum4 ) * ppbg2[i] + ( randnum4 ) * (int)( 255 * ( parag ));
                                    ppbg3[i] = ( 1 - randnum4 ) * ppbg3[i] + ( randnum4 ) * (int)( 255 * ( parab ));
                                }                                  

                                anpx1 = (int)( anpx[i] ); anpx2 = (int)( anpx[i] ) + 1; anpy1 = (int)( anpy[i] ); anpy2 = (int)( anpy[i] ) + 1; // nearest neighbour cells
                                anwht[0] = ( anpx2 - anpx[i] ) * ( anpy2 - anpy[i] ); anwht[1] = ( anpx2 - anpx[i] ) * ( anpy[i] - anpy1 );
                                anwht[2] = ( anpx[i] - anpx1 ) * ( anpy2 - anpy[i] ); anwht[3] = ( anpx[i] - anpx1 ) * ( anpy[i] - anpy1 );

                                andist = ffmax( abs(anpx[i] - px[i]) + 2, abs(anpy[i] - py[i]) + 2 );

                                iz = ffmax( 1, ( px[i] - andist ) * sico.N ) + ffmax( 1, py[i] - andist ); // loop over all relevant cells
                                for ( x=ffmax(1,px[i]-andist); x<ffmin((sico.M-1),px[i]+andist+1); x++ ) { 
                                    for ( y=ffmax(1,py[i]-andist); y<ffmin((sico.N-1),py[i]+andist+1); y++ ) {

                                        anid = -1;

                                        if ( px[iz] == anpx1 && py[iz] == anpy1 ) { anid = 0; anctr[iz] += anwht[0]; } // searching for target cells
                                        if ( px[iz] == anpx1 && py[iz] == anpy2 ) { anid = 1; anctr[iz] += anwht[1]; }
                                        if ( px[iz] == anpx2 && py[iz] == anpy1 ) { anid = 2; anctr[iz] += anwht[2]; }
                                        if ( px[iz] == anpx2 && py[iz] == anpy2 ) { anid = 3; anctr[iz] += anwht[3]; }

                                        if ( anid >= 0 ) {

                                            ppbgt1[iz] += (int)( ppbg1[i] * anwht[anid] ); // assigning r,g, and b values to target cells
                                            ppbgt2[iz] += (int)( ppbg2[i] * anwht[anid] ); 
                                            ppbgt3[iz] += (int)( ppbg3[i] * anwht[anid] );
 
                                        }
                                        iz += 1;
                                    }
                                    iz += ffmax( 2, sico.N - 2 * andist - 1 );
                                }
                            }
                        }

                        for ( ix=0; ix<iloop; ix++ ) {

                            i = ibasket[0][ix];
                        
                            if (( aw[i][0] + aw[i][3] + aw[i][6] ) > paramin ) {

                                if ( anctr[i] > 0 && ( ppbgt1[i] > 0 || ppbgt2[i] > 0 || ppbgt3[i] > 0 )) {

                                    ppbgx1[i] = (int)( ppbgt1[i] / anctr[i] ); // scaling r,g, and b values to the range 0-255
                                    ppbgx2[i] = (int)( ppbgt2[i] / anctr[i] ); 
                                    ppbgx3[i] = (int)( ppbgt3[i] / anctr[i] );
                                }
                              
                                ppbgt1[i] = 0; ppbgt2[i] = 0; ppbgt3[i] = 0; anctr[i] = 0; 
                                
                                if ( ctr_ansteps == ansteps ) { 
                                
                                    anpx[i] = px[i]; anpy[i] = py[i];
                                
                                    highpass_min = 255;
                                    highpass_max = 0;
                                    
                                    for ( k=0; k<9; k++ ) {
                                    
                                        highpass_sum = (int)((float)( ppbgx1[in[i][k]] + ppbgx2[in[i][k]] + ppbgx3[in[i][k]] ) / 3 + 0.5 );
                                        if ( highpass_sum < highpass_min ) highpass_min = highpass_sum;
                                        if ( highpass_sum > highpass_max ) highpass_max = highpass_sum;
                                    }
                                    
                                    highpass_halfdiff = ( highpass_max - highpass_min ) * 0.5;
                                    highpass_avg = ( highpass_max + highpass_min ) * 0.5;
                                
                                    ppbg1[i] = (int)( fmin( 255, fmax( 0, ppbgx1[i] + highpass * ( ppbgx1[i] - highpass_avg ) / highpass_halfdiff )));
                                    ppbg2[i] = (int)( fmin( 255, fmax( 0, ppbgx2[i] + highpass * ( ppbgx2[i] - highpass_avg ) / highpass_halfdiff ))); 
                                    ppbg3[i] = (int)( fmin( 255, fmax( 0, ppbgx3[i] + highpass * ( ppbgx3[i] - highpass_avg ) / highpass_halfdiff )));   
                                        // scaling r,g, and b values to the range 0-255                            
                                }
                            }
                        }
                    
                        if ( ctr_ansteps == ansteps ) ctr_ansteps = 1; else ctr_ansteps += 1;
                    }

                    if ( nout < 10 ) sprintf( madd, "000"); // fill string (for maintaining correct order in list of maps)
                    else if ( nout < 100 ) sprintf( madd, "00");
                    else if ( nout < 1000 ) sprintf( madd, "0");

                    sprintf(path, "%sdata/pv%s%i.csv", outvr, madd, nout+j-1);
                    f_vr[nout+j-1]=fopen(path, "w");

                    fprintf(f_vr[nout+j-1], "x,y,z,zb,z1,r,g,b,i\n");

                    if ( sico.PBG == 0 ) {

                        for ( i=0; i<sico.IMAX; i++ ) {

                            alpha = falpha( betax[i], betay[i], sico );
                            ppbg1[i] = (int)( 255.0 * ( cos( 45 * sico.PI / 180.0 ) * cos( betaxy[i] ) 
                                + sin( 45 * sico.PI / 180.0 ) * sin( betaxy[i] ) * cos( 135 * sico.PI / 180.0 - alpha )));
                            ppbg2[i] = ppbg1[i]; ppbg3[i] = ppbg1[i];
                        }
                    }

                    for ( i=0; i<sico.IMAX; i++ ) {
                
                        paraxmetric = ( float ) py[i] * sico.CSZ + sico.BDWEST;
                        paraymetric = ( float ) sico.BDNORTH - px[i] * sico.CSZ;

                        if ( j == 1 && sico.MODEL <= 3 ) {
                    
                            parahflow = aw[i][0];
                            parahrelease = phrelease[0];
                            parahs = pelev[i] + aw[i][0];
                            parahmax = awt[i][8];
                    
                        } else if ( j == 1 && sico.MODEL == 7 ) {
                    
                            parahflow = aw[i][0] + aw[i][3] + aw[i][6];
                            parahrelease = phrelease[i] + phrelease2[i] + phrelease3[i];
                            parahs = pelev[i] + aw[i][0] + aw[i][3];
                            parahmax = awt[i][32];
                            parahflow1 = aw[i][0];
                            parahflow2 = aw[i][3];
                            parahflow3 = aw[i][6];
                            
                        } else if ( j == 0 && sico.MODEL <= 3 ) {
                    
                            parahflow = phrelease[i];
                            parahrelease = phrelease[i];
                            parahs = pelev0[i] + phrelease[i];
                            parahmax = phrelease[i];
                    
                        } else if ( j == 0 && sico.MODEL == 7 ) {
                    
                            parahflow = phrelease[i] + phrelease2[i] + phrelease3[i];
                            parahrelease = phrelease[i] + phrelease2[i] + phrelease3[i];
                            parahs = pelev0[i] + phrelease[i] + phrelease2[i];
                            parahmax = phrelease[i] + phrelease2[i] + phrelease3[i];
                            parahflow1 = phrelease[i];
                            parahflow2 = phrelease2[i];
                            parahflow3 = phrelease3[i];                                                       
                        }

                        if ( j == 1 ) paraelev = pelev[i]; else paraelev = pelev0[i];
                        parah = paraelev + parahflow;                          

                        if ( parahflow >= paramin ) paraalpha = pow( fmin( 1, parahflow / pararef ), parad );
                        else paraalpha = 0;
        
                        if( parahmax >= paramin ) paraalphamax = fmin( 0.35, pow( fmin( 1, parahmax / pararef ), parad ));
                        else paraalphamax = 0;
        
                        if ( sico.LAYERS0 > 0 && sico.MODEL == 7 && sico.TSUNAMI == 0 ) {

                            if ( parahflow >= paramin ) {

                                if ( sico.PBG == 1 && paradistort >= 1 ) {

                                    parared = (float)ppbgx1[i] / 255;                       
                                    paragreen = (float)ppbgx2[i] / 255;                 
                                    parablue = (float)ppbgx3[i] / 255;

                                } else {
                       
                                    if( parahflow2 >= paramin || parahflow3 >= paramin ) parafactr = 0.25; else parafactr = 0.5;
                                    if( parahflow2 >= paramin || parahflow3 >= paramin ) parafactg = 0.25; else parafactg = 0.5;
                                    if( parahflow3 >= paramin ) parafactb = 0.5; else parafactb = 0.25;
                            
                                    parared = parafactr * paraalpha + (float)ppbg1[i] / 255 * ( 1 - paraalpha );
                                    paragreen = parafactg * paraalpha + (float)ppbg2[i] / 255 * ( 1 - paraalpha );
                                    parablue = parafactb * paraalpha + (float)ppbg3[i] / 255 * ( 1 - paraalpha );
                                }
                      
                            } else if ( parahmax > paramin ) {

                                if ( sico.PBG == 1 && paradistort >= 1 ) {
                        
                                    parared = parar;                       
                                    paragreen = parag;                 
                                    parablue = parab;
                                    
                                } else {
                       
                                    parared = 0.7 * paraalphamax + (float)ppbg1[i] / 255 * ( 1 - paraalphamax);                       
                                    paragreen = 0.3 * paraalphamax + (float)ppbg2[i] / 255 * ( 1 - paraalphamax);                 
                                    parablue = 0.0 * paraalphamax + (float)ppbg3[i] / 255 * ( 1 - paraalphamax);
                                }
                                                   
                            } else {

                                parared = (float)ppbg1[i] / 255;                       
                                paragreen = (float)ppbg2[i] / 255;                 
                                parablue = (float)ppbg3[i] / 255;
                            }

                        } else if ( sico.MODEL == 7 ) {
    
                            if ( sico.TSUNAMI == 1 && parahflow3 >= paramin ) {
                            
                                paracorrtsun = 0.5 * ( parahflow + paraelev - parahrelease - pelev0[i] ) / paratsunref;
                                paraaddtsun = 0.5;

                            } else {
        
                                paracorrtsun = 0;
                                paraaddtsun = 0;
                            }

                            if ( parahflow >= paramin ) {

                                if ( sico.PBG == 1 && aw[i][6] <= paramin ) {

                                    parared = (float)ppbgx1[i] / 255;                       
                                    paragreen = (float)ppbgx2[i] / 255;                 
                                    parablue = (float)ppbgx3[i] / 255;

                                } else if ( sico.PBG == 1 ) {
                        
                                    parared = fmax( 0.0, fmin(1.0, ( paraaddtsun + 0 + paracorrtsun ))) * paraalpha + (float)ppbg1[i] / 255 * ( 1 - paraalpha );
                                    paragreen = fmax( 0.0, fmin( 1.0, (paraaddtsun + 0 + paracorrtsun ))) * paraalpha + (float)ppbg2[i] / 255 * ( 1 - paraalpha );
                                    parablue = 1 * paraalpha + (float)ppbg3[i] / 255 * ( 1 - paraalpha );

                                } else {
                        
                                    parared = fmax( 0.0, fmin(1.0, ( paraaddtsun + parahflow1 / parahflow + paracorrtsun ))) * paraalpha 
                                        + (float)ppbg1[i] / 255 * ( 1 - paraalpha );
                                    paragreen = fmax( 0.0, fmin( 1.0, (paraaddtsun + parahflow2 / parahflow + paracorrtsun ))) * paraalpha 
                                        + (float)ppbg2[i] / 255 * ( 1 - paraalpha );
                                    parablue = parahflow3 / parahflow * paraalpha + (float)ppbg3[i] / 255 * ( 1 - paraalpha );
                                }
                      
                            } else if ( parahmax > paramin ) {
                        
                                if ( sico.PBG == 1 ) {
                        
                                    parared = parar;                       
                                    paragreen = parag;                 
                                    parablue = parab;
                                    
                                } else {
                        
                                    parared = 0.7 * paraalphamax + (float)ppbg1[i] / 255 * ( 1 - paraalphamax);                       
                                    paragreen = 0.3 * paraalphamax + (float)ppbg2[i] / 255 * ( 1 - paraalphamax);                 
                                    parablue = 0.0 * paraalphamax + (float)ppbg3[i] / 255 * ( 1 - paraalphamax);
                                }
                                                   
                            } else {

                                parared = (float)ppbg1[i] / 255;                       
                                paragreen = (float)ppbg2[i] / 255;                 
                                parablue = (float)ppbg3[i] / 255;
                            }

                        } else {
    
                            if ( sico.GLACIER == 1 ) {
            
                                paracolfactr = 1.00;
                                paracolfactg = 1.00;
                                paracolfactb = 1.00;
                
                            } else {
            
                                paracolfactr=parar;
                                paracolfactg=parag;
                                paracolfactb=parab;
                            }

                            if ( parahflow >= paramin ) {

                                if ( sico.PBG == 1 ) {

                                    parared = (float)ppbgx1[i] / 255;                       
                                    paragreen = (float)ppbgx2[i] / 255;                 
                                    parablue = (float)ppbgx3[i] / 255;

                                } else {

                                    parared = paracolfactr * paraalpha + (float)ppbg1[i] / 255 * ( 1 - paraalpha );
                                    paragreen = paracolfactg * paraalpha + (float)ppbg2[i] /255 * ( 1 - paraalpha );
                                    parablue = paracolfactb * paraalpha + (float)ppbg3[i] /255 * ( 1 - paraalpha );
                                }
                      
                            } else if ( parahmax > paramin ) {

                                if ( sico.PBG == 1 ) {

                                    parared = parar;                       
                                    paragreen = parag;                 
                                    parablue = parab;

                                } else {

                                    parared = 0.7 * paraalphamax + (float)ppbg1[i] / 255 * ( 1 - paraalphamax);                       
                                    paragreen = 0.3 * paraalphamax + (float)ppbg2[i] / 255 * ( 1 - paraalphamax);                 
                                    parablue = 0.0 * paraalphamax + (float)ppbg3[i] / 255 * ( 1 - paraalphamax);
                                }
                                                   
                            } else {

                                parared = (float)ppbg1[i] / 255;                       
                                paragreen = (float)ppbg2[i] / 255;                 
                                parablue = (float)ppbg3[i] / 255;
                            }                    
                        }

                        if ( parahmax >= sico.IMPTHR[0] ) paraimp = 1; else paraimp = 0;

                        fprintf(f_vr[nout+j-1], "%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%i\n", paraxmetric, paraymetric, parah, 
                            paraelev, parahs, parared, paragreen, parablue, paraimp);                     
                    }
                
                    fclose(f_vr[nout+j-1]);
                }
            }


// -- STOP --- Writing csv files for virtual reality ------------------------------------------------------------


// -- STOP -- Writing control point data ------------------------------------------------------------------------


// -- START -- Display and files of status of simulation --------------------------------------------------------


            if ( sico.MODEL <= 3 ) { // one-phase models

                if ( sico.MULT == 0 ) {

                    printf("   %i\t%i\t%.3f\t%.1f\t%.1f\t%.*f\t%.2f\t%.*f\t%.*f\n", nout, nsum, cflmax, 100*tlength, tsum, prec_hflow, hflow_max, vflow_max,
                        prec_vol, vol_flow/1000, prec_ekin, ekin_flow/1000000); // display

                    fflush(stdout);
                }

                fprintf(f_summary, "%i\t%i\t%.3f\t%.1f\t%.1f\t%.*f\t%.2f\t%.*f\t%.*f\n", nout, nsum, cflmax, 100*tlength, tsum, prec_hflow, hflow_max, vflow_max,
                    prec_vol, vol_flow, prec_ekin, ekin_flow); // summary file
                    
            } else if ( sico.MODEL == 7 ) { // multi-phase model

                if ( sico.MULT == 0 ) {

                    printf("   %i\t%i\t%.3f\t%.1f\t%.1f\t%.*f\t%.2f\t%.*f\t%.2f\t%.*f\t%.2f\t%.*f\t%.*f\t%.*f\t%.*f\n",
                        nout, nsum, cflmax, 100*tlength, tsum, prec_hflow, hflow_max, vflow_max, prec_hflow, hflow_max2, vflow_max2,
                        prec_hflow, hflow_max3, vflow_max3, prec_vol, vol_flow/1000, prec_vol, vol_flow2/1000, prec_vol, vol_flow3/1000, prec_ekin, ekin_flow/1000000); 
                            // display

                    fflush(stdout);
                }
                
                fprintf(f_summary, "%i\t%i\t%.3f\t%.1f\t%.1f\t%.*f\t%.2f\t%.*f\t%.2f\t%.*f\t%.2f\t%.*f\t%.*f\t%.*f\t%.*f\n",
                    nout, nsum, cflmax, 100*tlength, tsum, prec_hflow, hflow_max, vflow_max, prec_hflow, hflow_max2, vflow_max2,
                    prec_hflow, hflow_max3, vflow_max3, prec_vol, vol_flow, prec_vol, vol_flow2, prec_vol, vol_flow3, prec_ekin, ekin_flow); // summary file
            }

            fprintf(f_volumes, "%i\t%.1f", nout, tsum ); // volumes file
    
            for ( z=0; z<nzones; z++ ) {
    
                fprintf(f_volumes, "\t%.3f\t%.3f", vol_zone1[z], vol_czone1[z] );
                if ( sico.MODEL == 7 ) fprintf(f_volumes, "\t%.3f\t%.3f\t%.3f\t%.3f", vol_zone2[z], vol_czone2[z], vol_zone3[z], vol_czone3[z] );
            }
    
            fprintf(f_volumes, "\n");


// -- STOP --- Display and files of status of simulation --------------------------------------------------------


// -- START -- Preparing and writing output raster maps and velocity fields -------------------------------------


            #ifdef WITHGRASS


                for ( i=0; i<sico.IMAX; i++ ) {

                    if ( sico.TOPOFACT == 0 ) topofact = 1.0;
                    else topofact = cos( betaxy[i] );

                    v[0] = awt[i][0]; // mixture or PHASE 1 flow height
                    if ( awt[i][0] > sico.IMPTHR[3] ) v[1] = fdiv( aw[i][2], aw[i][0] * topofact, sico.IMPTHR[3] ); else v[1] = 0;
                    if ( awt[i][0] > sico.IMPTHR[3] ) v[2] = fdiv( -aw[i][1], aw[i][0] * topofact, sico.IMPTHR[3] ); else v[2] = 0; // PHASE 1 x and y velocities

                    if ( sico.MODEL <= 3 ) { // one-phase models

                        if ( awt[i][0] > sico.IMPTHR[3] ) v[3] = aw[i][3] / awt[i][0]; else v[3] = 0; // temperature
                        v[4] = awt[i][4]; // change of basal surface
                        v[5] = awt[i][5]; // velocity
                        v[6] = awt[i][6]; // flow kinetic energy
                        v[7] = awt[i][7]; // flow pressure
                        
                    } else if ( sico.MODEL == 7 ) { // multi-phase model

                        v[3] = awt[i][3]; // PHASE 2 flow height
                        if ( awt[i][3] > sico.IMPTHR[3] ) v[4] = fdiv( aw[i][5], aw[i][3] * topofact, sico.IMPTHR[3] ); else v[4] = 0;
                        if ( awt[i][3] > sico.IMPTHR[3] ) v[5] = fdiv( -aw[i][4], aw[i][3] * topofact, sico.IMPTHR[3] ); else v[5] = 0; // PHASE 2 x and y velocities
                        v[6] = awt[i][6]; // PHASE 3 flow height
                        if ( awt[i][6] > sico.IMPTHR[3] ) v[7] = fdiv( aw[i][5], aw[i][6] * topofact, sico.IMPTHR[3] ); else v[7] = 0;
                        if ( awt[i][6] > sico.IMPTHR[3] ) v[8] = fdiv( -aw[i][4], aw[i][6] * topofact, sico.IMPTHR[3] ); else v[8] = 0; // PHASE 3 x and y velocities
                        if ( awt[i][0] + awt[i][3] + awt[i][6] > sico.IMPTHR[0] ) v[9] = aw[i][9] / ( awt[i][0] + awt[i][3] + awt[i][6] ); else v[9] = 0; // temperature
                        v[10] = awt[i][10];
                        v[11] = awt[i][11];
                        v[12] = awt[i][12]; // change of basal surface
                        v[13] = awt[i][13];
                        v[14] = awt[i][14];
                        v[15] = awt[i][15]; // flow velocities
                        v[16] = awt[i][16]; // total flow height
                        v[17] = awt[i][17];
                        v[18] = awt[i][18];
                        v[19] = awt[i][19];
                        v[20] = awt[i][20]; // flow kinetic energies
                        v[21] = awt[i][21];
                        v[22] = awt[i][22];
                        v[23] = awt[i][23];
                        v[24] = awt[i][24]; // flow pressures
                        v[25] = awt[i][25]; //hentri + hentri2 + hentri3; // total change of basal surface !!!CHECK
                    }

                    for ( k=0; k<nvect_red; k++ ) outv[px[i]][py[i]][k] = v[k];
                }


            #endif


            if ( sico.MULT == 0 ) { // for single model run

                if ( nout < 10 ) sprintf( madd, "000"); // fill string (for maintaining correct order in list of maps)
                else if ( nout < 100 ) sprintf( madd, "00");
                else if ( nout < 1000 ) sprintf( madd, "0");

                for ( k=0; k<nvect_red; k++ ) {
                
                    sprintf( mv, "%s%s%s%i", prefix, mv0[k], madd, nout ); // names of output raster maps

                    if ( sico.AFLAG == 1 || ( sico.MODEL <= 3 && ( k==0 || k==3 || k==4 )) || ( sico.MODEL == 7 
                        && ( k==0 || k==3 || k==6 || k==9 || k==10 || k==11 || k==12 || k == 16 || k == 25 ))) {


                        #ifdef WITHGRASS


                            foutrast ( mv, outv, sico, k, tsum ); // if GRASS is used, writing GRASS raster maps


                        #endif


                        foutasc( awt, px, py, outmaps, mv, k, sico ); // writing ascii raster maps
                    }
                }

                if ( sico.MODEL <= 3 ) { // ascii raster maps of maximum flow height at time step
                    sprintf( mv, "%s%s%s%i", prefix, mv0[8], madd, nout );
                    foutasc( awt, px, py, outmaps, mv, 8, sico );
                    
                } else if ( sico.MODEL == 7 ) {
                    sprintf( mv, "%s%s%s%i", prefix, mv0[32], madd, nout );
                    foutasc( awt, px, py, outmaps, mv, 32, sico );
                }

                if ( nout == 1 && sico.MODEL <= 3 ) foutdircoord ( f_directions, 1, px, py, sico ); // file for display of flow vectors as arrows
                else if ( nout == 1 && sico.MODEL == 7 ) {
                    foutdircoord ( f_directions, 4, px, py, sico );
                    foutdircoord ( f_directions2, 5, px, py, sico );
                    foutdircoord ( f_directions3, 6, px, py, sico );
                } // writing coordinates to file

                if ( sico.MODEL <= 3 ) { for ( k=0; k<3; k++ ) foutdir ( f_directions, awt, k, 1, px, py, sico );

                } else if ( sico.MODEL == 7 ) {
                    for ( k=0; k<3; k++ ) foutdir ( f_directions, awt, k, 4, px, py, sico );
                    for ( k=3; k<6; k++ ) foutdir ( f_directions2, awt, k, 5, px, py, sico );
                    for ( k=6; k<9; k++ ) foutdir ( f_directions3, awt, k, 6, px, py, sico );
                } // writing parameters to file

                if ( sico.PBG == 0 ) {

                    for ( i=0; i<sico.IMAX; i++ ) {

                        alpha = falpha( betax[i], betay[i], sico );
                        ppbg1[i] = (int)( 255.0 * ( cos( 45 * sico.PI / 180.0 ) * cos( betaxy[i] ) + sin( 45 * sico.PI / 180.0 ) 
                            * sin( betaxy[i] ) * cos( 135 * sico.PI / 180.0 - alpha )));
                    }

                    sprintf( mv, "%shillshade%s%i", prefix, madd, nout );
                    foutascind( ppbg1, px, py, outmaps, mv, sico );
                }
                
                if ( sico.TSUNAMI != 0 ) {

                    for ( i=0; i<sico.IMAX; i++ ) {

                        if ( sico.TOPOFACT == 0 ) topofact = 1.0;
                        else topofact = cos( betaxy[i] );

                        if ( sico.MODEL <= 3 ) {
                    
                            hflow = awt[i][0];
                            hrelease = phrelease[0] * topofact;
                    
                        } else if ( sico.MODEL == 7 ) {
                    
                            hflow = awt[i][16];
                            hrelease = ( phrelease[i] + phrelease2[i] + phrelease3[i] ) * topofact;
                        }

                        htsun[i] = hflow + pelev[i] - hrelease - pelev0[i];
                        if ( htsun[i] > htsunmax[i] ) htsunmax[i] = htsun[i];
                        if ( htsun[i] > htsunmaxmax ) htsunmaxmax = htsun[i];
                    }

                    sprintf( mv, "%shtsun%s%i", prefix, madd, nout );
                    foutascindf( htsun, px, py, outmaps, mv, sico );
                }
            }

            if ( (int)( round( 1000 * tsum )) >= (int)( round( 1000 * tmax )) || ccontinue == 0 ) { // raster maps for last time step

                for ( k=0; k<nvect_red; k++ ) {
                
                    if (( sico.MODEL <= 3 && ( k == 0 || k == 4 )) || ( sico.MODEL == 7 && ( k==0 || k==3 || k==6 || k==10 || k==11 || k==12 || k==16 || k==25 ))) {

                        if ( sico.MULT == 0 ) sprintf( mv, "%s%s_fin", prefix, mv0[k] );
                        else sprintf( mv, "%s%s_fin%d", prefix, mv0[k], xint ); // names of output raster maps


                        #ifdef WITHGRASS


                            foutrast ( mv, outv, sico, k, tsum ); // if GRASS is used, writing GRASS raster maps


                        #endif


                        if ( sico.MODEL <= 3 ) foutasc( awt, px, py, outmaps, mv, k, sico ); // writing ascii raster maps (one-phase models)
                        else if ( sico.MODEL == 7 ) foutasc( awt, px, py, outmaps, mv, k, sico ); // writing ascii raster maps (multi-phase model)
                    }
                }
            }


// -- STOP --- Preparing and writing output raster maps and velocity fields -------------------------------------


            cflmax = 0; // resetting maximum cfl value
            nout += 1; // updating number of output time steps
        }

        nsum += 1; // updating total number of time steps


// *** End of loop over time steps ------------------------------------------------------------------------------


    }


// -- START -- Writing profile ----------------------------------------------------------------------------------


    if ( sico.PROFILE > 0 ) {

        if ( sico.MULT == 0 ) sprintf(path, "%s%sprofile.txt", outfiles, prefix); // profile file
        else sprintf(path, "%s%sprofile%d.txt", outfiles, prefix, xint);
        f_profile=fopen(path, "w");

        for ( profk = -1; profk < profctrl; profk++ ) {

            if ( sico.HDEPOSIT != 0 ) {

                if ( profk == -1 ) fprintf( f_profile, "dist\telev\thdeposit" );        
                else fprintf( f_profile, "%.3f\t%.3f\t%.3f", profdiffabs[profk], profelev[profk], profdeposit[profk] );
                
            } else {

                if ( profk == -1 ) fprintf( f_profile, "dist\telev" );        
                else fprintf( f_profile, "%.3f\t%.3f", profdiffabs[profk], profelev[profk] );
            }
        
            for ( profi = 0; profi < nout; profi++ ) {

                if ( profk == -1 && sico.MODEL <= 3 ) fprintf( f_profile, "\thflow_%i\tvflow_%i\ttflow_%i\tpflow_%i\tbasechange_%i", profi, profi, profi, profi, profi );
                else if ( profk == -1 && sico.MODEL == 7 ) fprintf( f_profile, "\thflow1_%i\thflow2_%i\thflow3_%i\tvflow1_%i\tvflow2_%i\tvflow3_%i\ttflow1_%i\ttflow2_%i\ttflow3_%i\tpflow1_%i\tpflow2_%i\tpflow3_%i\tbasechange1_%i\tbasechange2_%i\tbasechange3_%i", profi, profi, profi, profi, profi, profi, profi, profi, profi, profi, profi, profi, profi, profi, profi ); 

                for ( profj = 0; profj < profmax; profj++ ) {

                    if ( profk != -1 ) fprintf( f_profile, "\t%.3f", profdata[profk][profi][profj] );
                } 
            }
            
            fprintf( f_profile, "\n" );
        }
        
        fclose(f_profile);
        
        if ( sico.MULT == 1 && xint == 1 ) {
        
            sprintf(path, "%s%sprofile.xyz", outaimec, prefix); // aimec profile file
            f_profile_aimec=fopen(path, "w");

            for ( profk = -1; profk < profctrl; profk++ ) {

                if ( profk == -1 ) fprintf( f_profile_aimec, "x, y, z" );        
                else fprintf( f_profile_aimec, "%.3f, %.3f, %.3f", ( float ) profny[profk] * sico.CSZ + sico.BDWEST, 
                    ( float ) sico.BDNORTH - profnx[profk] * sico.CSZ, pelev[profk] );
                fprintf( f_profile_aimec, "\n" );
            }
        
            fclose(f_profile_aimec);
        }
    }                     


// -- STOP --- Writing profile ----------------------------------------------------------------------------------


// -- START -- Writing control point data -----------------------------------------------------------------------


    if ( sico.CTRLPOINTS > 0 ) {

        fprintf( f_ctrlpoints, "id\tt\tx\ty\thflow1\thflow2\thflow3\thentr1\thentr2\thentr3\n" );

        for ( j=0; j<sico.CTRLPOINTS; j++ ) {

            for ( k = 0; k < nout; k++ ) {

                if ( sico.MODEL <= 3 )

                    fprintf( f_ctrlpoints, "C%i\t%.0f\t%.0f\t%.0f\t%.3f\t0.000\t0.000\t%.3f\t0.000\t0.000\n", j+1, 
                    k*tout, ctrlpxm[j], ctrlpym[j], ctrlpdata[j][k][0], ctrlpdata[j][k][1]);
                    
                else if ( sico.MODEL == 7 )

                    fprintf( f_ctrlpoints, "C%i\t%.0f\t%.0f\t%.0f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", j+1, 
                    k*tout, ctrlpxm[j], ctrlpym[j], ctrlpdata[j][k][0], ctrlpdata[j][k][1], ctrlpdata[j][k][2], ctrlpdata[j][k][3], ctrlpdata[j][k][4], ctrlpdata[j][k][5]);
            }
            
            if ( sico.MODEL <= 3 )

                fprintf( f_ctrlpoints, "C%i\tMAX\t%.0f\t%.0f\t%.3f\t0.000\t0.000\t%.3f\t0.000\t0.000\n", j+1, 
                ctrlpxm[j], ctrlpym[j], awt[icor[ctrlpx[j]][ctrlpy[j]]][8], awt[icor[ctrlpx[j]][ctrlpy[j]]][4]);
                    
            else if ( sico.MODEL == 7 )

                fprintf( f_ctrlpoints, "C%i\tMAX\t%.0f\t%.0f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", j+1, 
                ctrlpxm[j], ctrlpym[j], awt[icor[ctrlpx[j]][ctrlpy[j]]][26], awt[icor[ctrlpx[j]][ctrlpy[j]]][28], awt[icor[ctrlpx[j]][ctrlpy[j]]][30], 
                awt[icor[ctrlpx[j]][ctrlpy[j]]][10], awt[icor[ctrlpx[j]][ctrlpy[j]]][11], awt[icor[ctrlpx[j]][ctrlpy[j]]][12]);            
        }
    }


// -- STOP -- Writing control point data ------------------------------------------------------------------------


// -- START -- Writing file with key output ---------------------------------------------------------------------


    ctrl_basechange = 0;
    for ( i=0; i<sico.IMAX; i++ ) {

        if ( cdomain[i] != 0 && (( sico.MODEL <= 3 && awt[i][4] != 0 ) || ( sico.MODEL == 7 && awt[i][25] != 0 ))) ctrl_basechange = 1;  
    }

    if (( sico.MODEL <= 3 && sico.RELM == 1 ) || ( sico.MODEL == 7 && ( sico.RELM == 1 || sico.RELM2 == 1 || sico.RELM3 == 1 ))) ctrl_release = 1;
    else ctrl_release = 0;

    sprintf(path, "%s%snout%d.txt", outfiles, prefix, xint); // file of number of time steps, control for success, and volumes
    f_nout=fopen(path, "w");
    fprintf(f_nout, 
        "%i\n%i\n%i\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%i\n%i\n%i\n%i\n%i\n%.0f\n%0f\n%0f\n%.0f\n%0f\n%0f\n", 
        nout-1, csuccess, ctrl_basechange, hflow_maxmax, tflow_maxmax1, tflow_maxmax2, tflow_maxmax3, tflow_maxmax, pflow_maxmax1, pflow_maxmax2, pflow_maxmax3, pflow_maxmax, 
        vflow_maxmax1, vflow_maxmax2, vflow_maxmax3, vflow_maxmax, basechange_max, basechange_min, ctrl_release, sico.IMPACTAREA, sico.HDEPOSIT, hydnin, hydnout, 
        vol_entr, vol_entr2, vol_entr3, vol_edge, vol_edge2, vol_edge3 );
    fclose(f_nout);


// -- STOP --- Writing file with key output ---------------------------------------------------------------------


// -- START -- Evaluating simulation ----------------------------------------------------------------------------


    if ( sico.IMPACTAREA != 0 || sico.HDEPOSIT != 0 || sico.ENTRAINMENT != 0 ) {

        if ( sico.MULT == 0 ) { 
        
            sprintf(path, "%s%sevaluation.txt", outfiles, prefix); // evaluation file
            f_evaluation=fopen(path, "w");
                    
        } else {
        
            sprintf(path, "%s%sevaluation%i.txt", outfiles, prefix, xint);
            f_evaluation=fopen(path, "w");
            
            sprintf(path, "%s%saimec%i.txt", outaimec, prefix, xint); // aimec file
            f_aimec=fopen(path, "w");
        }

        if ( sico.MULT == 0 ) fprintf( f_evaluation, "Criterion\tValue\n" );

        if ( sico.ENTRAINMENT != 0 && sico.MODEL <= 3 && sico.MULT == 0 ) fprintf( f_evaluation, "\nEntrained volume VE\t%.0f cubic metres\n", vol_entr );
        
        else if ( sico.ENTRAINMENT != 0 && sico.MODEL == 7 && sico.MULT == 0 ) fprintf( f_evaluation, 
            "\nEntrained PHASE 1 volume VE1\t%.0f cubic metres\nEntrained PHASE 2 volume VE2\t%.0f cubic metres\nEntrained PHASE 3 volume VE3\t%.0f cubic metres\n", 
            vol_entr, vol_entr2, vol_entr3 );
            
        else if ( sico.MULT == 1 ) { 
        
            if ( sico.MODEL <= 3 ) {
            
                fprintf( f_evaluation, "%i\t%.2f\t%.2f", xint, vhrelease, vhentrmax );
                fprintf( f_aimec, "%i, %.2f, %.2f", xint, vhrelease, vhentrmax );
                
            } else {
            
                fprintf( f_evaluation, "%i\t%.2f\t%.2f\t%.2f\t%.2f", xint, vhrelease, rhrelease1, vhentrmax, rhentrmax1 );
                fprintf( f_aimec, "%i, %.2f, %.2f, %.2f, %.2f", xint, vhrelease, rhrelease1, vhentrmax, rhentrmax1 );
            }
            
            for ( l=0; l<lmax; l++ ) {
            
                fprintf( f_evaluation, "\t%.4f", flowpar[l] );
                fprintf( f_aimec, ", %.4f", flowpar[l] );
            }

            fprintf( f_evaluation, "\t%.0f\t%.0f\t%.0f", vol_entr, vol_entr2, vol_entr3 );
            //fprintf( f_aimec, ", %.0f, %.0f, %.0f", vol_entr, vol_entr2, vol_entr3 );
        }
    }

    if ( sico.IMPACTAREA != 0 || sico.HDEPOSIT != 0 ) {
    
        if ( sico.PROFILE > 0 ) {
           
            if ( sico.MODEL <= 3 ) profk = 8; else if ( sico.MODEL == 7 ) profk = 32;

            for ( profi = 1; profi < profctrl; profi++ ) {

                profwhtx1 = 1 - ( profnx[profi] - (int)profnx[profi] );
                profwhtx2 = 1 - ((int)profnx[profi] + 1 - profnx[profi] );
                profwhty1 = 1 - ( profny[profi] - (int)profny[profi] );
                profwhty2 = 1 - ((int)profny[profi] + 1 - profny[profi] );
              
                profeval[profi][0] = awt[icor[(int)profnx[profi]][(int)profny[profi]]][profk] * profwhtx1 * profwhty1
                    + awt[icor[(int)profnx[profi]][(int)profny[profi]+1]][profk] * profwhtx1 * profwhty2
                    + awt[icor[(int)profnx[profi]+1][(int)profny[profi]]][profk] * profwhtx2 * profwhty1
                    + awt[icor[(int)profnx[profi]+1][(int)profny[profi]+1]][profk] * profwhtx2 * profwhty2;

                if ( sico.MODEL <= 3 ) profhrelease = profdata[profi][0][0];
                else if ( sico.MODEL == 7 ) profhrelease = profdata[profi][0][0] + profdata[profi][0][1] + profdata[profi][0][2];
             
                if ( profhrelease >= sico.IMPTHR[0] || ( profeval[profi][0] >= sico.IMPTHR[0] && proflength_hmax_ctrl > 0 )) 
                    proflength_hmax += ( profdiffabs[profi] - profdiffabs[profi-1] );
                else if ( proflength_hmax_ctrl == 0 ) proflength_hmax = profdiffabs[profi] - profdiffabs[profi+1];
                    //!!!CHECK may not give correct results if runup at flow terminus
                if ( profhrelease >= sico.IMPTHR[0] ) proflength_hmax_ctrl = 1;

                if ( profeval[profi][0] >= sico.IMPTHR[0] ) proflength_simulated = profdiffabs[profi];

                if ( sico.IMPACTAREA == 1 ) {
                
                    profeval[profi][1] = pimpactarea[icor[(int)profnx[profi]][(int)profny[profi]]] * profwhtx1 * profwhty1
                        + pimpactarea[icor[(int)profnx[profi]][(int)profny[profi]+1]] * profwhtx1 * profwhty2
                        + pimpactarea[icor[(int)profnx[profi]+1][(int)profny[profi]]] * profwhtx2 * profwhty1
                        + pimpactarea[icor[(int)profnx[profi]+1][(int)profny[profi]+1]] * profwhtx2 * profwhty2;
                    
                    if ( profeval[profi][1] > 0 ) proflength_impactarea = profdiffabs[profi];
                }
                
                if ( sico.HDEPOSIT == 1 ) { 
                
                    if ( profdeposit[profi] >= sico.IMPTHR[0] ) proflength_hdeposit = profdiffabs[profi];
                }
            }

            if ( sico.IMPACTAREA == 1 ) {
            
                profdiff_impactarea = proflength_simulated - proflength_impactarea;
                proflength_impactarea -= ( proflength_simulated - proflength_hmax );
                profratio_impactarea = proflength_hmax / proflength_impactarea - 1;
                
            } else {

                profdiff_impactarea = sico.UNDEF;
                profratio_impactarea = sico.UNDEF;
            }

            if ( sico.HDEPOSIT == 1 ) {
            
                profdiff_hdeposit = proflength_simulated - proflength_hdeposit;
                proflength_hdeposit -= ( proflength_simulated - proflength_hmax );
                profratio_hdeposit = proflength_hmax / proflength_hdeposit - 1;

            } else {

                profdiff_hdeposit = sico.UNDEF;
                profratio_hdeposit = sico.UNDEF;
            }

            if ( sico.MULT == 0 ) {
            
                fprintf( f_evaluation, "\nFlow length to terminus of observed impact area Lio\t%.0f\nFlow length to terminus of observed deposit Ldo\t%.0f\nFlow length to terminus of simulated impact area Lim\t%.0f\nDifference between simulated and observed flow lengths (impact area) dLi\t%.0f\nRatio between simulated and observed flow lengths (impact area) rLi\t%.2f\nFlow length to terminus of simulated deposit Ldm\t%.0f\nDifference between simulated and observed flow lengths (deposit) dLd\t%.0f\nRatio between simulated and observed flow lengths (deposit) rLd\t%.2f\n", proflength_impactarea, proflength_hdeposit, proflength_hmax, profdiff_impactarea, profratio_impactarea, proflength_hmax, profdiff_hdeposit, profratio_hdeposit );

            } else {

                fprintf( f_evaluation, "\t%.0f\t%.0f\t%.0f\t%.0f\t%.2f\t%.0f\t%.0f\t%.2f", 
                    proflength_impactarea, proflength_hdeposit, proflength_hmax, profdiff_impactarea, profratio_impactarea, proflength_hmax, profdiff_hdeposit, 
                    profratio_hdeposit );
            }
       
        } else {
        
            profratio_impactarea = -1;
            profratio_hdeposit = -1;
            
            if ( sico.MULT == 1 ) {
            
                fprintf( f_evaluation, "\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f", 
                    sico.UNDEF, sico.UNDEF, sico.UNDEF, sico.UNDEF, sico.UNDEF, sico.UNDEF, sico.UNDEF, sico.UNDEF );
            }
        }            

        if ( sico.MODEL <= 3 ) k = 8; else if ( sico.MODEL == 7 ) k = 32;
        if ( sico.MODEL <= 3 && ctrl_basechange == 1 ) l = 4;
        else if ( sico.MODEL <= 3 && ctrl_basechange == 0 ) l = 0;
        else if ( sico.MODEL == 7 && ctrl_basechange == 1 ) l = 25;
        else l = 16;

        for ( i=0; i<sico.IMAX; i++ ) {

            if ( cdomain[i] != 0 && sico.IMPACTAREA == 1 ) {
            
                if ( pimpactarea[i] > 0 && awt[i][k] >= sico.IMPTHR[0] ) evaltp_imp += 1;
                else if ( pimpactarea[i] > 0 && awt[i][k] < sico.IMPTHR[0] ) evalfn_imp += 1;
                else if ( pimpactarea[i] == 0 && awt[i][k] >= sico.IMPTHR[0] ) evalfp_imp += 1;
                else if ( pimpactarea[i] == 0 && awt[i][k] < 0 ) evaltn_imp += 1;                                            
            }
            
            if ( cdomain[i] != 0 && sico.HDEPOSIT == 1 ) {       

                if ( phdeposit[i] >= sico.IMPTHR[0] && awt[i][l] >= sico.IMPTHR[0] ) evaltp_dep += 1;
                else if ( phdeposit[i] >= sico.IMPTHR[0] && awt[i][l] < sico.IMPTHR[0] ) evalfn_dep += 1;
                else if ( phdeposit[i] < sico.IMPTHR[0] && awt[i][l] >= sico.IMPTHR[0] ) evalfp_dep += 1;
                else if ( phdeposit[i] < sico.IMPTHR[0] && awt[i][l] < sico.IMPTHR[0] ) evaltn_dep += 1;        
            }
        }

        if ( sico.IMPACTAREA == 1 ) {

            evaltn_imp = 5 * ( evaltp_imp + evalfn_imp ) - evalfp_imp;
            evalall = evaltp_imp + evaltn_imp + evalfp_imp + evalfn_imp;

            evalrp_imp = 100 * ( evaltp_imp + evalfn_imp ) / (float)evalall;
            evalrn_imp = 100 * ( evaltn_imp + evalfp_imp ) / (float)evalall;
            
            if ( evalrp_imp > 0 && evalrn_imp > 0 ) {
                       
                evalrtp_imp = 100 * evaltp_imp / (float)evalall;
                evalrtn_imp = 100 * evaltn_imp / (float)evalall;
                evalrfp_imp = 100 * evalfp_imp / (float)evalall;
                evalrfn_imp = 100 * evalfn_imp / (float)evalall;
                evalrtpp_imp = evalrtp_imp / evalrp_imp;
                evalrfpp_imp = evalrfp_imp / evalrn_imp;
                evalrtnn_imp = evalrtn_imp / evalrn_imp;
                evalrpp_imp = 100 * ( evaltp_imp + evalfp_imp ) / (float)evalall;
                evalrnp_imp = 100 * ( evaltn_imp + evalfn_imp ) / (float)evalall;            
                evalrt_imp = 100 * ( evaltp_imp + evaltn_imp ) / (float)evalall;
                evalrf_imp = 100 * ( evalfp_imp + evalfn_imp ) / (float)evalall;
            }

            if ( evalrtp_imp + evalrfp_imp + evalrfn_imp != 0 ) evalcsi_imp = evalrtp_imp / ( evalrtp_imp + evalrfp_imp + evalrfn_imp ); else evalcsi_imp = sico.UNDEF; 
                // critical success index

            if (( evalrtp_imp + evalrfn_imp ) * ( evalrfn_imp + evalrtn_imp ) + ( evalrfp_imp + evalrtp_imp ) * ( evalrfp_imp + evalrtn_imp ) != 0 )
                evalhss_imp = ( 2 * ( evalrtp_imp * evalrtn_imp ) - ( evalrfp_imp * evalrfn_imp )) 
                    / (( evalrtp_imp + evalrfn_imp ) * ( evalrfn_imp + evalrtn_imp ) + ( evalrfp_imp + evalrtp_imp ) * ( evalrfp_imp + evalrtn_imp ));
            else evalhss_imp = sico.UNDEF; // Heidke skill score
            
            evalauroc_imp = ( 1 - evalrfpp_imp ) * evalrtpp_imp + evalrfpp_imp * evalrtpp_imp * 0.5 + ( 1 - evalrfpp_imp ) * ( 1 - evalrtpp_imp ) * 0.5; 
                // area under the ROC curve
            evald2pc_imp = pow( pow( 1 - evalrtpp_imp, 2 ) + pow( evalrfpp_imp, 2 ), 0.5 ); // distance to perfect classification          
            if ( evalrp_imp > 0 ) evalfoc_imp = (float)( evaltp_imp + evalfp_imp ) / (float)( evaltp_imp + evalfn_imp ); else evalfoc_imp = 0; // factor of conservativeness

            if ( profratio_impactarea == -1 ) profratio_impactarea = -1.1;
            if ( evalfoc_imp == 0 ) evalfoc_imp = -1;
            
            evalspi_imp = ( fmax(0, fmin( 0.2, 0.2 * evalcsi_imp )) + fmax( 0, fmin( 0.2, 0.2 * evalhss_imp )) + fmax(0, fmin( 0.2, 0.2 * ( 1 - evald2pc_imp ))) 
                + fmax(0, fmin( 0.2, 0.2 * pow( sin( sico.PI * 0.5 * fmin( profratio_impactarea + 1, 1 / ( profratio_impactarea + 1 ))), 5 ))) 
                + fmax(0, fmin( 0.2, 0.2 * pow( sin( sico.PI * 0.5 * fmin( evalfoc_imp, 1 / evalfoc_imp )), 5 )))); // synthetic performance index

            if ( profratio_impactarea == -1.1 ) profratio_impactarea = sico.UNDEF;
            if ( evalfoc_imp == -1 ) evalfoc_imp = sico.UNDEF;

            if ( sico.MULT == 0 ) {

                fprintf( f_evaluation, "\nEvaluation scores for observed impact area\nRate of positive observations rOP\t%.1f%%\nRate of negative observations rON\t%.1f%%\nTrue positive rate rTP\t%.1f%%\nTrue negative rate rTN\t%.1f%%\nFalse positive rate rFP\t%.1f%%\nFalse negative rate rFN\t%.1f%%\nTrue positive rate out of all positive observations rTP/rOP\t%.1f%%\nTrue negative rate out of all negative observations rTN/rON\t%.1f%%\nRate of total positive predictions rMP\t%.1f%%\nRate of total negative predictions rMN\t%.1f%%\nRate of true predictions rT\t%.1f%%\nRate of false predictions rF\t%.1f%%\nCritical success index CSI\t%.3f\nHeidke Skill Score HSS\t%.3f\nArea under the ROC curve AUROC\t%.3f\nDistance to perfect classification D2PC\t%.3f\nFactor of conservativeness FoC\t%.3f\nSynthetic performance index SPI\t%.3f\n", evalrp_imp, evalrn_imp, evalrtp_imp, evalrtn_imp, evalrfp_imp, evalrfn_imp, 100*evalrtpp_imp, 100*evalrtnn_imp, evalrpp_imp, evalrnp_imp, evalrt_imp, evalrf_imp, evalcsi_imp, evalhss_imp, evalauroc_imp, evald2pc_imp, evalfoc_imp, evalspi_imp );

            } else {

                fprintf( f_evaluation, "\t%.1f\t%.1f\t%.1f\t%.1f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f", 
                    evalrp_imp, evalrn_imp, evalrtp_imp, evalrtn_imp, evalcsi_imp, evalhss_imp, evalauroc_imp, evald2pc_imp, evalfoc_imp, evalspi_imp );
            }
            
        }
        
        else if ( sico.MULT == 1 ) {
        
            fprintf( f_evaluation, "\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f", 
                sico.UNDEF, sico.UNDEF, sico.UNDEF, sico.UNDEF, sico.UNDEF, sico.UNDEF, sico.UNDEF, sico.UNDEF, sico.UNDEF, sico.UNDEF );
        }

        if ( sico.HDEPOSIT == 1 ) {

            evaltn_dep = 5 * ( evaltp_dep + evalfn_dep ) - evalfp_dep;
            evalall = evaltp_dep + evaltn_dep + evalfp_dep + evalfn_dep;

            evalrp_dep = 100 * ( evaltp_dep + evalfn_dep ) / (float)evalall;
            evalrn_dep = 100 * ( evaltn_dep + evalfp_dep ) / (float)evalall;
            
            if ( evalrp_dep > 0 && evalrn_dep > 0 ) {
                      
                evalrtp_dep = 100 * evaltp_dep / (float)evalall;
                evalrtn_dep = 100 * evaltn_dep / (float)evalall;
                evalrfp_dep = 100 * evalfp_dep / (float)evalall;
                evalrfn_dep = 100 * evalfn_dep / (float)evalall;
                evalrtpp_dep = evalrtp_dep / evalrp_dep;
                evalrfpp_dep = evalrfp_dep / evalrn_dep;
                evalrtnn_dep = evalrtn_dep / evalrn_dep;
                evalrpp_dep = 100 * ( evaltp_dep + evalfp_dep ) / (float)evalall;
                evalrnp_dep = 100 * ( evaltn_dep + evalfn_dep ) / (float)evalall;            
                evalrt_dep = 100 * ( evaltp_dep + evaltn_dep ) / (float)evalall;
                evalrf_dep = 100 * ( evalfp_dep + evalfn_dep ) / (float)evalall;
            }

            if ( evalrtp_dep + evalrfp_dep + evalrfn_dep != 0 ) evalcsi_dep = evalrtp_dep / ( evalrtp_dep + evalrfp_dep + evalrfn_dep ); else evalcsi_dep = sico.UNDEF; 
                // critical success index

            if (( evalrtp_dep + evalrfn_dep ) * ( evalrfn_dep + evalrtn_dep ) + ( evalrfp_dep + evalrtp_dep ) * ( evalrfp_dep + evalrtn_dep ) != 0 )
                evalhss_dep = ( 2 * ( evalrtp_dep * evalrtn_dep ) - ( evalrfp_dep * evalrfn_dep )) 
                    / (( evalrtp_dep + evalrfn_dep ) * ( evalrfn_dep + evalrtn_dep ) + ( evalrfp_dep + evalrtp_dep ) * ( evalrfp_dep + evalrtn_dep ));
            else evalhss_dep = sico.UNDEF; // Heidke skill score
            
            evalauroc_dep = ( 1 - evalrfpp_dep ) * evalrtpp_dep + evalrfpp_dep * evalrtpp_dep * 0.5 + ( 1 - evalrfpp_dep ) * ( 1 - evalrtpp_dep ) * 0.5; 
                // area under the ROC curve
            evald2pc_dep = pow( pow( 1 - evalrtpp_dep, 2 ) + pow( evalrfpp_dep, 2 ), 0.5 ); // distance to perfect classification          
            if ( evalrp_dep > 0 ) evalfoc_dep = (float)( evaltp_dep + evalfp_dep ) / (float)( evaltp_dep + evalfn_dep ); else evalfoc_dep = 0; // factor of conservativeness
            
            if ( profratio_hdeposit == -1 ) profratio_hdeposit = -1.1;
            if ( evalfoc_dep == 0 ) evalfoc_dep = -1;
            
            evalspi_dep = ( fmax(0, fmin( 0.2, 0.2 * evalcsi_dep )) + fmax( 0, fmin( 0.2, 0.2 * evalhss_dep )) + fmax(0, fmin( 0.2, 0.2 * ( 1 - evald2pc_dep ))) 
                + fmax(0, fmin( 0.2, 0.2 * pow( sin( sico.PI * 0.5 * fmin( profratio_hdeposit + 1, 1 / ( profratio_hdeposit + 1 ))), 5 ))) 
                + fmax(0, fmin( 0.2, 0.2 * pow( sin( sico.PI * 0.5 * fmin( evalfoc_dep, 1 / evalfoc_dep )), 5 )))); // synthetic performance index

            if ( profratio_hdeposit == -1.1 ) profratio_hdeposit = sico.UNDEF;
            if ( evalfoc_dep == -1 ) evalfoc_dep = sico.UNDEF;

            if ( sico.MULT == 0 ) {

                fprintf( f_evaluation, "\nEvaluation scores for observed deposit\nRate of positive observations rOP\t%.1f%%\nRate of negative observations rON\t%.1f%%\nTrue positive rate rTP\t%.1f%%\nTrue negative rate rTN\t%.1f%%\nFalse positive rate rFP\t%.1f%%\nFalse negative rate rFN\t%.1f%%\nTrue positive rate out of all positive observations rTP/rOP\t%.1f%%\nTrue negative rate out of all negative observations rTN/rON\t%.1f%%\nRate of total positive predictions rMP\t%.1f%%\nRate of total negative predictions rMN\t%.1f%%\nRate of true predictions rT\t%.1f%%\nRate of false predictions rF\t%.1f%%\nCritical success index CSI\t%.3f\nHeidke Skill Score HSS\t%.3f\nArea under the ROC curve AUROC\t%.3f\nDistance to perfect classification D2PC\t%.3f\nFactor of conservativeness FoC\t%.3f\nSynthetic performance index SPI\t%.3f\n", evalrp_dep, evalrn_dep, evalrtp_dep, evalrtn_dep, evalrfp_dep, evalrfn_dep, 100*evalrtpp_dep, 100*evalrtnn_dep, evalrpp_dep, evalrnp_dep, evalrt_dep, evalrf_dep, evalcsi_dep, evalhss_dep, evalauroc_dep, evald2pc_dep, evalfoc_dep, evalspi_dep ); 
                
            } else {

                fprintf( f_evaluation, "\t%.1f\t%.1f\t%.1f\t%.1f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f", 
                    evalrp_dep, evalrn_dep, evalrtp_dep, evalrtn_dep, evalcsi_dep, evalhss_dep, evalauroc_dep, evald2pc_dep, evalfoc_dep, evalspi_dep );
            }
        }

        else if ( sico.MULT == 1 ) {
        
            fprintf( f_evaluation, "\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f", 
                sico.UNDEF, sico.UNDEF, sico.UNDEF, sico.UNDEF, sico.UNDEF, sico.UNDEF, sico.UNDEF, sico.UNDEF, sico.UNDEF, sico.UNDEF );
        }
    }
    
    if ( sico.IMPACTAREA != 0 || sico.HDEPOSIT != 0 || sico.ENTRAINMENT != 0 ) { 
    
        fclose(f_evaluation);
        if ( sico.MULT == 1 ) fclose(f_aimec);
    }


// -- STOP --- Evaluating simulation ----------------------------------------------------------------------------


// -- START -- Preparing and writing output raster maps of cumulative and maximum values ------------------------


    for ( k=nvect_red; k<nvect_all; k++ ) {

        if ( sico.AFLAG == 1 || ( sico.MODEL <= 3 && ( k==8 || k == 12 || k==17 || k==18 )) 
            || ( sico.MODEL == 7 && ( k==26 || k==28 || k==30 || k==32 || k==41 || k==42 || k==43 || k==44 || k==45 || k==46 || k==47 ))) {

            if ( sico.MULT == 0 ) sprintf( mv, "%s%s", prefix, mv0[k] );
            else sprintf( mv, "%s%s%d", prefix, mv0[k], xint ); // names of output raster maps


            #ifdef WITHGRASS


                for ( i=0; i<sico.IMAX; i++ ) {

                    pout = awt[i][k]; // for edge cells, applying depths as heights
                    outv[px[i]][py[i]][k] = pout;
                }

                foutrast ( mv, outv, sico, k, tsum ); // writing GRASS raster maps


            #endif


            foutasc( awt, px, py, outmaps, mv, k, sico ); // writing ascii raster maps
        }
    }
    
    if ( sico.MULT == 0 ) sprintf( mv, "%shflow_dep", prefix );
    else sprintf( mv, "%shflow_dep%d", prefix, xint ); // names of output raster maps

    if ( sico.MODEL <= 3 && ctrl_basechange == 1 ) k = 4;
    else if ( sico.MODEL <= 3 && ctrl_basechange == 0 ) k = 0;
    else if ( sico.MODEL == 7 && ctrl_basechange == 1 ) k = 25;
    else k = 16;


    #ifdef WITHGRASS


        for ( i=0; i<sico.IMAX; i++ ) outv[px[i]][py[i]][k] = fmax( awt[i][k], 0 );
        foutrast ( mv, outv, sico, k, tsum ); // writing GRASS raster maps
        
        
    #endif
        
    for ( i=0; i<sico.IMAX; i++ ) awt[i][0] = fmax( awt[i][k], 0 ); // !!!CHECK     
    foutasc( awt, px, py, outmaps, mv, 0, sico ); // writing ascii raster maps    

    if ( sico.MULT == 0 ) sprintf( mv, "%selev_mod", prefix );
    else sprintf( mv, "%selev_mod%d", prefix, xint ); // names of output raster maps   

    if ( sico.MODEL <= 3  ) k = 4;
    else if ( sico.MODEL == 7 ) k = 25;


    #ifdef WITHGRASS


        for ( i=0; i<sico.IMAX; i++ ) outv[px[i]][py[i]][k] = pelev[i];
        foutrast ( mv, outv, sico, k, tsum ); // writing GRASS raster maps

    #endif


    for ( i=0; i<sico.IMAX; i++ ) awt[i][k] = pelev[i]; // !!!CHECK provisional
    foutasc( awt, px, py, outmaps, mv, k, sico ); // writing ascii raster maps

    if ( sico.TSUNAMI != 0 ) {

        sprintf( mv, "%shtsun_max", prefix );
        foutascindf( htsunmax, px, py, outmaps, mv, sico );
    }

    
// -- STOP --- Preparing and writing output raster maps of cumulative and maximum values ------------------------


// -- START -- Writing pvpython script for import to Paraview ---------------------------------------------------


    if ( sico.MULT == 0 ) {

        sprintf(path, "%spvimport.py", outvr); // pvpython script
        f_paraviewi=fopen(path, "w");

        fprintf(f_paraviewi, "import glob\n");
        fprintf(f_paraviewi, "import os\n");
        fprintf(f_paraviewi, "from paraview.simple import *\n");
        fprintf(f_paraviewi, "import shutil\n");
        fprintf(f_paraviewi, "\n");
    
        if ( sico.TSUNAMI == 1 ) {
    
            fprintf(f_paraviewi, "hcont = [%i]\n", paracontourshmin);
            fprintf(f_paraviewi, "hmin = %i\n", paracontourshint);
        
        } else {
    
            fprintf(f_paraviewi, "hcont = []\n");
            fprintf(f_paraviewi, "hmin = %i\n", paracontourshmin);
        }
    
        fprintf(f_paraviewi, "hmax = %i\n", paracontourshmax);
        fprintf(f_paraviewi, "hint = %i\n", paracontourshint);
        fprintf(f_paraviewi, "for i in range(hmin, hmax, hint): hcont.append(i)\n");
        fprintf(f_paraviewi, "\n");
        fprintf(f_paraviewi, "elevcont = []\n");
        fprintf(f_paraviewi, "elevmin = %i\n", paracontourszmin);
        fprintf(f_paraviewi, "elevmax = %i\n", paracontourszmax);
        fprintf(f_paraviewi, "elevint = %i\n", paracontourszint);
        fprintf(f_paraviewi, "for i in range(elevmin, elevmax, elevint): elevcont.append(i)\n");
        fprintf(f_paraviewi, "\n");
        fprintf(f_paraviewi, "if os.path.exists('surface'): shutil.rmtree('surface')\n");
        fprintf(f_paraviewi, "if os.path.exists('contoursh'): shutil.rmtree('contoursh')\n");
        fprintf(f_paraviewi, "if os.path.exists('contoursz'): shutil.rmtree('contoursz')\n");
        fprintf(f_paraviewi, "\n");
        fprintf(f_paraviewi, "intab = paraview.simple.CSVReader(FileName=glob.glob('data/pv*.csv'))\n");
        fprintf(f_paraviewi, "points = paraview.simple.TableToPoints(intab, XColumn='x', YColumn='y', ZColumn='z', KeepAllDataArrays=True)\n");
        fprintf(f_paraviewi, "surface = paraview.simple.Delaunay2D(points)\n");
        fprintf(f_paraviewi, "surfcalc = paraview.simple.Calculator(Input=surface, ResultArrayName='rgb', Function='r * iHat + g * jHat + b * kHat')\n");
        fprintf(f_paraviewi, "contcalc = paraview.simple.Calculator(Input=surface, ResultArrayName='h', Function='z - zb')\n");
        fprintf(f_paraviewi, "contoursh = paraview.simple.Contour(Input=contcalc, ContourBy='h', Isosurfaces=hcont)\n");
        if ( sico.TSUNAMI == 1 && sico.LAYERS0 > 0 ) fprintf(f_paraviewi, "contoursz = paraview.simple.Contour(Input=surface, ContourBy='z1', Isosurfaces=elevcont)\n");
        else fprintf(f_paraviewi, "contoursz = paraview.simple.Contour(Input=surface, ContourBy='z', Isosurfaces=elevcont)\n");
        fprintf(f_paraviewi, "\n");
        fprintf(f_paraviewi, "print('Writing surfaces ...')\n");
        fprintf(f_paraviewi, "paraview.simple.SaveData('surface.pvd', surfcalc, WriteTimeSteps=True)\n");
        fprintf(f_paraviewi, "\n");
        fprintf(f_paraviewi, "print('Writing flow contours ...')\n");
        fprintf(f_paraviewi, "paraview.simple.SaveData('contoursh.pvd', contoursh, WriteTimeSteps=True)\n");
        fprintf(f_paraviewi, "\n");
        fprintf(f_paraviewi, "print('Writing surface contours ...')\n");
        fprintf(f_paraviewi, "paraview.simple.SaveData('contoursz.pvd', contoursz, WriteTimeSteps=True)\n");
        fprintf(f_paraviewi, "\n");
        fprintf(f_paraviewi, "print ('completed.')\n");
        
        fclose(f_paraviewi);

        if ( strcmp( parapy, "None" ) != 0 ) {
        
            sprintf(path, "%spvimport.cmd", outvr); // cmd script to start pvpython script
            f_paraviewc=fopen(path, "w");

            fprintf(f_paraviewc, "\"%s/pvpython.exe\" pvimport.py\n", parapy );       
        
            fclose(f_paraviewc);
        }
    }


// -- STOP --- Writing pvpython script for import to Paraview ---------------------------------------------------


// -- START -- Writing python script for import to Blender ------------------------------------------------------


    if ( sico.MULT == 0 ) {

        sprintf(path, "%sblimport.py", outvr); // python script
        f_blenderi=fopen(path, "w");

        fprintf(f_blenderi, "import bpy\n");
        fprintf(f_blenderi, "import csv\n");
        fprintf(f_blenderi, "import math\n");
        fprintf(f_blenderi, "import os\n");
        fprintf(f_blenderi, "import shutil\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "# **********************************************************************\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "# These parameters are automatically set by r.avaflow\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "pf = '%s' # prefix for objects\n", prefix);
        fprintf(f_blenderi, "xdim = %i\n", sico.N);
        fprintf(f_blenderi, "ydim = %i\n", sico.M);
        fprintf(f_blenderi, "xmin = %.3f\n", sico.BDWEST);
        fprintf(f_blenderi, "xmax = %.3f\n", sico.BDEAST);
        fprintf(f_blenderi, "ymin = %.3f\n", sico.BDSOUTH);
        fprintf(f_blenderi, "ymax = %.3f\n", sico.BDNORTH);
        fprintf(f_blenderi, "nfiles = %i\n", nout-1);
        fprintf(f_blenderi, "undef = %.3f\n", sico.UNDEF);
        fprintf(f_blenderi, "layers = %i\n", sico.LAYERS0);
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "# **********************************************************************\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "# These parameters have to be adapted manually, if required\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "wkdir = '' # location of data directory in file system (full path with normal slashes, trailing slash required)\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "impfiles = 1 # control for csv file import (1 = activated, 0 = deactivated)\n");
        fprintf(f_blenderi, "standard = 0 # control for standard video (1 = activated, 0 = deactivated)\n");
        fprintf(f_blenderi, "anaglyph = 0 # control for anaglyph video (1 = activated, 0 = deactivated)\n");
        fprintf(f_blenderi, "stereo = 0 # control for stereo 3D video (1 = activated, 0 = deactivated)\n");
        fprintf(f_blenderi, "animation = 0 # control for animation\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "vbackground = -9999.0 # background surface elevation value (excluded from meshes)\n");
        fprintf(f_blenderi, "include_bg = 2 # control for including background (0 = no, 1 = only for impact area, 2 = yes)\n");
        fprintf(f_blenderi, "subsurface = 0 # control for including phase 1 surface (only useful for layered multi-phase model: 0 = no, 1 = yes)\n");
        fprintf(f_blenderi, "successive = 0 # control for successive simulation (0 = no, 1 = yes)\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "light = bpy.data.objects['Light'] # name of light object\n");
        fprintf(f_blenderi, "light.location = (3, 4, 5) # light location\n");
        fprintf(f_blenderi, "light.data.energy=200.0 # light energy\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "cam = bpy.data.objects['Camera'] # name of camera object\n");
        fprintf(f_blenderi, "camx = 0.0 # x position of camera\n");
        fprintf(f_blenderi, "camy = 0.0 # y position of camera\n");
        fprintf(f_blenderi, "camz = 8.0 # z position of camera\n");
        fprintf(f_blenderi, "cam.location = (camx, camy, camz) # camera location\n");
        fprintf(f_blenderi, "cam.rotation_euler = (0.0 * math.pi/180, 0.0 * math.pi/180, 0.0 * math.pi/180) # camera rotation\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "camspeedx = 0.0 # camera speed in x direction during video (km per frame)\n");
        fprintf(f_blenderi, "camspeedy = 0.0 # camera speed in y direction during video (km per frame)\n");
        fprintf(f_blenderi, "camspeedz = 0.0 # camera speed in z direction during video (km per frame)\n");
        fprintf(f_blenderi, "convdist = 3 # convergence distance of stereo camera\n");
        fprintf(f_blenderi, "resx = 11 # control for x dimension of stereo 3D video (bits)\n");
        fprintf(f_blenderi, "interval = 3 # interval between animation time steps (frames)\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "# **********************************************************************\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "# The following code should only be modified by experienced users\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "prefix = 'pv'\n");
        fprintf(f_blenderi, "zref = 0.0\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "scn = bpy.context.scene\n");
        fprintf(f_blenderi, "scn.camera = cam\n");
        fprintf(f_blenderi, "scn.render.film_transparent=False\n");
        fprintf(f_blenderi, "scn.render.image_settings.file_format='PNG'\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "collection = bpy.context.collection\n");
        fprintf(f_blenderi, "meshes = set()\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "bpy.types.View3DShading.color = 'ATTRIBUTE'\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "if standard == 1 and os.path.exists(wkdir + pf + 'standard') and os.path.isdir(wkdir + pf + 'standard'):\n");
        fprintf(f_blenderi, "    shutil.rmtree(wkdir + pf + 'standard')\n");
        fprintf(f_blenderi, "    os.makedirs(wkdir + pf + 'standard')\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "if anaglyph == 1 and os.path.exists(wkdir + pf + 'anaglyph') and os.path.isdir(wkdir + pf + 'anaglyph'):\n");
        fprintf(f_blenderi, "    shutil.rmtree(wkdir + pf + 'anaglyph')\n");
        fprintf(f_blenderi, "    os.makedirs(wkdir + pf + 'anaglyph')\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "if stereo == 1 and os.path.exists(wkdir + pf + 'stereo') and os.path.isdir(wkdir + pf + 'stereo'):\n");
        fprintf(f_blenderi, "    shutil.rmtree(wkdir + pf + 'stereo')\n");
        fprintf(f_blenderi, "    os.makedirs(wkdir + pf + 'stereo')\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "ob = []\n");
        fprintf(f_blenderi, "if layers > 0 and subsurface == 1: ob1 = []\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "if impfiles == 1:\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "    if successive == 0:\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        for obj in [o for o in collection.objects if o.type == 'MESH']:\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "            meshes.add( obj.data )\n");
        fprintf(f_blenderi, "            bpy.data.objects.remove( obj )\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        for mesh in [m for m in meshes if m.users == 0]:\n");
        fprintf(f_blenderi, "            bpy.data.meshes.remove( mesh )\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        for material in bpy.data.materials: bpy.data.materials.remove(material)\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        newmat = bpy.data.materials.new('VertCol')\n");
        fprintf(f_blenderi, "        newmat.use_nodes = True\n");
        fprintf(f_blenderi, "        \n");
        fprintf(f_blenderi, "    else:\n");
        fprintf(f_blenderi, "        \n");
        fprintf(f_blenderi, "        newmat = bpy.data.materials['VertCol']\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "    if include_bg == 1:\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        if nfiles < 10: filetext = '000' + str(nfiles)\n");
        fprintf(f_blenderi, "        elif nfiles < 100: filetext = '00' + str(nfiles)\n");
        fprintf(f_blenderi, "        elif nfiles < 1000: filetext = '0' + str(nfiles)\n");
        fprintf(f_blenderi, "        else: filetext = str(file)\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        path = wkdir + 'data/' + prefix + filetext + '.csv'\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        imp = []\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        with open(path, 'r') as csvfile:\n");
        fprintf(f_blenderi, "            datareader = csv.reader(csvfile)\n");
        fprintf(f_blenderi, "            next(datareader)\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "            i = 0\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "            for row in datareader:\n");
        fprintf(f_blenderi, "        \n");
        fprintf(f_blenderi, "                imp.append(int(row[8]))\n");
        fprintf(f_blenderi, "           \n");
        fprintf(f_blenderi, "                i += 1\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "    for file in range(0, nfiles + 1):\n");
        fprintf(f_blenderi, "    \n");
        fprintf(f_blenderi, "        if file < 10: filetext = '000' + str(file)\n");
        fprintf(f_blenderi, "        elif file < 100: filetext = '00' + str(file)\n");
        fprintf(f_blenderi, "        elif file < 1000: filetext = '0' + str(file)\n");
        fprintf(f_blenderi, "        else: filetext = str(file)\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        path = wkdir + 'data/' + prefix + filetext + '.csv'\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        vertices = []\n");
        fprintf(f_blenderi, "        edges = []\n");
        fprintf(f_blenderi, "        faces = []\n");
        fprintf(f_blenderi, "        colors = []\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        if subsurface == 1 and file == 0: verticesb = []\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        if layers > 0 and subsurface == 1:\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "            vertices1 = []\n");
        fprintf(f_blenderi, "            edges1 = []\n");
        fprintf(f_blenderi, "            faces1 = []\n");
        fprintf(f_blenderi, "            colors1 = []\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        if not include_bg == 1: imp = []\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        xcotr = 0\n");
        fprintf(f_blenderi, "        ycotr = 0\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        with open(path, 'r') as csvfile:\n");
        fprintf(f_blenderi, "            datareader = csv.reader(csvfile)\n");
        fprintf(f_blenderi, "            next(datareader)\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "            i = 0\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "            for row in datareader:\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "                x = float(row[0]) * 0.001 - ( xmax + xmin ) * 0.5 * 0.001\n");
        fprintf(f_blenderi, "                y = float(row[1]) * 0.001 - ( ymax + ymin ) * 0.5 * 0.001\n");
        fprintf(f_blenderi, "                z0 = float(row[2])\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "                if z0 == vbackground: z = undef\n");
        fprintf(f_blenderi, "                else: z = z0 * 0.001 - zref * 0.001\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "                if subsurface == 1 and file == 0:\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "                    zb0 = float(row[3])\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "                    if zb0 == vbackground: zb = undef\n");
        fprintf(f_blenderi, "                    else: zb = zb0 * 0.001 - zref * 0.001\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "                if layers > 0 and subsurface == 1:\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "                    z10 = float(row[4])\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "                    if z10 == vbackground: z1 = undef\n");
        fprintf(f_blenderi, "                    else: z1 = z10 * 0.001 - zref * 0.001\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "                r = float(row[5])\n");
        fprintf(f_blenderi, "                g = float(row[6])\n");
        fprintf(f_blenderi, "                b = float(row[7])\n");
        fprintf(f_blenderi, "                if include_bg == 0: imp.append(int(row[8]))\n");
        fprintf(f_blenderi, "                elif include_bg == 2: imp.append(1)\n");
        fprintf(f_blenderi, "           \n");
        fprintf(f_blenderi, "                vertices.append((x, y, z))\n");
        fprintf(f_blenderi, "                if subsurface == 1 and file == 0: verticesb.append((x, y, zb))\n");
        fprintf(f_blenderi, "                if layers > 0 and subsurface == 1: vertices1.append((x, y, z1))\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "                if not z0 == undef and xcotr > 0 and ycotr > 0:\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "                    i1 = i-1\n");
        fprintf(f_blenderi, "                    i2 = i\n");
        fprintf(f_blenderi, "                    i3 = i-xdim\n");
        fprintf(f_blenderi, "                    i4 = i-xdim-1\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "                    if ( not vertices[i1][2] == undef and not vertices[i2][2] == undef \n");
        fprintf(f_blenderi, "                        and not vertices[i3][2] == undef and not vertices[i4][2] == undef and not imp[i] == 0 ):\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "                        faces.append((i1, i2, i3, i4))\n");
        fprintf(f_blenderi, "                        colors.append((r, g, b, 1))\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "                    if ( layers > 0 and subsurface == 1 and not vertices1[i1][2] == undef and not vertices1[i2][2] == undef \n");
        fprintf(f_blenderi, "                        and not vertices1[i3][2] == undef and not vertices1[i4][2] == undef and not imp[i] == 0 \n");
        fprintf(f_blenderi, "                        and ( vertices1[i1][2] < vertices[i1][2] or vertices1[i2][2] < vertices[i2][2] \n");
        fprintf(f_blenderi, "                        or vertices1[i3][2] < vertices[i3][2] or vertices1[i4][2] < vertices[i4][2] )):\n");
        fprintf(f_blenderi, "                \n");
        fprintf(f_blenderi, "                        faces1.append((i1, i2, i3, i4))\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "                if xcotr < xdim-1:\n");
        fprintf(f_blenderi, "                    xcotr += 1\n");
        fprintf(f_blenderi, "                else:\n");
        fprintf(f_blenderi, "                    xcotr = 0\n");
        fprintf(f_blenderi, "                    ycotr += 1\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "                i += 1\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        mesh = bpy.data.meshes.new(pf + 'Mesh' + str(filetext))\n");
        fprintf(f_blenderi, "        mesh.from_pydata(vertices, edges, faces)\n");
        fprintf(f_blenderi, "        mesh.update()\n");
        fprintf(f_blenderi, "        ob.append(bpy.data.objects.new(pf + 'Surface' + str(filetext), mesh))\n");
        fprintf(f_blenderi, "        ob[file].data.update()\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        if subsurface == 1 and file == 0:\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "            meshb = bpy.data.meshes.new(pf + 'Meshb')\n");
        fprintf(f_blenderi, "            meshb.from_pydata(verticesb, edges, faces)\n");
        fprintf(f_blenderi, "            meshb.update()\n");
        fprintf(f_blenderi, "            obb = bpy.data.objects.new(pf + 'Surfaceb', meshb)\n");
        fprintf(f_blenderi, "            obb.data.update()\n");
        fprintf(f_blenderi, "            bpy.context.collection.objects.link(obb)\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        if layers > 0 and subsurface == 1:\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "            mesh1 = bpy.data.meshes.new(pf + 'Mesh1' + str(filetext))\n");
        fprintf(f_blenderi, "            mesh1.from_pydata(vertices1, edges1, faces1)\n");
        fprintf(f_blenderi, "            mesh1.update()\n");
        fprintf(f_blenderi, "            ob1.append(bpy.data.objects.new(pf + 'Surface1' + str(filetext), mesh1))\n");
        fprintf(f_blenderi, "            ob1[file].data.update()\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        if not mesh.vertex_colors:\n");
        fprintf(f_blenderi, "            mesh.vertex_colors.new()\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        color_layer = mesh.vertex_colors.active\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        j=0\n");
        fprintf(f_blenderi, "        k=0\n");
        fprintf(f_blenderi, "        for face in mesh.polygons:\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "            for idx in face.loop_indices:\n");
        fprintf(f_blenderi, "                color_layer.data[k].color = colors[j]\n");
        fprintf(f_blenderi, "                k += 1\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "            j += 1\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        if mesh.materials:\n");
        fprintf(f_blenderi, "            mesh.materials[0] = newmat\n");
        fprintf(f_blenderi, "        else:\n");
        fprintf(f_blenderi, "            mesh.materials.append(newmat) \n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        nodes = newmat.node_tree.nodes\n");
        fprintf(f_blenderi, "        bsdf = nodes.get('Principled BSDF') \n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        vertex_color_node = None\n");
        fprintf(f_blenderi, "        if not 'Vertex Color' in [node.type for node in nodes]:\n");
        fprintf(f_blenderi, "            vertex_color_node = nodes.new(type = 'ShaderNodeVertexColor')\n");
        fprintf(f_blenderi, "        else:\n");
        fprintf(f_blenderi, "            vertex_color_node = nodes.get('Vertex Color')\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        vertex_color_node.layer_name = 'Col'\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        links = newmat.node_tree.links\n");
        fprintf(f_blenderi, "        link = links.new(vertex_color_node.outputs[0], bsdf.inputs[0])\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        bpy.context.collection.objects.link(ob[file])\n");
        fprintf(f_blenderi, "        if layers > 0 and subsurface == 1: bpy.context.collection.objects.link(ob1[file])\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "if standard == 1 or anaglyph == 1 or stereo == 1:\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "    if impfiles == 0:\n");
        fprintf(f_blenderi, "        \n");
        fprintf(f_blenderi, "        for file in range(0, nfiles + 1):\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "            if file < 10: filetext = '000' + str(file)\n");
        fprintf(f_blenderi, "            elif file < 100: filetext = '00' + str(file)\n");
        fprintf(f_blenderi, "            elif file < 1000: filetext = '0' + str(file)\n");
        fprintf(f_blenderi, "            else: filetext = str(file)\n");
        fprintf(f_blenderi, "            \n");
        fprintf(f_blenderi, "            ob.append(bpy.data.objects[pf + 'Surface' + str(filetext)])\n");
        fprintf(f_blenderi, "            ob[file].hide_render = True\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "    for file in range(0, nfiles + 1):\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        if file < 10: filetext = '000' + str(file)\n");
        fprintf(f_blenderi, "        elif file < 100: filetext = '00' + str(file)\n");
        fprintf(f_blenderi, "        elif file < 1000: filetext = '0' + str(file)\n");
        fprintf(f_blenderi, "        else: filetext = str(file)\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        cam.location = (camx + float(file) * camspeedx, \n");
        fprintf(f_blenderi, "            camy + float(file) * camspeedy, \n");
        fprintf(f_blenderi, "            camz + float(file) * camspeedz)\n");
        fprintf(f_blenderi, "    \n");
        fprintf(f_blenderi, "        ob[file].hide_render = False\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        if standard == 1:\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "            scn.render.use_multiview = False\n");
        fprintf(f_blenderi, "            scn.render.resolution_x = 1920\n");
        fprintf(f_blenderi, "            scn.render.resolution_y = 1080\n");
        fprintf(f_blenderi, "            \n");
        fprintf(f_blenderi, "            scn.render.filepath = wkdir + pf + 'standard/standard' + filetext + '.png'\n");
        fprintf(f_blenderi, "            bpy.ops.render.render(write_still=1)\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        if anaglyph == 1:\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "            scn.render.image_settings.views_format = 'STEREO_3D'\n");
        fprintf(f_blenderi, "            scn.render.image_settings.stereo_3d_format.display_mode = 'ANAGLYPH'\n");
        fprintf(f_blenderi, "            scn.render.use_multiview = True\n");
        fprintf(f_blenderi, "            scn.render.resolution_x = 1920\n");
        fprintf(f_blenderi, "            scn.render.resolution_y = 1080\n");
        fprintf(f_blenderi, "            cam.data.stereo.convergence_distance = convdist\n");
        fprintf(f_blenderi, "            cam.data.stereo.interocular_distance = convdist / 30.0\n");
        fprintf(f_blenderi, "    \n");
        fprintf(f_blenderi, "            scn.render.filepath = wkdir + pf + 'anaglyph/anaglyph' + filetext + '.png'\n");
        fprintf(f_blenderi, "            bpy.ops.render.render(write_still=1)\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        if stereo == 1:\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "            scn.render.image_settings.views_format = 'STEREO_3D'\n");
        fprintf(f_blenderi, "            scn.render.image_settings.stereo_3d_format.display_mode = 'SIDEBYSIDE'\n");
        fprintf(f_blenderi, "            scn.render.use_multiview = True\n");
        fprintf(f_blenderi, "            scn.render.resolution_x = int(pow(2, resx))\n");
        fprintf(f_blenderi, "            scn.render.resolution_y = int(pow(2, resx)) #int(scn.render.resolution_x / 2)\n");
        fprintf(f_blenderi, "            cam.data.type = 'PANO'\n");
        fprintf(f_blenderi, "            cam.data.stereo.use_spherical_stereo = True\n");
        fprintf(f_blenderi, "            cam.data.stereo.convergence_mode = 'OFFAXIS'\n");
        fprintf(f_blenderi, "            cam.data.stereo.convergence_distance = convdist\n");
        fprintf(f_blenderi, "            cam.data.stereo.interocular_distance = convdist / 30.0\n");
        fprintf(f_blenderi, "    \n");
        fprintf(f_blenderi, "            scn.render.filepath = wkdir + pf + 'stereo/stereo' + filetext + '.png'\n");
        fprintf(f_blenderi, "            bpy.ops.render.render(write_still=1)\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        ob[file].hide_render = True\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "    for file in range(0, nfiles + 1):\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        ob[file].hide_render = False\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "if animation == 1:\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "    arm = bpy.data.objects['Armature']\n");
        fprintf(f_blenderi, "    arm.hide_render = True\n");
        fprintf(f_blenderi, "    arm.hide_viewport = False\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "    for action in bpy.data.actions: bpy.data.actions.remove(action)\n");
        fprintf(f_blenderi, "    \n");
        fprintf(f_blenderi, "    bpy.context.view_layer.objects.active = arm\n");
        fprintf(f_blenderi, "    bpy.ops.object.mode_set(mode='EDIT', toggle=False)\n");
        fprintf(f_blenderi, "    \n");
        fprintf(f_blenderi, "    edit_bones = arm.data.edit_bones\n");
        fprintf(f_blenderi, "    for bone in edit_bones:\n");
        fprintf(f_blenderi, "        edit_bones.remove(bone)\n");
        fprintf(f_blenderi, "    \n");
        fprintf(f_blenderi, "    armt = []\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "    for file in range(0, nfiles + 1):\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        armt.append(edit_bones.new('Bone' + str(file)))\n");
        fprintf(f_blenderi, "        armt[file].head = (0.0, -1.0, 0.0)\n");
        fprintf(f_blenderi, "        armt[file].tail = (0.0, 0.0, 0.0)\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "    bpy.ops.object.mode_set(mode='OBJECT')\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "    arm.animation_data_create()\n");
        fprintf(f_blenderi, "    arm.animation_data.action = bpy.data.actions.new(name='AnimationAction')\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "    current_frame = 0\n");
        fprintf(f_blenderi, "    pos = []\n");
        fprintf(f_blenderi, "    scn.frame_start = 0\n");
        fprintf(f_blenderi, "    scn.frame_end = ( nfiles + 1 ) * interval - 1\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "    for file in range(0, nfiles + 1):\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        if file < 10: filetext = '000' + str(file)\n");
        fprintf(f_blenderi, "        elif file < 100: filetext = '00' + str(file)\n");
        fprintf(f_blenderi, "        elif file < 1000: filetext = '0' + str(file)\n");
        fprintf(f_blenderi, "        else: filetext = str(file)\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        if impfiles == 0: ob.append(bpy.data.objects[pf + 'Surface' + str(filetext)])\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        ob[file].parent=arm\n");
        fprintf(f_blenderi, "        ob[file].parent_type = 'BONE'\n");
        fprintf(f_blenderi, "        ob[file].parent_bone = 'Bone' + str(file)\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        if layers > 0 and subsurface == 1:\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "            ob1[file].parent=arm\n");
        fprintf(f_blenderi, "            ob1[file].parent_type = 'BONE'\n");
        fprintf(f_blenderi, "            ob1[file].parent_bone = 'Bone' + str(file)\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        armf = arm.pose.bones['Bone' + str(file)]\n");   
        fprintf(f_blenderi, "\n");    
        fprintf(f_blenderi, "        pos.append(armf.location)\n");
        fprintf(f_blenderi, "        pos[file][0] += 1000\n");
        fprintf(f_blenderi, "        pos[file][1] += 1000\n");
        fprintf(f_blenderi, "        pos[file][2] += 1000\n");
        fprintf(f_blenderi, "        armf.location = pos[file]\n");
        fprintf(f_blenderi, "        ob[file].hide_render = False\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "    for file in range(0, nfiles + 1):\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        armf = arm.pose.bones['Bone' + str(file)]\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        fcurve = arm.animation_data.action.fcurves.new(data_path='location', index=file)\n");
        fprintf(f_blenderi, "        k1 = fcurve.keyframe_points.insert(frame=current_frame, value=0)\n");
        fprintf(f_blenderi, "        k1.interpolation = 'CONSTANT'\n");
        fprintf(f_blenderi, "        \n");
        fprintf(f_blenderi, "        if file > 0: armf.keyframe_insert(data_path = 'location', frame = current_frame - 1)\n");
        fprintf(f_blenderi, "\n");
        fprintf(f_blenderi, "        pos[file][0] -= 1000\n");
        fprintf(f_blenderi, "        pos[file][1] -= 1000\n");
        fprintf(f_blenderi, "        pos[file][2] -= 1000\n");
        fprintf(f_blenderi, "        armf.location = pos[file]\n");
        fprintf(f_blenderi, "        armf.keyframe_insert(data_path = 'location', frame = current_frame)\n");
        fprintf(f_blenderi, "        \n");
        fprintf(f_blenderi, "        if file < nfiles: \n");
        fprintf(f_blenderi, "            \n");
        fprintf(f_blenderi, "            armf.keyframe_insert(data_path = 'location', frame = current_frame + interval-1)\n");
        fprintf(f_blenderi, "        \n");
        fprintf(f_blenderi, "            pos[file][0] += 1000\n");
        fprintf(f_blenderi, "            pos[file][1] += 1000\n");
        fprintf(f_blenderi, "            pos[file][2] += 1000\n");
        fprintf(f_blenderi, "            armf.location = pos[file]\n");
        fprintf(f_blenderi, "            armf.keyframe_insert(data_path = 'location', frame = current_frame + interval)\n");
        fprintf(f_blenderi, "            \n");
        fprintf(f_blenderi, "        current_frame += interval\n");
        fprintf(f_blenderi, "    \n");
        fprintf(f_blenderi, "    bpy.ops.nla.bake(frame_start=0, frame_end=current_frame-interval, step=1, only_selected=False, \n");
        fprintf(f_blenderi, "        visual_keying=False, clear_constraints=False, clear_parents=False, \n");
        fprintf(f_blenderi, "        use_current_action=False, clean_curves=False, bake_types={'POSE'})\n");

        fclose(f_blenderi);
    }


// -- STOP --- Writing python script for import to Blender ------------------------------------------------------


// -- START -- Writing python script for animation in Unreal Engine ---------------------------------------------


    if ( sico.MULT == 0 ) {

        sprintf(path, "%sunimport.py", outvr); // python script
        f_unreali=fopen(path, "w");

        fprintf(f_unreali, "import unreal\n");
        fprintf(f_unreali, "\n");
        fprintf(f_unreali, "# **********************************************************************\n");
        fprintf(f_unreali, "\n");
        fprintf(f_unreali, "# The following parameters are partly set in r.avaflow, but should be revised\n");
        fprintf(f_unreali, "\n");
        fprintf(f_unreali, "nfiles = %i # number of time steps\n", nout-1);
        fprintf(f_unreali, "blendfile = '' # prefix from blender file\n");
        fprintf(f_unreali, "interval = 4 # length of time step in frames\n");
        fprintf(f_unreali, "framerate = 60 # frame rate\n");
        fprintf(f_unreali, "addactor = '' # additional actor to appear at end ('' to skip)\n");
        fprintf(f_unreali, "addtime = 100 # time for additional actor to appear\n");
        fprintf(f_unreali, "addx = 0.0 # x position of additional actor\n");
        fprintf(f_unreali, "addy = 0.0 # y position of additional actor\n");
        fprintf(f_unreali, "addstartz = 0.0 # start z position of additional actor\n");
        fprintf(f_unreali, "addendz = 0.0 # end z position of additional actor\n");
        fprintf(f_unreali, "addscale = 1.0 # scale of additional actor\n");
        fprintf(f_unreali, "\n");
        fprintf(f_unreali, "# **********************************************************************\n");
        fprintf(f_unreali, "\n");
        fprintf(f_unreali, "# The following code should only be modified by experienced users\n");
        fprintf(f_unreali, "\n");
        fprintf(f_unreali, "frate = unreal.FrameRate(numerator = framerate, denominator = 1)\n");
        fprintf(f_unreali, "\n");
        fprintf(f_unreali, "lseq = unreal.LevelSequenceEditorBlueprintLibrary.get_current_level_sequence()\n");
        fprintf(f_unreali, "lseq.set_display_rate(frate)\n");
        fprintf(f_unreali, "\n");
        fprintf(f_unreali, "lseq.set_playback_start(0)\n");
        fprintf(f_unreali, "lseq.set_playback_end(interval * nfiles + addtime)\n");
        fprintf(f_unreali, "\n");
        fprintf(f_unreali, "actor_system = unreal.get_editor_subsystem(unreal.EditorActorSubsystem)\n");
        fprintf(f_unreali, "actors = unreal.EditorLevelLibrary.get_all_level_actors()\n");
        fprintf(f_unreali, "\n");
        fprintf(f_unreali, "for i in range(0, nfiles+1):\n");
        fprintf(f_unreali, "\n");
        fprintf(f_unreali, "    for actor in actors:\n");
        fprintf(f_unreali, "    \n");
        fprintf(f_unreali, "        if i<10: itext = '000' + str(i)\n");
        fprintf(f_unreali, "        elif i<100: itext = '00' + str(i)\n");
        fprintf(f_unreali, "        elif i<1000: itext = '0' + str(i)\n");
        fprintf(f_unreali, "        else: itext = str(i)\n");
        fprintf(f_unreali, "        \n");
        fprintf(f_unreali, "        if actor.get_actor_label() == blendfile + '_Surface' + itext:\n");
        fprintf(f_unreali, "\n");
        fprintf(f_unreali, "            poss = lseq.add_possessable(actor)\n");
        fprintf(f_unreali, "            unreal.LevelSequenceEditorBlueprintLibrary.refresh_current_level_sequence()\n");
        fprintf(f_unreali, "\n");
        fprintf(f_unreali, "            vistrack = poss.add_track(unreal.MovieSceneVisibilityTrack)\n");
        fprintf(f_unreali, "            vistrack.set_property_name_and_path('bHiddenInGame', 'bHiddenInGame')\n");
        fprintf(f_unreali, "\n");
        fprintf(f_unreali, "            vissection = vistrack.add_section()\n");
        fprintf(f_unreali, "            vissection.set_start_frame_seconds(0)\n");
        fprintf(f_unreali, "            vissection.set_end_frame_seconds(interval * nfiles + addtime)\n");
        fprintf(f_unreali, "\n");
        fprintf(f_unreali, "            vischannel = vissection.get_channels()[0]\n");
        fprintf(f_unreali, "            vischannel.add_key(time=unreal.FrameNumber(0), new_value=False)\n");
        fprintf(f_unreali, "            vischannel.add_key(time=unreal.FrameNumber(interval * i), new_value=True)\n");
        fprintf(f_unreali, "            if i < nfiles: vischannel.add_key(time=unreal.FrameNumber(interval * i + interval), new_value=False)\n");
        fprintf(f_unreali, "\n");
        fprintf(f_unreali, "if not addactor == '':\n");
        fprintf(f_unreali, "\n");
        fprintf(f_unreali, "    for actor in actors:\n");
        fprintf(f_unreali, "\n");
        fprintf(f_unreali, "        if actor.get_actor_label() == addactor:\n");
        fprintf(f_unreali, "\n");
        fprintf(f_unreali, "            poss = lseq.add_possessable(actor)\n");
        fprintf(f_unreali, "            unreal.LevelSequenceEditorBlueprintLibrary.refresh_current_level_sequence()\n");
        fprintf(f_unreali, "\n");
        fprintf(f_unreali, "            vistrack = poss.add_track(unreal.MovieSceneVisibilityTrack)\n");
        fprintf(f_unreali, "            vistrack.set_property_name_and_path('bHiddenInGame', 'bHiddenInGame')\n");
        fprintf(f_unreali, "\n");
        fprintf(f_unreali, "            vissection = vistrack.add_section()\n");
        fprintf(f_unreali, "            vissection.set_start_frame_seconds(0)\n");
        fprintf(f_unreali, "            vissection.set_end_frame_seconds(interval * nfiles + addtime)\n");
        fprintf(f_unreali, "\n");
        fprintf(f_unreali, "            vischannel = vissection.get_channels()[0]\n");
        fprintf(f_unreali, "            vischannel.add_key(time=unreal.FrameNumber(interval * nfiles -1), new_value=False)\n");
        fprintf(f_unreali, "            vischannel.add_key(time=unreal.FrameNumber(interval * nfiles), new_value=True)\n");
        fprintf(f_unreali, "\n");
        fprintf(f_unreali, "            transtrack = poss.add_track(unreal.MovieScene3DTransformTrack)\n");
        fprintf(f_unreali, "\n");
        fprintf(f_unreali, "            transsection = transtrack.add_section()\n");
        fprintf(f_unreali, "            transsection.set_start_frame_seconds(0)\n");
        fprintf(f_unreali, "            transsection.set_end_frame_seconds(interval * nfiles + addtime)\n");
        fprintf(f_unreali, "\n");
        fprintf(f_unreali, "            transchannel = transsection.get_channels()[0]\n");
        fprintf(f_unreali, "            transchannel.add_key(time=unreal.FrameNumber(0), new_value=addx)\n");
        fprintf(f_unreali, "\n");
        fprintf(f_unreali, "            transchannel = transsection.get_channels()[1]\n");
        fprintf(f_unreali, "            transchannel.add_key(time=unreal.FrameNumber(0), new_value=addy)\n");
        fprintf(f_unreali, "\n");
        fprintf(f_unreali, "            transchannel = transsection.get_channels()[2]\n");
        fprintf(f_unreali, "            transchannel.add_key(time=unreal.FrameNumber(interval * nfiles), new_value=addstartz)\n");
        fprintf(f_unreali, "            transchannel.add_key(time=unreal.FrameNumber(int(interval * nfiles + addtime * 0.9)), new_value=addendz)\n");
        fprintf(f_unreali, "\n");
        fprintf(f_unreali, "            transchannel = transsection.get_channels()[6]\n");
        fprintf(f_unreali, "            transchannel.add_key(time=unreal.FrameNumber(0), new_value=addscale)\n");
        fprintf(f_unreali, "            \n");
        fprintf(f_unreali, "            transchannel = transsection.get_channels()[7]\n");
        fprintf(f_unreali, "            transchannel.add_key(time=unreal.FrameNumber(0), new_value=addscale)\n");
        fprintf(f_unreali, "\n");
        fprintf(f_unreali, "            transchannel = transsection.get_channels()[8]\n");
        fprintf(f_unreali, "            transchannel.add_key(time=unreal.FrameNumber(0), new_value=addscale)\n");

        fclose(f_unreali);
    }


// -- STOP --- Writing python script for animation in Unreal Engine ---------------------------------------------


// -- START -- Writing R script for hydrograph plots ------------------------------------------------------------


#ifdef WITHGRASS


    if ( sico.MULT == 1 ) {

        sprintf(path, "%sr.avaflow.multval.R", outrplots); // R script for evaluation of multiple model runs
        f_rmultval=fopen(path, "w\n");

        fprintf(f_rmultval, "#!/usr/bin/R\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "##############################################################################\n");
        fprintf(f_rmultval, "#\n");
        fprintf(f_rmultval, "# MODULE:       r.avaflow.multval.R\n");
        fprintf(f_rmultval, "# AUTHOR:       Martin Mergili\n");
        fprintf(f_rmultval, "#\n");
        fprintf(f_rmultval, "# PURPOSE:      Script for systematic analysis\n");
        fprintf(f_rmultval, "#               of multiple evaluation parameters\n");
        fprintf(f_rmultval, "#\n");
        fprintf(f_rmultval, "# COPYRIGHT:    (c) 2016 - 2024 by the author\n");
        fprintf(f_rmultval, "#               (c) 2020 - 2024 by the University of Graz\n");
        fprintf(f_rmultval, "#               (c) 2016 - 2021 by the BOKU University, Vienna\n");
        fprintf(f_rmultval, "#               (c) 2016 - 2020 by the University of Vienna\n");
        fprintf(f_rmultval, "#               (c) 1993 - 2024 by the R Development Core Team\n");
        fprintf(f_rmultval, "#\n");
        fprintf(f_rmultval, "#               This program is free software under the GNU General Public\n");
        fprintf(f_rmultval, "#               License (>=v2).\n");
        fprintf(f_rmultval, "#\n");
        fprintf(f_rmultval, "##############################################################################\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "# Loading library\n");
        fprintf(f_rmultval, "options(warn = -1)\n");
        fprintf(f_rmultval, "library(stats)\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "# Importing arguments\n");
        fprintf(f_rmultval, "args <- commandArgs(trailingOnly = TRUE)\n");
        fprintf(f_rmultval, "prefix <- args[1]\n");
        fprintf(f_rmultval, "phases <- as.integer(args[2])\n");
        fprintf(f_rmultval, "nruns <-as.integer(args[3])\n");
        fprintf(f_rmultval, "ipar1 <- as.integer(args[4])\n");
        fprintf(f_rmultval, "ipar2 <- as.integer(args[5])\n");
        fprintf(f_rmultval, "obstype <- as.character(args[6])\n");
        fprintf(f_rmultval, "ictrlpoint <- as.character(args[7])\n");
        fprintf(f_rmultval, "cat(paste('Creating performance heat map plots for ', obstype, ' ... ', sep=''))\n");  
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "# Merging evaluation files from different model runs\n");
        fprintf(f_rmultval, "idata = paste(prefix, '_results/', prefix, '_files/', prefix, '_evaluationh.txt', sep= '')\n");
        fprintf(f_rmultval, "xdata <- read.table(idata, sep='\\t', colClasses = 'character')\n");
        fprintf(f_rmultval, "for ( i in 1:nruns ) {\n");
        fprintf(f_rmultval, "    jdata = paste(prefix, '_results/', prefix, '_files/', prefix, '_evaluation', as.character(i), '.txt', sep= '')\n");
        fprintf(f_rmultval, "    ydata <- read.table(jdata, sep='\\t', colClasses = 'character')\n");
        fprintf(f_rmultval, "    xdata <- rbind(xdata, ydata)\n");
        fprintf(f_rmultval, "}\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "# Merging aimec files from different model runs\n");
        fprintf(f_rmultval, "iadata = paste(prefix, '_results/', prefix, '_aimec/', prefix, '_aimech.txt', sep= '')\n");
        fprintf(f_rmultval, "xadata <- read.table(iadata, sep=' ', colClasses = 'character')\n");
        fprintf(f_rmultval, "for ( i in 1:nruns ) {\n");
        fprintf(f_rmultval, "    jadata = paste(prefix, '_results/', prefix, '_aimec/', prefix, '_aimec', as.character(i), '.txt', sep= '')\n");
        fprintf(f_rmultval, "    yadata <- read.table(jadata, sep=' ', colClasses = 'character')\n");
        fprintf(f_rmultval, "    xadata <- rbind(xadata, yadata)\n");
        fprintf(f_rmultval, "}\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "# Defining evaluation parameter file name\n");
        fprintf(f_rmultval, "invalid = paste(prefix, '_results/', prefix, '_files/', prefix, '_evaluation.txt', sep= '')\n");
        fprintf(f_rmultval, "write.table(xdata, file = invalid, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "# Defining aimec file name\n");
        fprintf(f_rmultval, "invalida = paste(prefix, '_results/', prefix, '_aimec/', prefix, '_aimec.txt', sep= '')\n");
        fprintf(f_rmultval, "write.table(xadata, file = invalida, sep=' ', quote=FALSE, row.names=FALSE, col.names=FALSE)\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "# Creating vectors from input parameters\n");
        fprintf(f_rmultval, "par1 <- read.table(invalid, skip = 1)[, ipar1 + 1]  #1st parameter\n");
        fprintf(f_rmultval, "par2 <- read.table(invalid, skip = 1)[, ipar2 + 1]  #2nd parameter\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "# Building geometry of profile plot\n");
        fprintf(f_rmultval, "maxx <- max(par1, na.rm = TRUE)  #maximum of horizontal coordinate\n");
        fprintf(f_rmultval, "minx <- min(par1, na.rm = TRUE)  #minimum of horizontal coordinate\n");
        fprintf(f_rmultval, "maxy <- max(par2, na.rm = TRUE)  #maximum of vertical coordinate\n");
        fprintf(f_rmultval, "miny <- min(par2, na.rm = TRUE)  #minimum of vertical coordinate\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "# Defining vectors characterizing each evaluation parameter\n");
        fprintf(f_rmultval, "if (phases==1) {\n");
        fprintf(f_rmultval, "    colparams <- 18\n");
        fprintf(f_rmultval, "} else if (phases==3) {\n");
        fprintf(f_rmultval, "    colparams <- 55\n");
        fprintf(f_rmultval, "}\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "if (obstype == 'i') {\n");
        fprintf(f_rmultval, "    coldata <- c(8, 16, 17, 19, 20, 21)\n");
        fprintf(f_rmultval, "} else if (obstype == 'd') {\n");
        fprintf(f_rmultval, "    coldata <- c(11, 26, 27, 29, 30, 31)\n");
        fprintf(f_rmultval, "} else if (obstype == 't') {\n");
        fprintf(f_rmultval, "    coldata <- c(31 + 3 * as.numeric(ictrlpoint))\n");
        fprintf(f_rmultval, "}\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "if (obstype == 't') {\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "    colname <- c(paste('treachrat', ictrlpoint, sep = ''))\n");
        fprintf(f_rmultval, "    coltext <- c('Time of reach ratio')\n");
        fprintf(f_rmultval, "    mdig <- c(2)\n");
        fprintf(f_rmultval, "    unit <- c('')\n");
        fprintf(f_rmultval, "    ngraphs <- 1\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "} else {\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "    colname <- c('Lratio', 'CSI', 'HSS', 'D2PC', 'FoC', 'SPI')\n");
        fprintf(f_rmultval, "    coltext <- c('Excess travel distance ratio', 'Critical success index', 'Heidke skill score', \n");
        fprintf(f_rmultval, "        'Distance to perfect classification', 'Factor of conservativeness', 'Synthetic Performance Index')\n");
        fprintf(f_rmultval, "    mdig <- c(2, 2, 2, 2, 2, 2)\n");
        fprintf(f_rmultval, "    unit <- c('', '', '', '', '', '')\n");
        fprintf(f_rmultval, "    ngraphs <- 6\n");
        fprintf(f_rmultval, "}\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "# Starting loop over all evaluation parameters\n");
        fprintf(f_rmultval, "for (i in 1:ngraphs) {\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "    # Reading evaluation parameter-specific data from file\n");
        fprintf(f_rmultval, "    data <- read.table(invalid, skip = 1)[, colparams+coldata[i]]\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "    # Creating plot file\n");
        fprintf(f_rmultval, "    profileplot = paste(prefix, '_results/', prefix, '_plots/', prefix, '_multval_', \n");
        fprintf(f_rmultval, "        obstype, '_', colname[i], '.png', sep = '')\n");
        fprintf(f_rmultval, "    png(filename = profileplot, width = 10, height = 10, units = 'cm', res = 300)\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "    # Defining margins\n");
        fprintf(f_rmultval, "    par(mar = c(3.2, 3.2, 3.2, 3.2))\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "    # Defining thresholds and colours\n");
        fprintf(f_rmultval, "    if (obstype == 't') {\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "        data1 <- 0.4\n");
        fprintf(f_rmultval, "        data2 <- 0.7\n");
        fprintf(f_rmultval, "        data3 <- 1\n");
        fprintf(f_rmultval, "        data4 <- 1.4\n");
        fprintf(f_rmultval, "        data5 <- 2.5\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "        plotblue <- 2 * log10(data)\n");
        fprintf(f_rmultval, "        plotred <- -2 * log10(data)\n");
        fprintf(f_rmultval, "        plotblue[plotblue < 0] <- 0\n");
        fprintf(f_rmultval, "        plotblue[plotblue > 1] <- 1\n");
        fprintf(f_rmultval, "        plotred[plotred < 0] <- 0\n");
        fprintf(f_rmultval, "        plotred[plotred > 1] <- 1\n");
        fprintf(f_rmultval, "        plotgreen <- 1 - plotred - plotblue\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "        colleg <- c(rgb(0.796, 0.204, 0), rgb(0.31, 0.69, 0), rgb(0, 1, 0), rgb(0, \n");
        fprintf(f_rmultval, "            0.708, 0.292), rgb(0, 0.204, 0.796))\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "    } else if (i == 1) {\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "        maxdata <- 0.5\n");
        fprintf(f_rmultval, "        mindata <- -0.5\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "        data1 <- mindata\n");
        fprintf(f_rmultval, "        data2 <- mindata * 0.5\n");
        fprintf(f_rmultval, "        data3 <- 0\n");
        fprintf(f_rmultval, "        data4 <- maxdata * 0.5\n");
        fprintf(f_rmultval, "        data5 <- maxdata\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "        plotblue <- data/maxdata\n");
        fprintf(f_rmultval, "        plotred <- data/mindata\n");
        fprintf(f_rmultval, "        plotblue[plotblue > 1] <- 1\n");
        fprintf(f_rmultval, "        plotblue[plotblue < 0] <- 0\n");
        fprintf(f_rmultval, "        plotred[plotred > 1] <- 1\n");
        fprintf(f_rmultval, "        plotred[plotred < 0] <- 0\n");
        fprintf(f_rmultval, "        plotgreen <- 1 - plotblue - plotred\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "        colleg <- c(rgb(1, 0, 0), rgb(0.5, 0.5, 0), rgb(0, 1, 0), rgb(0, 0.5, 0.5), \n");
        fprintf(f_rmultval, "            rgb(0, 0, 1))\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "    } else if (i < 4 || i == 6) {\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "        maxdata = 1\n");
        fprintf(f_rmultval, "        mindata = 0\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "        data1 <- 0\n");
        fprintf(f_rmultval, "        data2 <- 0.25\n");
        fprintf(f_rmultval, "        data3 <- 0.5\n");
        fprintf(f_rmultval, "        data4 <- 0.75\n");
        fprintf(f_rmultval, "        data5 <- 1\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "        plotred <- (maxdata - data)/(maxdata - mindata)\n");
        fprintf(f_rmultval, "        plotred[plotred > 1] <- 1\n");
        fprintf(f_rmultval, "        plotred[plotred < 0] <- 0\n");
        fprintf(f_rmultval, "        plotgreen <- (data - mindata)/(maxdata - mindata)\n");
        fprintf(f_rmultval, "        plotgreen[plotgreen > 1] <- 1\n");
        fprintf(f_rmultval, "        plotgreen[plotgreen < 0] <- 0\n");
        fprintf(f_rmultval, "        plotblue <- 0 * data\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "        colleg <- c(rgb(1, 0, 0), rgb(0.75, 0.25, 0), rgb(0.5, 0.5, 0), rgb(0.25, \n");
        fprintf(f_rmultval, "            0.75, 0), rgb(0, 1, 0))\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "    } else if (i < 5) {\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "        maxdata = 1\n");
        fprintf(f_rmultval, "        mindata = 0\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "        data1 <- 0\n");
        fprintf(f_rmultval, "        data2 <- 0.25\n");
        fprintf(f_rmultval, "        data3 <- 0.5\n");
        fprintf(f_rmultval, "        data4 <- 0.75\n");
        fprintf(f_rmultval, "        data5 <- 1\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "        plotred <- (data - mindata)/(maxdata - mindata)\n");
        fprintf(f_rmultval, "        plotred[plotred > 1] <- 1\n");
        fprintf(f_rmultval, "        plotred[plotred < 0] <- 0\n");
        fprintf(f_rmultval, "        plotgreen <- (maxdata - data)/(maxdata - mindata)\n");
        fprintf(f_rmultval, "        plotgreen[plotgreen > 1] <- 1\n");
        fprintf(f_rmultval, "        plotgreen[plotgreen < 0] <- 0\n");
        fprintf(f_rmultval, "        plotblue <- 0 * data\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "        colleg <- c(rgb(0, 1, 0), rgb(0.25, 0.75, 0), rgb(0.5, 0.5, 0), rgb(0.75, \n");
        fprintf(f_rmultval, "            0.25, 0), rgb(1, 0, 0))\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "    } else {\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "        data1 <- 0.1\n");
        fprintf(f_rmultval, "        data2 <- 0.32\n");
        fprintf(f_rmultval, "        data3 <- 1\n");
        fprintf(f_rmultval, "        data4 <- 3.16\n");
        fprintf(f_rmultval, "        data5 <- 10\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "        plotblue <- log10(data)\n");
        fprintf(f_rmultval, "        plotred <- -log10(data)\n");
        fprintf(f_rmultval, "        plotblue[plotblue < 0] <- 0\n");
        fprintf(f_rmultval, "        plotblue[plotblue > 1] <- 1\n");
        fprintf(f_rmultval, "        plotred[plotred < 0] <- 0\n");
        fprintf(f_rmultval, "        plotred[plotred > 1] <- 1\n");
        fprintf(f_rmultval, "        plotgreen <- 1 - plotred - plotblue\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "        colleg <- c(rgb(1, 0, 0), rgb(0.5, 0.5, 0), rgb(0, 1, 0), rgb(0, 0.5, 0.5), \n");
        fprintf(f_rmultval, "            rgb(0, 0, 1))\n");
        fprintf(f_rmultval, "    }\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "    plotblue[data == -9999] <- 0.9\n");
        fprintf(f_rmultval, "    plotred[data == -9999] <- 0.9\n");
        fprintf(f_rmultval, "    plotgreen[data == -9999] <- 0.9\n");
        fprintf(f_rmultval, "    plotcol <- rgb(plotred, plotgreen, plotblue)\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "    # Plotting points\n");
        fprintf(f_rmultval, "    plot(x = par1, y = par2, axes = F, xlab = NA, ylab = NA, type = 'p', xlim = c(minx, \n");
        fprintf(f_rmultval, "        maxx), ylim = c(miny, maxy), pch = 15, cex = 3, col = plotcol)\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "    # Plotting header\n");
        fprintf(f_rmultval, "    par(xpd = TRUE)\n");
        fprintf(f_rmultval, "    text(x = minx + (maxx - minx)/2, y = maxy + (maxy - miny)/4, labels = coltext[i], \n");
        fprintf(f_rmultval, "        cex = 1.2, col = 'black')\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "    # Plotting legend\n");
        fprintf(f_rmultval, "    ltext = vector('expression', 5)\n");
        fprintf(f_rmultval, "    ltext[1] <- substitute(expression(text ~ funit), list(text = format(round(as.numeric(data1), \n");
        fprintf(f_rmultval, "        mdig[i]), nsmall = mdig[i]), funit = format(unit[i])))[2]\n");
        fprintf(f_rmultval, "    ltext[2] <- substitute(expression(text ~ funit), list(text = format(round(as.numeric(data2), \n");
        fprintf(f_rmultval, "        mdig[i]), nsmall = mdig[i]), funit = format(unit[i])))[2]\n");
        fprintf(f_rmultval, "    ltext[3] <- substitute(expression(text ~ funit), list(text = format(round(as.numeric(data3), \n");
        fprintf(f_rmultval, "        mdig[i]), nsmall = mdig[i]), funit = format(unit[i])))[2]\n");
        fprintf(f_rmultval, "    ltext[4] <- substitute(expression(text ~ funit), list(text = format(round(as.numeric(data4), \n");
        fprintf(f_rmultval, "        mdig[i]), nsmall = mdig[i]), funit = format(unit[i])))[2]\n");
        fprintf(f_rmultval, "    ltext[5] <- substitute(expression(text ~ funit), list(text = format(round(as.numeric(data5), \n");
        fprintf(f_rmultval, "        mdig[i]), nsmall = mdig[i]), funit = format(unit[i])))[2]\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "    legend('topleft', legend = ltext, pch = 15, pt.cex = 3, col = colleg, text.col = 'black', \n");
        fprintf(f_rmultval, "        bty = 'n', inset = c(-0.2, -0.15), horiz = TRUE, x.intersp = 1)\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "    # Plotting bounding box and axes\n");
        fprintf(f_rmultval, "    box()\n");
        fprintf(f_rmultval, "    axis(side = 1, tck = -0.02, labels = NA)  #x axis\n");
        fprintf(f_rmultval, "    axis(side = 1, lwd = 0, line = -0.4)\n");
        fprintf(f_rmultval, "    axis(side = 2, tck = -0.02, labels = NA)  #y axis\n");
        fprintf(f_rmultval, "    axis(side = 2, lwd = 0, line = -0.4)\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "    # Plotting axis labels\n");
        fprintf(f_rmultval, "    mtext(side = 1, paste('Parameter', as.character(ipar1), sep = ' '), line = 1.8)  #x axis\n");
        fprintf(f_rmultval, "    mtext(side = 2, paste('Parameter', as.character(ipar2), sep = ' '), line = 1.8)  #y axis\n");
        fprintf(f_rmultval, "\n");
        fprintf(f_rmultval, "    # Closing plot file\n");
        fprintf(f_rmultval, "    dev.off()\n");
        fprintf(f_rmultval, "}\n");
        fprintf(f_rmultval, "cat('completed.\\n')\n"); 

        fclose(f_rmultval);

        sprintf(path, "%sr.avaflow.roc.R", outrplots); // R script for ROC plots
        f_rroc=fopen(path, "w\n");

        fprintf(f_rroc, "#!/usr/bin/R\n");
        fprintf(f_rroc, "\n");
        fprintf(f_rroc, "##############################################################################\n");
        fprintf(f_rroc, "#\n");
        fprintf(f_rroc, "# MODULE:       r.avaflow.roc.R\n");
        fprintf(f_rroc, "#\n");
        fprintf(f_rroc, "# AUTHOR:       Martin Mergili\n");
        fprintf(f_rroc, "#\n");
        fprintf(f_rroc, "# PURPOSE:      The simulation model for avalanche and debris flows\n");
        fprintf(f_rroc, "#               Script for ROC plot\n");
        fprintf(f_rroc, "#\n");
        fprintf(f_rroc, "# COPYRIGHT:    (c) 2013 - 2024 by the author\n");
        fprintf(f_rroc, "#               (c) 2020 - 2024 by the University of Graz\n");
        fprintf(f_rroc, "#               (c) 2013 - 2021 by the BOKU University, Vienna\n");
        fprintf(f_rroc, "#               (c) 2015 - 2020 by the University of Vienna\n");
        fprintf(f_rroc, "#               (c) 1993 - 2024 by the R Development Core Team\n");
        fprintf(f_rroc, "#\n");
        fprintf(f_rroc, "#               This program is free software under the GNU General Public\n");
        fprintf(f_rroc, "#               License (>=v2).\n");
        fprintf(f_rroc, "#\n");
        fprintf(f_rroc, "##############################################################################\n");
        fprintf(f_rroc, "\n");
        fprintf(f_rroc, "# Loading libraries\n");
        fprintf(f_rroc, "options(warn = -1)\n");
        fprintf(f_rroc, "suppressPackageStartupMessages(library('raster'))\n");
        fprintf(f_rroc, "suppressPackageStartupMessages(library('ROCR'))\n");
        fprintf(f_rroc, "suppressPackageStartupMessages(library('sp'))\n");
        fprintf(f_rroc, "suppressPackageStartupMessages(library('terra'))\n");
        fprintf(f_rroc, "library(stats, quietly = T)\n");
        fprintf(f_rroc, "library(raster)\n");
        fprintf(f_rroc, "library(ROCR)\n");
        fprintf(f_rroc, "library(sp)\n");
        fprintf(f_rroc, "library(terra)\n");
        fprintf(f_rroc, "\n");
        fprintf(f_rroc, "# Importing arguments (defined in r.avaflow.py)\n");
        fprintf(f_rroc, "args <- commandArgs(trailingOnly = TRUE)\n");
        fprintf(f_rroc, "temppath <- args[1]\n");
        fprintf(f_rroc, "prefix <- args[2]\n");
        fprintf(f_rroc, "mode <- as.integer(args[3])\n");
        fprintf(f_rroc, "normalization <- as.integer(args[4])\n");
        fprintf(f_rroc, "type <- as.integer(args[5])\n");
        fprintf(f_rroc, "\n");
        fprintf(f_rroc, "# Creating plot file\n");
        fprintf(f_rroc, "cat(paste('Creating ROC plot ', prefix, '_roc_iii ... ', sep=''))\n"); 
        fprintf(f_rroc, "if (normalization == 1 && type == 1) rocplot = paste(prefix, '_results/', prefix, \n");
        fprintf(f_rroc, "    '_plots/', prefix, '_roc_iii.png', sep= '')\n");
        fprintf(f_rroc, "if (normalization == 2 && type == 1) rocplot = paste(prefix, '_results/', prefix, \n");
        fprintf(f_rroc, "    '_plots/', prefix, '_roc_iii_n.png', sep= '')\n");
        fprintf(f_rroc, "if (normalization == 1 && type == 2) rocplot = paste(prefix, '_results/', prefix, \n");
        fprintf(f_rroc, "    '_plots/', prefix, '_roc_dii.png', sep= '')\n");
        fprintf(f_rroc, "if (normalization == 2 && type == 2) rocplot = paste(prefix, '_results/', prefix, \n");
        fprintf(f_rroc, "    '_plots/', prefix, '_roc_dii_n.png', sep= '')\n");
        fprintf(f_rroc, "png(filename = rocplot, width = 12, height = 10, units = 'cm', res = 300)\n");
        fprintf(f_rroc, "\n");
        fprintf(f_rroc, "# Defining margins\n");
        fprintf(f_rroc, "par(mar = c(3.2, 3.2, 0.5, 0.5))\n");
        fprintf(f_rroc, "\n");
        fprintf(f_rroc, "# Initializing ROC plot\n");
        fprintf(f_rroc, "plot(x = c(-2, 2), y = c(-2, 2), col = 'gray', lty = 3, lwd = 1, xlim = c(-0.05, \n");
        fprintf(f_rroc, "    1), ylim = c(0, 1.02), axes = FALSE, xlab = NA, ylab = NA)\n");
        fprintf(f_rroc, "clip(0, 1, 0, 1)\n");
        fprintf(f_rroc, "abline(a = 0, b = 1, col = 'gray', lty = 2, lwd = 1)\n");
        fprintf(f_rroc, "clip(-1, 2, -1, 2)\n");
        fprintf(f_rroc, "\n");
        fprintf(f_rroc, "# Initializing vectors\n");
        fprintf(f_rroc, "arocb <- vector('numeric', 1)\n");
        fprintf(f_rroc, "aroc <- vector('expression', 1)\n");
        fprintf(f_rroc, "lwidth <- vector('character', 1)\n");
        fprintf(f_rroc, "lltp <- vector('integer', 1)\n");
        fprintf(f_rroc, "lcol <- vector('character', 1)\n");
        fprintf(f_rroc, "\n");
        fprintf(f_rroc, "# Defining formatting of ROC plot\n");
        fprintf(f_rroc, "i = 1\n");
        fprintf(f_rroc, "mval = 0\n");
        fprintf(f_rroc, "ia = ''\n");
        fprintf(f_rroc, "lintsp = 1\n");
        fprintf(f_rroc, "mstring = ''\n");
        fprintf(f_rroc, "lwidth[i] = '2'\n");
        fprintf(f_rroc, "lltp[i] = 1\n");
        fprintf(f_rroc, "if (mode == 1) lcol[i] = 'aquamarine4'\n");
        fprintf(f_rroc, "if (mode == 2) lcol[i] = 'cornflowerblue'\n");
        fprintf(f_rroc, "if (mode == 3) lcol[i] = 'darkorange'\n");
        fprintf(f_rroc, "if (mode == 4) lcol[i] = 'brown1'\n");
        fprintf(f_rroc, "\n");
        fprintf(f_rroc, "# Generating data frames from raster maps\n");
        fprintf(f_rroc, "observ <- paste(temppath, '/observed', mstring, '.asc', sep= '')\n");
        fprintf(f_rroc, "ind <- paste(temppath, '/index', mstring, '.asc', sep= '')\n");
        fprintf(f_rroc, "observedr <- raster(observ)  #observation\n");
        fprintf(f_rroc, "indexr <- raster(ind)  #index\n");
        fprintf(f_rroc, "observed <- t(as.matrix(observedr))\n");
        fprintf(f_rroc, "index <- t(as.matrix(indexr))\n");
        fprintf(f_rroc, "\n");
        fprintf(f_rroc, "# Preprocessing data frames\n");
        fprintf(f_rroc, "index[index < 0] <- NaN\n");
        fprintf(f_rroc, "if (normalization == 2) index[index == 0 & observed == 0] <- NaN\n");
        fprintf(f_rroc, "observed[observed > 0] <- 1\n");
        fprintf(f_rroc, "observed[observed < 0] <- NaN\n");
        fprintf(f_rroc, "index[is.na(observed) == TRUE] <- NaN\n");
        fprintf(f_rroc, "observed[is.na(index) == TRUE] <- NaN\n");
        fprintf(f_rroc, "\n");
        fprintf(f_rroc, "indexf <- index[!is.na(index)]\n");
        fprintf(f_rroc, "observedf <- observed[!is.na(observed)]\n");
        fprintf(f_rroc, "\n");
        fprintf(f_rroc, "if (normalization == 2) {\n");
        fprintf(f_rroc, "\n");
        fprintf(f_rroc, "    npos <- sum(observedf == 1, na.rm = TRUE)\n");
        fprintf(f_rroc, "    nneg <- sum(observedf == 0, na.rm = TRUE)\n");
        fprintf(f_rroc, "    nadd <- 5 * npos - nneg\n");
        fprintf(f_rroc, "\n");
        fprintf(f_rroc, "    if (nadd >= 0) {\n");
        fprintf(f_rroc, "        vectadd <- rep(0, nadd)\n");
        fprintf(f_rroc, "        obsroc <- c(observedf, vectadd)\n");
        fprintf(f_rroc, "        indroc <- c(indexf, vectadd)\n");
        fprintf(f_rroc, "        ctrlval = 1\n");
        fprintf(f_rroc, "\n");
        fprintf(f_rroc, "    } else {\n");
        fprintf(f_rroc, "        ctrlval = 0\n");
        fprintf(f_rroc, "    }\n");
        fprintf(f_rroc, "} else {\n");
        fprintf(f_rroc, "    obsroc <- observedf\n");
        fprintf(f_rroc, "    indroc <- indexf\n");
        fprintf(f_rroc, "    ctrlval = 1\n");
        fprintf(f_rroc, "}\n");
        fprintf(f_rroc, "\n");
        fprintf(f_rroc, "# Writing averages of observation and model to output file\n");
        fprintf(f_rroc, "if (normalization == 2 && ctrlval == 1) {\n");
        fprintf(f_rroc, "    write(c(paste('AVGobs', as.character(round(mean(obsroc, na.rm = TRUE), 4)), sep = '\\t'), \n");
        fprintf(f_rroc, "        paste('AVGmodel', as.character(round(mean(indroc, na.rm = TRUE), 4)), sep = '\\t')), \n");
        fprintf(f_rroc, "        file = paste(prefix, '_results/', prefix, '_files/', prefix, '_averages.txt', \n");
        fprintf(f_rroc, "            sep= ''), append = FALSE)\n");
        fprintf(f_rroc, "}\n");
        fprintf(f_rroc, "\n");
        fprintf(f_rroc, "if (ctrlval == 1) {\n");
        fprintf(f_rroc, "\n");
        fprintf(f_rroc, "    # Performing ROC analysis\n");
        fprintf(f_rroc, "    roc <- prediction(indexf, observedf)  #ROC prediction\n");
        fprintf(f_rroc, "    proc <- performance(roc, 'tpr', 'fpr')  #performance\n");
        fprintf(f_rroc, "    arocv <- performance(roc, 'auc')  #area under curve\n");
        fprintf(f_rroc, "\n");
        fprintf(f_rroc, "    # Formatting area under curve value for display\n");
        fprintf(f_rroc, "    arocb[i] = as.numeric(arocv@y.values)\n");
        fprintf(f_rroc, "    aroc[i] <- substitute(expression(iaa ~ aroca), list(iaa = format(ia), aroca = format(round(as.numeric(arocv@y.values), \n");
        fprintf(f_rroc, "        digits = 3), nsmall = 3)))[2]\n");
        fprintf(f_rroc, "\n");
        fprintf(f_rroc, "    # Plotting ROC curve\n");
        fprintf(f_rroc, "    plot(proc, add = TRUE, axes = F, xlab = NA, ylab = NA, lwd = lwidth[i], col = lcol[i], \n");
        fprintf(f_rroc, "        lty = lltp[i])  #ROC curve\n");
        fprintf(f_rroc, "    mval = mval + 1\n");
        fprintf(f_rroc, "\n");
        fprintf(f_rroc, "} else {\n");
        fprintf(f_rroc, "    arocb[i] <- NA\n");
        fprintf(f_rroc, "    aroc[i] <- NA\n");
        fprintf(f_rroc, "}\n");
        fprintf(f_rroc, "\n");
        fprintf(f_rroc, "# Plotting bounding box and axes\n");
        fprintf(f_rroc, "box()\n");
        fprintf(f_rroc, "axis(side = 1, tck = -0.02, labels = NA)\n");
        fprintf(f_rroc, "axis(side = 2, tck = -0.02, labels = NA)\n");
        fprintf(f_rroc, "axis(side = 1, lwd = 0, line = -0.4)\n");
        fprintf(f_rroc, "axis(side = 2, lwd = 0, line = -0.4)\n");
        fprintf(f_rroc, "\n");
        fprintf(f_rroc, "# Plotting axis labels\n");
        fprintf(f_rroc, "mtext(side = 1, 'rFP/rON', line = 1.9)\n");
        fprintf(f_rroc, "mtext(side = 2, 'rTP/rOP', line = 2)\n");
        fprintf(f_rroc, "\n");
        fprintf(f_rroc, "# Plotting legend\n");
        fprintf(f_rroc, "aroch <- vector('expression', 1)\n");
        fprintf(f_rroc, "aroch[1] <- substitute(expression(AUROC))[2]\n");
        fprintf(f_rroc, "legend('bottomright', pch = NA, lty = c(0, lltp), lwd = c(0, lwidth), col = c('black', \n");
        fprintf(f_rroc, "    lcol), legend = c(aroch[1], aroc), y.intersp = c(1, lintsp), text.col = c('black', \n");
        fprintf(f_rroc, "    lcol), bty = 'n', adj = c(0.1, 0.5))\n");
        fprintf(f_rroc, "\n");
        fprintf(f_rroc, "# Writing area under curve to text file if ( ctrlval==1 ) { if ( normalization==1\n");
        fprintf(f_rroc, "# ) write(paste(prefix, '\t', round(as.numeric(arocv@y.values),3), sep=''),\n");
        fprintf(f_rroc, "# file='aucroc.txt', append=TRUE) if ( normalization==2 ) write(paste(prefix,\n");
        fprintf(f_rroc, "# '\t', round(as.numeric(arocv@y.values),3), sep=''), file='aucrocn.txt',\n");
        fprintf(f_rroc, "# append=TRUE) }\n");
        fprintf(f_rroc, "\n");
        fprintf(f_rroc, "# Closing plot file\n");
        fprintf(f_rroc, "invisible(dev.off())\n");
        fprintf(f_rroc, "cat('completed.\\n')\n");       

        fclose(f_rroc);
    }


#endif


    if ( hydrograph == 1 ) {

        sprintf(path, "%sr.avaflow.hydrograph.R", outrplots); // R script for hydrograph plots
        f_rhydrograph=fopen(path, "w\n");

        fprintf(f_rhydrograph, "#!/usr/bin/R\n");
        fprintf(f_rhydrograph, "\n");
        fprintf(f_rhydrograph, "##############################################################################\n");
        fprintf(f_rhydrograph, "#\n");
        fprintf(f_rhydrograph, "# MODULE:       r.avaflow.hydrograph.R\n");
        fprintf(f_rhydrograph, "# AUTHOR:       Martin Mergili\n");
        fprintf(f_rhydrograph, "#\n");
        fprintf(f_rhydrograph, "# PURPOSE:      The simulation model for avalanche and debris flows\n");
        fprintf(f_rhydrograph, "#               Script for the creation of hydrograph plots\n");
        fprintf(f_rhydrograph, "#\n");
        fprintf(f_rhydrograph, "# COPYRIGHT:    (c) 2013 - 2024 by the author\n");
        fprintf(f_rhydrograph, "#               (c) 2020 - 2024 by the University of Graz\n");
        fprintf(f_rhydrograph, "#               (c) 2013 - 2021 by the BOKU University, Vienna\n");
        fprintf(f_rhydrograph, "#               (c) 2015 - 2020 by the University of Vienna\n");
        fprintf(f_rhydrograph, "#               (c) 1993 - 2024 by the R Development Core Team\n");
        fprintf(f_rhydrograph, "#\n");
        fprintf(f_rhydrograph, "#               This program is free software under the GNU General Public\n");
        fprintf(f_rhydrograph, "#               License (>=v2).\n");
        fprintf(f_rhydrograph, "#\n");
        fprintf(f_rhydrograph, "##############################################################################\n");
        fprintf(f_rhydrograph, "\n");
        fprintf(f_rhydrograph, "# Loading library\n");
        fprintf(f_rhydrograph, "library(stats, quietly = T)\n");
        fprintf(f_rhydrograph, "cat('Creating hydrograph plots ... ')\n"); 
        fprintf(f_rhydrograph, "\n");
        fprintf(f_rhydrograph, "# Defining arguments\n");
        
    #ifdef WITHGRASS

        fprintf(f_rhydrograph, "wkdir <- ''  #working directory\n" );

    #else

        fprintf(f_rhydrograph, "wkdir <- '%s/'  #working directory\n", wkdir );

    #endif

        fprintf(f_rhydrograph, "prefix <- '%s'  #prefix for file names\n", prefix );
        fprintf(f_rhydrograph, "model <- %i  #type of model\n", sico.MODEL );
        fprintf(f_rhydrograph, "nhydin <- %i  #total number of input hydrographs\n", hydnin );
        fprintf(f_rhydrograph, "nhydout <- %i  #total number of output hydrographs\n", hydnout );
        fprintf(f_rhydrograph, "tint <- %.5f  #length of output time step\n", tout );
        fprintf(f_rhydrograph, "\n");
        fprintf(f_rhydrograph, "for (nhyd in 1:(nhydin+nhydout)) { #loop over all hydrographs\n");
        fprintf(f_rhydrograph, "\n");
        fprintf(f_rhydrograph, "    # Defining data file name\n");
        fprintf(f_rhydrograph, "    hydtable = paste(wkdir, prefix, 'results/', prefix, 'files/', prefix, 'hydinfo', as.character(nhyd),\n");
        fprintf(f_rhydrograph, "        '.txt', sep = '')\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "    # Creating vectors from data\n");
        fprintf(f_rhydrograph, "    thyd <- read.table(hydtable, skip = 1)[, 1]  #time passed\n");
        fprintf(f_rhydrograph, "    hhyd <- read.table(hydtable, skip = 1)[, 2]  #(PHASE 1) flow height\n");
        fprintf(f_rhydrograph, "    vhyd <- read.table(hydtable, skip = 1)[, 3]  #(PHASE 1) velocity\n");
        fprintf(f_rhydrograph, "    ehyd <- read.table(hydtable, skip = 1)[, 4]  #depth of (PHASE 1) entrainment/deposition\n");
        fprintf(f_rhydrograph, "    qhyd <- read.table(hydtable, skip = 1)[, 5]  #(PHASE 1) discharge\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "    tmax <- round(max(thyd, na.rm = TRUE), 0)  #maximum value of time passed\n");
        fprintf(f_rhydrograph, "    hmax <- max(hhyd, na.rm = TRUE)  #maximum value of (PHASE 1) flow height\n");
        fprintf(f_rhydrograph, "    vmax <- max(abs(vhyd), na.rm = TRUE)  #maximum value of (PHASE 1) velocity\n");
        fprintf(f_rhydrograph, "    emax <- max(ehyd, na.rm = TRUE)  #maximum value of depth of (PHASE 1) entrainment/deposition\n");
        fprintf(f_rhydrograph, "    qmax <- max(abs(qhyd), na.rm = TRUE)  #maximum value of (PHASE 1) discharge\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "    if (model > 3) {\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "        hhyd2 <- read.table(hydtable, skip = 1)[, 6]  #PHASE 2 flow height\n");
        fprintf(f_rhydrograph, "        vhyd2 <- read.table(hydtable, skip = 1)[, 7]  #PHASE 2 velocity\n");
        fprintf(f_rhydrograph, "        ehyd2 <- read.table(hydtable, skip = 1)[, 8]  #depth of PHASE 2 entrainment/deposition\n");
        fprintf(f_rhydrograph, "        qhyd2 <- read.table(hydtable, skip = 1)[, 9]  #PHASE 2 discharge\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "        hmax2 <- max(hhyd2, na.rm = TRUE)  #maximum value of PHASE 2 flow height\n");
        fprintf(f_rhydrograph, "        vmax2 <- max(abs(vhyd2), na.rm = TRUE)  #maximum value of PHASE 2 velocity\n");
        fprintf(f_rhydrograph, "        emax2 <- max(ehyd2, na.rm = TRUE)  #maximum value of depth of PHASE 2 entrainment/deposition\n");
        fprintf(f_rhydrograph, "        qmax2 <- max(abs(qhyd2), na.rm = TRUE)  #maximum value of PHASE 2 discharge\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "    } else {\n");
        fprintf(f_rhydrograph, "        hmax2 <- 0\n");
        fprintf(f_rhydrograph, "        qmax2 <- 0\n");
        fprintf(f_rhydrograph, "    }\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "    if (model == 7) {\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "        hhyd3 <- read.table(hydtable, skip = 1)[, 10]  #PHASE 3 flow height\n");
        fprintf(f_rhydrograph, "        vhyd3 <- read.table(hydtable, skip = 1)[, 11]  #PHASE 3 velocity\n");
        fprintf(f_rhydrograph, "        ehyd3 <- read.table(hydtable, skip = 1)[, 12]  #depth of PHASE 3 entrainment/deposition\n");
        fprintf(f_rhydrograph, "        qhyd3 <- read.table(hydtable, skip = 1)[, 13]  #PHASE 3 discharge\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "        hmax3 <- max(hhyd3, na.rm = TRUE)  #maximum value of PHASE 3 flow height\n");
        fprintf(f_rhydrograph, "        vmax3 <- max(abs(vhyd3), na.rm = TRUE)  #maximum value of PHASE 3 velocity\n");
        fprintf(f_rhydrograph, "        emax3 <- max(ehyd3, na.rm = TRUE)  #maximum value of depth of PHASE 3 entrainment/deposition\n");
        fprintf(f_rhydrograph, "        qmax3 <- max(abs(qhyd3), na.rm = TRUE)  #maximum value of PHASE 3 discharge\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "    } else {\n");
        fprintf(f_rhydrograph, "        hmax3 <- 0\n");
        fprintf(f_rhydrograph, "        qmax3 <- 0\n");
        fprintf(f_rhydrograph, "    }\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "    # Defining labels and scaling\n");
        fprintf(f_rhydrograph, "    ltext = vector('expression', 10)  #initializing label vector\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "    if (nhyd <= nhydin) {\n");
        fprintf(f_rhydrograph, "        textcol = 'green'\n");
        fprintf(f_rhydrograph, "        ltext[1] = substitute(expression(Input ~ hydrograph ~ IH ~ fnhyd), list(fnhyd = format(nhyd)))[2]  #title text for input hydrograph\n");
        fprintf(f_rhydrograph, "    } else {\n");
        fprintf(f_rhydrograph, "        textcol = 'purple'\n");
        fprintf(f_rhydrograph, "        ltext[1] = substitute(expression(Output ~ hydrograph ~ OH ~ fnhyd), list(fnhyd = format(nhyd -\n");
        fprintf(f_rhydrograph, "            nhydin)))[2]  #title text for output hydrograph\n");
        fprintf(f_rhydrograph, "    }\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "    ltext[2] = substitute(expression(Time ~ passed ~ italic(T) ~ (s)))[2]  #x axis label\n");
        fprintf(f_rhydrograph, "    ltext[3] = substitute(expression(Flow ~ height ~ italic(H) ~ (m)))[2]  #left y axis label\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "    if (max(qmax, qmax2 + qmax3) < 1000) {\n");
        fprintf(f_rhydrograph, "        # right y axis label:\n");
        fprintf(f_rhydrograph, "        mconv = '1'\n");
        fprintf(f_rhydrograph, "        ltext[4] = substitute(expression(Discharge ~ italic(Q) ~ (m^3 ~ s^-1)))[2]\n");
        fprintf(f_rhydrograph, "    } else if (max(qmax, qmax2 + qmax3) < 1e+06) {\n");
        fprintf(f_rhydrograph, "        mconv = '0.001'\n");
        fprintf(f_rhydrograph, "        ltext[4] = substitute(expression(Discharge ~ italic(Q) ~ (10^3 ~ m^3 ~ s^-1)))[2]\n");
        fprintf(f_rhydrograph, "    } else {\n");
        fprintf(f_rhydrograph, "        mconv = '0.000001'\n");
        fprintf(f_rhydrograph, "        ltext[4] = substitute(expression(Discharge ~ italic(Q) ~ (10^6 ~ m^3 ~ s^-1)))[2]\n");
        fprintf(f_rhydrograph, "    }\n");
        fprintf(f_rhydrograph, "    qmax <- qmax * as.numeric(mconv)  #scaling of right y axis\n");
        fprintf(f_rhydrograph, "    qmax2 <- qmax2 * as.numeric(mconv)\n");
        fprintf(f_rhydrograph, "    qmax3 <- qmax3 * as.numeric(mconv)\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "    # Creating plot file\n");
        fprintf(f_rhydrograph, "    hydplot = paste(wkdir, prefix, 'results/', prefix, 'plots/', prefix, 'hydrograph', nhyd,\n");
        fprintf(f_rhydrograph, "        '.png', sep = '')\n");
        fprintf(f_rhydrograph, "    png(filename = hydplot, width = 15, height = 10, units = 'cm', res = 300)\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "    # Defining external margins\n");
        fprintf(f_rhydrograph, "    par(mar = c(3.2, 3.2, 1, 3.2))\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "    # Defining internal margins\n");
        fprintf(f_rhydrograph, "    if (((hmax + hmax2) > (2 * hmax3)) && ((qmax + qmax2) > (2 * qmax3))) {\n");
        fprintf(f_rhydrograph, "        marg <- 1.1\n");
        fprintf(f_rhydrograph, "        tpos <- -(hmax + hmax2) * 0.85\n");
        fprintf(f_rhydrograph, "        lpos <- 'bottomright'\n");
        fprintf(f_rhydrograph, "    } else if ((hmax2 + hmax3) > 2 * hmax && (qmax2 + qmax3) > 2 * qmax) {\n");
        fprintf(f_rhydrograph, "        marg <- 1.1\n");
        fprintf(f_rhydrograph, "        tpos <- hmax3 * 0.85\n");
        fprintf(f_rhydrograph, "        lpos <- 'topright'\n");
        fprintf(f_rhydrograph, "    } else {\n");
        fprintf(f_rhydrograph, "        marg <- 1.3\n");
        fprintf(f_rhydrograph, "        lpos <- 'bottomright'\n");
        fprintf(f_rhydrograph, "        tpos <- max(hmax + hmax2, hmax3) * 1.25\n");
        fprintf(f_rhydrograph, "    }\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "    # Plotting bars\n");
        fprintf(f_rhydrograph, "    if (model <= 3) {\n");
        fprintf(f_rhydrograph, "        barplot(qhyd * as.numeric(mconv), axes = F, xlab = NA, ylab = NA, ylim = c(-qmax *\n");
        fprintf(f_rhydrograph, "            marg, qmax * marg), border = NA, col = rgb(0.5, 0.85, 0.5))  #discharge\n");
        fprintf(f_rhydrograph, "    } else if (model > 3 && model < 7) {\n");
        fprintf(f_rhydrograph, "        barplot(qhyd * as.numeric(mconv), axes = F, xlab = NA, ylab = NA, ylim = c(-max(qmax2 *\n");
        fprintf(f_rhydrograph, "            marg, qmax * marg), max(qmax2 * marg, qmax * marg)), border = NA, col = rgb(0.45,\n");
        fprintf(f_rhydrograph, "            0.3, 0, 0.5))  #PHASE 1 discharge\n");
        fprintf(f_rhydrograph, "    } else if (model == 7) {\n");
        fprintf(f_rhydrograph, "        barplot(qhyd * as.numeric(mconv) + qhyd2 * as.numeric(mconv), axes = F, xlab = NA,\n");
        fprintf(f_rhydrograph, "            ylab = NA, ylim = c(-max(qmax3 * marg, qmax2 * marg + qmax * marg), max(qmax3 *\n");
        fprintf(f_rhydrograph, "                marg, qmax2 * marg + qmax * marg)), border = NA, col = rgb(0.6, 0, 0,\n");
        fprintf(f_rhydrograph, "                0.5))  #PHASE 1 discharge\n");
        fprintf(f_rhydrograph, "    }\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "    par(new = TRUE)\n");
        fprintf(f_rhydrograph, "    if (model == 7) {\n");
        fprintf(f_rhydrograph, "        barplot(-qhyd3 * as.numeric(mconv), axes = F, xlab = NA, ylab = NA, ylim = c(-max(qmax3 *\n");
        fprintf(f_rhydrograph, "            marg, qmax2 * marg + qmax * marg), max(qmax3 * marg, qmax2 * marg + qmax *\n");
        fprintf(f_rhydrograph, "            marg)), border = NA, col = rgb(0, 0, 0.6, 0.5))  #PHASE 3 discharge\n");
        fprintf(f_rhydrograph, "    }\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "    par(new = TRUE)\n");
        fprintf(f_rhydrograph, "    if (model > 3 && model < 7) {\n");
        fprintf(f_rhydrograph, "        barplot(-qhyd2 * as.numeric(mconv), axes = F, xlab = NA, ylab = NA, ylim = c(-max(qmax2 *\n");
        fprintf(f_rhydrograph, "            marg, qmax * marg), max(qmax2 * marg, qmax * marg)), border = NA, col = rgb(0,\n");
        fprintf(f_rhydrograph, "            0, 0.6, 0.5))  #PHASE 2 discharge\n");
        fprintf(f_rhydrograph, "    } else if (model == 7) {\n");
        fprintf(f_rhydrograph, "        barplot(qhyd2 * as.numeric(mconv), axes = F, xlab = NA, ylab = NA, ylim = c(-max(qmax3 *\n");
        fprintf(f_rhydrograph, "            marg, qmax2 * marg + qmax * marg), max(qmax3 * marg, qmax2 * marg + qmax *\n");
        fprintf(f_rhydrograph, "            marg)), border = NA, col = rgb(0, 0.6, 0, 0.5))  #PHASE 2 discharge\n");
        fprintf(f_rhydrograph, "    }\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "    axis(side = 4, tck = -0.02, labels = NA)  #y axis for velocity and random kinetic energy\n");
        fprintf(f_rhydrograph, "    axis(side = 4, lwd = 0, line = -0.4)\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "    # Plotting lines\n");
        fprintf(f_rhydrograph, "    par(new = TRUE)\n");
        fprintf(f_rhydrograph, "    if (model <= 3) {\n");
        fprintf(f_rhydrograph, "        plot(x = thyd, y = hhyd, lty = 1, lwd = 2, col = rgb(0.45, 0.3, 0), axes = F,\n");
        fprintf(f_rhydrograph, "            xlab = NA, ylab = NA, type = 'l', ylim = c(-hmax * marg, hmax * marg))  #flow height\n");
        fprintf(f_rhydrograph, "    } else if (model > 3 && model < 7) {\n");
        fprintf(f_rhydrograph, "        plot(x = thyd, y = hhyd, lty = 1, lwd = 2, col = rgb(0.45, 0.3, 0), axes = F,\n");
        fprintf(f_rhydrograph, "            xlab = NA, ylab = NA, type = 'l', ylim = c(-max(hmax * marg, hmax2 * marg,\n");
        fprintf(f_rhydrograph, "                hmax3 * marg), max(hmax * marg, hmax2 * marg, hmax3 * marg)), xlim = c(-tint/2,\n");
        fprintf(f_rhydrograph, "                tmax + tint/2))  #PHASE 1 flow height\n");
        fprintf(f_rhydrograph, "        lines(x = thyd, y = -hhyd2, lty = 1, lwd = 2, col = rgb(0, 0, 0.8), xlim = c(-tint/2,\n");
        fprintf(f_rhydrograph, "            tmax + tint/2))  #PHASE 2 flow height\n");
        fprintf(f_rhydrograph, "    } else if (model == 7) {\n");
        fprintf(f_rhydrograph, "        plot(x = thyd, y = hhyd, lty = 1, lwd = 2, col = rgb(0.6, 0, 0), axes = F, xlab = NA,\n");
        fprintf(f_rhydrograph, "            ylab = NA, type = 'l', ylim = c(-max(hmax * marg, hmax2 * marg, hmax3 * marg),\n");
        fprintf(f_rhydrograph, "                max(hmax * marg, hmax2 * marg, hmax3 * marg)), xlim = c(-tint/2, tmax +\n");
        fprintf(f_rhydrograph, "                tint/2))  #PHASE 1 flow height\n");
        fprintf(f_rhydrograph, "        lines(x = thyd, y = hhyd2, lty = 1, lwd = 2, col = rgb(0, 0.6, 0), xlim = c(-tint/2,\n");
        fprintf(f_rhydrograph, "            tmax + tint/2))  #PHASE 2 flow height\n");
        fprintf(f_rhydrograph, "    }\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "    if (model == 7) {\n");
        fprintf(f_rhydrograph, "        lines(x = thyd, y = -hhyd3, lty = 1, lwd = 2, col = rgb(0, 0, 0.6), xlim = c(-tint/2,\n");
        fprintf(f_rhydrograph, "            tmax + tint/2))  #PHASE 3 flow height\n");
        fprintf(f_rhydrograph, "    }\n");
        fprintf(f_rhydrograph, "    abline(0, 0, col = 'black', lwd = 1)  #zero line\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "    # Text and legends\n");
        fprintf(f_rhydrograph, "    text(x = tmax/2, y = tpos, labels = ltext[1], font = 2, col = textcol, cex = 1.1)  #title text\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "    if (model <= 3) {\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "        ltext[5] = substitute(expression(italic(H)))[2]  #flow height\n");
        fprintf(f_rhydrograph, "        ltext[6] = substitute(expression(italic(Q)))[2]  #discharge\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "        legend(lpos, legend = ltext[5], lty = 1, lwd = 2, col = 'brown', text.col = 'black',\n");
        fprintf(f_rhydrograph, "            bty = 'n', inset = c(1.025, -0.025), horiz = TRUE)\n");
        fprintf(f_rhydrograph, "        # flow height\n");
        fprintf(f_rhydrograph, "        legend(lpos, legend = ltext[6], border = NA, fill = rgb(0.5, 0.85, 0.5), text.col = 'black',\n");
        fprintf(f_rhydrograph, "            bty = 'n', inset = c(0.225, -0.025), horiz = TRUE)\n");
        fprintf(f_rhydrograph, "        # discharge\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "    } else if (model > 3 && model < 7) {\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "        ltext[5] = substitute(expression(italic(H)[s]))[2]  #PHASE 1 flow height\n");
        fprintf(f_rhydrograph, "        ltext[6] = substitute(expression(italic(Q)[s]))[2]  #PHASE 1 discharge\n");
        fprintf(f_rhydrograph, "        ltext[7] = substitute(expression(italic(-H)[f]))[2]  #PHASE 2 flow height\n");
        fprintf(f_rhydrograph, "        ltext[8] = substitute(expression(italic(-Q)[f]))[2]  #PHASE 2 discharge\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "        legend(lpos, legend = ltext[5], lty = 1, lwd = 2, col = rgb(0.45, 0.3, 0), text.col = 'black',\n");
        fprintf(f_rhydrograph, "            bty = 'n', inset = c(0.825, -0.025), horiz = TRUE)\n");
        fprintf(f_rhydrograph, "        # PHASE 1 flow height\n");
        fprintf(f_rhydrograph, "        legend(lpos, legend = ltext[7], lty = 1, lwd = 2, col = rgb(0, 0, 0.8), text.col = 'black',\n");
        fprintf(f_rhydrograph, "            bty = 'n', inset = c(0.625, -0.025), horiz = TRUE)\n");
        fprintf(f_rhydrograph, "        # PHASE 2 flow height\n");
        fprintf(f_rhydrograph, "        legend(lpos, legend = ltext[6], border = NA, fill = rgb(0.45, 0.3, 0, 0.5), text.col = 'black',\n");
        fprintf(f_rhydrograph, "            bty = 'n', inset = c(0.225, -0.025), horiz = TRUE)\n");
        fprintf(f_rhydrograph, "        # PHASE 1 discharge\n");
        fprintf(f_rhydrograph, "        legend(lpos, legend = ltext[8], border = NA, fill = rgb(0, 0, 0.6, 0.5), text.col = 'black',\n");
        fprintf(f_rhydrograph, "            bty = 'n', inset = c(0.025, -0.025), horiz = TRUE)\n");
        fprintf(f_rhydrograph, "        # PHASE 2 discharge\n");
        fprintf(f_rhydrograph, "    } else if (model == 7) {\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "        par(cex = 0.75)\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "        ltext[5] = substitute(expression(italic(H)[P1]))[2]  #PHASE 1 flow height\n");
        fprintf(f_rhydrograph, "        ltext[6] = substitute(expression(italic(Q)[P1]))[2]  #PHASE 1 discharge\n");
        fprintf(f_rhydrograph, "        ltext[7] = substitute(expression(italic(H)[P2]))[2]  #PHASE 2 flow height\n");
        fprintf(f_rhydrograph, "        ltext[8] = substitute(expression(italic(Q)[P2]))[2]  #PHASE 2 discharge\n");
        fprintf(f_rhydrograph, "        ltext[9] = substitute(expression(italic(-H)[P3]))[2]  #PHASE 3 flow height\n");
        fprintf(f_rhydrograph, "        ltext[10] = substitute(expression(italic(-Q)[P3]))[2]  #PHASE 3 discharge\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "        legend(lpos, legend = ltext[5], lty = 1, lwd = 2, col = rgb(0.6, 0, 0), text.col = 'black',\n");
        fprintf(f_rhydrograph, "            bty = 'n', inset = c(0.85, -0.025), horiz = TRUE)\n");
        fprintf(f_rhydrograph, "        # PHASE 1 flow height\n");
        fprintf(f_rhydrograph, "        legend(lpos, legend = ltext[7], lty = 1, lwd = 2, col = rgb(0, 0.6, 0), text.col = 'black',\n");
        fprintf(f_rhydrograph, "            bty = 'n', inset = c(0.7, -0.025), horiz = TRUE)\n");
        fprintf(f_rhydrograph, "        # PHASE 2 flow height\n");
        fprintf(f_rhydrograph, "        legend(lpos, legend = ltext[9], lty = 1, lwd = 2, col = rgb(0, 0, 0.6), text.col = 'black',\n");
        fprintf(f_rhydrograph, "            bty = 'n', inset = c(0.55, -0.025), horiz = TRUE)\n");
        fprintf(f_rhydrograph, "        # PHASE 3 flow height\n");
        fprintf(f_rhydrograph, "        legend(lpos, legend = ltext[6], border = NA, fill = rgb(0.6, 0, 0, 0.5), text.col = 'black',\n");
        fprintf(f_rhydrograph, "            bty = 'n', inset = c(0.35, -0.025), horiz = TRUE)\n");
        fprintf(f_rhydrograph, "        # PHASE 1 discharge\n");
        fprintf(f_rhydrograph, "        legend(lpos, legend = ltext[8], border = NA, fill = rgb(0, 0.6, 0, 0.5), text.col = 'black',\n");
        fprintf(f_rhydrograph, "            bty = 'n', inset = c(0.2, -0.025), horiz = TRUE)\n");
        fprintf(f_rhydrograph, "        # PHASE 2 discharge\n");
        fprintf(f_rhydrograph, "        legend(lpos, legend = ltext[10], border = NA, fill = rgb(0, 0, 0.6, 0.5), text.col = 'black',\n");
        fprintf(f_rhydrograph, "            bty = 'n', inset = c(0.05, -0.025), horiz = TRUE)\n");
        fprintf(f_rhydrograph, "        # PHASE 3 discharge\n");
        fprintf(f_rhydrograph, "        par(cex = 1)\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "    }\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "    # Plotting bounding box and axes\n");
        fprintf(f_rhydrograph, "    box()\n");
        fprintf(f_rhydrograph, "    axis(side = 1, tck = -0.02, labels = NA, at = c(0, tmax/10, 2 * tmax/10, 3 * tmax/10,\n");
        fprintf(f_rhydrograph, "        4 * tmax/10, 5 * tmax/10, 6 * tmax/10, 7 * tmax/10, 8 * tmax/10, 9 * tmax/10,\n");
        fprintf(f_rhydrograph, "        tmax))\n");
        fprintf(f_rhydrograph, "    # x axis\n");
        fprintf(f_rhydrograph, "    axis(side = 1, lwd = 0, line = -0.4, at = c(0, tmax/10, 2 * tmax/10, 3 * tmax/10,\n");
        fprintf(f_rhydrograph, "        4 * tmax/10, 5 * tmax/10, 6 * tmax/10, 7 * tmax/10, 8 * tmax/10, 9 * tmax/10,\n");
        fprintf(f_rhydrograph, "        tmax))\n");
        fprintf(f_rhydrograph, "    axis(side = 2, tck = -0.02, labels = NA)  #y axis for lines\n");
        fprintf(f_rhydrograph, "    axis(side = 2, lwd = 0, line = -0.4)\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "    # Plotting axis labels\n");
        fprintf(f_rhydrograph, "    mtext(side = 1, ltext[2], line = 2)  #x axis\n");
        fprintf(f_rhydrograph, "    mtext(side = 2, ltext[3], line = 2)  #y axis for lines\n");
        fprintf(f_rhydrograph, "    if (model <= 3) {\n");
        fprintf(f_rhydrograph, "        mtext(side = 4, ltext[4], line = 2)  #y axis for bars\n");
        fprintf(f_rhydrograph, "    } else if (model > 3 && model < 7) {\n");
        fprintf(f_rhydrograph, "        mtext(side = 4, ltext[4], line = 2)  #y axis for bars\n");
        fprintf(f_rhydrograph, "    } else if (model == 7) {\n");
        fprintf(f_rhydrograph, "        mtext(side = 4, ltext[4], line = 2)  #y axis for bars\n");
        fprintf(f_rhydrograph, "    }\n");
        fprintf(f_rhydrograph, "    \n");
        fprintf(f_rhydrograph, "    # Closing plot file\n");
        fprintf(f_rhydrograph, "    dev.off()\n");
        fprintf(f_rhydrograph, "}\n");
        fprintf(f_rhydrograph, "cat('completed.\\n')\n");

        fclose(f_rhydrograph);

        if ( strcmp( rscript, "None" ) != 0 ) {
        
            sprintf(path, "%sr.avaflow.hydrograph.cmd", outrplots); // cmd script to start R script
            f_chydrograph=fopen(path, "w\n");

            fprintf(f_chydrograph, "\"%s\" \"%s/r.avaflow.hydrograph.R\"\n", rscript, outrplots );       
        
            fclose(f_chydrograph);
        }
    }


// -- STOP --- Writing R script for hydrograph plots ------------------------------------------------------------


// -- START -- Writing R script for profile plots ---------------------------------------------------------------


    if ( sico.PROFILE > 0 ) {

        sprintf(path, "%sr.avaflow.profile.R", outrplots); // R script for profile plots
        f_rprofile=fopen(path, "w\n");

        fprintf(f_rprofile, "#!/usr/bin/R\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "##############################################################################\n");
        fprintf(f_rprofile, "#\n");
        fprintf(f_rprofile, "# MODULE:       r.avaflow.profile.R\n");
        fprintf(f_rprofile, "# AUTHOR:       Martin Mergili\n");
        fprintf(f_rprofile, "#\n");
        fprintf(f_rprofile, "# PURPOSE:      The simulation model for avalanche and debris flows\n");
        fprintf(f_rprofile, "#               Script for the creation of profile plots of model results\n");
        fprintf(f_rprofile, "#\n");
        fprintf(f_rprofile, "# COPYRIGHT:    (c) 2013 - 2024 by the author\n");
        fprintf(f_rprofile, "#               (c) 2020 - 2024 by the University of Graz\n");
        fprintf(f_rprofile, "#               (c) 2013 - 2020 by the BOKU University, Vienna\n");
        fprintf(f_rprofile, "#               (c) 2015 - 2020 by the University of Vienna\n");
        fprintf(f_rprofile, "#               (c) 1993 - 2024 by the R Development Core Team\n");
        fprintf(f_rprofile, "#\n");
        fprintf(f_rprofile, "#               This program is free software under the GNU General Public\n");
        fprintf(f_rprofile, "#               License (>=v2).\n");
        fprintf(f_rprofile, "#\n");
        fprintf(f_rprofile, "##############################################################################\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "# Loading library\n");
        fprintf(f_rprofile, "library(stats, quietly = T)\n");
        fprintf(f_rprofile, "cat('Creating profile plots ... ')\n");  
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "# Defining arguments\n");
        
    #ifdef WITHGRASS

        fprintf(f_rprofile, "wkdir <- ''  #working directory\n" );

    #else

        fprintf(f_rprofile, "wkdir <- '%s/'  #working directory\n", wkdir );

    #endif
        
        fprintf(f_rprofile, "prefix <- '%s' #prefix for file names\n", prefix );
        fprintf(f_rprofile, "model <- %i  #prefix for file names\n", sico.MODEL );
        fprintf(f_rprofile, "maxtimesteps <- %i  #number of time steps\n", nout-1 );
        fprintf(f_rprofile, "cellsize <- %.5f  #pixel size\n", sico.CSZ );
        fprintf(f_rprofile, "hmax <- %.5f  #maximum flow height for display\n", hflow_maxmax );
        fprintf(f_rprofile, "hexagg <- %.5f  #exaggeration of flow height for display\n", phexagg );
        fprintf(f_rprofile, "tint <- %.5f  #length of time step\n", tout );
        fprintf(f_rprofile, "tstop <- %.5f  #total time of simulation\n", tmax );
        fprintf(f_rprofile, "depdef <- %i  #control for observed deposit map\n", sico.HDEPOSIT );
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "# Defining data file name\n");
        fprintf(f_rprofile, "intable = paste(wkdir, prefix, 'results/', prefix, 'files/', prefix, 'profile.txt', sep = '')\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "# Creating vectors from data\n");
        fprintf(f_rprofile, "horpos <- read.table(intable, skip = 1)[, 1]  #horizontal position\n");
        fprintf(f_rprofile, "elev <- read.table(intable, skip = 1)[, 2]  #elevation\n");
        fprintf(f_rprofile, "if (depdef == 1) {\n");
        fprintf(f_rprofile, "    hdeposit <- read.table(intable, skip = 1)[, 3]  #height of observed deposit\n");
        fprintf(f_rprofile, "    hdeposit[hdeposit < 0] <- 0\n");
        fprintf(f_rprofile, "}\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "stext = vector('expression', 7)  #initializing right y axis labels\n");
        fprintf(f_rprofile, "ltext = vector('expression', 8)  #initializing legend labels\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "for (sparam in 1:3) {\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "    if (sparam == 1) {\n");
        fprintf(f_rprofile, "        if (model <= 3) {\n");
        fprintf(f_rprofile, "            smax0 <- %.5f  #maximum flow kinetic energy for display\n", tflow_maxmax );
        fprintf(f_rprofile, "            smaxf0 <- 0.0\n");
        fprintf(f_rprofile, "            smaxw0 <- 0.0\n");
        fprintf(f_rprofile, "        } else {\n");
        fprintf(f_rprofile, "            smax0 <- %.5f  #maximum PHASE 1 flow kinetic energy for display\n", tflow_maxmax1 );       
        fprintf(f_rprofile, "            smaxf0 <- %.5f  #maximum PHASE 2 flow kinetic energy for display\n", tflow_maxmax2 );
        fprintf(f_rprofile, "            smaxw0 <- %.5f  #maximum PHASE 3 flow kinetic energy for display\n", tflow_maxmax3 );
        fprintf(f_rprofile, "        }\n");
        fprintf(f_rprofile, "    } else if (sparam == 2) {\n");
        fprintf(f_rprofile, "        if (model <= 3) {\n");
        fprintf(f_rprofile, "            smax0 <- %.5f  #maximum flow pressure for display\n", pflow_maxmax );
        fprintf(f_rprofile, "            smaxf0 <- 0.0\n");
        fprintf(f_rprofile, "            smaxw0 <- 0.0\n");
        fprintf(f_rprofile, "        } else {\n");
        fprintf(f_rprofile, "            smax0 <- %.5f  #maximum (PHASE 1) flow pressure for display\n", pflow_maxmax1 );
        fprintf(f_rprofile, "            smaxf0 <- %.5f  #maximum PHASE 2 flow pressure for display\n", pflow_maxmax2 );
        fprintf(f_rprofile, "            smaxw0 <- %.5f  #maximum PHASE 3 flow pressure for display\n", pflow_maxmax3 );
        fprintf(f_rprofile, "        }\n");
        fprintf(f_rprofile, "    } else {\n");
        fprintf(f_rprofile, "        if (model <= 3) {\n");
        fprintf(f_rprofile, "            smax0 <- %.5f  #maximum flow velocity for display\n", vflow_maxmax );
        fprintf(f_rprofile, "            smaxf0 <- 0.0\n");
        fprintf(f_rprofile, "            smaxw0 <- 0.0\n");
        fprintf(f_rprofile, "        } else {\n");
        fprintf(f_rprofile, "            smax0 <- %.5f  #maximum (PHASE 1) flow velocity for display\n", vflow_maxmax1 );
        fprintf(f_rprofile, "            smaxf0 <- %.5f  #maximum PHASE 2 flow velocity for display\n", vflow_maxmax2 );
        fprintf(f_rprofile, "            smaxw0 <- %.5f  #maximum PHASE 3 flow velocity for display\n", vflow_maxmax3 );
        fprintf(f_rprofile, "        }\n");
        fprintf(f_rprofile, "    }\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "    hdisp = max(smax0, smaxf0, smaxw0) # unit conversion factor for scaling of plot\n");
        fprintf(f_rprofile, "    if ( hdisp != 0 ) {\n");
        fprintf(f_rprofile, "        mconv0 <- log10(hdisp)\n");
        fprintf(f_rprofile, "    } else {\n");
        fprintf(f_rprofile, "        mconv0 <- 3\n");
        fprintf(f_rprofile, "    }\n");
        fprintf(f_rprofile, "    if ( mconv0 < 3 ) {\n");
        fprintf(f_rprofile, "        mconv <- '1'\n");
        fprintf(f_rprofile, "    } else if ( mconv0 < 6 ) {\n");
        fprintf(f_rprofile, "        mconv <- '0.001'\n");
        fprintf(f_rprofile, "    } else if ( mconv0 < 9 ) {\n");
        fprintf(f_rprofile, "        mconv <- '0.000001'\n");
        fprintf(f_rprofile, "    } else if ( mconv0 < 12 ) {\n");
        fprintf(f_rprofile, "        mconv <- '0.000000001'\n");
        fprintf(f_rprofile, "    } else if ( mconv0 < 15 ) {\n");
        fprintf(f_rprofile, "        mconv <- '0.000000000001'\n");
        fprintf(f_rprofile, "    } else {\n");
        fprintf(f_rprofile, "        mconv <- '0.000000000000001'\n");
        fprintf(f_rprofile, "    }\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "    smax <- smax0 * as.numeric(mconv)  # converting units of maxima for second y axis\n");
        fprintf(f_rprofile, "    smaxf <- smaxf0 * as.numeric(mconv)\n");
        fprintf(f_rprofile, "    smaxw <- smaxw0 * as.numeric(mconv)\n");
        fprintf(f_rprofile, "\n");      
        fprintf(f_rprofile, "    if (model <= 3) {\n");
        fprintf(f_rprofile, "        # parameter to be displayed as bar diagram:\n");
        fprintf(f_rprofile, "        ltext[2] = substitute(expression(italic(H)))[2]\n");
        fprintf(f_rprofile, "        ltext[4] = substitute(expression(italic(H)[d]))[2]\n");
        fprintf(f_rprofile, "        if (sparam == 1) {\n");
        fprintf(f_rprofile, "            scol <- 6\n");
        fprintf(f_rprofile, "            sstring <- 'tflow'\n");
        fprintf(f_rprofile, "            if (mconv == '1')\n");
        fprintf(f_rprofile, "                munit <- 'J'\n");
        fprintf(f_rprofile, "            if (mconv == '0.001')\n");
        fprintf(f_rprofile, "                munit <- 'kJ'\n");
        fprintf(f_rprofile, "            if (mconv == '0.000001')\n");
        fprintf(f_rprofile, "                munit <- 'MJ'\n");
        fprintf(f_rprofile, "            if (mconv == '0.000000001')\n");
        fprintf(f_rprofile, "                munit <- 'GJ'\n");
        fprintf(f_rprofile, "            if (mconv == '0.000000000001')\n");
        fprintf(f_rprofile, "                munit <- 'TJ'\n");
        fprintf(f_rprofile, "            stext[1] = substitute(expression(Flow ~ kinetic ~ energy ~ italic(T) ~ (funit)),\n");
        fprintf(f_rprofile, "                list(funit = format(munit)))[2]\n");
        fprintf(f_rprofile, "            ltext[5] = substitute(expression(italic(T)))[2]\n");
        fprintf(f_rprofile, "            col <- rgb(0.6, 0.3, 0.1, 0.5)\n");
        fprintf(f_rprofile, "        } else if (sparam == 2) {\n");
        fprintf(f_rprofile, "            scol <- 7\n");
        fprintf(f_rprofile, "            sstring <- 'pflow'\n");
        fprintf(f_rprofile, "            if (mconv == '1')\n");
        fprintf(f_rprofile, "                munit <- 'Pa'\n");
        fprintf(f_rprofile, "            if (mconv == '0.001')\n");
        fprintf(f_rprofile, "                munit <- 'kPa'\n");
        fprintf(f_rprofile, "            if (mconv == '0.000001')\n");
        fprintf(f_rprofile, "                munit <- 'MPa'\n");
        fprintf(f_rprofile, "            if (mconv == '0.000000001')\n");
        fprintf(f_rprofile, "                munit <- 'GPa'\n");
        fprintf(f_rprofile, "            if (mconv == '0.000000000001')\n");
        fprintf(f_rprofile, "                munit <- 'TPa'\n");
        fprintf(f_rprofile, "            stext[1] = substitute(expression(Flow ~ pressure ~ italic(P) ~ (funit)),\n");
        fprintf(f_rprofile, "                list(funit = format(munit)))[2]\n");
        fprintf(f_rprofile, "            ltext[5] = substitute(expression(italic(P)))[2]\n");
        fprintf(f_rprofile, "            col <- rgb(0.6, 0.6, 0, 0.5)\n");
        fprintf(f_rprofile, "        } else {\n");
        fprintf(f_rprofile, "            scol <- 5\n");
        fprintf(f_rprofile, "            sstring <- 'vflow'\n");
        fprintf(f_rprofile, "            stext[1] = substitute(expression(Flow ~ velocity ~ italic(V) ~ (m/s)))[2]\n");
        fprintf(f_rprofile, "            ltext[5] = substitute(expression(italic(T)))[2]\n");
        fprintf(f_rprofile, "            col <- rgb(0, 0, 0.8, 0.5)\n");
        fprintf(f_rprofile, "        }\n");
        fprintf(f_rprofile, "    } else if (model > 3 && model < 7) {\n");
        fprintf(f_rprofile, "        ltext[1] = substitute(expression(italic(H)[P1]))[2]\n");
        fprintf(f_rprofile, "        ltext[2] = substitute(expression(italic(H)[P2]))[2]\n");
        fprintf(f_rprofile, "        ltext[5] = substitute(expression(italic(H)[P3]))[2]\n");
        fprintf(f_rprofile, "        if (sparam == 1) {\n");
        fprintf(f_rprofile, "            scol <- 8\n");
        fprintf(f_rprofile, "            sstring <- 'tflow'\n");
        fprintf(f_rprofile, "            if (mconv == '1')\n");
        fprintf(f_rprofile, "                munit <- 'J'\n");
        fprintf(f_rprofile, "            if (mconv == '0.001')\n");
        fprintf(f_rprofile, "                munit <- 'kJ'\n");
        fprintf(f_rprofile, "            if (mconv == '0.000001')\n");
        fprintf(f_rprofile, "                munit <- 'MJ'\n");
        fprintf(f_rprofile, "            if (mconv == '0.000000001')\n");
        fprintf(f_rprofile, "                munit <- 'GJ'\n");
        fprintf(f_rprofile, "            if (mconv == '0.000000000001')\n");
        fprintf(f_rprofile, "                munit <- 'TJ'\n");
        fprintf(f_rprofile, "            stext[1] = substitute(expression(Flow ~ kinetic ~ energy ~ italic(T) ~ (funit)),\n");
        fprintf(f_rprofile, "                list(funit = format(munit)))[2]\n");
        fprintf(f_rprofile, "            ltext[6] = substitute(expression(italic(T)[s]))[2]\n");
        fprintf(f_rprofile, "            ltext[7] = substitute(expression(-italic(T)[f]))[2]\n");
        fprintf(f_rprofile, "            cols <- rgb(0.6, 0.3, 0.1, 0.5)\n");
        fprintf(f_rprofile, "            colf <- rgb(0, 0.5, 0.5, 0.5)\n");
        fprintf(f_rprofile, "        } else if (sparam == 2) {\n");
        fprintf(f_rprofile, "            scol <- 10\n");
        fprintf(f_rprofile, "            sstring <- 'pflow'\n");
        fprintf(f_rprofile, "            if (mconv == '1')\n");
        fprintf(f_rprofile, "                munit <- 'Pa'\n");
        fprintf(f_rprofile, "            if (mconv == '0.001')\n");
        fprintf(f_rprofile, "                munit <- 'kPa'\n");
        fprintf(f_rprofile, "            if (mconv == '0.000001')\n");
        fprintf(f_rprofile, "                munit <- 'MPa'\n");
        fprintf(f_rprofile, "            if (mconv == '0.000000001')\n");
        fprintf(f_rprofile, "                munit <- 'GPa'\n");
        fprintf(f_rprofile, "            if (mconv == '0.000000000001')\n");
        fprintf(f_rprofile, "                munit <- 'TPa'\n");
        fprintf(f_rprofile, "            stext[1] = substitute(expression(Flow ~ pressure ~ italic(P) ~ (funit)),\n");
        fprintf(f_rprofile, "                list(funit = format(munit)))[2]\n");
        fprintf(f_rprofile, "            ltext[6] = substitute(expression(italic(P)[s]))[2]\n");
        fprintf(f_rprofile, "            ltext[7] = substitute(expression(-italic(P)[f]))[2]\n");
        fprintf(f_rprofile, "            cols <- rgb(0.6, 0.6, 0, 0.5)\n");
        fprintf(f_rprofile, "            colf <- rgb(0.2, 0, 0.6, 0.5)\n");
        fprintf(f_rprofile, "        } else {\n");
        fprintf(f_rprofile, "            scol <- 6\n");
        fprintf(f_rprofile, "            sstring <- 'vflow'\n");
        fprintf(f_rprofile, "            stext[1] = substitute(expression(Flow ~ velocity ~ italic(V) ~ (m/s)))[2]\n");
        fprintf(f_rprofile, "            ltext[6] = substitute(expression(italic(V)[s]))[2]\n");
        fprintf(f_rprofile, "            ltext[7] = substitute(expression(-italic(V)[f]))[2]\n");
        fprintf(f_rprofile, "            colf <- rgb(0, 0, 0.8, 0.5)\n");
        fprintf(f_rprofile, "            cols <- rgb(0.6, 0.4, 0, 0.5)\n");
        fprintf(f_rprofile, "        }\n");
        fprintf(f_rprofile, "    } else if (model == 7) {\n");
        fprintf(f_rprofile, "        ltext[1] = substitute(expression(italic(H)[P1]))[2]\n");
        fprintf(f_rprofile, "        ltext[2] = substitute(expression(italic(H)[P2]))[2]\n");
        fprintf(f_rprofile, "        ltext[3] = substitute(expression(-italic(H)[P3]))[2]\n");
        fprintf(f_rprofile, "        ltext[5] = substitute(expression(italic(H)[d]))[2]\n");
        fprintf(f_rprofile, "        if (sparam == 1) {\n");
        fprintf(f_rprofile, "            scol <- 10\n");
        fprintf(f_rprofile, "            sstring <- 'tflow'\n");
        fprintf(f_rprofile, "            if (mconv == '1')\n");
        fprintf(f_rprofile, "                munit <- 'J'\n");
        fprintf(f_rprofile, "            if (mconv == '0.001')\n");
        fprintf(f_rprofile, "                munit <- 'kJ'\n");
        fprintf(f_rprofile, "            if (mconv == '0.000001')\n");
        fprintf(f_rprofile, "                munit <- 'MJ'\n");
        fprintf(f_rprofile, "            if (mconv == '0.000000001')\n");
        fprintf(f_rprofile, "                munit <- 'GJ'\n");
        fprintf(f_rprofile, "            if (mconv == '0.000000000001')\n");
        fprintf(f_rprofile, "                munit <- 'TJ'\n");
        fprintf(f_rprofile, "            stext[1] = substitute(expression(Flow ~ kinetic ~ energy ~ italic(T) ~ (funit)),\n");
        fprintf(f_rprofile, "                list(funit = format(munit)))[2]\n");
        fprintf(f_rprofile, "            ltext[6] = substitute(expression(italic(T)[P1]))[2]\n");
        fprintf(f_rprofile, "            ltext[7] = substitute(expression(italic(T)[P2]))[2]\n");
        fprintf(f_rprofile, "            ltext[8] = substitute(expression(-italic(T)[P3]))[2]\n");
        fprintf(f_rprofile, "            cols <- rgb(0.8, 0, 0, 0.5)\n");
        fprintf(f_rprofile, "            colf <- rgb(0, 0.8, 0, 0.5)\n");
        fprintf(f_rprofile, "            colw <- rgb(0, 0, 0.8, 0.5)\n");
        fprintf(f_rprofile, "        } else if (sparam == 2) {\n");
        fprintf(f_rprofile, "            scol <- 13\n");
        fprintf(f_rprofile, "            sstring <- 'pflow'\n");
        fprintf(f_rprofile, "            if (mconv == '1')\n");
        fprintf(f_rprofile, "                munit <- 'Pa'\n");
        fprintf(f_rprofile, "            if (mconv == '0.001')\n");
        fprintf(f_rprofile, "                munit <- 'kPa'\n");
        fprintf(f_rprofile, "            if (mconv == '0.000001')\n");
        fprintf(f_rprofile, "                munit <- 'MPa'\n");
        fprintf(f_rprofile, "            if (mconv == '0.000000001')\n");
        fprintf(f_rprofile, "                munit <- 'GPa'\n");
        fprintf(f_rprofile, "            if (mconv == '0.000000000001')\n");
        fprintf(f_rprofile, "                munit <- 'TPa'\n");
        fprintf(f_rprofile, "            stext[1] = substitute(expression(Flow ~ pressure ~ italic(P) ~ (funit)),\n");
        fprintf(f_rprofile, "                list(funit = format(munit)))[2]\n");
        fprintf(f_rprofile, "            ltext[6] = substitute(expression(italic(P)[P1]))[2]\n");
        fprintf(f_rprofile, "            ltext[7] = substitute(expression(italic(P)[P2]))[2]\n");
        fprintf(f_rprofile, "            ltext[8] = substitute(expression(-italic(P)[P3]))[2]\n");
        fprintf(f_rprofile, "            cols <- rgb(0.8, 0, 0, 0.5)\n");
        fprintf(f_rprofile, "            colf <- rgb(0, 0.8, 0, 0.5)\n");
        fprintf(f_rprofile, "            colw <- rgb(0, 0, 0.8, 0.5)\n");
        fprintf(f_rprofile, "        } else {\n");
        fprintf(f_rprofile, "            scol <- 7\n");
        fprintf(f_rprofile, "            sstring <- 'vflow'\n");
        fprintf(f_rprofile, "            stext[1] = substitute(expression(Flow ~ velocity ~ italic(V) ~ (m/s)))[2]\n");
        fprintf(f_rprofile, "            ltext[6] = substitute(expression(italic(V)[P1]))[2]\n");
        fprintf(f_rprofile, "            ltext[7] = substitute(expression(italic(V)[P2]))[2]\n");
        fprintf(f_rprofile, "            ltext[8] = substitute(expression(-italic(V)[P3]))[2]\n");
        fprintf(f_rprofile, "            cols <- rgb(0.8, 0, 0, 0.5)\n");
        fprintf(f_rprofile, "            colf <- rgb(0, 0.8, 0, 0.5)\n");
        fprintf(f_rprofile, "            colw <- rgb(0, 0, 0.8, 0.5)\n");
        fprintf(f_rprofile, "        }\n");
        fprintf(f_rprofile, "    }\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "    for (ntimesteps in 0:maxtimesteps) {\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "        if (ntimesteps<10) {\n");       
        fprintf(f_rprofile, "            fill = paste('000', ntimesteps, sep='')\n");      
        fprintf(f_rprofile, "        } else if (ntimesteps<100) {\n");       
        fprintf(f_rprofile, "            fill = paste('00', ntimesteps, sep='')\n");        
        fprintf(f_rprofile, "        } else if (ntimesteps<1000) {\n");       
        fprintf(f_rprofile, "            fill = paste('0', ntimesteps, sep='')\n");        
        fprintf(f_rprofile, "        } else {\n");       
        fprintf(f_rprofile, "            fill = paste('', ntimesteps, sep='')\n");
        fprintf(f_rprofile, "        }\n");
        fprintf(f_rprofile, "\n");        
        fprintf(f_rprofile, "        if (model <= 3) {\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "            hrelease <- read.table(intable, skip = 1)[, 3 + depdef]  #release height\n");
        fprintf(f_rprofile, "            hflow <- read.table(intable, skip = 1)[, 5 * ntimesteps + 3 + depdef]  #flow height\n");
        fprintf(f_rprofile, "            sflow <- read.table(intable, skip = 1)[, 5 * ntimesteps + scol - 1 + depdef]  #flow parameter\n");
        fprintf(f_rprofile, "            hentr <- read.table(intable, skip = 1)[, 5 * ntimesteps + 7 + depdef]  #entrained/deposited depth\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "        } else if (model > 3 && model < 7) {\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "            hreleases <- read.table(intable, skip = 1)[, 3 + depdef]  #PHASE 1 release height\n");
        fprintf(f_rprofile, "            hreleasef <- read.table(intable, skip = 1)[, 4 + depdef]  #fluid release height\n");
        fprintf(f_rprofile, "            hflows <- read.table(intable, skip = 1)[, 10 * ntimesteps + 3 + depdef]  #PHASE 1 flow height\n");
        fprintf(f_rprofile, "            hflowf <- read.table(intable, skip = 1)[, 10 * ntimesteps + 4 + depdef]  #fluid flow height\n");
        fprintf(f_rprofile, "            sflows <- read.table(intable, skip = 1)[, 10 * ntimesteps + scol - 1 + depdef]  #PHASE 1 flow parameter\n");
        fprintf(f_rprofile, "            sflowf <- read.table(intable, skip = 1)[, 10 * ntimesteps + scol + depdef]  #fluid flow parameter\n");
        fprintf(f_rprofile, "            hentrs <- read.table(intable, skip = 1)[, 10 * ntimesteps + 11 + depdef]  #PHASE 1 entrained/deposited depth\n");
        fprintf(f_rprofile, "            hentrf <- read.table(intable, skip = 1)[, 10 * ntimesteps + 12 + depdef]  #fluid entrained/deposited depth\n");
        fprintf(f_rprofile, "        } else if (model == 7) {\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "            hreleases <- read.table(intable, skip = 1)[, 3 + depdef]  #PHASE 1 release height\n");
        fprintf(f_rprofile, "            hreleasef <- read.table(intable, skip = 1)[, 4 + depdef]  #PHASE 2 release height\n");
        fprintf(f_rprofile, "            hreleasew <- read.table(intable, skip = 1)[, 5 + depdef]  #PHASE 3 release height\n");
        fprintf(f_rprofile, "            hflows <- read.table(intable, skip = 1)[, 15 * ntimesteps + 3 + depdef]  #PHASE 1 flow height\n");
        fprintf(f_rprofile, "            hflowf <- read.table(intable, skip = 1)[, 15 * ntimesteps + 4 + depdef]  #PHASE 2 flow height\n");
        fprintf(f_rprofile, "            hfloww <- read.table(intable, skip = 1)[, 15 * ntimesteps + 5 + depdef]  #PHASE 3 flow height\n");
        fprintf(f_rprofile, "            sflows <- read.table(intable, skip = 1)[, 15 * ntimesteps + scol - 1 + depdef]  #PHASE 1 flow parameter\n");
        fprintf(f_rprofile, "            sflowf <- read.table(intable, skip = 1)[, 15 * ntimesteps + scol + depdef]  #PHASE 2 flow parameter\n");
        fprintf(f_rprofile, "            sfloww <- read.table(intable, skip = 1)[, 15 * ntimesteps + scol + 1 + depdef]  #PHASE 3 flow parameter\n");
        fprintf(f_rprofile, "            hentrs <- read.table(intable, skip = 1)[, 15 * ntimesteps + 15 + depdef]  #PHASE 1 entrained/deposited depth\n");
        fprintf(f_rprofile, "            hentrf <- read.table(intable, skip = 1)[, 15 * ntimesteps + 16 + depdef]  #PHASE 2 entrained/deposited depth\n");
        fprintf(f_rprofile, "            hentrw <- read.table(intable, skip = 1)[, 15 * ntimesteps + 17 + depdef]  #PHASE 3 entrained/deposited depth\n");
        fprintf(f_rprofile, "        }\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "        # Building geometry of profile plot\n");
        fprintf(f_rprofile, "        maxhorpos <- max(horpos, na.rm = TRUE)  #maximum of horizontal coordinate\n");
        fprintf(f_rprofile, "        minhorpos <- min(horpos, na.rm = TRUE)  #minimum of horizontal coordinate\n");
        fprintf(f_rprofile, "        topelev <- elev + hexagg * hmax  #maximum elevation of each pixel\n");
        fprintf(f_rprofile, "        maxelev <- max(topelev, na.rm = TRUE)  #maximum elevation over entire area\n");
        fprintf(f_rprofile, "        minelev <- min(elev, na.rm = TRUE)  #minimum elevation over entire area\n");
        fprintf(f_rprofile, "        maxelev <- maxelev + (maxelev - minelev) * 0.1  #updating maximum elevation over entire area (for legend)\n");
        fprintf(f_rprofile, "        minelev <- minelev - (maxelev - minelev) * 0.1  #updating maximum elevation over entire area (for texts)\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "        # Creating plot file\n");
        fprintf(f_rprofile, "        profileplot = paste(wkdir, prefix, 'results/', prefix, 'plots/', prefix, 'profiles_timesteps/',\n");
        fprintf(f_rprofile, "            prefix, sstring, fill, '.png', sep = '')\n");
        fprintf(f_rprofile, "        png(filename = profileplot, width = 15, height = 10, units = 'cm', res = 300)\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "        # Defining margins\n");
        fprintf(f_rprofile, "        par(mar = c(3.2, 3.2, 1, 3.2))\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "        # Plotting bars for all time steps except initial condition:\n");
        fprintf(f_rprofile, "        if (ntimesteps > 0) {\n");
        fprintf(f_rprofile, "            if (model <= 3) {\n");
        fprintf(f_rprofile, "                barplot(sflow * as.numeric(mconv), axes = F, xlab = NA, ylab = NA, ylim = c(-smax * 1.4, \n");
        fprintf(f_rprofile, "                    smax * 1.4), border = NA, col = col) # flow parameter\n");
        fprintf(f_rprofile, "            } else if (model > 3 && model < 7) {\n");
        fprintf(f_rprofile, "                barplot(sflows * as.numeric(mconv), axes = F, xlab = NA, ylab = NA, ylim = c(min(-smax * 1.4, \n");
        fprintf(f_rprofile, "                    -smaxf * 1.4), max(smax * 1.4, smaxf * 1.4)), border = NA, col = cols)  #PHASE 1 flow parameter\n");
        fprintf(f_rprofile, "                par(new = TRUE)\n");
        fprintf(f_rprofile, "                barplot(-sflowf * as.numeric(mconv), axes = F, xlab = NA, ylab = NA, ylim = c(min(-smax * 1.4, \n");
        fprintf(f_rprofile, "                    -smaxf * 1.4), max(smax * 1.4, smaxf * 1.4)), border = NA, col = colf)  #PHASE 2 flow parameter\n");
        fprintf(f_rprofile, "            } else if (model == 7) {\n");
        fprintf(f_rprofile, "                barplot((sflows + sflowf) * as.numeric(mconv), axes = F, xlab = NA, ylab = NA,\n");
        fprintf(f_rprofile, "                    ylim = c(min(-smax * 1.4 - smaxf * 1.4, -smaxw * 1.4), max(smax * 1.4 + smaxf * 1.4, smaxw * 1.4)), \n");
        fprintf(f_rprofile, "                    border = NA, col = cols)  #PHASE 1 flow parameter\n");
        fprintf(f_rprofile, "                par(new = TRUE)\n");
        fprintf(f_rprofile, "                barplot(sflowf * as.numeric(mconv), axes = F, xlab = NA, ylab = NA, ylim = c(min(-smax * 1.4 - smaxf * 1.4, \n");
        fprintf(f_rprofile, "                    -smaxw * 1.4), max(smax * 1.4 + smaxf * 1.4, smaxw * 1.4)), border = NA, col = colf)  #PHASE 2 flow parameter\n");
        fprintf(f_rprofile, "                par(new = TRUE)\n");
        fprintf(f_rprofile, "                barplot(-sfloww * as.numeric(mconv), axes = F, xlab = NA, ylab = NA, ylim = c(min(-smax * 1.4 - smaxf * 1.4, \n");
        fprintf(f_rprofile, "                    -smaxw * 1.4), max(smax * 1.4 + smaxf * 1.4, smaxw * 1.4)), border = NA, col = colw)  #PHASE 3 flow parameter\n");
        fprintf(f_rprofile, "            }\n");
        fprintf(f_rprofile, "            axis(side = 4, tck = -0.02, labels = NA)  #y axis for velocity\n");
        fprintf(f_rprofile, "            axis(side = 4, lwd = 0, line = -0.4)\n");
        fprintf(f_rprofile, "            par(new = TRUE)\n");
        fprintf(f_rprofile, "        }\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "        # Plotting lines\n");
        fprintf(f_rprofile, "        plot(x = horpos, y = elev, axes = F, xlab = NA, ylab = NA, type = 'l', ylim = c(minelev,\n");
        fprintf(f_rprofile, "            maxelev), lty = 1, lwd = 1, col = 'gray')  #elevation\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "        if (depdef == 1) lines(x = horpos, y = elev + hexagg * hdeposit, lty = 3, lwd = 1,\n");
        fprintf(f_rprofile, "            col = 'black')  #height of deposit\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "        if (model <= 3) {\n");
        fprintf(f_rprofile, "            lines(x = horpos, y = elev + hexagg * hrelease, lty = 3, lwd = 1, col = 'red')  #initial flow height\n");
        fprintf(f_rprofile, "            lines(x = horpos, y = elev + hexagg * hentr, lty = 1, lwd = 1, col = 'black')  #entrained/deposited depth\n");
        fprintf(f_rprofile, "            lines(x = horpos[2:(length(horpos) - 1)], y = elev[2:(length(horpos) - 1)] +\n");
        fprintf(f_rprofile, "                hexagg * (hflow[2:(length(horpos) - 1)] + hentr[2:(length(horpos) - 1)]) +\n");
        fprintf(f_rprofile, "                0 * (1/(hflow[2:(length(horpos) - 1)] + hflow[1:(length(horpos) - 2)] + hflow[3:length(horpos)])),\n");
        fprintf(f_rprofile, "                lty = 1, lwd = 2, col = 'red')  #flow height\n");
        fprintf(f_rprofile, "        } else if (model > 3 && model < 7) {\n");
        fprintf(f_rprofile, "            lines(x = horpos, y = elev + hexagg * (hreleases + hreleasef), lty = 3, lwd = 1,\n");
        fprintf(f_rprofile, "                col = rgb(0, 0, 0.8))  #initial PHASE 2 flow height\n");
        fprintf(f_rprofile, "            lines(x = horpos, y = elev + hexagg * hreleases, lty = 3, lwd = 1, col = rgb(0.45,\n");
        fprintf(f_rprofile, "                0.3, 0))  #initial PHASE 1 flow height\n");
        fprintf(f_rprofile, "            lines(x = horpos, y = elev + hexagg * (hentrs + hentrf), lty = 1, lwd = 1, col = 'black')  #entrained/deposited PHASE 1 depth\n");
        fprintf(f_rprofile, "            lines(x = horpos, y = elev + hexagg * hentrf, lty = 1, lwd = 1, col = 'black')  #entrained/deposited PHASE 2 depth\n");
        fprintf(f_rprofile, "            lines(x = horpos[2:(length(horpos) - 1)], y = elev[2:(length(horpos) - 1)] +\n");
        fprintf(f_rprofile, "                hexagg * (hflows[2:(length(horpos) - 1)] + hflowf[2:(length(horpos) - 1)] +\n");
        fprintf(f_rprofile, "                    hentrs[2:(length(horpos) - 1)] + hentrf[2:(length(horpos) - 1)]) + 0 *\n");
        fprintf(f_rprofile, "                (1/(hflowf[2:(length(horpos) - 1)] + hflowf[1:(length(horpos) - 2)] + hflowf[3:length(horpos)])),\n");
        fprintf(f_rprofile, "                lty = 1, lwd = 2, col = rgb(0, 0, 0.8))  #PHASE 2 flow height\n");
        fprintf(f_rprofile, "            lines(x = horpos[2:(length(horpos) - 1)], y = elev[2:(length(horpos) - 1)] +\n");
        fprintf(f_rprofile, "                hexagg * (hflows[2:(length(horpos) - 1)] + hentrs[2:(length(horpos) - 1)] +\n");
        fprintf(f_rprofile, "                    hentrf[2:(length(horpos) - 1)]) + 0 * (1/(hflows[2:(length(horpos) -\n");
        fprintf(f_rprofile, "                1)] + hflows[1:(length(horpos) - 2)] + hflows[3:length(horpos)])), lty = 1,\n");
        fprintf(f_rprofile, "                lwd = 2, col = rgb(0.45, 0.3, 0))  #PHASE 1 flow height\n");
        fprintf(f_rprofile, "        } else if (model == 7) {\n");
        fprintf(f_rprofile, "            lines(x = horpos, y = elev + hexagg * (hreleases + hreleasef + hreleasew), lty = 3,\n");
        fprintf(f_rprofile, "                lwd = 1, col = rgb(0, 0, 0.8))  #initial PHASE 3 flow height\n");
        fprintf(f_rprofile, "            lines(x = horpos, y = elev + hexagg * (hreleases + hreleasef), lty = 3, lwd = 1,\n");
        fprintf(f_rprofile, "                col = rgb(0, 0.8, 0))  #initial PHASE 2 flow height\n");
        fprintf(f_rprofile, "            lines(x = horpos, y = elev + hexagg * hreleases, lty = 3, lwd = 1, col = rgb(0.8,\n");
        fprintf(f_rprofile, "                0, 0))  #initial PHASE 1 flow height\n");
        fprintf(f_rprofile, "            lines(x = horpos, y = elev + hexagg * (hentrs + hentrf + hentrw), lty = 1, lwd = 1,\n");
        fprintf(f_rprofile, "                col = 'black')  #entrained/deposited PHASE 3 depth\n");
        fprintf(f_rprofile, "            lines(x = horpos, y = elev + hexagg * (hentrs + hentrf), lty = 1, lwd = 1, col = 'black')  #entrained/deposited PHASE 1 depth\n");
        fprintf(f_rprofile, "            lines(x = horpos, y = elev + hexagg * hentrf, lty = 1, lwd = 1, col = 'black')  #entrained/deposited PHASE 2 depth\n");
        fprintf(f_rprofile, "            lines(x = horpos[2:(length(horpos) - 1)], y = elev[2:(length(horpos) - 1)] +\n");
        fprintf(f_rprofile, "                hexagg * (hflows[2:(length(horpos) - 1)] + hflowf[2:(length(horpos) - 1)] +\n");
        fprintf(f_rprofile, "                    hfloww[2:(length(horpos) - 1)] + hentrs[2:(length(horpos) - 1)] + hentrf[2:(length(horpos) -\n");
        fprintf(f_rprofile, "                    1)] + hentrw[2:(length(horpos) - 1)]) + 0 * (1/(hfloww[2:(length(horpos) -\n");
        fprintf(f_rprofile, "                1)] + hfloww[1:(length(horpos) - 2)] + hfloww[3:length(horpos)])), lty = 1,\n");
        fprintf(f_rprofile, "                lwd = 2, col = rgb(0, 0, 0.8))  #PHASE 3 flow height\n");
        fprintf(f_rprofile, "            lines(x = horpos[2:(length(horpos) - 1)], y = elev[2:(length(horpos) - 1)] +\n");
        fprintf(f_rprofile, "                hexagg * (hflows[2:(length(horpos) - 1)] + hflowf[2:(length(horpos) - 1)] +\n");
        fprintf(f_rprofile, "                    hentrs[2:(length(horpos) - 1)] + hentrf[2:(length(horpos) - 1)] + hentrw[2:(length(horpos) -\n");
        fprintf(f_rprofile, "                    1)]) + 0 * (1/(hflowf[2:(length(horpos) - 1)] + hflowf[1:(length(horpos) -\n");
        fprintf(f_rprofile, "                2)] + hflowf[3:length(horpos)])), lty = 1, lwd = 2, col = rgb(0, 0.8, 0))  #PHASE 2 flow height\n");
        fprintf(f_rprofile, "            lines(x = horpos[2:(length(horpos) - 1)], y = elev[2:(length(horpos) - 1)] +\n");
        fprintf(f_rprofile, "                hexagg * (hflows[2:(length(horpos) - 1)] + hentrs[2:(length(horpos) - 1)] +\n");
        fprintf(f_rprofile, "                    hentrf[2:(length(horpos) - 1)] + hentrw[2:(length(horpos) - 1)]) + 0 *\n");
        fprintf(f_rprofile, "                (1/(hflows[2:(length(horpos) - 1)] + hflows[1:(length(horpos) - 2)] + hflows[3:length(horpos)])),\n");
        fprintf(f_rprofile, "                lty = 1, lwd = 2, col = rgb(0.8, 0, 0))  #PHASE 1 flow height\n");
        fprintf(f_rprofile, "        }\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "        # Legends\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "        if (model <= 3) {\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "            legend('topright', legend = ltext[2], lty = 1, lwd = 2, col = 'red', text.col = 'red',\n");
        fprintf(f_rprofile, "                bty = 'n', inset = c(0.825, -0.025), horiz = TRUE, x.intersp = 0.5)  #flow height\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "            if (depdef == 1)\n");
        fprintf(f_rprofile, "                legend('topright', legend = ltext[4], lty = 3, lwd = 2, col = 'black', text.col = 'black',\n");
        fprintf(f_rprofile, "                    bty = 'n', inset = c(0.425, -0.025), horiz = TRUE, x.intersp = 0.5)  #height of deposit\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "            if (ntimesteps > 0) {\n");
        fprintf(f_rprofile, "                legend('topright', legend = ltext[5], border = NA, fill = col, text.col = col,\n");
        fprintf(f_rprofile, "                    bty = 'n', inset = c(0.225, -0.025), horiz = TRUE)  #bar plot parameter\n");
        fprintf(f_rprofile, "            }\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "        } else if (model > 3 && model < 7) {\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "            legend('topright', legend = ltext[1], lty = 1, lwd = 2, col = rgb(0.45, 0.3,\n");
        fprintf(f_rprofile, "                0), text.col = rgb(0.45, 0.3, 0), bty = 'n', inset = c(0.825, -0.025), horiz = TRUE,\n");
        fprintf(f_rprofile, "                x.intersp = 0.5)  #PHASE 1 flow height\n");
        fprintf(f_rprofile, "            legend('topright', legend = ltext[2], lty = 1, lwd = 2, col = rgb(0, 0, 0.8),\n");
        fprintf(f_rprofile, "                text.col = rgb(0, 0, 0.8), bty = 'n', inset = c(0.625, -0.025), horiz = TRUE,\n");
        fprintf(f_rprofile, "                x.intersp = 0.5)  #PHASE 2 flow height\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "            if (depdef == 1)\n");
        fprintf(f_rprofile, "                legend('topright', legend = ltext[5], lty = 3, lwd = 2, col = 'black', text.col = 'black',\n");
        fprintf(f_rprofile, "                    bty = 'n', inset = c(0.425, -0.025), horiz = TRUE, x.intersp = 0.5)  #height of deposit\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "            if (ntimesteps > 0) {\n");
        fprintf(f_rprofile, "                legend('topright', legend = ltext[6], border = NA, fill = cols, text.col = cols,\n");
        fprintf(f_rprofile, "                    bty = 'n', inset = c(0.225, -0.025), horiz = TRUE)  #PHASE 1 bar plot parameter\n");
        fprintf(f_rprofile, "                legend('topright', legend = ltext[7], border = NA, fill = colf, text.col = colf,\n");
        fprintf(f_rprofile, "                    bty = 'n', inset = c(0.025, -0.025), horiz = TRUE)  #PHASE 2 bar plot parameter\n");
        fprintf(f_rprofile, "            }\n");
        fprintf(f_rprofile, "        } else if (model == 7) {\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "            par(cex = 0.75)\n");
        fprintf(f_rprofile, "            legend('topright', legend = ltext[1], lty = 1, lwd = 2, col = rgb(0.8, 0, 0),\n");
        fprintf(f_rprofile, "                text.col = rgb(0.8, 0, 0), bty = 'n', inset = c(0.85, -0.025), horiz = TRUE,\n");
        fprintf(f_rprofile, "                x.intersp = 0.5)  #PHASE 1 flow height\n");
        fprintf(f_rprofile, "            legend('topright', legend = ltext[2], lty = 1, lwd = 2, col = rgb(0, 0.8, 0),\n");
        fprintf(f_rprofile, "                text.col = rgb(0, 0.8, 0), bty = 'n', inset = c(0.72, -0.025), horiz = TRUE,\n");
        fprintf(f_rprofile, "                x.intersp = 0.5)  #PHASE 2 flow height\n");
        fprintf(f_rprofile, "            legend('topright', legend = ltext[3], lty = 1, lwd = 2, col = rgb(0, 0, 0.8),\n");
        fprintf(f_rprofile, "                text.col = rgb(0, 0, 0.8), bty = 'n', inset = c(0.575, -0.025), horiz = TRUE,\n");
        fprintf(f_rprofile, "                x.intersp = 0.5)  #PHASE 3 flow height\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "            if (depdef == 1)\n");
        fprintf(f_rprofile, "                legend('topright', legend = ltext[5], lty = 3, lwd = 2, col = 'black', text.col = 'black',\n");
        fprintf(f_rprofile, "                    bty = 'n', inset = c(0.425, -0.025), horiz = TRUE, x.intersp = 0.5)  #height of deposit\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "            if (ntimesteps > 0) {\n");
        fprintf(f_rprofile, "                legend('topright', legend = ltext[6], border = NA, fill = cols, text.col = cols,\n");
        fprintf(f_rprofile, "                    bty = 'n', inset = c(0.275, -0.025), horiz = TRUE)  #PHASE 1 bar plot parameter\n");
        fprintf(f_rprofile, "                legend('topright', legend = ltext[7], border = NA, fill = colf, text.col = colf,\n");
        fprintf(f_rprofile, "                    bty = 'n', inset = c(0.165, -0.025), horiz = TRUE)  #PHASE 2 bar plot parameter\n");
        fprintf(f_rprofile, "                legend('topright', legend = ltext[8], border = NA, fill = colw, text.col = colw,\n");
        fprintf(f_rprofile, "                    bty = 'n', inset = c(0.025, -0.025), horiz = TRUE)  #PHASE 3 bar plot parameter\n");
        fprintf(f_rprofile, "            }\n");
        fprintf(f_rprofile, "            par(cex = 1)\n");
        fprintf(f_rprofile, "        }\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "        # Plotting bounding box and axes\n");
        fprintf(f_rprofile, "        box()\n");
        fprintf(f_rprofile, "        axis(side = 1, tck = -0.02, labels = NA)  #x axis\n");
        fprintf(f_rprofile, "        axis(side = 1, lwd = 0, line = -0.4)\n");
        fprintf(f_rprofile, "        axis(side = 2, tck = -0.02, labels = NA)  #y axis for lines\n");
        fprintf(f_rprofile, "        axis(side = 2, lwd = 0, line = -0.4)\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "        # Plotting axis labels\n");
        fprintf(f_rprofile, "        mtext(side = 1, 'Horizontal distance (m)', line = 2)  #x axis\n");
        fprintf(f_rprofile, "        mtext(side = 2, 'Elevation (m)', line = 2)  #y axis for lines\n");
        fprintf(f_rprofile, "        if (ntimesteps > 0) mtext(side = 4, stext[1], line = 2)  #y axis for bars\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "        # Texts showing time passed since start of flow and exaggeration of flow height\n");
        fprintf(f_rprofile, "        secpassed <- tint * ntimesteps  #time passed since start of flow\n");
        fprintf(f_rprofile, "        if (secpassed > tstop) secpassed <- tstop\n");
        fprintf(f_rprofile, "        textpassed = vector('expression', 1)  #initializing text label\n");
        fprintf(f_rprofile, "        textpassed[1] = substitute(expression(italic(t) == secformat ~ s), list(secformat = format(round(secpassed,\n");
        fprintf(f_rprofile, "            1), nsmall = 1)))[2]  #defining text label\n");
        fprintf(f_rprofile, "        text(x = minhorpos + (maxhorpos - minhorpos)/6, y = minelev + (maxelev - minelev)/25,\n");
        fprintf(f_rprofile, "            labels = textpassed[1], cex = 1.4, col = 'black', font = 2)  #printing text for time passed since start of flow\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "        if (hexagg != 1) {\n");
        fprintf(f_rprofile, "            texthexagg <- paste(hexagg, '-fold exaggeration of flow height', sep = '')  #text label for exaggeration of flow height\n");
        fprintf(f_rprofile, "            text(x = minhorpos + (maxhorpos - minhorpos) * 4/6, y = minelev + (maxelev -\n");
        fprintf(f_rprofile, "                minelev)/35, labels = texthexagg, col = 'black')  #printing text for exaggeration of flow height\n");
        fprintf(f_rprofile, "        }\n");
        fprintf(f_rprofile, "\n");
        fprintf(f_rprofile, "        # Closing plot file\n");
        fprintf(f_rprofile, "        dev.off()\n");
        fprintf(f_rprofile, "    }\n");
        fprintf(f_rprofile, "}\n");
        fprintf(f_rprofile, "cat('completed.\\n')\n"); 

        fclose(f_rprofile);

        if ( strcmp( rscript, "None" ) != 0 ) {
        
            sprintf(path, "%sr.avaflow.profile.cmd", outrplots); // cmd script to start R script
            f_cprofile=fopen(path, "w\n");

            fprintf(f_cprofile, "\"%s\" \"%s/r.avaflow.profile.R\"\n", rscript, outrplots );       
        
            fclose(f_cprofile);
        }
    }


// -- STOP --- Writing R script for profile plots ---------------------------------------------------------------


// -- START -- Writing R script for map plots -------------------------------------------------------------------


    if ( sico.MULT == 0 || xint == 1 ) {

        sprintf(path, "%sr.avaflow.map.R", outrplots); // R script for map plots
        f_rmap=fopen(path, "w\n");

        fprintf(f_rmap, "#!/usr/bin/R\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "##############################################################################\n");
        fprintf(f_rmap, "#\n");
        fprintf(f_rmap, "# MODULE:       r.avaflow.map.R\n");
        fprintf(f_rmap, "# AUTHOR:       Martin Mergili\n");
        fprintf(f_rmap, "#\n");
        fprintf(f_rmap, "# PURPOSE:      The simulation model for avalanche and debris flows\n");
        fprintf(f_rmap, "#               Script for the creation of map plots of simulation results\n");
        fprintf(f_rmap, "#\n");
        fprintf(f_rmap, "# COPYRIGHT:    (c) 2013 - 2024 by the author\n");
        fprintf(f_rmap, "#               (c) 2020 - 2024 by the University of Graz\n");
        fprintf(f_rmap, "#               (c) 2013 - 2021 by the BOKU University, Vienna\n");
        fprintf(f_rmap, "#               (c) 2015 - 2020 by the University of Vienna\n");
        fprintf(f_rmap, "#               (c) 1993 - 2024 by the R Development Core Team\n");
        fprintf(f_rmap, "#\n");
        fprintf(f_rmap, "#               This program is free software under the GNU General Public\n");
        fprintf(f_rmap, "#               License (>=v2).\n");
        fprintf(f_rmap, "#\n");
        fprintf(f_rmap, "##############################################################################\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "# Loading libraries\n");
        fprintf(f_rmap, "options(warn = -1)\n");
        fprintf(f_rmap, "suppressPackageStartupMessages(library('sp'))\n");
        fprintf(f_rmap, "suppressPackageStartupMessages(library('terra'))\n");
        fprintf(f_rmap, "suppressPackageStartupMessages(library('raster'))\n");
        
        if ( strcmp( rlibs, "None" ) ==  0 ) {

            fprintf(f_rmap, "library(sp, quietly = T)\n");
            fprintf(f_rmap, "library(terra, quietly = T)\n");
            fprintf(f_rmap, "library(raster, quietly = T)\n");
            
        } else {

            fprintf(f_rmap, "library(sp, lib.loc = '%s', quietly = T)\n", rlibs );
            fprintf(f_rmap, "library(terra, lib.loc = '%s', quietly = T)\n", rlibs );
            fprintf(f_rmap, "library(raster, lib.loc = '%s', quietly = T)\n", rlibs );       
        }
        
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "# Defining arguments\n");

    #ifdef WITHGRASS

        fprintf(f_rmap, "wkdir <- ''  #working directory\n" );

    #else

        fprintf(f_rmap, "wkdir <- '%s/'  #working directory\n", wkdir );

    #endif

        fprintf(f_rmap, "prefix <- '%s'  #prefix for file names\n", prefix );
        fprintf(f_rmap, "aflag <- %i  #control for additional output\n", sico.AFLAG );
        fprintf(f_rmap, "tflag <- %i  #control for tsunami\n", sico.TSUNAMI );
        fprintf(f_rmap, "model <- %i  #model\n", sico.MODEL );
        fprintf(f_rmap, "ymax0 <- %.6f  #coordinates defining map boundaries\n", sico.BDNORTH );
        fprintf(f_rmap, "ymin0 <- %.6f\n", sico.BDSOUTH );
        fprintf(f_rmap, "xmin0 <- %.6f\n", sico.BDWEST );
        fprintf(f_rmap, "xmax0 <- %.6f\n", sico.BDWEST + sico.N * sico.CSZ );
        fprintf(f_rmap, "cellsize <- %.6f\n", sico.CSZ );
        fprintf(f_rmap, "tint <- %.6f\n", tout );
        fprintf(f_rmap, "tstop <- %.6f\n", tmax );
        fprintf(f_rmap, "impdef <- %i\n", sico.IMPACTAREA );
        fprintf(f_rmap, "depdef <- %i\n", sico.HDEPOSIT );
        fprintf(f_rmap, "thresholdsh <- %.6f\n", sico.IMPTHR[0] );
        fprintf(f_rmap, "thresholdst <- %.6f\n", sico.IMPTHR[1] );
        fprintf(f_rmap, "thresholdsp <- %.6f\n", sico.IMPTHR[2] );        
        fprintf(f_rmap, "maxhflow <- %.6f\n", hflow_maxmax );
        fprintf(f_rmap, "maxtflow <- %.6f\n", tflow_maxmax );
        fprintf(f_rmap, "maxpflow <- %.6f\n", pflow_maxmax );
        fprintf(f_rmap, "maxvflow <- %.6f\n", vflow_maxmax );
        fprintf(f_rmap, "maxtreach <- %.6f\n", treachmaxmax );
        fprintf(f_rmap, "maxhtsun <- %.6f\n", htsunmaxmax );
        fprintf(f_rmap, "maxhentr <- %.6f\n", basechange_min );
        fprintf(f_rmap, "maxhdep <- %.6f\n", basechange_max );
        fprintf(f_rmap, "mapprof <- %i\n", sico.PROFILE );
        fprintf(f_rmap, "hydrograph <- %i\n", hydrograph );
        fprintf(f_rmap, "releasemass <- %i\n", ctrl_release );
        fprintf(f_rmap, "ninhyd <- %i\n", hydnin );
        fprintf(f_rmap, "nouthyd <- %i\n", hydnout );
        fprintf(f_rmap, "ntimemax <- %i\n", nout - 1 );
        fprintf(f_rmap, "ctrlpts <- %i\n", sico.CTRLPOINTS );
        fprintf(f_rmap, "ctrl_basechange <- %i\n", ctrl_basechange );
        
    #ifdef WITHGRASS

        if ( strcmp( pbg1name, "None" ) != 0 ) {

            fprintf(f_rmap, "ortho1 <- '%sresults/%sascii/%s.asc'  #path to orthophoto red channel\n", prefix, prefix, pbg1name );
            fprintf(f_rmap, "ortho2 <- '%sresults/%sascii/%s.asc'  #path to orthophoto green channel\n", prefix, prefix, pbg2name );
            fprintf(f_rmap, "ortho3 <- '%sresults/%sascii/%s.asc'  #path to orthophoto blue channel\n", prefix, prefix, pbg3name );
            
        } else {

            fprintf(f_rmap, "ortho1 <- 'None'\n" );
            fprintf(f_rmap, "ortho2 <- 'None'\n" );
            fprintf(f_rmap, "ortho3 <- 'None'\n" );
        }

    #else

        if ( strcmp( pbg1name, "None" ) != 0 ) {

            fprintf(f_rmap, "ortho1 <- '%s/%s.asc'  #path to orthophoto red channel\n", wkdir, pbg1name );
            fprintf(f_rmap, "ortho2 <- '%s/%s.asc'  #path to orthophoto green channel\n", wkdir, pbg2name );
            fprintf(f_rmap, "ortho3 <- '%s/%s.asc'  #path to orthophoto blue channel\n", wkdir, pbg3name );
            
        } else {

            fprintf(f_rmap, "ortho1 <- 'None'\n" );
            fprintf(f_rmap, "ortho2 <- 'None'\n" );
            fprintf(f_rmap, "ortho3 <- 'None'\n" );
        }

    #endif
        
        fprintf(f_rmap, "nphases <- %i  #number of phases\n", sico.PMAX );
        fprintf(f_rmap, "slomo <- %.6f  #factor for slow motion\n", sico.SLOMO );
        fprintf(f_rmap, "mult <- %i  #control for multiple simulations\n", sico.MULT );
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "if ( maxhentr == 0 ) {\n");
        fprintf(f_rmap, "  maxhentr <- maxhdep * -1\n");
        fprintf(f_rmap, "} else if ( maxhdep == 0 ) {\n");
        fprintf(f_rmap, "  maxhdep = maxhentr * -1\n");
        fprintf(f_rmap, "}\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "# Preparing vectors for looping over different sets of maps\n");
        fprintf(f_rmap, "if ( mult == 0 ) {\n");
        fprintf(f_rmap, "  jrange <- c(0)\n");
        fprintf(f_rmap, "  mstringlist <- c('hflow', 'None', 'None', 'None', 'None', 'None', 'None')\n");
        fprintf(f_rmap, "  if ( aflag == 1 ) {\n");
        fprintf(f_rmap, "      jrange <- append(jrange, 1)\n");
        fprintf(f_rmap, "      jrange <- append(jrange, 2)\n");
        fprintf(f_rmap, "      mstringlist[2] <- 'tflow'\n");
        fprintf(f_rmap, "      mstringlist[3] <- 'pflow'\n");
        fprintf(f_rmap, "  }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "  if ( ctrl_basechange == 1 ) {\n");
        fprintf(f_rmap, "    jrange <- append(jrange, 3)\n");
        fprintf(f_rmap, "    mstringlist[4] <- 'basechange'\n");
        fprintf(f_rmap, "  }\n");
        fprintf(f_rmap, "  if ( model == 7 && tflag == 1 ) {\n");
        fprintf(f_rmap, "    jrange <- append(jrange, 5)\n");
        fprintf(f_rmap, "    mstringlist[6] <- 'htsun'\n");
        fprintf(f_rmap, "  }\n");
        fprintf(f_rmap, "  jrange <- append(jrange, 6)\n");
        fprintf(f_rmap, "  mstringlist[7] <- 'treach'\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "} else {\n");
        fprintf(f_rmap, "  jrange <- c(0, 1, 2, 3)\n");        
        fprintf(f_rmap, "  mstringlist <- c('iii_hflow', 'iii_tflow', 'iii_pflow', 'dii')\n");
        fprintf(f_rmap, "}\n");
        fprintf(f_rmap, "  for ( j in jrange ) {  # loop over all sets of maps to be created\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    mstring <- mstringlist[j+1]\n");
        fprintf(f_rmap, "    if ( mult != 0 ) {\n");
        fprintf(f_rmap, "      d_hdisp <- 1  # absolute maximum raster value\n");
        fprintf(f_rmap, "      d_hdispn <- 0  # absolute minimum raster value\n");
        fprintf(f_rmap, "      thresholds <- 0.0\n");
        fprintf(f_rmap, "    } else if ( j == 6 ) {\n");
        fprintf(f_rmap, "      d_hdisp <- maxtreach  # absolute maximum raster value\n");
        fprintf(f_rmap, "      d_hdispn <- 0  # absolute minimum raster value\n");
        fprintf(f_rmap, "      thresholds <- 0.0\n"); 
        fprintf(f_rmap, "    } else if ( j == 0 ) {\n");
        fprintf(f_rmap, "      d_hdisp <- maxhflow  # absolute maximum raster value\n");
        fprintf(f_rmap, "      d_hdispn <- 0  # absolute minimum raster value\n");
        fprintf(f_rmap, "      thresholds <- thresholdsh\n"); 
        fprintf(f_rmap, "    } else if ( j == 1 ) {\n");
        fprintf(f_rmap, "      d_hdisp <- maxtflow  # absolute maximum raster value\n");
        fprintf(f_rmap, "      d_hdispn <- 0  # absolute minimum raster value\n");
        fprintf(f_rmap, "      thresholds <- thresholdst\n"); 
        fprintf(f_rmap, "    } else if ( j == 2 ) {\n");
        fprintf(f_rmap, "      d_hdisp <- maxpflow  # absolute maximum raster value\n");
        fprintf(f_rmap, "      d_hdispn <- 0  # absolute minimum raster value\n");
        fprintf(f_rmap, "      thresholds <- thresholdsp\n");         
        fprintf(f_rmap, "    } else if ( j == 3 ) {\n");
        fprintf(f_rmap, "      d_hdisp <- maxhdep\n");
        fprintf(f_rmap, "      d_hdispn <- maxhentr\n");
        fprintf(f_rmap, "      thresholds <- thresholdsh\n");         
        fprintf(f_rmap, "    } else {\n");
        fprintf(f_rmap, "      d_hdisp <- maxhtsun\n");
        fprintf(f_rmap, "      d_hdispn <- maxhtsun * -1\n");
        fprintf(f_rmap, "      thresholds <- thresholdsh\n");         
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    # Unit conversion and number of digits for legend according to maximum value:\n");
        fprintf(f_rmap, "    if ( j == 6 ) {\n");
        fprintf(f_rmap, "      mconv <- 1\n");
        fprintf(f_rmap, "    } else if ( max(d_hdisp, -d_hdispn) != 0 ) {\n");
        fprintf(f_rmap, "      mconv <- log10(max(d_hdisp, -d_hdispn))\n");
        fprintf(f_rmap, "    } else {\n");
        fprintf(f_rmap, "      mconv <- 3\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    if ( mconv < 3 ) {\n");
        fprintf(f_rmap, "        mconv2 <- '1'\n");
        fprintf(f_rmap, "        if ( mconv < 0 ) {\n");
        fprintf(f_rmap, "            mdig <- 3\n");
        fprintf(f_rmap, "        } else if ( mconv < 1 ) {\n");
        fprintf(f_rmap, "            mdig <- 2\n");
        fprintf(f_rmap, "        } else if ( mconv < 2 ) {\n");
        fprintf(f_rmap, "            mdig <- 1\n");
        fprintf(f_rmap, "        } else {\n");
        fprintf(f_rmap, "            mdig <- 0\n");
        fprintf(f_rmap, "        }\n");
        fprintf(f_rmap, "    } else if ( mconv < 6 ) {\n");
        fprintf(f_rmap, "        mconv2 <- '0.001'\n");
        fprintf(f_rmap, "        if ( mconv < 4 ) {\n");
        fprintf(f_rmap, "            mdig <- 2\n");
        fprintf(f_rmap, "        } else if ( mconv < 5 ) {\n");
        fprintf(f_rmap, "            mdig <- 1\n");
        fprintf(f_rmap, "        } else {\n");
        fprintf(f_rmap, "            mdig <- 0\n");
        fprintf(f_rmap, "        }\n");
        fprintf(f_rmap, "    } else if ( mconv < 9 ) {\n");
        fprintf(f_rmap, "        mconv2 <- '0.000001'\n");
        fprintf(f_rmap, "        if (mconv < 7) {\n");
        fprintf(f_rmap, "            mdig <- 2\n");
        fprintf(f_rmap, "        } else if (mconv < 8) {\n");
        fprintf(f_rmap, "            mdig <- 1\n");
        fprintf(f_rmap, "        } else {\n");
        fprintf(f_rmap, "            mdig <- 0\n");
        fprintf(f_rmap, "        }\n");
        fprintf(f_rmap, "    } else if ( mconv < 12 ) {\n");
        fprintf(f_rmap, "        mconv2 <- '0.000000001'\n");
        fprintf(f_rmap, "        if ( mconv < 10 ) {\n");
        fprintf(f_rmap, "            mdig <- 2\n");
        fprintf(f_rmap, "        } else if ( mconv < 11 ) {\n");
        fprintf(f_rmap, "            mdig <- 1\n");
        fprintf(f_rmap, "        } else {\n");
        fprintf(f_rmap, "            mdig <- 0\n");
        fprintf(f_rmap, "        }\n");
        fprintf(f_rmap, "    } else if ( mconv < 15 ) {\n");
        fprintf(f_rmap, "        mconv2 <- '0.000000000001'\n");
        fprintf(f_rmap, "        if ( mconv < 13 ) {\n");
        fprintf(f_rmap, "            mdig <- 2\n");
        fprintf(f_rmap, "        } else if ( mconv < 14 ) {\n");
        fprintf(f_rmap, "            mdig <- 1\n");
        fprintf(f_rmap, "        } else {\n");
        fprintf(f_rmap, "            mdig <- 0\n");
        fprintf(f_rmap, "        }\n");
        fprintf(f_rmap, "    } else {\n");
        fprintf(f_rmap, "        mconv2 <- '0.000000000000001'\n");
        fprintf(f_rmap, "        if ( mconv < 16 ) {\n");
        fprintf(f_rmap, "            mdig <- 2\n");
        fprintf(f_rmap, "        } else if ( mconv < 17 ) {\n");
        fprintf(f_rmap, "            mdig <- 1\n");
        fprintf(f_rmap, "        } else {\n");
        fprintf(f_rmap, "            mdig <- 0\n");
        fprintf(f_rmap, "        }\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "    d_hdisp <- d_hdisp * as.numeric(mconv2)\n");
        fprintf(f_rmap, "  \n");
        fprintf(f_rmap, "    if ( mult != 0 ) {\n");
        fprintf(f_rmap, "        thrs5 <- 1.00\n");
        fprintf(f_rmap, "        thrs4 <- 0.80\n");
        fprintf(f_rmap, "        thrs3 <- 0.60\n");
        fprintf(f_rmap, "        thrs2 <- 0.40\n");
        fprintf(f_rmap, "        thrs1 <- 0.20\n");
        fprintf(f_rmap, "        thrsm <- 0.01\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    } else if ( j == 3 ) {  # for entrainment and deposition map\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        d_hdispn <- d_hdispn * as.numeric(mconv2)\n");
        fprintf(f_rmap, "        jimp <- 0  # index of impact parameter\n");
        fprintf(f_rmap, "        thrs5 <- 0\n");
        fprintf(f_rmap, "        thrs4 <- round( d_hdisp, digits <- mdig) # breaks for raster maps\n");
        fprintf(f_rmap, "        thrs3 <- round( d_hdisp * 10 / 20, digits <- mdig )\n");
        fprintf(f_rmap, "        thrs2 <- round( d_hdisp * 6 / 20, digits <- mdig )\n");
        fprintf(f_rmap, "        thrs1 <- round( d_hdisp * 3 / 20, digits <- mdig )\n");
        fprintf(f_rmap, "        if (d_hdisp == 0) {\n");
        fprintf(f_rmap, "            thrsm <- 0\n");
        fprintf(f_rmap, "        } else {\n");
        fprintf(f_rmap, "            thrsm <- min( as.numeric(thresholds) * as.numeric(mconv2), as.numeric(thrs1) * 0.75) # minimum value displayed\n");
        fprintf(f_rmap, "        }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        thrs5n <- 0\n");
        fprintf(f_rmap, "        thrs4n <- round( d_hdispn, digits <- mdig )\n");
        fprintf(f_rmap, "        thrs3n <- round( d_hdispn * 10 / 20, digits <- mdig )\n");
        fprintf(f_rmap, "        thrs2n <- round( d_hdispn * 6 / 20, digits <- mdig )\n");
        fprintf(f_rmap, "        thrs1n <- round( d_hdispn * 3 / 20, digits <- mdig )\n");
        fprintf(f_rmap, "        if ( d_hdispn == 0 ) {\n");
        fprintf(f_rmap, "            thrsmn <- 0\n");
        fprintf(f_rmap, "        } else {\n");
        fprintf(f_rmap, "            thrsmn <- max(-1 * as.numeric(thresholds) * as.numeric(mconv2), as.numeric(thrs1n) * 0.75) # minimum value displayed\n");
        fprintf(f_rmap, "        }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    } else if ( model == 7 && j == 5 ) {  # for tsunami map\n");
        fprintf(f_rmap, "        d_hdispn <- d_hdispn * as.numeric(mconv2)\n");
        fprintf(f_rmap, "        thrs5 <- 0\n");
        fprintf(f_rmap, "        thrs4 <- round( d_hdisp, digits <- mdig ) # breaks for raster maps\n");
        fprintf(f_rmap, "        thrs3 <- round( d_hdisp * 15 / 20, digits <- mdig )\n");
        fprintf(f_rmap, "        thrs2 <- round( d_hdisp * 10 / 20, digits <- mdig )\n");
        fprintf(f_rmap, "        thrs1 <- round( d_hdisp * 5 / 20, digits <- mdig )\n");
        fprintf(f_rmap, "        thrsm <- round( d_hdisp * 1 / 25, digits <- mdig )\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        thrs5n <- 0\n");
        fprintf(f_rmap, "        thrs4n <- round( d_hdispn, digits <- mdig )\n");
        fprintf(f_rmap, "        thrs3n <- round( d_hdispn * 15 / 20, digits <- mdig )\n");
        fprintf(f_rmap, "        thrs2n <- round( d_hdispn * 10 / 20, digits <- mdig )\n");
        fprintf(f_rmap, "        thrs1n <- round( d_hdispn * 5 / 20, digits <- mdig )\n");
        fprintf(f_rmap, "        thrsmn <- round( d_hdispn * 1 / 25, digits <- mdig )\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    } else if ( j == 6 ) {  # for time of reach map\n");
        fprintf(f_rmap, "        jimp <- j\n");
        fprintf(f_rmap, "        thrs5 <- round( d_hdisp, digits <- mdig ) # breaks for raster maps\n");
        fprintf(f_rmap, "        thrs4 <- round( d_hdisp * 10 / 20, digits <- mdig )\n");
        fprintf(f_rmap, "        thrs3 <- round( d_hdisp * 7 / 20, digits <- mdig )\n");
        fprintf(f_rmap, "        thrs2 <- round( d_hdisp * 4 / 20, digits <- mdig )\n");
        fprintf(f_rmap, "        thrs1 <- round( d_hdisp * 2 / 20, digits <- mdig )\n");
        fprintf(f_rmap, "        thrsm <- 0\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        thrs5n <- 0\n");
        fprintf(f_rmap, "        thrs4n <- 0\n");
        fprintf(f_rmap, "        thrs3n <- 0\n");
        fprintf(f_rmap, "        thrs2n <- 0\n");
        fprintf(f_rmap, "        thrs1n <- 0\n");
        fprintf(f_rmap, "        thrsmn <- 0\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    } else {  # for all other maps\n");
        fprintf(f_rmap, "        jimp <- j\n");
        fprintf(f_rmap, "        thrs5 <- round( d_hdisp, digits <- mdig )  # breaks\n");
        fprintf(f_rmap, "        thrs4 <- round( d_hdisp * 12 / 20, digits <- mdig )\n");
        fprintf(f_rmap, "        thrs3 <- round( d_hdisp * 8 / 20, digits <- mdig )\n");
        fprintf(f_rmap, "        thrs2 <- round( d_hdisp * 5 / 20, digits <- mdig )\n");
        fprintf(f_rmap, "        thrs1 <- round( d_hdisp * 2 / 20, digits <- mdig )\n");
        fprintf(f_rmap, "        if ( d_hdisp == 0) {\n");
        fprintf(f_rmap, "            thrsm <- 0\n");
        fprintf(f_rmap, "        } else {\n");
        fprintf(f_rmap, "            thrsm <- min(as.numeric(thresholds) * as.numeric(mconv2), as.numeric(thrs1) * 0.75) # minimum value displayed\n");
        fprintf(f_rmap, "        }\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    if (j == 0 && impdef == 1 ) {\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        impactareaname <- paste(wkdir, prefix, 'results/', prefix, 'ascii/', prefix, 'impactarea.asc', sep = '')\n");
        fprintf(f_rmap, "        impactarea0 <- raster(impactareaname)\n");
        fprintf(f_rmap, "        impactarea <- t(as.matrix(impactarea0))\n");
        fprintf(f_rmap, "        impactarea <- impactarea[, ncol(impactarea):1]\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    if ( j == 0 && depdef == 1 ) {\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        hdepositname <- paste(wkdir, prefix, 'results/', prefix, 'ascii/', prefix, 'hdeposit.asc', sep = '')\n");
        fprintf(f_rmap, "        hdeposit0 <- raster(hdepositname)\n");
        fprintf(f_rmap, "        hdeposit1 <- reclassify(hdeposit0,c(-Inf,as.numeric(thresholdsh),0,as.numeric(thresholdsh),Inf,1))\n");
        fprintf(f_rmap, "        hdeposit <- t(as.matrix(hdeposit1))\n");
        fprintf(f_rmap, "        hdeposit <- hdeposit[, ncol(hdeposit):1]\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "  if ( j == 6) {\n");
        fprintf(f_rmap, "      ntimemin <- ntimemax+1\n");
        fprintf(f_rmap, "  } else {\n");
        fprintf(f_rmap, "      ntimemin <- 0\n");
        fprintf(f_rmap, "  }  \n");
        fprintf(f_rmap, "  if ( mult != 0 ) {\n");
        fprintf(f_rmap, "      ntimemax <- -1\n");
        fprintf(f_rmap, "  }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "  for ( ntimesteps in ntimemin:( ntimemax+1 )) { # loop over all time steps plus one for maps of maximum values\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    if ( mult != 0 ) {\n");
        fprintf(f_rmap, "        fill <- 'iscore'\n");
        fprintf(f_rmap, "    } else if ( ntimesteps < 10 ) {\n");
        fprintf(f_rmap, "        fill <- paste('000', as.character(ntimesteps), sep='')  # formatting time step string\n");
        fprintf(f_rmap, "    } else if ( ntimesteps < 100 ) {\n");
        fprintf(f_rmap, "        fill <- paste('00', as.character(ntimesteps), sep='')\n");
        fprintf(f_rmap, "    } else if ( ntimesteps < 1000 ) {\n");
        fprintf(f_rmap, "        fill <- paste('0', as.character(ntimesteps), sep='')\n");
        fprintf(f_rmap, "    } else {\n");
        fprintf(f_rmap, "        fill <- as.character(ntimesteps)\n");
        fprintf(f_rmap, "    }  \n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    if ( mult != 0 ) {  # strings for names of maps:\n");
        fprintf(f_rmap, "        mstringt <- paste(prefix, mstring, sep='')\n");
        fprintf(f_rmap, "    } else if ( ntimesteps <= ntimemax && model == 7 && j == 5 ) {\n");
        fprintf(f_rmap, "        mstringt <- paste(prefix, 'htsun', fill, sep='')\n");
        fprintf(f_rmap, "    } else if ( ntimesteps <= ntimemax && j != 6 ) {\n");
        fprintf(f_rmap, "        mstringt <- paste( prefix, mstring, fill, sep='')\n");
        fprintf(f_rmap, "        mstrings <- paste( prefix, mstring, '1', fill, sep='')\n");
        fprintf(f_rmap, "        if ( model == 7 ) {\n");
        fprintf(f_rmap, "            mstringf <- paste( prefix, mstring, '2', fill, sep='')\n");
        fprintf(f_rmap, "            mstringw <- paste( prefix, mstring, '3', fill, sep='')\n");
        fprintf(f_rmap, "        }\n");
        fprintf(f_rmap, "    } else if ( j < 3 ) {\n");
        fprintf(f_rmap, "        mstringt <- paste( prefix, mstring, '_max', sep='')\n");
        fprintf(f_rmap, "        mstrings <- paste( prefix, mstring, '1_max', sep='')\n");
        fprintf(f_rmap, "        if ( model == 7 ) {\n");
        fprintf(f_rmap, "            mstringf <- paste( prefix, mstring, '2_max', sep='')\n");
        fprintf(f_rmap, "            mstringw <- paste( prefix, mstring, '3_max', sep='')\n");
        fprintf(f_rmap, "        }\n");
        fprintf(f_rmap, "    } else if ( j == 3 ) {\n");
        fprintf(f_rmap, "        mstringt <- paste( prefix, mstring, '_fin', sep='')\n");
        fprintf(f_rmap, "        mstrings <- paste( prefix, mstring, '1_fin', sep='')\n");
        fprintf(f_rmap, "        if ( model == 7 ) {\n");
        fprintf(f_rmap, "            mstringf <- paste( prefix, mstring, '2_fin', sep='')\n");
        fprintf(f_rmap, "            mstringw <- paste( prefix, mstring, '3_fin', sep='')\n");
        fprintf(f_rmap, "        }\n");
        fprintf(f_rmap, "    } else if ( model == 7 && j == 5 ) {\n");
        fprintf(f_rmap, "        mstringt <- paste( prefix, 'htsun_max', sep='')\n");
        fprintf(f_rmap, "    } else if ( ntimesteps == ntimemax + 1 && j == 6 ) {\n");
        fprintf(f_rmap, "        mstringt <- paste( prefix, 'treach', sep='')\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    # Reading flow velocity and direction data from file defining control variable:\n");
        fprintf(f_rmap, "    if (ntimesteps > 0 && ntimesteps <= ntimemax) {\n");
        fprintf(f_rmap, "        velocity <- 1\n");
        fprintf(f_rmap, "    } else {\n");
        fprintf(f_rmap, "        velocity <- 0\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "    if (velocity == 1) {\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        dirline <- ntimesteps - 1\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        intable = paste(wkdir, prefix, 'results/', prefix, 'files/', prefix, 'directions1.txt', \n");
        fprintf(f_rmap, "            sep = '')\n");
        fprintf(f_rmap, "        dirx <- scan(intable, nlines = 1, quiet = TRUE)  #(PHASE 1) x coordinate\n");
        fprintf(f_rmap, "        diry <- scan(intable, skip = 1, nlines = 1, quiet = TRUE)  #(PHASE 1) y coordinate\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        if (model <= 3) {\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "            dirh <- scan(intable, skip = 2 + 3 * dirline, nlines = 1, quiet = TRUE)  #flow height\n");
        fprintf(f_rmap, "            dirvx <- scan(intable, skip = 3 + 3 * dirline, nlines = 1, quiet = TRUE)  #flow velocity in x direction\n");
        fprintf(f_rmap, "            dirvy <- scan(intable, skip = 4 + 3 * dirline, nlines = 1, quiet = TRUE)  #flow velocity in y direction\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        } else if (model == 7) {\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "            intable2 = paste(wkdir, prefix, 'results/', prefix, 'files/', prefix, 'directions2.txt', \n");
        fprintf(f_rmap, "                sep = '')\n");
        fprintf(f_rmap, "            dirxf <- scan(intable2, nlines = 1, quiet = TRUE)  #PHASE 2 x coordinate\n");
        fprintf(f_rmap, "            diryf <- scan(intable2, skip = 1, nlines = 1, quiet = TRUE)  #PHASE 2 y coordinate\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "            dirh <- scan(intable, skip = 2 + 3 * dirline, nlines = 1, quiet = TRUE)  #PHASE 1 flow height\n");
        fprintf(f_rmap, "            dirvx <- scan(intable, skip = 3 + 3 * dirline, nlines = 1, quiet = TRUE)  #PHASE 1 flow velocity in x direction\n");
        fprintf(f_rmap, "            dirvy <- scan(intable, skip = 4 + 3 * dirline, nlines = 1, quiet = TRUE)  #PHASE 1 flow velocity in y direction\n");
        fprintf(f_rmap, "            dirhf <- scan(intable2, skip = 2 + 3 * dirline, nlines = 1, quiet = TRUE)  #PHASE 2 flow height\n");
        fprintf(f_rmap, "            dirvxf <- scan(intable2, skip = 3 + 3 * dirline, nlines = 1, quiet = TRUE)  #PHASE 2 flow velocity in x direction\n");
        fprintf(f_rmap, "            dirvyf <- scan(intable2, skip = 4 + 3 * dirline, nlines = 1, quiet = TRUE)  #PHASE 2 flow velocity in y direction\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "            dirhf[dirhf > 0] <- 1\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "            intable3 = paste(wkdir, prefix, 'results/', prefix, 'files/', prefix, 'directions3.txt', \n");
        fprintf(f_rmap, "                sep = '')\n");
        fprintf(f_rmap, "            dirxw <- scan(intable3, nlines = 1, quiet = TRUE)  #PHASE 3 x coordinate\n");
        fprintf(f_rmap, "            diryw <- scan(intable3, skip = 1, nlines = 1, quiet = TRUE)  #PHASE 3 y coordinate\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "            dirhw <- scan(intable3, skip = 2 + 3 * dirline, nlines = 1, quiet = TRUE)  #PHASE 3 flow height\n");
        fprintf(f_rmap, "            dirvxw <- scan(intable3, skip = 3 + 3 * dirline, nlines = 1, quiet = TRUE)  #PHASE 3 flow velocity in x direction\n");
        fprintf(f_rmap, "            dirvyw <- scan(intable3, skip = 4 + 3 * dirline, nlines = 1, quiet = TRUE)  #PHASE 3 flow velocity in y direction\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "            dirhw[dirhw > 0] <- 1\n");
        fprintf(f_rmap, "        }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        dirh[dirh > 0] <- 1\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    # Computing map extent\n");
        fprintf(f_rmap, "    xmin <- xmin0\n");
        fprintf(f_rmap, "    xmax <- xmax0\n");
        fprintf(f_rmap, "    ymin <- ymin0\n");
        fprintf(f_rmap, "    ymax <- ymax0\n");        
        fprintf(f_rmap, "    xdiff <- xmax - xmin  #extent in x direction\n");
        fprintf(f_rmap, "    ydiff <- ymax - ymin  #extent in y direction\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    # if ( ydiff > xdiff ) { #ensuring minimum extent in x direction\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    # xmin<-xmin-(ydiff-xdiff)/2 xmax<-xmax+(ydiff-xdiff)/2 xdiff<-xmax-xmin\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    # } else\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    if (xdiff > 1.6 * ydiff) {\n");
        fprintf(f_rmap, "        # ensuring minimum extent in y direction\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        ymin <- ymin - (xdiff/1.6 - ydiff)/2\n");
        fprintf(f_rmap, "        ymax <- ymax + (xdiff/1.6 - ydiff)/2\n");
        fprintf(f_rmap, "        ydiff <- ymax - ymin\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    vunit <- 0.05 * min(xdiff, ydiff)/maxvflow\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    # Building map geometry\n");
        fprintf(f_rmap, "    asprat <- xdiff/(1.4 * ydiff)  #x/y ratio\n");
        fprintf(f_rmap, "    margx <- 6  #margin in x direction\n");
        fprintf(f_rmap, "    margy <- 2.5  #margin in y direction\n");
        fprintf(f_rmap, "    dispx <- 15  #total width of output image\n");
        fprintf(f_rmap, "    dispy <- (dispx - margx)/asprat + margy  #total height of output image\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    # Adapting boundaries for vector layers\n");
        fprintf(f_rmap, "    xsmin <- xmin + 0.035 * xdiff\n");
        fprintf(f_rmap, "    xsmax <- xmax - 0.035 * xdiff\n");
        fprintf(f_rmap, "    ysmin <- ymin + 0.035 * ydiff\n");
        fprintf(f_rmap, "    ysmax <- ymax - 0.035 * ydiff\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    # Defining graphic parameters for raster plot\n");
        fprintf(f_rmap, "    if (model <= 3 || fill == 'iscore' || mstring == 'basechange' || mstring == 'htsun' || mstring == 'treach') {\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        if (mstring == 'hflow' || mstring == 'treach' || mstring == 'iii_hflow' || \n");
        fprintf(f_rmap, "            mstring == 'dii') {\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "            if (fill == 'iscore') {\n");
        fprintf(f_rmap, "                munit <- ''\n");
        fprintf(f_rmap, "                if (mstring == 'iii_hflow') \n");
        fprintf(f_rmap, "                    mlabel <- 'Impact indicator index'\n");
        fprintf(f_rmap, "                if (mstring == 'dii') \n");
        fprintf(f_rmap, "                    mlabel <- 'Deposition indicator index'\n");
        fprintf(f_rmap, "                rcol1 <- rgb(0, 0, 0.8, 0.4)\n");
        fprintf(f_rmap, "                rcol7 <- rgb(0.15, 0.1, 0.6, 0.5)\n");
        fprintf(f_rmap, "                rcol13 <- rgb(0.3, 0.2, 0.4, 0.6)\n");
        fprintf(f_rmap, "                rcol19 <- rgb(0.45, 0.3, 0.2, 0.7)\n");
        fprintf(f_rmap, "                rcol25 <- rgb(0.6, 0.4, 0, 0.8)\n");
        fprintf(f_rmap, "            } else if (mstring == 'hflow') {\n");
        fprintf(f_rmap, "                if (mconv2 == '1') \n");
        fprintf(f_rmap, "                    munit <- 'm'\n");
        fprintf(f_rmap, "                if (mconv2 == '0.001') \n");
        fprintf(f_rmap, "                    munit <- 'km'\n");
        fprintf(f_rmap, "                mlabel <- 'Flow height'\n");
        fprintf(f_rmap, "                rcol1 <- rgb(0.5, 0.5, 0, 0.2)\n");
        fprintf(f_rmap, "                rcol7 <- rgb(0.75, 0.25, 0, 0.4)\n");
        fprintf(f_rmap, "                rcol13 <- rgb(1, 0, 0, 0.6)\n");
        fprintf(f_rmap, "                rcol19 <- rgb(0.8, 0, 0.2, 0.8)\n");
        fprintf(f_rmap, "                rcol25 <- rgb(0.6, 0, 0.4, 1)\n");
        fprintf(f_rmap, "            } else if (mstring == 'treach') {\n");
        fprintf(f_rmap, "                munit <- 's'\n");
        fprintf(f_rmap, "                mlabel <- 'Time of reach'\n");
        fprintf(f_rmap, "                rcol1 <- rgb(0.7, 0.3, 0.1, 1)\n");
        fprintf(f_rmap, "                rcol7 <- rgb(0.55, 0.275, 0.2, 0.8)\n");
        fprintf(f_rmap, "                rcol13 <- rgb(0.4, 0.25, 0.3, 0.6)\n");
        fprintf(f_rmap, "                rcol19 <- rgb(0.25, 0.225, 0.4, 0.4)\n");
        fprintf(f_rmap, "                rcol25 <- rgb(0.1, 0.2, 0.5, 0.2)\n");
        fprintf(f_rmap, "            }\n");
        fprintf(f_rmap, "        } else if (mstring == 'tflow' || mstring == 'iii_tflow') {\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "            if (fill != 'iscore') {\n");
        fprintf(f_rmap, "                if (mconv2 == '1') \n");
        fprintf(f_rmap, "                    munit <- 'J'\n");
        fprintf(f_rmap, "                if (mconv2 == '0.001') \n");
        fprintf(f_rmap, "                    munit <- 'kJ'\n");
        fprintf(f_rmap, "                if (mconv2 == '0.000001') \n");
        fprintf(f_rmap, "                    munit <- 'MJ'\n");
        fprintf(f_rmap, "                if (mconv2 == '0.000000001') \n");
        fprintf(f_rmap, "                    munit <- 'GJ'\n");
        fprintf(f_rmap, "                if (mconv2 == '0.000000000001') \n");
        fprintf(f_rmap, "                    munit <- 'TJ'\n");
        fprintf(f_rmap, "                mlabel <- 'Flow kinetic energy'\n");
        fprintf(f_rmap, "                rcol1 <- rgb(0, 0.5, 0.5, 0.2)\n");
        fprintf(f_rmap, "                rcol7 <- rgb(0.15, 0.45, 0.4, 0.4)\n");
        fprintf(f_rmap, "                rcol13 <- rgb(0.3, 0.4, 0.3, 0.6)\n");
        fprintf(f_rmap, "                rcol19 <- rgb(0.45, 0.35, 0.2, 0.8)\n");
        fprintf(f_rmap, "                rcol25 <- rgb(0.6, 0.3, 0.1, 1)\n");
        fprintf(f_rmap, "            } else {\n");
        fprintf(f_rmap, "                munit <- ''\n");
        fprintf(f_rmap, "                mlabel <- 'Impact indicator index'\n");
        fprintf(f_rmap, "                rcol1 <- rgb(0, 0.5, 0.5, 0.4)\n");
        fprintf(f_rmap, "                rcol7 <- rgb(0.15, 0.45, 0.4, 0.5)\n");
        fprintf(f_rmap, "                rcol13 <- rgb(0.3, 0.4, 0.3, 0.6)\n");
        fprintf(f_rmap, "                rcol19 <- rgb(0.45, 0.35, 0.2, 0.7)\n");
        fprintf(f_rmap, "                rcol25 <- rgb(0.6, 0.3, 0.1, 0.8)\n");
        fprintf(f_rmap, "            }\n");
        fprintf(f_rmap, "        } else if (mstring == 'pflow' || mstring == 'iii_pflow') {\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "            if (fill != 'iscore') {\n");
        fprintf(f_rmap, "                if (mconv2 == '1') \n");
        fprintf(f_rmap, "                    munit <- 'Pa'\n");
        fprintf(f_rmap, "                if (mconv2 == '0.001') \n");
        fprintf(f_rmap, "                    munit <- 'kPa'\n");
        fprintf(f_rmap, "                if (mconv2 == '0.000001') \n");
        fprintf(f_rmap, "                    munit <- 'MPa'\n");
        fprintf(f_rmap, "                if (mconv2 == '0.000000001') \n");
        fprintf(f_rmap, "                    munit <- 'GPa'\n");
        fprintf(f_rmap, "                if (mconv2 == '0.000000000001') \n");
        fprintf(f_rmap, "                    munit <- 'TPa'\n");
        fprintf(f_rmap, "                mlabel <- 'Flow pressure'\n");
        fprintf(f_rmap, "                rcol1 <- rgb(0.2, 0, 0.6, 0.2)\n");
        fprintf(f_rmap, "                rcol7 <- rgb(0.4, 0.3, 0.35, 0.4)\n");
        fprintf(f_rmap, "                rcol13 <- rgb(0.4, 0.3, 0.3, 0.6)\n");
        fprintf(f_rmap, "                rcol19 <- rgb(0.5, 0.45, 0.15, 0.8)\n");
        fprintf(f_rmap, "                rcol25 <- rgb(0.6, 0.6, 0, 1)\n");
        fprintf(f_rmap, "            } else {\n");
        fprintf(f_rmap, "                munit <- ''\n");
        fprintf(f_rmap, "                mlabel <- 'Impact indicator index'\n");
        fprintf(f_rmap, "                rcol1 <- rgb(0.2, 0, 0.6, 0.4)\n");
        fprintf(f_rmap, "                rcol7 <- rgb(0.4, 0.3, 0.35, 0.5)\n");
        fprintf(f_rmap, "                rcol13 <- rgb(0.4, 0.3, 0.3, 0.6)\n");
        fprintf(f_rmap, "                rcol19 <- rgb(0.5, 0.45, 0.15, 0.7)\n");
        fprintf(f_rmap, "                rcol25 <- rgb(0.6, 0.6, 0, 0.8)\n");
        fprintf(f_rmap, "            }\n");
        fprintf(f_rmap, "        } else if (mstring == 'basechange') {\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "            if (mconv2 == '1') \n");
        fprintf(f_rmap, "                munit <- 'm'\n");
        fprintf(f_rmap, "            if (mconv2 == '0.001') \n");
        fprintf(f_rmap, "                munit <- 'km'\n");
        fprintf(f_rmap, "            mlabel <- 'Change of basal topography'\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "            rcol3 <- rgb(0.1, 0.2, 0.6, 1)\n");
        fprintf(f_rmap, "            rcol7 <- rgb(0.1, 0.2, 0.6, 0.75)\n");
        fprintf(f_rmap, "            rcol11 <- rgb(0.1, 0.2, 0.6, 0.5)\n");
        fprintf(f_rmap, "            rcol15 <- rgb(0.1, 0.2, 0.6, 0.25)\n");
        fprintf(f_rmap, "            rcol17 <- rgb(1, 1, 1, 0)\n");
        fprintf(f_rmap, "            rcol20 <- rgb(0.2, 0.6, 0.1, 0.25)\n");
        fprintf(f_rmap, "            rcol24 <- rgb(0.2, 0.6, 0.1, 0.5)\n");
        fprintf(f_rmap, "            rcol28 <- rgb(0.2, 0.6, 0.1, 0.75)\n");
        fprintf(f_rmap, "            rcol32 <- rgb(0.2, 0.6, 0.1, 1)\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        } else if (mstring == 'htsun') {\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "            if (mconv2 == '1') \n");
        fprintf(f_rmap, "                munit <- 'm'\n");
        fprintf(f_rmap, "            if (mconv2 == '0.001') \n");
        fprintf(f_rmap, "                munit <- 'km'\n");
        fprintf(f_rmap, "            mlabel <- 'Tsunami height'\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "            rcol3 <- rgb(0, 0.6, 0.1, 1)\n");
        fprintf(f_rmap, "            rcol7 <- rgb(0, 0.6, 0.1, 0.75)\n");
        fprintf(f_rmap, "            rcol11 <- rgb(0, 0.6, 0.1, 0.5)\n");
        fprintf(f_rmap, "            rcol15 <- rgb(0, 0.6, 0.1, 0.25)\n");
        fprintf(f_rmap, "            rcol17 <- rgb(1, 1, 1, 0)\n");
        fprintf(f_rmap, "            rcol20 <- rgb(0, 0.1, 0.6, 0.25)\n");
        fprintf(f_rmap, "            rcol24 <- rgb(0, 0.1, 0.6, 0.5)\n");
        fprintf(f_rmap, "            rcol28 <- rgb(0, 0.1, 0.6, 0.75)\n");
        fprintf(f_rmap, "            rcol32 <- rgb(0, 0.1, 0.6, 1)\n");
        fprintf(f_rmap, "        }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        if (mstring == 'basechange' || mstring == 'htsun') {\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "            ctext = vector('expression', 10)  #vector for colour bar labels\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "            if (mstring == 'basechange') {\n");
        fprintf(f_rmap, "                ctext[1] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs4n), \n");
        fprintf(f_rmap, "                    mdig), nsmall = mdig), fmunit = format(munit)))[2]\n");
        fprintf(f_rmap, "            } else {\n");
        fprintf(f_rmap, "                ctext[1] <- ''\n");
        fprintf(f_rmap, "            }\n");
        fprintf(f_rmap, "            ctext[2] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs3n), \n");
        fprintf(f_rmap, "                mdig), nsmall = mdig), fmunit = format(munit)))[2]\n");
        fprintf(f_rmap, "            ctext[3] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs2n), \n");
        fprintf(f_rmap, "                mdig), nsmall = mdig), fmunit = format(munit)))[2]\n");
        fprintf(f_rmap, "            ctext[4] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs1n), \n");
        fprintf(f_rmap, "                mdig), nsmall = mdig), fmunit = format(munit)))[2]\n");
        fprintf(f_rmap, "            ctext[5] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(as.numeric(thrsmn)), \n");
        fprintf(f_rmap, "                fmunit = format(munit)))[2]\n");
        fprintf(f_rmap, "            ctext[6] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(as.numeric(thrsm)), \n");
        fprintf(f_rmap, "                fmunit = format(munit)))[2]\n");
        fprintf(f_rmap, "            ctext[7] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs1), \n");
        fprintf(f_rmap, "                mdig), nsmall = mdig), fmunit = format(munit)))[2]\n");
        fprintf(f_rmap, "            ctext[8] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs2), \n");
        fprintf(f_rmap, "                mdig), nsmall = mdig), fmunit = format(munit)))[2]\n");
        fprintf(f_rmap, "            ctext[9] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs3), \n");
        fprintf(f_rmap, "                mdig), nsmall = mdig), fmunit = format(munit)))[2]\n");
        fprintf(f_rmap, "            ctext[10] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs4), \n");
        fprintf(f_rmap, "                mdig), nsmall = mdig), fmunit = format(munit)))[2]\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "            rthresholds <- c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5)\n");
        fprintf(f_rmap, "            rlabels <- ctext\n");
        fprintf(f_rmap, "            rcolours <- c(rcol3, rcol7, rcol11, rcol15, rcol17, rcol20, rcol24, rcol28, \n");
        fprintf(f_rmap, "                rcol32)\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        } else {\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "            ctext = vector('expression', 6)  #vector for colour bar labels\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "            ctext[1] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(as.numeric(thrsm)), \n");
        fprintf(f_rmap, "                fmunit = format(munit)))[2]\n");
        fprintf(f_rmap, "            ctext[2] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs1), \n");
        fprintf(f_rmap, "                mdig), nsmall = mdig), fmunit = format(munit)))[2]\n");
        fprintf(f_rmap, "            ctext[3] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs2), \n");
        fprintf(f_rmap, "                mdig), nsmall = mdig), fmunit = format(munit)))[2]\n");
        fprintf(f_rmap, "            ctext[4] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs3), \n");
        fprintf(f_rmap, "                mdig), nsmall = mdig), fmunit = format(munit)))[2]\n");
        fprintf(f_rmap, "            ctext[5] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs4), \n");
        fprintf(f_rmap, "                mdig), nsmall = mdig), fmunit = format(munit)))[2]\n");
        fprintf(f_rmap, "            ctext[6] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs5), \n");
        fprintf(f_rmap, "                mdig), nsmall = mdig), fmunit = format(munit)))[2]\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "            rthresholds <- c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5)\n");
        fprintf(f_rmap, "            rlabels <- ctext\n");
        fprintf(f_rmap, "            rcolours <- c(rcol1, rcol7, rcol13, rcol19, rcol25)\n");
        fprintf(f_rmap, "        }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    } else if (model == 7) {\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        if (mstring == 'hflow') {\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "            if (mconv2 == '1') \n");
        fprintf(f_rmap, "                munit <- 'm'\n");
        fprintf(f_rmap, "            if (mconv2 == '0.001') \n");
        fprintf(f_rmap, "                munit <- 'km'\n");
        fprintf(f_rmap, "            mlabel <- 'Flow height'\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "            rcolours <- rgb(rep(0:4, 25, each = 25) * 0.25, rep(0:4, 125, each = 5) * \n");
        fprintf(f_rmap, "                0.25, rep(0:4, 625, each = 1) * 0.25, rep(1:5, 5, each = 125) * 0.2)\n");
        fprintf(f_rmap, "            rthresholds <- rep(0:625) + 0.5\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        } else if (mstring == 'tflow') {\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "            if (mconv2 == '1') \n");
        fprintf(f_rmap, "                munit <- 'J'\n");
        fprintf(f_rmap, "            if (mconv2 == '0.001') \n");
        fprintf(f_rmap, "                munit <- 'kJ'\n");
        fprintf(f_rmap, "            if (mconv2 == '0.000001') \n");
        fprintf(f_rmap, "                munit <- 'MJ'\n");
        fprintf(f_rmap, "            if (mconv2 == '0.000000001') \n");
        fprintf(f_rmap, "                munit <- 'GJ'\n");
        fprintf(f_rmap, "            if (mconv2 == '0.000000000001') \n");
        fprintf(f_rmap, "                munit <- 'TJ'\n");
        fprintf(f_rmap, "            mlabel <- 'Flow kinetic\\nenergy'\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "            rcolours <- rgb(rep(0:4, 25, each = 25) * 0.25, rep(0:4, 125, each = 5) * \n");
        fprintf(f_rmap, "                0.25, rep(0:4, 625, each = 1) * 0.25, rep(1:5, 5, each = 125) * 0.2)\n");
        fprintf(f_rmap, "            rthresholds <- rep(0:625) + 0.5\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        } else if (mstring == 'pflow') {\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "            if (mconv2 == '1') \n");
        fprintf(f_rmap, "                munit <- 'Pa'\n");
        fprintf(f_rmap, "            if (mconv2 == '0.001') \n");
        fprintf(f_rmap, "                munit <- 'kPa'\n");
        fprintf(f_rmap, "            if (mconv2 == '0.000001') \n");
        fprintf(f_rmap, "                munit <- 'MPa'\n");
        fprintf(f_rmap, "            if (mconv2 == '0.000000001') \n");
        fprintf(f_rmap, "                munit <- 'GPa'\n");
        fprintf(f_rmap, "            if (mconv2 == '0.000000000001') \n");
        fprintf(f_rmap, "                munit <- 'TPa'\n");
        fprintf(f_rmap, "            mlabel <- 'Flow pressure'\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "            rcolours <- rgb(rep(0:4, 25, each = 25) * 0.25, rep(0:4, 125, each = 5) * \n");
        fprintf(f_rmap, "                0.25, rep(0:4, 625, each = 1) * 0.25, rep(1:5, 5, each = 125) * 0.2)\n");
        fprintf(f_rmap, "            rthresholds <- rep(0:625) + 0.5\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        } else if (mstring == 'htsun') {\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "            if (mconv2 == '1') \n");
        fprintf(f_rmap, "                munit <- 'm'\n");
        fprintf(f_rmap, "            if (mconv2 == '0.001') \n");
        fprintf(f_rmap, "                munit <- 'km'\n");
        fprintf(f_rmap, "            mlabel <- 'Height of tsunami'\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "            rcolours <- rgb(rep(0:4, 25, each = 25) * 0.25, rep(0:4, 125, each = 5) * \n");
        fprintf(f_rmap, "                0.25, rep(0:4, 625, each = 1) * 0.25, rep(1:5, 5, each = 125) * 0.2)\n");
        fprintf(f_rmap, "            rthresholds <- rep(0:625) + 0.5\n");
        fprintf(f_rmap, "        }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        ctext = vector('expression', 6)  #vector for contour line legend labels\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "         ctext[6] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(as.numeric(thrsm)), \n");
        fprintf(f_rmap, "            fmunit = format(munit)))[2]\n");
        fprintf(f_rmap, "        ctext[5] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs1), \n");
        fprintf(f_rmap, "            mdig), nsmall = mdig), fmunit = format(munit)))[2]\n");
        fprintf(f_rmap, "        ctext[4] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs2), \n");
        fprintf(f_rmap, "            mdig), nsmall = mdig), fmunit = format(munit)))[2]\n");
        fprintf(f_rmap, "        ctext[3] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs3), \n");
        fprintf(f_rmap, "            mdig), nsmall = mdig), fmunit = format(munit)))[2]\n");
        fprintf(f_rmap, "        ctext[2] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs4), \n");
        fprintf(f_rmap, "            mdig), nsmall = mdig), fmunit = format(munit)))[2]\n");
        fprintf(f_rmap, "        ctext[1] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs5), \n");
        fprintf(f_rmap, "            mdig), nsmall = mdig), fmunit = format(munit)))[2]\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    # Importing raster layers\n");
        fprintf(f_rmap, "    cat(paste('Plotting ', mstringt, ' ...', sep=''))\n"); 
        fprintf(f_rmap, "    mstringtname <- paste(wkdir, prefix, 'results/', prefix, 'ascii/', mstringt, '.asc', sep = '')\n");
        fprintf(f_rmap, "    mstringtr <- raster(mstringtname)  #raster of total value\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    if (ortho1 != 'None') {\n");
        fprintf(f_rmap, "        ortho1r <- raster(ortho1)\n");
        fprintf(f_rmap, "        ortho2r <- raster(ortho2)\n");
        fprintf(f_rmap, "        ortho3r <- raster(ortho3)\n");
        fprintf(f_rmap, "        orthogray <- 153+(ortho1r+ortho2r+ortho3r)/7.5 #grayscale orthophoto\n");
        fprintf(f_rmap, "        orthogray[orthogray<153] <- 153\n");
        fprintf(f_rmap, "        orthogray[is.na(orthogray)==TRUE] <- 153\n");
        fprintf(f_rmap, "    } else if ( mult != 0 ) {\n");
        fprintf(f_rmap, "        hillshadename <- paste(wkdir, prefix, 'results/', prefix, 'ascii/', prefix, 'hillshade0000.asc', sep = '')\n");
        fprintf(f_rmap, "        hillshade <- raster(hillshadename)  #raster of hillshade\n");
        fprintf(f_rmap, "    } else if ( ntimesteps <= ntimemax ) {\n");
        fprintf(f_rmap, "        hillshadename <- paste(wkdir, prefix, 'results/', prefix, 'ascii/', prefix, 'hillshade', fill, '.asc', sep = '')\n");
        fprintf(f_rmap, "        hillshade <- raster(hillshadename)  #raster of hillshade\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    # Computing vectors\n");
        fprintf(f_rmap, "    if ( mult == 0 ) { # for single model run\n");
        fprintf(f_rmap, "        mstringrname <- paste(wkdir, prefix, 'results/', prefix, 'ascii/', prefix, 'hflow0000.asc', sep = '')\n");
        fprintf(f_rmap, "        startshape0 <- raster(mstringrname)  #raster of release area\n");
        fprintf(f_rmap, "        startshape <- t(as.matrix(startshape0)) # raster map for contour creation\n");
        fprintf(f_rmap, "        startshape <- startshape[, ncol(startshape):1]\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    # Computing composite raster and defining contour line parameters\n");
        fprintf(f_rmap, "    if ( mult != 0 ) { # for multiple model runs\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "      a2 <- reclassify(mstringtr, c(-Inf,thrsm,0,thrsm,thrs1,1,thrs1,thrs2,2,thrs2,thrs3,3,thrs3,thrs4,4,thrs4,thrs5,5))\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "      a2c <- t(as.matrix(a2)) # raster map for contour creation\n");
        fprintf(f_rmap, "      a2c <- a2c[, ncol(a2c):1]\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "      contstep = 1  # contour interval\n");
        fprintf(f_rmap, "      contmin = 1  # minimum contour\n");
        fprintf(f_rmap, "      contmax = 5  # maximum contour  \n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    } else if ( model == 7 && j == 5 ) { # for tsunami map\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "      a2 <- reclassify(mstringtr * as.numeric(mconv2), c(-Inf,thrs3n,1,thrs3n,thrs2n,2,thrs2n,thrs1n,3,thrs1n,thrsmn,4,thrsmn,thrsm,5,\n");
        fprintf(f_rmap, "          thrsm,thrs1,6,thrs1,thrs2,7,thrs2,thrs3,8,thrs3,Inf,9))\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "      a2c <- t(as.matrix(a2)) # raster map for contour creation\n");
        fprintf(f_rmap, "      a2c <- a2c[, ncol(a2c):1]\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "      contstep = 1  # contour interval\n");
        fprintf(f_rmap, "      contmin = 0  # minimum contour\n");
        fprintf(f_rmap, "      contmax = 9  # maximum contour  \n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    } else if ( model <= 3 || j == 3 || (ntimesteps == ntimemax + 1 && j == 6 )) {\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "      if ( j == 3 ) { # for entrainment and deposition map:\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        a2 <- reclassify((mstringtr) * as.numeric(mconv2), c(-Inf,thrs3n,1,thrs3n,thrs2n,2,thrs2n,thrs1n,3,thrs1n,thrsmn,4,thrsmn,thrsm,5,\n");
        fprintf(f_rmap, "            thrsm,thrs1,6,thrs1,thrs2,7,thrs2,thrs3,8,thrs3,Inf,9))\n");
        fprintf(f_rmap, "    \n");
        fprintf(f_rmap, "        a2c <- t(as.matrix(a2)) # raster map for contour creation\n");
        fprintf(f_rmap, "        a2c <- a2c[, ncol(a2c):1]\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        contstep = 1  # contour interval\n");
        fprintf(f_rmap, "        contmin = 0  # minimum contour\n");
        fprintf(f_rmap, "        contmax = 9  # maximum contour\n");
        fprintf(f_rmap, "  \n");
        fprintf(f_rmap, "      } else if ( j == 6 ) { # for time of reach map:\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        a2 <- reclassify(mstringtr * as.numeric(mconv2), c(-Inf,thrsm,0,thrsm,thrs1,1,thrs1,thrs2,2,thrs2,thrs3,3,thrs3,thrs4,4,thrs4,thrs5,5,thrs5,Inf,6))\n");
        fprintf(f_rmap, "    \n");
        fprintf(f_rmap, "        a2c <- t(as.matrix(a2)) # raster map for contour creation\n");
        fprintf(f_rmap, "        a2c <- a2c[, ncol(a2c):1]\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        contstep = 1  # contour interval\n");
        fprintf(f_rmap, "        contmin = 1  # minimum contour\n");
        fprintf(f_rmap, "        contmax = 6  # maximum contour\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "      } else { # for all other maps:\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        a2 <- reclassify(mstringtr * as.numeric(mconv2), c(-Inf,thrsm,0,thrsm,thrs1,1,thrs1,thrs2,2,thrs2,thrs3,3,thrs3,thrs4,4,thrs4,Inf,5))\n");
        fprintf(f_rmap, "    \n");
        fprintf(f_rmap, "        a2c <- t(as.matrix(a2)) # raster map for contour creation\n");
        fprintf(f_rmap, "        a2c <- a2c[, ncol(a2c):1]\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        contstep = 1  # contour interval\n");
        fprintf(f_rmap, "        contmin = 1  # minimum contour\n");
        fprintf(f_rmap, "        contmax = 5  # maximum contour\n");
        fprintf(f_rmap, "      }\n");
        fprintf(f_rmap, "    } else if ( j != 6 ) {\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "      mstringsname <- paste(wkdir, prefix, 'results/', prefix, 'ascii/', mstrings, '.asc', sep = '')\n");
        fprintf(f_rmap, "      mstringsr <- raster(mstringsname)  #raster of phase 1 value\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "      mstringfname <- paste(wkdir, prefix, 'results/', prefix, 'ascii/', mstringf, '.asc', sep = '')\n");
        fprintf(f_rmap, "      mstringfr <- raster(mstringfname)  #raster of phase 2 value\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "      mstringwname <- paste(wkdir, prefix, 'results/', prefix, 'ascii/', mstringw, '.asc', sep = '')\n");
        fprintf(f_rmap, "      mstringwr <- raster(mstringwname)  #raster of phase 3 value\n");
        fprintf(f_rmap, "\n");        
        fprintf(f_rmap, "      a21 <- mstringtr * as.numeric(mconv2)\n");
        fprintf(f_rmap, "  \n");
        fprintf(f_rmap, "      a22r <- abs(mstringsr/mstringtr)\n");
        fprintf(f_rmap, "      a23r <- abs(mstringfr/mstringtr)\n");
        fprintf(f_rmap, "      a24r <- abs(mstringwr/mstringtr)\n");
        fprintf(f_rmap, "      a22 <- mask(a22r, mstringsr, maskvalue=0, updatevalue=0) # ratio of PHASE 1 component\n");
        fprintf(f_rmap, "      a23 <- mask(a23r, mstringfr, maskvalue=0, updatevalue=0) # ratio of PHASE 2 component\n");
        fprintf(f_rmap, "      a24 <- mask(a24r, mstringwr, maskvalue=0, updatevalue=0) # ratio of PHASE 3 component\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "      a2 <- reclassify(a21, c(-Inf,thrsm,0,thrsm,thrs1,1,thrs1,thrs2,2,thrs2,thrs3,3,thrs3,thrs4,4,thrs4,Inf,5)) # reclass for magnitude\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "      a2sr <- reclassify(a22, c(-Inf,0.2,1,0.2,0.4,2,0.4,0.6,3,0.6,0.8,4,0.8,Inf,5))\n");
        fprintf(f_rmap, "      a2fr <- reclassify(a23, c(-Inf,0.2,1,0.2,0.4,2,0.4,0.6,3,0.6,0.8,4,0.8,Inf,5))\n");
        fprintf(f_rmap, "      a2wr <- reclassify(a24, c(-Inf,0.2,1,0.2,0.4,2,0.4,0.6,3,0.6,0.8,4,0.8,Inf,5))\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "      a2s <- mask(a2sr, a2, maskvalue=0, updatevalue=0) # reclass for PHASE 1 ratio\n");
        fprintf(f_rmap, "      a2f <- mask(a2fr, a2, maskvalue=0, updatevalue=0) # reclass for PHASE 2 ratio\n");
        fprintf(f_rmap, "      a2w <- mask(a2wr, a2, maskvalue=0, updatevalue=0) # reclass for PHASE 3 ratio\n");
        fprintf(f_rmap, "      a2 <- (a2-1)*125+(a2s-1)*25+(a2f-1)*5+a2w #combined ratio\n");
        fprintf(f_rmap, "      a2c <- t(as.matrix(a2+124)) # raster map for contour creation\n");
        fprintf(f_rmap, "      a2c <- a2c[, ncol(a2c):1]\n");
        fprintf(f_rmap, "      a2c[is.na(a2c) == TRUE] <- 124\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "      contstep = 125  # contour interval\n");
        fprintf(f_rmap, "      contmin = 125  # minimum contour\n");
        fprintf(f_rmap, "      contmax = 625  # maximum contour\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    rastplot <- a2\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    # Creating plot file\n");
        fprintf(f_rmap, "    if (fill == 'iscore') {\n");
        fprintf(f_rmap, "        mapplot = paste(wkdir, prefix, 'results/', prefix, 'plots/', prefix, mstring, '.png', \n");
        fprintf(f_rmap, "            sep = '')\n");
        fprintf(f_rmap, "    } else if (ntimesteps <= ntimemax) {\n");
        fprintf(f_rmap, "        mapplot = paste(wkdir, prefix, 'results/', prefix, 'plots/', prefix, 'maps_timesteps/', \n");
        fprintf(f_rmap, "            prefix, mstring, fill, '.png', sep = '')\n");
        fprintf(f_rmap, "    } else if (mstring == 'basechange') {\n");
        fprintf(f_rmap, "        mapplot = paste(wkdir, prefix, 'results/', prefix, 'plots/', prefix, mstring, '_fin.png', \n");
        fprintf(f_rmap, "            sep = '')\n");
        fprintf(f_rmap, "    } else if (mstring == 'htsun') {\n");
        fprintf(f_rmap, "        mapplot = paste(wkdir, prefix, 'results/', prefix, 'plots/', prefix, mstring, '_max.png', \n");
        fprintf(f_rmap, "            sep = '')\n");
        fprintf(f_rmap, "    } else if (mstring == 'treach') {\n");
        fprintf(f_rmap, "        mapplot = paste(wkdir, prefix, 'results/', prefix, 'plots/', prefix, mstring, '.png', \n");
        fprintf(f_rmap, "            sep = '')\n");
        fprintf(f_rmap, "    } else {\n");
        fprintf(f_rmap, "        mapplot = paste(wkdir, prefix, 'results/', prefix, 'plots/', prefix, mstring, '_max.png', \n");
        fprintf(f_rmap, "            sep = '')\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "    png(filename = mapplot, width = dispx, height = dispy, units = 'cm', res = 300)\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    # Defining margins\n");
        fprintf(f_rmap, "    par(mar = c(2, 1.3, 0.5, 2.5))\n");
        fprintf(f_rmap, "    par(oma = c(0, 0, 0, 1.5))\n");
        fprintf(f_rmap, "    clip(xmin, xmax, ymin, ymax)  #constraining drawing area\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    if (ydiff > xdiff * 1.2) {\n");
        fprintf(f_rmap, "        # shrink factor and scaling of labels for colour bar legend\n");
        fprintf(f_rmap, "        lshrink <- 0.7\n");
        fprintf(f_rmap, "        lcex <- 1\n");
        fprintf(f_rmap, "    } else if (ydiff > xdiff/1.2) {\n");
        fprintf(f_rmap, "        lshrink <- 0.85\n");
        fprintf(f_rmap, "        lcex <- 1\n");
        fprintf(f_rmap, "    } else {\n");
        fprintf(f_rmap, "        lshrink <- 1\n");
        fprintf(f_rmap, "        lcex <- 0.8\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    # Plotting raster layers\n");
        fprintf(f_rmap, "    if (ortho1 != 'None') {\n");
        fprintf(f_rmap, "        plot(orthogray, legend.width = 1, col = gray((max(0, cellStats(orthogray, 'min')):min(255, \n");
        fprintf(f_rmap, "            cellStats(orthogray, 'max')))/255), legend = FALSE, axes = FALSE, box = FALSE, \n");
        fprintf(f_rmap, "            xlab = NA, ylab = NA, xlim = c(xmin, xmax), ylim = c(ymin, ymax), useRaster = TRUE)  #grayscale orthophoto\n");
        fprintf(f_rmap, "    } else {\n");
        fprintf(f_rmap, "        plot(hillshade, legend.width = 1, col = gray((max(0, cellStats(hillshade, 'min')):min(255, \n");
        fprintf(f_rmap, "            cellStats(hillshade, 'max')))/255), legend = FALSE, axes = FALSE, box = FALSE, \n");
        fprintf(f_rmap, "            xlab = NA, ylab = NA, xlim = c(xmin, xmax), ylim = c(ymin, ymax), useRaster = TRUE)  #hillshade\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    par(new = TRUE)\n");
        fprintf(f_rmap, "    par(cex = lcex)\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    if (model < 7 || fill == 'iscore' || mstring == 'basechange' || mstring == 'htsun' || mstring == 'treach' ) {\n");
        fprintf(f_rmap, "        plot(rastplot, legend.width = 1, useRaster = TRUE, col = rcolours, breaks = rthresholds, \n");
        fprintf(f_rmap, "            xlim = c(xmin, xmax), ylim = c(ymin, ymax), axes = FALSE, xlab = NA, ylab = NA, \n");
        fprintf(f_rmap, "            box = FALSE, legend.args = list(text = mlabel, side = 4, line = -2.1, cex = lcex), \n");
        fprintf(f_rmap, "            axis.args = list(labels = rlabels), legend.shrink = lshrink)  #flow parameter\n");
        fprintf(f_rmap, "    } else {\n");
        fprintf(f_rmap, "        plot(rastplot, legend.width = 1, useRaster = TRUE, col = rcolours, breaks = rthresholds, \n");
        fprintf(f_rmap, "            xlim = c(xmin, xmax), ylim = c(ymin, ymax), axes = FALSE, xlab = NA, ylab = NA, \n");
        fprintf(f_rmap, "            box = FALSE, legend = FALSE)  #flow parameter\n");     
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    # Plotting vector layers\n");
        fprintf(f_rmap, "    if (model <= 3 && releasemass > 0 && mult == 0 ) {\n");
        fprintf(f_rmap, "        contour(x = seq(xmin0, xmax0, length.out = nrow(startshape)), y = seq(ymin0, ymax0, length.out = ncol(startshape)), z = startshape, \n");
        fprintf(f_rmap, "          levels=c(thrsm), drawlabels=FALSE, col = 'red', lty = 3, lwd = 1.5, xlim = c(xmin, xmax), \n");
        fprintf(f_rmap, "          ylim = c(ymin, ymax), axes = F, xlab = NA, ylab = NA, add = TRUE) # release area\n");
        fprintf(f_rmap, "    } else if (model == 7 && mult == 0 ) {\n");
        fprintf(f_rmap, "        if (releasemass > 3) {\n");
        fprintf(f_rmap, "            contour(x = seq(xmin0, xmax0, length.out = nrow(startshape3)), y = seq(ymin0, ymax0, length.out = ncol(startshape3)), z = startshape3, \n");
        fprintf(f_rmap, "              levels=c(thrsm), drawlabels=FALSE, col = 'blue', lty = 3, lwd = 1.5, xlim = c(xmin, xmax), \n");
        fprintf(f_rmap, "              ylim = c(ymin, ymax), axes = F, xlab = NA, ylab = NA, add = TRUE) # PHASE 3 release area\n");
        fprintf(f_rmap, "        }\n");
        fprintf(f_rmap, "        if (releasemass == 2 || releasemass == 3 || releasemass == 7) {\n");
        fprintf(f_rmap, "            contour(x = seq(xmin0, xmax0, length.out = nrow(startshape2)), y = seq(ymin0, ymax0, length.out = ncol(startshape2)), z = startshape2, \n");
        fprintf(f_rmap, "              levels=c(thrsm), drawlabels=FALSE, col = 'green', lty = 3, lwd = 1.5, xlim = c(xmin, xmax), \n");
        fprintf(f_rmap, "              ylim = c(ymin, ymax), axes = F, xlab = NA, ylab = NA, add = TRUE) # PHASE 2 release area\n");
        fprintf(f_rmap, "        }\n");
        fprintf(f_rmap, "        if (releasemass == 1 || releasemass == 3 || releasemass == 5 || releasemass == 7) {\n");
        fprintf(f_rmap, "            contour(x = seq(xmin0, xmax0, length.out = nrow(startshape)), y = seq(ymin0, ymax0, length.out = ncol(startshape)), z = startshape, \n");
        fprintf(f_rmap, "              levels=c(thrsm), drawlabels=FALSE, col = 'red', lty = 3, lwd = 1.5, xlim = c(xmin, xmax), \n");
        fprintf(f_rmap, "              ylim = c(ymin, ymax), axes = F, xlab = NA, ylab = NA, add = TRUE) # PHASE 1 release area\n");
        fprintf(f_rmap, "        }\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    if (j == 0 && impdef == 1 && (fill != 'iscore' || mlabel == 'Impact indicator index')) {\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        contour(x = seq(xmin0, xmax0, length.out = nrow(impactarea)), y = seq(ymin0, ymax0, length.out = ncol(impactarea)), z = impactarea, \n");
        fprintf(f_rmap, "          levels=c(1), drawlabels=FALSE, col = 'red', lty = 1, lwd = 1.5, xlim = c(xmin, xmax), \n");
        fprintf(f_rmap, "          ylim = c(ymin, ymax), axes = F, xlab = NA, ylab = NA, add = TRUE) #observed impact area\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    if ( j == 0 && depdef == 1 && (fill != 'iscore' || mlabel == 'Deposition indicator index')) {\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        contour(x = seq(xmin0, xmax0, length.out = nrow(hdeposit)), y = seq(ymin0, ymax0, length.out = ncol(hdeposit)), z = hdeposit, \n");
        fprintf(f_rmap, "          levels=c(1), drawlabels=FALSE, col = 'orange', lty = 1, lwd = 1.5, xlim = c(xmin, xmax), \n");
        fprintf(f_rmap, "          ylim = c(ymin, ymax), axes = F, xlab = NA, ylab = NA, add = TRUE) #height of observed deposit\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    if (mapprof != 0) {\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        par(new = TRUE)\n");
        fprintf(f_rmap, "        intable <- paste(wkdir, prefix, 'results/', prefix, 'files/', prefix, 'mapprof.txt', sep='')\n");
        fprintf(f_rmap, "        profilex <- scan(intable, nlines = 1, quiet = TRUE)  #profile x coordinate\n");
        fprintf(f_rmap, "        profiley <- scan(intable, skip = 1, nlines = 1, quiet = TRUE)  #profile y coordinate\n");
        fprintf(f_rmap, "        lines(x = profilex, y = profiley, col = 'yellow', lwd = 1.5, lty = 2)  #profile\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    if (max(replace(a2c, is.na(a2c), 0)) != min(replace(a2c, is.na(a2c), 0))) {\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        contour(x = seq(xmin0, xmax0, length.out = nrow(a2c)), y = seq(ymin0, ymax0, length.out = ncol(a2c)), z = a2c, \n");
        fprintf(f_rmap, "          levels=seq(contmin, contmax, by=contstep), drawlabels=FALSE, col = 'black', lty = 1, lwd = 0.5, xlim = c(xmin, xmax), \n");
        fprintf(f_rmap, "          ylim = c(ymin, ymax), axes = F, xlab = NA, ylab = NA, add = TRUE) #contour lines\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    if (velocity == 1) {\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        if (model <= 3) {\n");
        fprintf(f_rmap, "            par(new = TRUE)\n");
        fprintf(f_rmap, "            arrows(x0 = dirx, y0 = diry, x1 = dirx + vunit * dirvy, y1 = diry - vunit * \n");
        fprintf(f_rmap, "                dirvx, length = 0.025, code = 2, angle = 45, col = 'red', lwd = 0.5 * \n");
        fprintf(f_rmap, "                dirh)  #directions\n");
        fprintf(f_rmap, "        } else if (model == 7) {\n");
        fprintf(f_rmap, "            par(new = TRUE)\n");
        fprintf(f_rmap, "            arrows(x0 = dirxw, y0 = diryw, x1 = dirxw + vunit * dirvyw, y1 = diryw - \n");
        fprintf(f_rmap, "                vunit * dirvxw, length = 0.025, code = 2, angle = 45, col = rgb(0, 0, \n");
        fprintf(f_rmap, "                1), lwd = 0.5 * dirhw)  #PHASE 3 directions\n");
        fprintf(f_rmap, "            par(new = TRUE)\n");
        fprintf(f_rmap, "            arrows(x0 = dirxf, y0 = diryf, x1 = dirxf + vunit * dirvyf, y1 = diryf - \n");
        fprintf(f_rmap, "                vunit * dirvxf, length = 0.025, code = 2, angle = 45, col = rgb(0, 1, \n");
        fprintf(f_rmap, "                0), lwd = 0.5 * dirhf)  #PHASE 2 directions\n");
        fprintf(f_rmap, "            par(new = TRUE)\n");
        fprintf(f_rmap, "            arrows(x0 = dirx, y0 = diry, x1 = dirx + vunit * dirvy, y1 = diry - vunit * \n");
        fprintf(f_rmap, "                dirvx, length = 0.025, code = 2, angle = 45, col = rgb(1, 0, 0), lwd = 0.5 * \n");
        fprintf(f_rmap, "                dirh)  #PHASE 1 directions\n");
        fprintf(f_rmap, "        }\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    # Plotting control points\n");
        fprintf(f_rmap, "    if (ctrlpts > 0) {\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        par(new = TRUE)\n");
        fprintf(f_rmap, "        if ( fill == 'iscore' ) {\n");
        fprintf(f_rmap, "            inctrl <- paste(wkdir, prefix, 'results/', prefix, 'files/', prefix, 'ctrlpoints1.txt', \n");
        fprintf(f_rmap, "                sep = '')\n");
        fprintf(f_rmap, "        } else {\n");
        fprintf(f_rmap, "            inctrl <- paste(wkdir, prefix, 'results/', prefix, 'files/', prefix, 'ctrlpoints.txt', \n");
        fprintf(f_rmap, "                sep = '')\n");
        fprintf(f_rmap, "        }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        jz <- 1  #initializing counter for loop over all control points\n");
        fprintf(f_rmap, "        repeat {\n");
        fprintf(f_rmap, "            # starting loop over all control points\n");
        fprintf(f_rmap, "            if (jz > ctrlpts) {\n");
        fprintf(f_rmap, "                break  #break if last control point was reached\n");
        fprintf(f_rmap, "            } else {\n");
        fprintf(f_rmap, "                ctrlx <- read.table(inctrl, header = FALSE)[(jz - 1) * (ntimemax + 2) + \n");
        fprintf(f_rmap, "                    2, 3]  #control point x coordinate\n");
        fprintf(f_rmap, "                ctrly <- read.table(inctrl, header = FALSE)[(jz - 1) * (ntimemax + 2) + \n");
        fprintf(f_rmap, "                    2, 4]  #control point y coordinate\n");
        fprintf(f_rmap, "                points(x = as.vector(t(ctrlx)), y = as.vector(t(ctrly)), col = 'yellow', \n");
        fprintf(f_rmap, "                    pch = 3)  #plotting control point\n");
        fprintf(f_rmap, "                ctrltext <- as.character(read.table(inctrl, header = FALSE)[(jz - 1) * \n");
        fprintf(f_rmap, "                    (ntimemax + 2) + 2, 1])\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "                text(x = as.numeric(as.vector(t(ctrlx))) + 0.035 * (xmax - xmin), y = as.numeric(as.vector(t(ctrly))) + \n");
        fprintf(f_rmap, "                    0.035 * (xmax - xmin), labels = ctrltext, cex = 1, col = 'yellow', \n");
        fprintf(f_rmap, "                    font = 1)  #control point label\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "                jz <- jz + 1\n");
        fprintf(f_rmap, "            }\n");
        fprintf(f_rmap, "        }\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    # Plotting hydrograph profiles\n");
        fprintf(f_rmap, "    if (hydrograph > 0) {\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        par(new = TRUE)\n");
        fprintf(f_rmap, "        inhyd <- paste(wkdir, prefix, 'results/', prefix, 'files/', prefix, 'hydprofiles.txt', \n");
        fprintf(f_rmap, "            sep = '')\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        jz <- 1  #initializing counter for loop over all hydrographs\n");
        fprintf(f_rmap, "        repeat {\n");
        fprintf(f_rmap, "            # starting loop over all hydrographs\n");
        fprintf(f_rmap, "            if (jz > (ninhyd + nouthyd)) {\n");
        fprintf(f_rmap, "                break  #break if last hydrograph was reached\n");
        fprintf(f_rmap, "            } else {\n");
        fprintf(f_rmap, "                if (jz <= ninhyd) {\n");
        fprintf(f_rmap, "                    # color of hydrograph profile depending on type:\n");
        fprintf(f_rmap, "                    hydcol <- 'green'\n");
        fprintf(f_rmap, "                } else {\n");
        fprintf(f_rmap, "                    hydcol <- 'purple'\n");
        fprintf(f_rmap, "                }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "                hydx <- read.table(inhyd, header = FALSE)[jz + 1, 2:4]  #hydrograph x coordinates\n");
        fprintf(f_rmap, "                hydy <- read.table(inhyd, header = FALSE)[jz + 1, 5:7]  #hydrograph y coordinates\n");
        fprintf(f_rmap, "                lines(x = as.vector(t(hydx)), y = as.vector(t(hydy)), col = hydcol, lwd = 1.5, \n");
        fprintf(f_rmap, "                    lty = 3)  #profile\n");
        fprintf(f_rmap, "                points(x = as.vector(t(hydx))[2], y = as.vector(t(hydy))[2], col = hydcol, \n");
        fprintf(f_rmap, "                    pch = 19, cex = 0.5)  #centre point\n");
        fprintf(f_rmap, "                hydtext <- as.character(read.table(inhyd, header = FALSE)[jz + 1, 1])\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "                hyddx1 <- as.numeric(as.vector(t(hydx))[1])  #coordinates for hydrograph label:\n");
        fprintf(f_rmap, "                hyddx3 <- as.numeric(as.vector(t(hydx))[3])\n");
        fprintf(f_rmap, "                if (hyddx1 > hyddx3) {\n");
        fprintf(f_rmap, "                    hyddx <- hyddx1\n");
        fprintf(f_rmap, "                    hyddy <- as.numeric(as.vector(t(hydy))[1])\n");
        fprintf(f_rmap, "                } else {\n");
        fprintf(f_rmap, "                    hyddx <- hyddx3\n");
        fprintf(f_rmap, "                    hyddy <- as.numeric(as.vector(t(hydy))[3])\n");
        fprintf(f_rmap, "                }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "                text(x = hyddx, y = hyddy, labels = hydtext, cex = 1, col = hydcol, font = 1, \n");
        fprintf(f_rmap, "                    pos = 4, offset = 0.3)  #hydrograph label\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "                jz <- jz + 1\n");
        fprintf(f_rmap, "            }\n");
        fprintf(f_rmap, "        }\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    # Plotting bounding box and axes\n");
        fprintf(f_rmap, "    box()\n");
        fprintf(f_rmap, "    par(cex = 0.8)\n");
        fprintf(f_rmap, "    axis(side = 1, tck = -0.01, labels = NA)  #x axis\n");
        fprintf(f_rmap, "    axis(side = 2, tck = -0.01, labels = NA)  #y axis\n");
        fprintf(f_rmap, "    axis(side = 1, lwd = 0, line = -0.6)  #x axis\n");
        fprintf(f_rmap, "    axis(side = 2, lwd = 0, line = -0.6)  #y axis\n");
        fprintf(f_rmap, "    par(cex = 1)\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    # Plotting header text\n");
        fprintf(f_rmap, "    htext = vector('expression', 1)\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    #Defining string for time unit\n");
        fprintf(f_rmap, "    if ( slomo == 1 ) {\n");
        fprintf(f_rmap, "        secunit0 <- 's'\n");
        fprintf(f_rmap, "    } else if ( slomo == 60 ) {\n");
        fprintf(f_rmap, "        secunit0 <- 'min'\n");
        fprintf(f_rmap, "    } else if ( slomo == 3600 ) {\n");
        fprintf(f_rmap, "        secunit0 <- 'h'\n");
        fprintf(f_rmap, "    } else if ( slomo == 86400 ) {\n");
        fprintf(f_rmap, "        secunit0 <- 'days'\n");
        fprintf(f_rmap, "    } else if ( slomo == 604800 ) {\n");
        fprintf(f_rmap, "        secunit0 <- 'weeks'\n");
        fprintf(f_rmap, "    } else if ( slomo == 2592000 ) {\n");
        fprintf(f_rmap, "        secunit0 <- 'months'\n");
        fprintf(f_rmap, "    } else if ( slomo == 31536000 ) {\n");
        fprintf(f_rmap, "        secunit0 <- 'years'\n");
        fprintf(f_rmap, "    } else if ( slomo == 315360000 ) {\n");
        fprintf(f_rmap, "        secunit0 <- 'dec'\n");
        fprintf(f_rmap, "    } else {\n");
        fprintf(f_rmap, "        secunit0 <- paste('x', as.character(slomo), 's')\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    if (fill != 'iscore') {\n");
        fprintf(f_rmap, "        # defining x position and label:\n");
        fprintf(f_rmap, "        if (ntimesteps <= ntimemax) {\n");
        fprintf(f_rmap, "            secpassed <- tint * ntimesteps  #time passed since start of simulation\n");
        fprintf(f_rmap, "            if (secpassed > tstop) \n");
        fprintf(f_rmap, "                secpassed <- tstop\n");
        fprintf(f_rmap, "            htext[1] <- substitute(expression(italic(t) == secformat ~ secunit), list(secformat = format(round(secpassed, \n");
        fprintf(f_rmap, "                1), nsmall = 1), secunit = format(secunit0)))[2]\n");
        fprintf(f_rmap, "            posx <- xmin + xdiff/5\n");
        fprintf(f_rmap, "        } else if (mstring == 'hflow') {\n");
        fprintf(f_rmap, "            htext[1] <- 'Maximum flow height'\n");
        fprintf(f_rmap, "            posx <- xmin + xdiff/2\n");
        fprintf(f_rmap, "        } else if (mstring == 'tflow') {\n");
        fprintf(f_rmap, "            htext[1] <- 'Maximum flow kinetic energy'\n");
        fprintf(f_rmap, "            posx <- xmin + xdiff/2\n");
        fprintf(f_rmap, "        } else if (mstring == 'pflow') {\n");
        fprintf(f_rmap, "            htext[1] <- 'Maximum flow pressure'\n");
        fprintf(f_rmap, "            posx <- xmin + xdiff/2\n");
        fprintf(f_rmap, "        } else if (mstring == 'basechange') {\n");
        fprintf(f_rmap, "            htext[1] <- 'Final change of basal topography'\n");
        fprintf(f_rmap, "            posx <- xmin + xdiff/2\n");
        fprintf(f_rmap, "        } else if (mstring == 'htsun') {\n");
        fprintf(f_rmap, "            htext[1] <- 'Maximum tsunami height'\n");
        fprintf(f_rmap, "            posx <- xmin + xdiff/2\n");
        fprintf(f_rmap, "        } else if (mstring == 'treach') {\n");
        fprintf(f_rmap, "            htext[1] <- 'Time of reach'\n");
        fprintf(f_rmap, "            posx <- xmin + xdiff/2\n");
        fprintf(f_rmap, "        }\n");
        fprintf(f_rmap, "    } else if (mstring == 'iii_hflow') {\n");
        fprintf(f_rmap, "        htext[1] <- 'Impact indicator index'\n");
        fprintf(f_rmap, "        posx <- xmin + xdiff/2\n");
        fprintf(f_rmap, "    } else if (mstring == 'iii_tflow') {\n");
        fprintf(f_rmap, "        htext[1] <- 'Impact indicator index'\n");
        fprintf(f_rmap, "        posx <- xmin + xdiff/2\n");
        fprintf(f_rmap, "    } else if (mstring == 'iii_pflow') {\n");
        fprintf(f_rmap, "        htext[1] <- 'Impact indicator index'\n");
        fprintf(f_rmap, "        posx <- xmin + xdiff/2\n");
        fprintf(f_rmap, "    } else if (mstring == 'dii') {\n");
        fprintf(f_rmap, "        htext[1] <- 'Deposition indicator index'\n");
        fprintf(f_rmap, "        posx <- xmin + xdiff/2\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    if (ydiff > xdiff * 1.2) {\n");
        fprintf(f_rmap, "        # y positions and scaling of header text and header legend\n");
        fprintf(f_rmap, "        posy1 <- ymax + ydiff/14\n");
        fprintf(f_rmap, "        posy2 <- ymax + ydiff/10\n");
        fprintf(f_rmap, "        posy3 <- ymax + ydiff/14\n");
        fprintf(f_rmap, "        posy4 <- ymax + ydiff/26\n");
        fprintf(f_rmap, "        lcex <- 1\n");
        fprintf(f_rmap, "    } else if (ydiff > xdiff/1.2) {\n");
        fprintf(f_rmap, "        posy1 <- ymax + ydiff/12\n");
        fprintf(f_rmap, "        posy2 <- ymax + ydiff/8\n");
        fprintf(f_rmap, "        posy3 <- ymax + ydiff/13\n");
        fprintf(f_rmap, "        posy4 <- ymax + ydiff/34\n");
        fprintf(f_rmap, "        lcex <- 1\n");
        fprintf(f_rmap, "    } else {\n");
        fprintf(f_rmap, "        posy1 <- ymax + ydiff/11\n");
        fprintf(f_rmap, "        posy2 <- ymax + ydiff/7\n");
        fprintf(f_rmap, "        posy3 <- ymax + ydiff/11.5\n");
        fprintf(f_rmap, "        posy4 <- ymax + ydiff/27\n");
        fprintf(f_rmap, "        lcex <- 0.8\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    text(x = posx, y = posy1, labels = htext[1], cex = 1.4 * lcex, col = 'black')  #printing text\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    #Defining string for velocity unit\n");
        fprintf(f_rmap, "    if ( slomo == 1 ) {\n");
        fprintf(f_rmap, "        velunit0 <- 'm/s'\n");
        fprintf(f_rmap, "    } else if ( slomo == 60 ) {\n");
        fprintf(f_rmap, "        velunit0 <- 'm/min'\n");
        fprintf(f_rmap, "    } else if ( slomo == 3600 ) {\n");
        fprintf(f_rmap, "        velunit0 <- 'm/h'\n");
        fprintf(f_rmap, "    } else if ( slomo == 86400 ) {\n");
        fprintf(f_rmap, "        velunit0 <- 'm/day'\n");
        fprintf(f_rmap, "    } else if ( slomo == 604800 ) {\n");
        fprintf(f_rmap, "        velunit0 <- 'm/week'\n");
        fprintf(f_rmap, "    } else if ( slomo == 2592000 ) {\n");
        fprintf(f_rmap, "        velunit0 <- 'm/month'\n");
        fprintf(f_rmap, "    } else if ( slomo == 31536000 ) {\n");
        fprintf(f_rmap, "        velunit0 <- 'm/year'\n");
        fprintf(f_rmap, "    } else if ( slomo == 315360000 ) {\n");
        fprintf(f_rmap, "        velunit0 <- 'm/dec'\n");
        fprintf(f_rmap, "    } else {\n");
        fprintf(f_rmap, "        velunit0 <- paste('m/', as.character(slomo), 's')\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    # Plotting header legend (flow velocity)\n");
        fprintf(f_rmap, "    if (fill != 'iscore' && ntimesteps <= ntimemax) {\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        par(new = TRUE)\n");
        fprintf(f_rmap, "        if (model <= 3) {\n");
        fprintf(f_rmap, "            arrows(x0 = xmin + 0.62 * xdiff, y0 = posy3, x1 = xmin + 0.62 * xdiff + 0.05 * \n");
        fprintf(f_rmap, "                min(xdiff, ydiff), y1 = posy3, length = 0.025, code = 2, angle = 45, \n");
        fprintf(f_rmap, "                col = 'red', lwd = 0.5)  #legend for velocity\n");
        fprintf(f_rmap, "        } else if (model > 3 && model < 7) {\n");
        fprintf(f_rmap, "            arrows(x0 = xmin + 0.62 * xdiff, y0 = posy3, x1 = xmin + 0.62 * xdiff + 0.05 * \n");
        fprintf(f_rmap, "                min(xdiff, ydiff), y1 = posy3, length = 0.025, code = 2, angle = 45, \n");
        fprintf(f_rmap, "                col = 'brown', lwd = 0.5)  #legend for PHASE 1 velocity\n");
        fprintf(f_rmap, "            arrows(x0 = xmin + 0.62 * xdiff + 0.05 * min(xdiff, ydiff), y0 = posy4, x1 = xmin + \n");
        fprintf(f_rmap, "                0.62 * xdiff, y1 = posy4, length = 0.025, code = 2, angle = 45, col = 'blue', \n");
        fprintf(f_rmap, "                lwd = 0.5)  #legend for PHASE 2 velocity\n");
        fprintf(f_rmap, "        } else if (model == 7) {\n");
        fprintf(f_rmap, "            arrows(x0 = xmin + 0.58 * xdiff, y0 = posy3, x1 = xmin + 0.58 * xdiff + 0.05 * \n");
        fprintf(f_rmap, "                min(xdiff, ydiff), y1 = posy3, length = 0.025, code = 2, angle = 45, \n");
        fprintf(f_rmap, "                col = rgb(1, 0, 0), lwd = 0.5)  #legend for PHASE 1 velocity\n");
        fprintf(f_rmap, "            arrows(x0 = xmin + 0.58 * xdiff, y0 = posy4, x1 = xmin + 0.58 * xdiff + 0.05 * \n");
        fprintf(f_rmap, "                min(xdiff, ydiff), y1 = posy4, length = 0.025, code = 2, angle = 45, \n");
        fprintf(f_rmap, "                col = rgb(0, 1, 0), lwd = 0.5)  #legend for PHASE 2 velocity\n");
        fprintf(f_rmap, "            arrows(x0 = xmin + 0.68 * xdiff, y0 = (posy3 + posy4)/2, x1 = xmin + 0.68 * \n");
        fprintf(f_rmap, "                xdiff + 0.05 * min(xdiff, ydiff), y1 = (posy3 + posy4)/2, length = 0.025, \n");
        fprintf(f_rmap, "                code = 2, angle = 45, col = rgb(0, 0, 1), lwd = 0.5)  #legend for PHASE 3 velocity\n");
        fprintf(f_rmap, "        }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        textarrow = vector('expression', 1)\n");
        fprintf(f_rmap, "        textarrow[1] = substitute(expression(italic(v)[max] == vformat ~ velunit), list(vformat = format(round(maxvflow, \n");
        fprintf(f_rmap, "            1), nsmall = 1), velunit = format(velunit0)))[2]  #defining label for flow velocity legend\n");
        fprintf(f_rmap, "        text(x = xmax - xdiff/2.15, y = posy2, labels = textarrow[1], col = 'black', \n");
        fprintf(f_rmap, "            pos = 4, cex = lcex)\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        if (model == 7) {\n");
        fprintf(f_rmap, "            text(x = xmin + 0.61 * xdiff, y = posy3, labels = 'P1', col = rgb(1, 0, 0), \n");
        fprintf(f_rmap, "                pos = 4, cex = lcex)  #printing legend text for PHASE 1 velocity\n");
        fprintf(f_rmap, "            text(x = xmin + 0.61 * xdiff, y = posy4, labels = 'P2', col = rgb(0, 1, 0), \n");
        fprintf(f_rmap, "                pos = 4, cex = lcex)  #printing legend text for PHASE 2 velocity\n");
        fprintf(f_rmap, "            text(x = xmin + 0.71 * xdiff, y = (posy3 + posy4)/2, labels = 'P3', col = rgb(0, \n");
        fprintf(f_rmap, "                0, 1), pos = 4, cex = lcex)  #printing legend text for PHASE 3 velocity\n");
        fprintf(f_rmap, "        }\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    # Plotting footer legend (release heights, hydrographs and observed impact area\n");
        fprintf(f_rmap, "    # and deposit) x positions:\n");
        fprintf(f_rmap, "    if (releasemass > 0) {\n");
        fprintf(f_rmap, "        if (model < 7) {\n");
        fprintf(f_rmap, "            pos1 <- 0.035\n");
        fprintf(f_rmap, "        } else {\n");
        fprintf(f_rmap, "            pos1 <- 0.01\n");
        fprintf(f_rmap, "            pos1w <- 0.125\n");
        fprintf(f_rmap, "        }\n");
        fprintf(f_rmap, "        pos2 <- 0.335\n");
        fprintf(f_rmap, "    } else {\n");
        fprintf(f_rmap, "        pos2 <- 0.035\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "    pos3 <- 0.635\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    if (ydiff > xdiff * 1.2) {\n");
        fprintf(f_rmap, "        # y positions and scaling of labels\n");
        fprintf(f_rmap, "        posy1 <- 0.05\n");
        fprintf(f_rmap, "        posy2 <- 0.025\n");
        fprintf(f_rmap, "        posy3 <- 0\n");
        fprintf(f_rmap, "        lcex <- 1\n");
        fprintf(f_rmap, "    } else if (ydiff > xdiff/1.2) {\n");
        fprintf(f_rmap, "        posy1 <- 0.066\n");
        fprintf(f_rmap, "        posy2 <- 0.026\n");
        fprintf(f_rmap, "        posy3 <- -0.014\n");
        fprintf(f_rmap, "        lcex <- 1\n");
        fprintf(f_rmap, "    } else {\n");
        fprintf(f_rmap, "        posy1 <- 0.071\n");
        fprintf(f_rmap, "        posy2 <- 0.029\n");
        fprintf(f_rmap, "        posy3 <- -0.013\n");
        fprintf(f_rmap, "        lcey <- 0.8\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    if (releasemass > 0) {\n");
        fprintf(f_rmap, "        if (model <= 3) {\n");
        fprintf(f_rmap, "            legend('bottomleft', lwd = NA, legend = 'Release area', text.col = 'black', \n");
        fprintf(f_rmap, "                bty = 'n', horiz = TRUE, x.intersp = 0.25, inset = c(pos1 - 0.1, posy1), \n");
        fprintf(f_rmap, "                cex = lcex)\n");
        fprintf(f_rmap, "            legend('bottomleft', lty = 3, lwd = 1.5, col = 'red', legend = '', text.col = 'black', \n");
        fprintf(f_rmap, "                bty = 'n', horiz = TRUE, x.intersp = 0.25, inset = c(pos1, posy2), cex = lcex)  #legend for release area\n");
        fprintf(f_rmap, "        } else if (model == 7) {\n");
        fprintf(f_rmap, "            legend('bottomleft', lwd = NA, legend = 'Release areas', text.col = 'black', \n");
        fprintf(f_rmap, "                bty = 'n', horiz = TRUE, x.intersp = 0.25, inset = c(pos1 - 0.075, posy1), \n");
        fprintf(f_rmap, "                cex = lcex)\n");
        fprintf(f_rmap, "            legend('bottomleft', lty = 3, lwd = 1.5, col = rgb(1, 0, 0), legend = 'P1', \n");
        fprintf(f_rmap, "                text.col = 'red', bty = 'n', horiz = TRUE, x.intersp = 0.25, inset = c(pos1, \n");
        fprintf(f_rmap, "                    posy2), cex = lcex)  #legend for PHASE 1 release area\n");
        fprintf(f_rmap, "            legend('bottomleft', lty = 3, lwd = 1.5, col = rgb(0, 1, 0), legend = 'P2', \n");
        fprintf(f_rmap, "                text.col = 'green', bty = 'n', horiz = TRUE, x.intersp = 0.25, inset = c(pos1w, \n");
        fprintf(f_rmap, "                    (posy2 + posy3)/2), cex = lcex)  #legend for PHASE 2 release area\n");
        fprintf(f_rmap, "            legend('bottomleft', lty = 3, lwd = 1.5, col = rgb(0, 0, 1), legend = 'P3', \n");
        fprintf(f_rmap, "                text.col = 'blue', bty = 'n', horiz = TRUE, x.intersp = 0.25, inset = c(pos1, \n");
        fprintf(f_rmap, "                    posy3), cex = lcex)  #legend for PHASE 3 release area\n");
        fprintf(f_rmap, "        }\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    if (model == 7 && mstring != 'basechange' && mstring != 'htsun' && mstring != 'treach' && fill != 'iscore') {\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        par(xpd = TRUE)\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        legend(x = xmax + xdiff/50, y = ymax - ydiff/25, xjust = 0, title = mlabel, legend = NA, \n");
        fprintf(f_rmap, "            bty = 'n', horiz = FALSE, cex = lcex)  #legend title\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        legend(x = xmax, y = ymax - ydiff/10, xjust = 0, legend = c(ctext), lty = c(1, \n");
        fprintf(f_rmap, "            1, 1, 1, 1, 1), lwd = 0.5, col = c('white', 'black', 'black', 'black', 'black', \n");
        fprintf(f_rmap, "            'black'), text.col = 'black', bty = 'n', horiz = FALSE, cex = lcex)  #main legend\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        legend(x = xmax, y = ymax - ydiff/1.75, xjust = 0, legend = c('Phase 1', 'Phase 2', 'Phase 3'), \n");
        fprintf(f_rmap, "            fill = c(rgb(1, 0, 0, 0.5), rgb(0, 1, 0, 0.5), rgb(0, 0, 1, 0.5)), border = 'black', \n");
        fprintf(f_rmap, "            text.col = 'black', bty = 'n', horiz = FALSE, cex = lcex)  #legend for fraction\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "        legend(x = xmax + xdiff/50, y = ymax - ydiff/9, xjust = 0, title = 'Max', legend = NA, \n");
        fprintf(f_rmap, "            bty = 'n', horiz = FALSE, cex = lcex)  #maximum\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    if (hydrograph > 0) {\n");
        fprintf(f_rmap, "        legend('bottomleft', lwd = NA, legend = 'Hydrographs', text.col = 'black', bty = 'n', \n");
        fprintf(f_rmap, "            horiz = TRUE, x.intersp = 0.25, inset = c(pos2 - 0.1, posy1), cex = lcex)\n");
        fprintf(f_rmap, "        legend('bottomleft', lty = 3, lwd = 1.5, pch = 19, pt.cex = 0.5, col = 'green', \n");
        fprintf(f_rmap, "            legend = 'Input', text.col = 'green', bty = 'n', horiz = TRUE, x.intersp = 0.25, \n");
        fprintf(f_rmap, "            inset = c(pos2, posy2), cex = lcex)  #legend for input hydrograph\n");
        fprintf(f_rmap, "        legend('bottomleft', lty = 3, lwd = 1.5, pch = 19, pt.cex = 0.5, col = 'purple', \n");
        fprintf(f_rmap, "            legend = 'Output', text.col = 'purple', bty = 'n', horiz = TRUE, x.intersp = 0.25, \n");
        fprintf(f_rmap, "            inset = c(pos2, posy3), cex = lcex)  #legend for output hydrograph\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "    if (depdef == 1 || impdef == 1) {\n");
        fprintf(f_rmap, "        legend('bottomleft', lwd = NA, legend = 'Observation', text.col = 'black', bty = 'n', \n");
        fprintf(f_rmap, "            horiz = TRUE, x.intersp = 0.25, inset = c(pos3 - 0.1, posy1), cex = lcex)\n");
        fprintf(f_rmap, "        legend('bottomleft', lty = 1, lwd = 1.5, col = 'red', legend = 'Impact area', \n");
        fprintf(f_rmap, "            text.col = 'red', bty = 'n', horiz = TRUE, x.intersp = 0.25, inset = c(pos3, \n");
        fprintf(f_rmap, "                posy2), cex = lcex)  #legend for observed impact area\n");
        fprintf(f_rmap, "        legend('bottomleft', lty = 1, lwd = 1.5, col = 'orange', legend = 'Deposit', \n");
        fprintf(f_rmap, "            text.col = 'orange', bty = 'n', horiz = TRUE, x.intersp = 0.25, inset = c(pos3, \n");
        fprintf(f_rmap, "                posy3), cex = lcex)  #legend for observed deposit\n");
        fprintf(f_rmap, "    }\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    # Closing plot file\n");
        fprintf(f_rmap, "    dev.off()\n");
        fprintf(f_rmap, "\n");
        fprintf(f_rmap, "    cat(' completed.\\n')\n"); 
        fprintf(f_rmap, "  }\n");
        fprintf(f_rmap, "}\n");
        
        fclose(f_rmap);

        if ( strcmp( rscript, "None" ) != 0 ) {
        
            sprintf(path, "%sr.avaflow.map.cmd", outrplots); // cmd script to start R script
            f_cmap=fopen(path, "w\n");

            fprintf(f_cmap, "\"%s\" \"%s/r.avaflow.map.R\"\n", rscript, outrplots );       
        
            fclose(f_cmap);
        }
    }


// -- STOP --- Writing R script for map plots -------------------------------------------------------------------


// -- START -- Cleaning system ----------------------------------------------------------------------------------


    #ifdef WITHGRASS


        Segment_release(&seg_elev); // releasing segment data (if GRASS is used)
        if ( sico.RELM == 1 ) Segment_release(&seg_hrelease);
        if ( sico.RELV == 1 ) { Segment_release(&seg_vinx); Segment_release(&seg_viny); }
        if ( sico.ENTR == 1 ) Segment_release(&seg_hentrmax);
        if ( sico.MODEL == 7 && sico.RELM2 == 1 ) Segment_release(&seg_hrelease2);
        if ( sico.MODEL == 7 && sico.RELM3 == 1 ) Segment_release(&seg_hrelease3);
        if ( sico.MODEL == 7 && sico.RELV2 == 1 ) { Segment_release(&seg_vinx2); Segment_release(&seg_viny2); }
        if ( sico.MODEL == 7 && sico.RELV3 == 1 ) { Segment_release(&seg_vinx3); Segment_release(&seg_viny3); }
        if ( sico.MODEL == 7 && sico.ENTR2 == 1 ) Segment_release(&seg_hentrmax2);
        if ( sico.MODEL == 7 && sico.ENTR3 == 1 ) Segment_release(&seg_hentrmax3);
        if ( sico.ZONES == 1 ) Segment_release(&seg_zones);
        if ( sico.CENTR == 1 ) Segment_release(&seg_centr);
        if ( sico.PHI == 1 ) Segment_release(&seg_phi);
        if ( sico.PHI2 == 1 ) Segment_release(&seg_phi2);
        if ( sico.PHI3 == 1 ) Segment_release(&seg_phi3);
        if ( sico.DELTA == 1 ) Segment_release(&seg_delta);
        if ( sico.DELTA2 == 1 ) Segment_release(&seg_delta2);
        if ( sico.DELTA3 == 1 ) Segment_release(&seg_delta3);
        if ( sico.ADDFRI == 1 ) Segment_release(&seg_addfri);
        if ( sico.ADDFRI2 == 1 ) Segment_release(&seg_addfri2);
        if ( sico.ADDFRI3 == 1 ) Segment_release(&seg_addfri3);
        if ( sico.COHESION == 1 ) Segment_release(&seg_cohesion);
        if ( sico.COHESION2 == 1 ) Segment_release(&seg_cohesion2);
        if ( sico.COHESION3 == 1 ) Segment_release(&seg_cohesion3);
        if ( sico.NY == 1 ) Segment_release(&seg_ny);
        if ( sico.NY2 == 1 ) Segment_release(&seg_ny2);
        if ( sico.NY3 == 1 ) Segment_release(&seg_ny3);
        if ( sico.CDEFORM == 1 ) Segment_release(&seg_cdeform);
        if ( sico.AMBDRAG == 1 ) Segment_release(&seg_ambdrag);
        if ( sico.TRANS12 == 1 ) Segment_release(&seg_trans12);
        if ( sico.TRANS13 == 1 ) Segment_release(&seg_trans13);
        if ( sico.TRANS23 == 1 ) Segment_release(&seg_trans23);
        if ( sico.ZFRAG == 1 ) Segment_release(&seg_zfrag);
        if ( sico.TRELEASE == 1 ) Segment_release(&seg_trelease);
        if ( sico.TRELSTOP == 1 ) Segment_release(&seg_trelstop);
        if ( sico.STOPTIME == 1 ) Segment_release(&seg_stoptime);
        if ( sico.TSLIDE == 1 ) Segment_release(&seg_tslide);
        if ( sico.IMPACTAREA == 1 ) Segment_release(&seg_impactarea);
        if ( sico.HDEPOSIT == 1 ) Segment_release(&seg_hdeposit);
        if ( sico.PBG == 1 ) { 
            Segment_release(&seg_pbg1);
            Segment_release(&seg_pbg2);
            Segment_release(&seg_pbg3);
        }
                
        free( v ); free_dmatrix3(outv, sico.M, sico.N); // freeing memory


    #endif


    fclose ( f_summary );
    if ( sico.CTRLPOINTS > 0 ) fclose(f_ctrlpoints);
    fclose ( f_volumes );
    if ( hydrograph == 1 ) {

        fclose ( f_hydout );
        for ( i = 0; i < hydnin + hydnout; i++ ) {
            if ( sico.MULT == 0 ) { 

                fprintf(f_hydtrans[i], "%.1f\t0.00\t0.00\t0.00\t0.00\t0.00\t0.00\n", tsum - qtinit[i] + tout );
                fprintf(f_hydtrans[i], "100000.0\t0.00\t0.00\t0.00\t0.00\t0.00\t0.00\n"); // finalizing hydrograph output files in input format
                fclose( f_hydinfo[i] ); fclose( f_hydtrans[i] );
            }
        }
    }

    if ( sico.MULT == 0 ) fclose ( f_directions );
    if ( sico.MULT == 0 && sico.MODEL == 7 ) { fclose ( f_directions2 ); fclose ( f_directions3 ); } // closing files

    free(sico.MAINMAPSET); free(mv); free(mv2); free(mv0); free(madd); free(prefix); free(outmaps); free(outfiles); free(outaimec); free(outvr); free(parapy); free(rscript); free(rlibs);
        free(outrplots); free(proflist0); free(ctrlplist0); free(hydrographslist0), free(hydrocoordslist0); free(in[0]); free(in); // freeing operational arrays

    free(flowpar); free(ib); free(ibasket[0]); free(ibasket); free(icheck[0]); free(icheck); free(cplain); free(cdomain); free(cdomain2); free(cstopped);
        // freeing flow parameter and control arrays

    free(elevname); free(pelev); free(pelev0); free(qelev); free(relev); free(px); free(py); free(dx); free(dy); free(betax); free(betay); free(betaxh); free(betayh); free(betaxy);
    if ( sico.LAYERS0 > 0 ) { free(betax2); free(betay2); free(betaxh2); free(betayh2); free(betax3); free(betay3); free(betaxy2); free(betaxy3); }
    if ( sico.LAYERS0 > 0 || sico.SURFACE > 1 ) { free(betaxh3); free(betayh3); }
        // freeing terrain arrays
    
    free(aw[0]); free(awt[0]); free(af[0]); free(ag[0]); free(as[0]); free(ad[0]); free(asigma_x[0]); free(asigma_y[0]); 
    free(asigma_f[0]); free(asigma_g[0]); free(wintd[0]); free(wintdtest[0]); free(f[0]); free(g[0]); free(s[0]); free(wintelev[0]);
    free(aw); free(awt); free(af); free(ag); free(as); free(ad); free(asigma_x); free(asigma_y); 
    free(asigma_f); free(asigma_g); free(asigma_xelev); free(asigma_yelev); free(wintd); free(wintdtest); free(f); free(g); free(s); free(wintelev); free(wintelevd);
    free_dmatrix3(winta, sico.IMAX, 4); free_dmatrix3(wintb, sico.IMAX, 6); free_dmatrix3(wintc, sico.IMAX, 4); free_dmatrix3(d, sico.IMAX, 6);
        // freeing state variable and numerical scheme arrays

    if ( hydrograph == 1 ) {
        free_dmatrix3(hydhyd, hydtmaxx+1, 5);
        free(hydp[0]), free(hydp);
        free(hydelev); free(hydalpha); free(hydx_metric); free(hydy_metric); free(hydl); free(hydtmax); free(hydi); free(hydx); free(hydy);
        
        for ( hydj=0; hydj<hydnin; hydj++ ) free(hydreleasename[hydj]);
        for ( hydj=0; hydj<hydnin; hydj++ ) free(hydrographslist[hydj]); // freeing hydrograph arrays
        free(hydrocoordslist);
    }
    
    free(hreleasename); 
    free(vinxname); free(vinyname);
    if ( sico.TRELEASE == 1 ) { free( ptrelease ); } 
    if ( sico.TRELSTOP == 1 ) { free( ptrelstop ); }
    free(phrelease); 
    if ( sico.TRELEASE == 1 ) free(qhrelease); 
    free(pvinx); 
    free(pviny);
    
    free(hreleasename2); 
    free(vinxname2); free(vinyname2);
    free(hreleasename3); 
    free(vinxname3); free(vinyname3);
    
    if ( sico.MODEL == 7 ) {
    
        free(phrelease2); 
        if ( sico.TRELEASE == 1 ) free(qhrelease2); 
        free(pvinx2); 
        free(pviny2);
        free(phrelease3); 
        if ( sico.TRELEASE == 1 ) free(qhrelease3); 
        free(pvinx3); 
        free(pviny3); // freeing release arrays
    }

    free(phentrmax);
    free(hentrmaxname);
    free(zonesname);
    free(pzones);
    if ( sico.CENTR == 1 ) { free(pcentr); } 
    if ( sico.MODEL == 7 ) { free(phentrmax2); free(phentrmax3); }
    free(hentrmaxname2); 
    free(hentrmaxname3); // freeing arrays for entrainment
    
    free( pstoptime );
    free( ptslide );
    free( pxslide );
    free(anx);
    free(anu); // freeing arrays for stopping and initial sliding

    if ( sico.IMPACTAREA != 0 ) free(pimpactarea);
    if ( sico.HDEPOSIT != 0 ) free(phdeposit); // freeing arrays for reference data

    free(ppbg1); free(ppbg2); free(ppbg3);
    if ( sico.TSUNAMI != 0 ) { free(htsun); free(htsunmax); }

    free(phiname); free(phi2name); free(phi3name); free(deltaname); free(delta2name); free(delta3name); free(addfriname); free(addfri2name); free(addfri3name); 
    free(cohesionname); free(cohesion2name); free(cohesion3name); free(nyname); free(ny2name); free(ny3name); free(ambdragname); free(trans12name); free(trans13name); free(trans23name); 
    free(zfragname); free(impactareaname); free(hdepositname); free(treleasename); free(trelstopname); free(centrname); free(cdeformname); free(tslidename); free(stoptimename); 
    free(pbg1name); free(pbg2name); free(pbg3name);

    if ( sico.PHI == 1 ) { free(pphi); }  
    if ( sico.PHI2 == 1 ) { free(pphi2); } 
    if ( sico.PHI3 == 1 ) { free(pphi3); } 
    if ( sico.DELTA == 1 ) { free(pdelta); }
    if ( sico.DELTA2 == 1 ) { free(pdelta2); } 
    if ( sico.DELTA3 == 1 ) { free(pdelta3); }
    if ( sico.ADDFRI == 1 ) { free(paddfri); }
    if ( sico.ADDFRI2 == 1 ) { free(paddfri2); } 
    if ( sico.ADDFRI3 == 1 ) { free(paddfri3); }
    if ( sico.COHESION == 1 ) { free(pcohesion); }
    if ( sico.COHESION2 == 1 ) { free(pcohesion2); } 
    if ( sico.COHESION3 == 1 ) { free(pcohesion3); } 
    if ( sico.CDEFORM == 1 ) { free(pcdeform); } // freeing friction arrays
    
    if ( sico.NY == 1 ) { free(pny); } 
    if ( sico.NY2 == 1 ) { free(pny2); } 
    if ( sico.NY3 == 1 ) { free(pny3); } 
    if ( sico.AMBDRAG == 1 ) { free(pambdrag); } // freeing viscosity and ambient drag arrays

    if ( sico.TRANS12 == 1 ) { free(ptrans12); } 
    if ( sico.TRANS13 == 1 ) { free(ptrans13); } 
    if ( sico.TRANS23 == 1 ) { free(ptrans23); } 
    if ( sico.ZFRAG == 1 ) { free(pzfrag); } // freeing transformation and fragmentation arrays
    
    if ( adaptograph == 1 ) { free(adaptoname); free(adaada[0]); free(adaada); } // freeing adaptograph arrays
    if ( frictiograph == 1 ) { free(frictioname); free(frifri[0]); free(frifri); } // freeing frictiograph arrays
    if ( transformograph == 1 ) { free(transformoname); free(tratra[0]); free(tratra); } // freeing transformograph arrays
  
    if ( sico.PROFILE > 0 ) {
    
        free(proflist); free(profelev); free(profdeposit); free_dmatrix3(profdata, sico.M+sico.N, (int)(tmax/tout+3)); free(profeval[0]); free(profeval); 
        free(profx); free(profy); free(profnx); free(profny); free(profdiffabs);
    }

    if ( sico.CTRLPOINTS > 0 ) {
    
        free(ctrlplist); free(ctrlpx); free(ctrlpy); free(ctrlpxm); free(ctrlpym); free_dmatrix3(ctrlpdata, sico.CTRLPOINTS, (int)(tmax/tout+3));
    }
    
    if ( sico.PBG == 1 ) { 

        free(ppbgt1); free(ppbgt2); free(ppbgt3); free(anpx); free(anpy); free(anctr);
    }
    
    free( kappau ); free( path ); free( wkdir ); free( indir );


// -- STOP --- Cleaning system ----------------------------------------------------------------------------------


    time_stop = clock(); // time at end of model execution
    time_elapsed = ( ( float ) ( time_stop - time_start ) ) / CLOCKS_PER_SEC; // time needed for model execution
    
    if ( sico.MULT == 0 ) {
    
        printf("\nModel execution completed in %.2f seconds.\n", time_elapsed);
        fflush(stdout);
    }    

    return 0;
}
