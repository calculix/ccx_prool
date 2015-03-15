
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <ctype.h>

#define     MAX_LINE_LENGTH 256
#define     MAX_INTEGER 2147483647
#define     MAX_FLOAT   1.e32

typedef struct {
  char  model[MAX_LINE_LENGTH]; /* model-name header*/
  char  **uheader; /* user header */
  char  **pheader; /* project header (remark: leading dataset-project-headers are stored in the related dataset!) */
  ITG   u;         /* number of user headers */
  ITG   p;         /* number of project headers */
  ITG   n;         /* number of nodes */
  ITG   e;         /* number of elements  */
  ITG   f;         /* number of faces */
  ITG   g;         /* number of edges */
  ITG   t;         /* number of texts */
  ITG   sets;      /* sets (groups) of entities */
  ITG   mats;      /* materials   */
  ITG   amps;      /* amplitudes  */
  ITG   l;         /* number of loadcases (Datasets) */
  ITG   b;         /* number of nodeBlocks */
  ITG   c;         /* number of 'cuts' over all nodeBlocks (block-to-block interfaces for isaac) */
  ITG   etype[100];/* number of elements of a certain type */
  ITG   nmax;      /* maximum nodenumber */
  ITG   nmin;      /* minimum nodenumber */
  ITG   emax;      /* maximum elemnumber */
  ITG   emin;      /* minimum elemnumber */
  ITG   orignmax;  /* max-node-nr-of-original-nodes (w/o nodes for drawing purposes) */
  ITG   orign;     /* nr-of-original-nodes (w/o nodes for drawing purposes) */
  ITG   olc;       /* nr-of-original-loadcases (w/o cgx generated datasets (lc)) */
  ITG   noffs;     /* node-nr offset, defined with asgn */
  ITG   eoffs;     /* elem-nr offset, defined with asgn */
} Summen;


typedef struct {
  ITG   nr;              /*   external node-nr (node[node-indx].nr) */
  ITG   indx;            /*   node-index (node[ext-node-nr].indx)   */
  char  pflag;           /*   1 if used for display purposes    */
                         /*  -1 if the node is deleted          */
                         /*   0 default                         */
  double nx;             /*   coordinates  node[ext-node-nr].nx */
  double ny;
  double nz;
} Nodes;

typedef struct {
  ITG nr;                /* external element-nr */
  //  ITG indx;              /* -index (elem[external elem-nr].indx)   */
  ITG type;              /* element type (1:Hexa8)  */
  ITG group;
  ITG mat;
  ITG attr;              /* -1: unstructured mesh tr3u (-2 for mesh with libGLu tr3g ) */
                         /*  0: default           */
                         /*  1: reduced integration he8r he20r */
                         /*  2: incompatible modes he8i */
                         /*  3: DASHPOTA be2d */
                         /*  4: plane strain (CPE) tr3e tr6e qu4e qu8e */
                         /*  5: plane stress (CPS) tr3s */
                         /*  6: axisymmetric  (CAX) tr3c */
                         /*  7: fluid he8f */
                         /*  8: tet10m */
  ITG nod[27];
  double **side;         /* side[Nr.][x|y|z]== normal vector */
} Elements;

typedef struct {
  char  **pheader;    /* project header */
  ITG   npheader;              /* number of headers */
  char  **compName;
  char  **icname;
  char  name[MAX_LINE_LENGTH];
  char  dataset_name[MAX_LINE_LENGTH];
  char  dataset_text[MAX_LINE_LENGTH];
  char  analysis_name[MAX_LINE_LENGTH];
  float value;
  char  filename[MAX_LINE_LENGTH];
  FILE *handle;
  fpos_t *fileptr;
  ITG   loaded;       /* if data are stored:1 else: 0 */
  ITG format_flag;
  ITG analysis_type;
  ITG step_number;
  ITG ncomps;
  ITG irtype;
  ITG *menu;
  ITG *ictype;
  ITG *icind1;
  ITG *icind2;
  ITG *iexist;
  float **dat;        /* node related data */
  float ***edat;      /* element related data */
  float *max;         /* maximum datum */
  float *min;         /* minimum datum */
  ITG *nmax;          /* node with maximum datum */
  ITG *nmin;          /* node with minimum datum */
} Datasets;

void freeDatasets(Datasets *lcase, ITG nr);
ITG readfrd(char *datin, Summen *anz, Nodes **nptr, Elements **eptr, Datasets **lptr, ITG read_mode );
ITG readfrdblock( ITG lc, Summen *anz,   Nodes     *node, Datasets *lcase );
ITG stoi(char *string, ITG a, ITG b);
double stof(char *string, ITG a, ITG b);
void stos(char *string, ITG a, ITG b, char *puffer);
ITG compare (char *str1, char *str2, ITG length);
ITG frecord( FILE *handle1,  char *string);
ITG elemChecker(ITG sum_e, Nodes *node, Elements *elem);
void v_result( double *A, double *B, double *C );
void v_prod( double *A, double *B, double *C );
double v_betrag(double *a);
ITG strsplt( char *rec_str, char breakchar, char ***ptr);
