/**
	@file
	mbc.allpole~ - an MSP object shell
	mark cartwright - mcartwright@gmail.com
 
	@ingroup	lpcToolkit
 */


#include "ext.h"
#include "ext_obex.h"
#include "z_dsp.h"

#include <Accelerate/Accelerate.h>

#define MAX_ORDER 200
#define DEFAULT_INTERP 1.0
#define DEFAULT_DEEMPH 0
#define t_floatarg double
#define t_mbcfloat double
#define MBC_VDSP

////////////////////////// object struct
typedef struct {
	t_pxobject x_obj;

    t_mbcfloat*					a_a;			//filer coefficients
    t_mbcfloat* 				a_aBuff;		//filter coeff input buffer
    t_mbcfloat* 				a_y;			//filter memory
    t_mbcfloat* 				a_tempVec;		//temporary buff for vector math
    float 						a_a1;			//deemph coeff
    float 						a_y1;			//deemph filter memory
    t_mbcfloat* 				a_Ar;			//tube areas
    t_mbcfloat**				a_A;			//for parcor to coeff conversion
    t_mbcfloat*					a_K;			//parcor coefficients
    t_mbcfloat*					a_interpCoeff;	//interpolated coefficients
    t_mbcfloat*					a_interpInc;	//interpolation increment
    t_mbcfloat 					a_G;
    double 						a_interp;		//interpolation time in ms
    int 						a_order;
    int 						a_deemph;
    int 						a_coeffType;	//type of coeffecients: CT_PARCOR,CT_FILTER,CT_AREA
    double 						a_fs;			//sampling rate
    int 						a_vsize;		//vector size

} t_allpole;


enum coeffTypes
{
    CT_PARCOR,
    CT_FILTER,
    CT_AREA
};


///////////////////////// function prototypes
//// standard set
void *allpole_new(t_symbol *s, long argc, t_atom *argv);
void allpole_free(t_allpole *x);
void allpole_assist(t_allpole *x, void *b, long m, long a, char *s);

void allpole_float(t_allpole *x, double f);

void allpole_dsp(t_allpole *x, t_signal **sp, short *count);
//t_int *allpole_perf_filter(t_int *w);
//t_int *allpole_perf_area(t_int *w);
t_int *allpole_perf_parcor(t_int *w);
t_int *allpole_perf_parcorI(t_int *w);
//t_int *allpole_perf_filterI(t_int *w);
//t_int *allpole_perf_areaI(t_int *w);

void allpole_dsp64(t_allpole *x, t_object *dsp64, short *count, double samplerate,
                   long maxvectorsize, long flags);
void allpole_perf64_parcor(t_allpole *x, t_object *dsp64, double **ins, long numins,
                       double **outs, long numouts, long sampleframes, long flags, void *userparam);
void allpole_perf64_parcorI(t_allpole *x, t_object *dsp64, double **ins, long numins,
                           double **outs, long numouts, long sampleframes, long flags, void *userparam);

void allpole_interp(t_allpole *x, t_floatarg interp);
void allpole_order(t_allpole *x, int order);
void allpole_deemph(t_allpole *x, int deemph);
void allpole_init(t_allpole *x);
void allpole_free(t_allpole *x);
void allpole_free_arrays(t_allpole *x);
void allpole_clear(t_allpole *x);
//inline
void allpole_highOrdFilter(t_allpole* x, int N, int order, t_double* in, t_double* out);
// inline
void allpole_solveForFiltCoefs(t_allpole* x, int order);
//inline
void allpole_deemphFilter(t_allpole *x, int N, t_double* vec);

//////////////////////// global class pointer variable
void *allpole_class;


int C74_EXPORT main(void) {
	t_class *c;
	

	
    c = class_new("mbc.allpole~", (method)allpole_new, (method)allpole_free, (long)sizeof(t_allpole), 0L, A_GIMME, 0);
    
    //class_addmethod(c, (method)allpole_dsp, "dsp", A_CANT, 0);
    class_addmethod(c, (method)allpole_dsp64, "dsp64", A_CANT, 0);
    class_addmethod(c, (method)allpole_interp,"interp",A_DEFFLOAT,0);
    class_addmethod(c, (method)allpole_order,"order",A_DEFLONG,0);
    class_addmethod(c, (method)allpole_clear,"clear",0);
    class_addmethod(c, (method)allpole_deemph,"deemph",A_DEFLONG,0);
    class_addmethod(c, (method)allpole_assist,"assist",A_CANT,0);
    
    class_dspinit(c);				// new style object version of dsp_initclass();
    class_register(CLASS_BOX, c);	// register class as a box class
    allpole_class = c;
    
    return 0;
}


/*
//perform for PARCOR COEFFICIENTS, INTERPOLATION OFF
t_int *allpole_perf_parcor(t_int *w)
{
    t_float *in = (t_float *)(w[1]);
    t_float *coeffIn = (t_float *)(w[2]);
    t_float *coeffIdxIn = (t_float *)(w[3]);
    t_float *G_n = (t_float *)(w[4]);
    t_allpole *x = (t_allpole *)(w[5]);
    t_float *out = (t_float *)(w[6]);
    int N = (int)(w[7]);
    int order = x->a_order;
    int i, n;
    
    for (n=0; n < N; n++)
    {
        if(IS_NAN_FLOAT(G_n[n])) G_n[n] = 0.0;
        if(IS_NAN_FLOAT(coeffIn[n])) coeffIn[n] = 0.0;
    }
    
    //look at coefficient index, if not zeros, buffer in coefficients
    for (n=0; n < N; n++)
    {
        if ((int)(coeffIdxIn[n]) > 0)
        {
            x->a_aBuff[(int)(coeffIdxIn[n])-1] = (t_mbcfloat)coeffIn[n];
            if ((int)(coeffIdxIn[n]) == order)
            {
                for(i = 0; i < order; i++)
                {
                    x->a_K[i+1] = x->a_aBuff[i];
                }
                x->a_G = (t_mbcfloat)(G_n[n]);
                
                allpole_solveForFiltCoefs(x, order);
            }
        }
    }
    
    allpole_highOrdFilter(x, N, order, in, out);
    
    //deemphasis?
    if (x->a_deemph) 
    {
        allpole_deemphFilter(x, N, out);
    }
    
    return (w+8);
}
*/

//perform for PARCOR COEFFICIENTS, INTERPOLATION OFF
void allpole_perf64_parcor(t_allpole *x, t_object *dsp64, double **ins, long numins,
                           double **outs, long numouts, long sampleframes, long flags, void *userparam)
{
    t_double *in = ins[0];
    t_double *coeffIn = ins[1];
    t_double *coeffIdxIn = ins[2];
    t_double *G_n = ins[3];
    t_double *out = outs[0];
    
    int N = (int)sampleframes;
    int order = x->a_order;
    int i, n;
    
    for (n=0; n < N; n++)
    {
        if(IS_NAN_DOUBLE(G_n[n])) G_n[n] = 0.0;
        if(IS_NAN_DOUBLE(coeffIn[n])) coeffIn[n] = 0.0;
    }
    
    //look at coefficient index, if not zeros, buffer in coefficients
    for (n=0; n < N; n++)
    {
        if ((int)(coeffIdxIn[n]) > 0)
        {
            x->a_aBuff[(int)(coeffIdxIn[n])-1] = (t_mbcfloat)coeffIn[n];
            if ((int)(coeffIdxIn[n]) == order)
            {
                for(i = 0; i < order; i++)
                {
                    x->a_K[i+1] = x->a_aBuff[i];
                }
                x->a_G = (t_mbcfloat)(G_n[n]);
                
                allpole_solveForFiltCoefs(x, order);
            }
        }
    }
    
    allpole_highOrdFilter(x, N, order, in, out);
    
    //deemphasis?
    if (x->a_deemph)
    {
        allpole_deemphFilter(x, N, out);
    }

}


void allpole_perf64_parcorI(t_allpole *x, t_object *dsp64, double **ins, long numins,
                           double **outs, long numouts, long sampleframes, long flags, void *userparam)
{
    t_double *in = ins[0];
    t_double *coeffIn = ins[1];
    t_double *coeffIdxIn = ins[2];
    t_double *G_n = ins[3];

    t_double *out = outs[0];
    int N = (int)sampleframes;
    int order = x->a_order;
    int i, n;
    float interpDiv; //interpolation denominator
    
    for (n=0; n < N; n++)
    {
        if(IS_NAN_DOUBLE(G_n[n])) G_n[n] = 0.0;
        if(IS_NAN_DOUBLE(coeffIn[n])) coeffIn[n] = 0.0;
    }
    
    //look at coefficient index, if not zeros, buffer in coefficients
    for (n=0; n < N; n++)
    {
        if ((int)(coeffIdxIn[n]) > 0)
        {
            x->a_aBuff[(int)(coeffIdxIn[n])-1] = (t_mbcfloat)coeffIn[n];
            if ((int)(coeffIdxIn[n]) == order)
            {
                interpDiv = x->a_interp * 0.001 * x->a_fs / x->a_vsize;
                for(i = 0; i < order; i++)
                {
                    x->a_K[i+1] = x->a_aBuff[i];
                    //determine interpolation parameters
                    x->a_interpInc[i+1] = (x->a_K[i+1] - x->a_interpCoeff[i+1]) / interpDiv;
                }
                
                x->a_G = (t_mbcfloat)(G_n[n]);
            }
        }
    }
    
    //interpolate coefficients
    for (i=1; i <= order; i++)
    {
        if (fabs(x->a_interpCoeff[i] - x->a_K[i]) > fabs(x->a_interpInc[i]))
        {
            x->a_interpCoeff[i] += x->a_interpInc[i];
        } 
        else 
        {
            x->a_interpCoeff[i] = x->a_K[i];
        }
    }
    
    allpole_solveForFiltCoefs(x, order);
    
    allpole_highOrdFilter(x, N, order, in, out);
    
    //deemphasis?
    if (x->a_deemph) 
    {
        allpole_deemphFilter(x, N, out);
    }	
    
}

/*
t_int *allpole_perf_parcorI(t_int *w)
{
    t_float *in = (t_float *)(w[1]);
    t_float *coeffIn = (t_float *)(w[2]);
    t_float *coeffIdxIn = (t_float *)(w[3]);
    t_float *G_n = (t_float *)(w[4]);
    t_allpole *x = (t_allpole *)(w[5]);
    t_float *out = (t_float *)(w[6]);
    int N = (int)(w[7]);
    int order = x->a_order;
    int i, n;
    float interpDiv; //interpolation denominator
    
    for (n=0; n < N; n++)
    {
        if(IS_NAN_FLOAT(G_n[n])) G_n[n] = 0.0;
        if(IS_NAN_FLOAT(coeffIn[n])) coeffIn[n] = 0.0;
    }
    
    //look at coefficient index, if not zeros, buffer in coefficients
    for (n=0; n < N; n++)
    {
        if ((int)(coeffIdxIn[n]) > 0)
        {
            x->a_aBuff[(int)(coeffIdxIn[n])-1] = (t_mbcfloat)coeffIn[n];
            if ((int)(coeffIdxIn[n]) == order)
            {
                interpDiv = x->a_interp * 0.001 * x->a_fs / x->a_vsize;
                for(i = 0; i < order; i++)
                {
                    x->a_K[i+1] = x->a_aBuff[i];
                    //determine interpolation parameters
                    x->a_interpInc[i+1] = (x->a_K[i+1] - x->a_interpCoeff[i+1]) / interpDiv;
                }
                
                x->a_G = (t_mbcfloat)(G_n[n]);
            }
        }
    }
    
    //interpolate coefficients
    for (i=1; i <= order; i++)
    {
        if (fabs(x->a_interpCoeff[i] - x->a_K[i]) > fabs(x->a_interpInc[i]))
        {
            x->a_interpCoeff[i] += x->a_interpInc[i];
        }
        else
        {
            x->a_interpCoeff[i] = x->a_K[i];
        }
    }
    
    allpole_solveForFiltCoefs(x, order);
    
    allpole_highOrdFilter(x, N, order, in, out);
    
    //deemphasis?
    if (x->a_deemph)
    {
        allpole_deemphFilter(x, N, out);
    }
    
    
    return (w+8);
}
*/


void allpole_init(t_allpole *x) {
    int i;
    /*
    x->a_a = (t_mbcfloat *) getbytes16( MAX_ORDER * sizeof(t_mbcfloat));
    x->a_aBuff = (t_mbcfloat *) getbytes( MAX_ORDER * sizeof(t_mbcfloat));
    x->a_y = (t_mbcfloat *) getbytes16( MAX_ORDER * sizeof(t_mbcfloat));
    x->a_tempVec = (t_mbcfloat *) getbytes16( MAX_ORDER * sizeof(t_mbcfloat));
    x->a_Ar = (t_mbcfloat *) getbytes( MAX_ORDER * sizeof(t_mbcfloat));
    x->a_K = (t_mbcfloat *) getbytes( (MAX_ORDER + 1) * sizeof(t_mbcfloat));
    x->a_interpCoeff = (t_mbcfloat *) getbytes( (MAX_ORDER + 1) * sizeof(t_mbcfloat));
    x->a_interpInc = (t_mbcfloat *) getbytes( (MAX_ORDER + 1) * sizeof(t_mbcfloat));
    x->a_A = (t_mbcfloat **) getbytes( (MAX_ORDER + 1) * sizeof(t_mbcfloat*));
    for(i=0; i<MAX_ORDER; i++) {
        x->a_A[i] = (t_mbcfloat *)getbytes( (MAX_ORDER + 1) * sizeof(t_mbcfloat));
    }
    */
    // TODO: memory alignment still an issue?
    x->a_a = (t_mbcfloat *) sysmem_newptr( MAX_ORDER * sizeof(t_mbcfloat));
    x->a_aBuff = (t_mbcfloat *) sysmem_newptr( MAX_ORDER * sizeof(t_mbcfloat));
    x->a_y = (t_mbcfloat *) sysmem_newptr( MAX_ORDER * sizeof(t_mbcfloat));
    x->a_tempVec = (t_mbcfloat *) sysmem_newptr( MAX_ORDER * sizeof(t_mbcfloat));
    x->a_Ar = (t_mbcfloat *) sysmem_newptr( MAX_ORDER * sizeof(t_mbcfloat));
    x->a_K = (t_mbcfloat *) sysmem_newptr( (MAX_ORDER + 1) * sizeof(t_mbcfloat));
    x->a_interpCoeff = (t_mbcfloat *) sysmem_newptr( (MAX_ORDER + 1) * sizeof(t_mbcfloat));
    x->a_interpInc = (t_mbcfloat *) sysmem_newptr( (MAX_ORDER + 1) * sizeof(t_mbcfloat));
    x->a_A = (t_mbcfloat **) sysmem_newptr( (MAX_ORDER + 1) * sizeof(t_mbcfloat*));
    for(i=0; i<MAX_ORDER; i++) {
        x->a_A[i] = (t_mbcfloat *)sysmem_newptr( (MAX_ORDER + 1) * sizeof(t_mbcfloat));
    }
    
    //allpole_clear(x);
    
    x->a_a1 = 0.98;
}

void allpole_free(t_allpole *x) {
    dsp_free((t_pxobject *) x);
    allpole_free_arrays(x);
}

void allpole_free_arrays(t_allpole *x)
{
    int i;
    if (x->a_a) {
        sysmem_freeptr((char *)x->a_a);
        x->a_a = NULL;
    }
    if (x->a_aBuff) {
        sysmem_freeptr(x->a_aBuff);
        x->a_aBuff = NULL;
    }
    if (x->a_y) {
        sysmem_freeptr((char *)x->a_y);
        x->a_y = NULL;
    }
    if (x->a_tempVec) {
        sysmem_freeptr((char *)x->a_tempVec);
        x->a_tempVec = NULL;
    }
    if (x->a_Ar) {
        sysmem_freeptr(x->a_Ar);
        x->a_Ar = NULL;
    }
    if (x->a_K) {
        sysmem_freeptr(x->a_K);
        x->a_K = NULL;
    }
    if (x->a_interpCoeff) {
        sysmem_freeptr(x->a_interpCoeff);
        x->a_interpCoeff = NULL;
    }
    if (x->a_interpInc) {
        sysmem_freeptr(x->a_interpInc);
        x->a_interpInc = NULL;
    }
    if (x->a_A) {
        for(i=0; i<MAX_ORDER; i++) {
            if (x->a_A[i]) {
                sysmem_freeptr(x->a_A[i]);
                x->a_A[i] = NULL;
            }
        }
        sysmem_freeptr(x->a_A);
        x->a_A = NULL;
    }
    /*
    if (x->a_a) {
        freebytes16((char *)x->a_a, MAX_ORDER * sizeof(t_mbcfloat));
        x->a_a = NULL;
    }
    if (x->a_aBuff) {
        freebytes(x->a_aBuff, MAX_ORDER * sizeof(t_mbcfloat));
        x->a_aBuff = NULL;
    }
    if (x->a_y) {
        freebytes16((char *)x->a_y, MAX_ORDER * sizeof(t_mbcfloat));
        x->a_y = NULL;
    }
    if (x->a_tempVec) {
        freebytes16((char *)x->a_tempVec, MAX_ORDER * sizeof(t_mbcfloat));
        x->a_tempVec = NULL;
    }
    if (x->a_Ar) {
        freebytes(x->a_Ar, MAX_ORDER * sizeof(t_mbcfloat));
        x->a_Ar = NULL;
    }
    if (x->a_K) {
        freebytes(x->a_K, (MAX_ORDER + 1) * sizeof(t_mbcfloat));
        x->a_K = NULL;
    }
    if (x->a_interpCoeff) {
        freebytes(x->a_interpCoeff, (MAX_ORDER + 1) * sizeof(t_mbcfloat));
        x->a_interpCoeff = NULL;
    }
    if (x->a_interpInc) {
        freebytes(x->a_interpInc, (MAX_ORDER + 1) * sizeof(t_mbcfloat));
        x->a_interpInc = NULL;
    }
    if (x->a_A) {
        for(i=0; i<MAX_ORDER; i++) {
            if (x->a_A[i]) {
                freebytes(x->a_A[i], (MAX_ORDER + 1) * sizeof(t_mbcfloat));
                x->a_A[i] = NULL;
            }
        }
        freebytes(x->a_A, MAX_ORDER * sizeof(t_mbcfloat*));
        x->a_A = NULL;
    }*/
}

void allpole_clear(t_allpole *x) {
    
    int i, j;
    for(i = 0; i < MAX_ORDER; i++) {
        x->a_a[i] = 0.0;
        x->a_aBuff[i] = 0.0;
        x->a_y[i] = 0.0;
        x->a_Ar[i] = 0.0;
        x->a_tempVec[i] = 0.0;
        for (j = 0; j < (MAX_ORDER+1); j++) {
            x->a_A[i][j] = 0.0;
        }
    }
    
    for(i = 0; i < (MAX_ORDER+1); i++) {
        x->a_K[i] = 0.0;
        x->a_interpCoeff[i] = 0.0;
        x->a_interpInc[i] = 0.0;
    }
    
    x->a_G = 0.0;
    x->a_y1 = 0.0;
}


// this function is called when the DAC is enabled, and "registers" a function
// for the signal chain. in this case, "allpole_perform"
/*
// 32-bit dsp method
void allpole_dsp(t_allpole *x, t_signal **sp, short *count) {
    
    //NOTE: need to specify parcor, filter, or area perform string !!!!!!
    x->a_fs = sp[0]->s_sr;
    x->a_vsize = sp[0]->s_n;
    
    allpole_clear(x);
    
    if (x->a_interp > 0.0) {
        switch (x->a_coeffType) {
            case CT_FILTER:
                //dsp_add(allpole_perf_filterI, 7, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, x, sp[4]->s_vec, sp[0]->s_n);
                //break;
                
            case CT_PARCOR:
                dsp_add(allpole_perf_parcorI, 7, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, x, sp[4]->s_vec, sp[0]->s_n);
                break;
                
            case CT_AREA:
                //dsp_add(allpole_perf_areaI, 7, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, x, sp[4]->s_vec, sp[0]->s_n);
                //break;
                
            default:
                error("mbc.allpole~: no coefficient type selected");
                break;
        }
    } else {
        switch (x->a_coeffType) {
            case CT_FILTER:
                //dsp_add(allpole_perf_filter, 7, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, x, sp[4]->s_vec, sp[0]->s_n);
                //break;
                
            case CT_PARCOR:
                dsp_add(allpole_perf_parcor, 7, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, x, sp[4]->s_vec, sp[0]->s_n);
                break;
                
            case CT_AREA:
                //dsp_add(allpole_perf_area, 7, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, x, sp[4]->s_vec, sp[0]->s_n);
                //break;
                
            default:
                error("mbc.allpole~: no coefficient type selected");
                break;
        }
    }	

}
*/


//64-bit dsp method
void allpole_dsp64(t_allpole *x, t_object *dsp64, short *count, double samplerate, 
				 long maxvectorsize, long flags) {
    
    //NOTE: need to specify parcor, filter, or area perform string !!!!!!
    x->a_fs = samplerate;
    x->a_vsize = (int)maxvectorsize;
    
    allpole_clear(x);
    
//	object_method(dsp64, gensym("dsp_add64"), x, allpole_perform64, 0, NULL);
	
    if (x->a_interp > 0.0) {
        switch (x->a_coeffType) {
            case CT_FILTER:
                //dsp_add(allpole_perf_filterI, 7, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, x, sp[4]->s_vec, sp[0]->s_n);
                //break;
                
            case CT_PARCOR:
                //dsp_add(allpole_perf_parcorI, 7, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, x, sp[4]->s_vec, sp[0]->s_n);
                object_method(dsp64, gensym("dsp_add64"), x, allpole_perf64_parcorI, 0, NULL);
                break;
                
            case CT_AREA:
                //dsp_add(allpole_perf_areaI, 7, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, x, sp[4]->s_vec, sp[0]->s_n);
                //break;
                
            default:
                error("mbc.allpole~: no coefficient type selected");
                break;
        }
    } else {
        switch (x->a_coeffType) {
            case CT_FILTER:
                //dsp_add(allpole_perf_filter, 7, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, x, sp[4]->s_vec, sp[0]->s_n);
                //break;
                
            case CT_PARCOR:
                //dsp_add(allpole_perf_parcor, 7, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, x, sp[4]->s_vec, sp[0]->s_n);
                object_method(dsp64, gensym("dsp_add64"), x, allpole_perf64_parcor, 0, NULL);
                break;
                
            case CT_AREA:
                //dsp_add(allpole_perf_area, 7, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, x, sp[4]->s_vec, sp[0]->s_n);
                //break;
                
            default:
                error("mbc.allpole~: no coefficient type selected");
                break;
        }
    }


}




// 64 bit signal input version
void allpole_perform64(t_allpole *x, t_object *dsp64, double **ins, long numins, 
					   double **outs, long numouts, long sampleframes, long flags, void *userparam){
	
	t_double *in = ins[0];
	t_double *out = outs[0];	
	long vs = sampleframes;

	
	if (x->x_obj.z_disabled)
		return;


	while(vs--) {
		*out++ = (*in++);
	}

	return;
	
}



void allpole_interp(t_allpole *x, t_floatarg interp) {
    
    x->a_interp = (float)interp;
    int i;
    int order = x->a_order;
    float interpDiv;
    
    switch (x->a_coeffType) {
        case CT_FILTER:
            error("mbc.allpole~: interpolation is inactive with coefficient type: filter");
            break;
            
        case CT_PARCOR:
            interpDiv = x->a_interp * 0.001 * x->a_fs / x->a_vsize;
            for(i = 0; i < order; i++) {
                x->a_interpInc[i+1] = (x->a_K[i+1] - x->a_interpCoeff[i+1]) / interpDiv;
            }
            break;
            
        case CT_AREA:
            interpDiv = x->a_interp * 0.001 * x->a_fs / x->a_vsize;
            for(i = 0; i < order; i++) {
                x->a_interpInc[i+1] = (x->a_Ar[i+1] - x->a_interpCoeff[i+1]) / interpDiv;
            }
            break;
            
        default:
            break;
    }
}

void allpole_order(t_allpole *x, int order) {
    if (order < 1) {
        error("mbc.allpole~: mbc.allpole~ needs a positive integer value for order");
        order = 1;
    } else if (order > 199) {
        error("mbc.allpole~: max order is 199");
        order = 199;
    }
    
    int i = x->a_order;
    x->a_order = order;
    for ( i = i+1; i < order; i++) {
        x->a_a[i] = 0.0;
        x->a_aBuff[i] = 0.0;
        x->a_y[i] = 0.0;
        x->a_Ar[i] = 0.0;
    }
}

void allpole_deemph(t_allpole *x, int deemph) {
    x->a_deemph = deemph;
}



void *allpole_new(t_symbol *s, long argc, t_atom *argv)
{
    t_allpole *x = NULL;
    //object_alloc(allpole_class);
    
    if (argc < 1)
    {
        error("mbc.allpole~: must specify filter order (this should match mbc.lpc~ and mbc.errfilt~)");
        return NULL;
    }
    
    x = object_alloc(allpole_class);
    
	if(x) {
        dsp_setup((t_pxobject *)x, 1);	// MSP inlets: arg is # of inlets and is REQUIRED!
        // use 0 if you don't need inlets
        
        dsp_setup((t_pxobject *)x,4);
        outlet_new(x, "signal");
        
        //get arguments out of gimme list
        long order = 0;
        double interp = DEFAULT_INTERP;
        long deemph = DEFAULT_DEEMPH;
        //t_symbol* coeffType = atom_getsymarg(3,argc,argv);
        
        switch (argc)
        {
            case 0:
                break;
                
            case 1:
                order = atom_getintarg(0,argc,argv);
                //order = atom_getlong(argv);
                break;
                
            case 2:
                order = atom_getintarg(0,argc,argv);
                interp = atom_getfloatarg(1,argc,argv);
                break;
                
            case 3:
                order = atom_getintarg(0,argc,argv);
                interp = atom_getfloatarg(1,argc,argv);
                deemph = atom_getintarg(2,argc,argv);
                break;
                
            default:
                order = atom_getintarg(0,argc,argv);
                interp = atom_getfloatarg(1,argc,argv);
                deemph = atom_getintarg(2,argc,argv);
                error("mbc.allpole~: too many arguments");
        }
        
        //order bounds
        if (order < 1) {
            error("mbc.allpole~: mbc.allpole~ needs a positive integer value for order");
            order = 1;
        } else if (order > 199) {
            error("mbc.allpole~: max order is 199");
            order = 199;
        }
        
        //parse filter coeff type
#if 0
        if (coeffType == gensym("parcor")) {
            x->a_coeffType = CT_PARCOR;
        } else if (coeffType == gensym("filter")) {
            x->a_coeffType = CT_FILTER;
        } else if (coeffType == gensym("area")) {
            x->a_coeffType = CT_AREA;
        } else {
            //error("mbc.allpole~: coefficient must be specified as 'parcor', 'filter', or 'area'");
            x->a_coeffType = CT_FILTER;
        }
#else
        // FORCE TO PARCOR FOR NOW (FOR SIMPLICITY)
        x->a_coeffType = CT_PARCOR;
#endif
        
        //assign locals to globals and init
        x->a_order = (int)order;
        x->a_interp = interp;
        x->a_deemph = (int)deemph;
        allpole_init(x);

	}
	
	else {
		object_free(x);
		x = NULL;
	}
		
	
	return x;
}


void allpole_assist(t_allpole *x, void *b, long m, long a, char *s) {
    if (m==ASSIST_INLET) {
        switch (a) {
            case 0: sprintf(s,"(signal) Filter Input"); break;
            case 1: sprintf(s,"(signal) PARCOR Coefficients"); break;
            case 2: sprintf(s,"(signal) Coeff Index"); break;
            case 3: sprintf(s,"(signal) Filter Gain"); break;
        }
    }
    else {
        sprintf(s,"(signal) Filter Output");
    }
}


inline void allpole_deemphFilter(t_allpole *x, int N, t_double* vec)
{
    int n;
    
    float a1 = x->a_a1;
    float y1 = x->a_y1;
    
    for (n=N; n!=0; n--)
    {
        *vec = *vec + (a1 * y1);
        y1 = *vec++;
    }
    
    x->a_y1 = y1;
}


inline void allpole_solveForFiltCoefs(t_allpole* x, int order)
{
    int i, i1, ji, j;
    
    for (i = 1; i <= order; i++)
    {
        x->a_A[i][i] = x->a_K[i];
        i1 = i - 1;
        
        if (i1 >= 1)
        {
            for (j = 1; j <=i1; j++)
            {
                ji = i - j;
                x->a_A[j][i] = x->a_A[j][i1] - x->a_K[i] * x->a_A[ji][i1];
            }
        }
    }
    
    for (j = 1; j <=order; j++)
    {
        x->a_a[j-1] = x->a_A[j][order];
    }
}

#if defined(MBC_VDSP)
//inline
void allpole_highOrdFilter(t_allpole* x, int N, int order, t_double* in, t_double* out)
{
    int n, i;
    t_mbcfloat val;
    
    for (n=N; n!=0; n--)
    {
        vDSP_vmulD(x->a_a,1,x->a_y,1,x->a_tempVec,1,order);
        vDSP_sveD(x->a_tempVec,1,&val,order);
        val =  x->a_G * (t_mbcfloat)(*(in++)) + val;
        for (i=order-1; i>0; i--) x->a_y[i] = x->a_y[i-1];
        x->a_y[0] = val;
        
        *(out++) = val;
    }
}
#else
//inline
void allpole_highOrdFilter(t_allpole* x, int N, int order, t_double* in, t_double* out)
{
    int n, i;
    t_mbcfloat val;
    
    for (n=N; n!=0; n--)
    {
        val = 0;
        t_mbcfoat* pA = x->a_a;
        t_mbcfloat* pY = x->a_y;
        
        for (i=0; i < order; i++)
        {
            val += *pY *
        }
        vDSP_vmulD(x->a_a,1,x->a_y,1,x->a_tempVec,1,order);
        vDSP_sveD(x->a_tempVec,1,&val,order);
        val =  x->a_G * (t_mbcfloat)(*(in++)) + val;
        for (i=order-1; i>0; i--) x->a_y[i] = x->a_y[i-1];
        x->a_y[0] = val;
        
        *(out++) = val;
    }
}
#endif // MBC_VDSP

