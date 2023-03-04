/* BEEBS edn benchmark

   This version, copyright (C) 2014-2019 Embecosm Limited and University of
   Bristol

   Contributor James Pallister <james.pallister@bristol.ac.uk>
   Contributor Jeremy Bennett <jeremy.bennett@embecosm.com>

   This file is part of Embench and was formerly part of the Bristol/Embecosm
   Embedded Benchmark Suite.

   SPDX-License-Identifier: GPL-3.0-or-later

   Original code from: WCET Benchmarks,
http://www.mrtc.mdh.se/projects/wcet/benchmarks.html

Permission to license under GPL obtained by email from Björn Lisper */

/*
 * MDH WCET BENCHMARK SUITE.
 */

/************************************************************************
 *	Simple vector multiply				*
 ************************************************************************/

/*
 * Changes: JG 2005/12/22: Inserted prototypes, changed type of main to int
 * etc. Added parenthesis in expressions in jpegdct. Removed unused variable
 * dx. Changed int to long to avoid problems when compiling to 16 bit target
 * Indented program.
 * JG 2006-01-27: Removed code in codebook
 */

#include <string.h>
#include "support.h"

#if USE_VECTOR==1
#include <riscv_vector.h>
#endif


/* This scale factor will be changed to equalise the runtime of the
   benchmarks. */
#define LOCAL_SCALE_FACTOR 87


#define N 100
#define ORDER 50

#if USE_VECTOR==0
void vec_mpy1 (short y[], const short x[], short scaler);
long int mac (const short *a, const short *b, long int sqr, long int *sum);
void fir (const short array1[], const short coeff[], long int output[]);
void fir_no_red_ld (const short x[], const short h[], long int y[]);
long int latsynth (short b[], const short k[], long int n, long int f);
void iir1 (const short *coefs, const short *input, long int *optr,
        long int *state);
void jpegdct (short *d, short *r);

void vec_mpy1 (short y[], const short x[], short scaler)
{
    long int i;

    for (i = 0; i < 150; i++)
        y[i] += ((scaler * x[i]) >> 15);
}


/*****************************************************
 *			Dot Product	      *
 *****************************************************/
long int mac (const short *a, const short *b, long int sqr, long int *sum)
{
    long int i;
    long int dotp = *sum;

    for (i = 0; i < 150; i++)
    {
        dotp += b[i] * a[i];
        sqr += b[i] * b[i];
    }

    *sum = dotp;
    return sqr;
}


/*****************************************************
 *		FIR Filter		     *
 *****************************************************/
void fir (const short array1[], const short coeff[], long int output[])
{
    long int i, j, sum;

    for (i = 0; i < N - ORDER; i++)
    {
        sum = 0;
        for (j = 0; j < ORDER; j++)
        {
            sum += array1[i + j] * coeff[j];
        }
        output[i] = sum >> 15;
    }
}

/****************************************************
 *	FIR Filter with Redundant Load Elimination

 By doing two outer loops simultaneously, you can potentially  reuse data (depending on the DSP architecture).
 x and h  only  need to be loaded once, therefore reducing redundant loads.
 This reduces memory bandwidth and power.
 *****************************************************/
void fir_no_red_ld (const short x[], const short h[], long int y[])
{
    long int i, j;
    long int sum0, sum1;
    short x0, x1, h0, h1;
    for (j = 0; j < 100; j += 2)
    {
        sum0 = 0;
        sum1 = 0;
        x0 = x[j];
        for (i = 0; i < 32; i += 2)
        {
            x1 = x[j + i + 1];
            h0 = h[i];
            sum0 += x0 * h0;
            sum1 += x1 * h0;
            x0 = x[j + i + 2];
            h1 = h[i + 1];
            sum0 += x1 * h1;
            sum1 += x0 * h1;
        }
        y[j] = sum0 >> 15;
        y[j + 1] = sum1 >> 15;
    }
}

/*******************************************************
 *	Lattice Synthesis	           *
 * This function doesn't follow the typical DSP multiply two  vector operation, but it will point out the compiler's flexibility   ********************************************************/
long int latsynth (short b[], const short k[], long int n, long int f)
{
    long int i;

    f -= b[n - 1] * k[n - 1];
    for (i = n - 2; i >= 0; i--)
    {
        f -= b[i] * k[i];
        b[i + 1] = b[i] + ((k[i] * (f >> 16)) >> 16);
    }
    b[0] = f >> 16;
    return f;
}

/*****************************************************
 *			IIR Filter		     *
 *****************************************************/
void iir1 (const short *coefs, const short *input, long int *optr, long int *state)
{
    long int x;
    long int t;
    long int n;

    x = input[0];
    for (n = 0; n < 50; n++)
    {
        t = x + ((coefs[2] * state[0] + coefs[3] * state[1]) >> 15);
        x = t + ((coefs[0] * state[0] + coefs[1] * state[1]) >> 15);
        state[1] = state[0];
        state[0] = t;
        coefs += 4;		/* point to next filter coefs  */
        state += 2;		/* point to next filter states */
    }
    *optr++ = x;
}

/*****************************************************
 *		JPEG Discrete Cosine Transform 		     *
 *****************************************************/
void jpegdct (short *d, short *r)
{
    long int t[12];
    short i, j, k, m, n, p;
    for (k = 1, m = 0, n = 13, p = 8; k <= 8;
            k += 7, m += 3, n += 3, p -= 7, d -= 64)
    {
        for (i = 0; i < 8; i++, d += p)
        {
            for (j = 0; j < 4; j++)
            {
                t[j] = d[k * j] + d[k * (7 - j)];
                t[7 - j] = d[k * j] - d[k * (7 - j)];
            }
            t[8] = t[0] + t[3];
            t[9] = t[0] - t[3];
            t[10] = t[1] + t[2];
            t[11] = t[1] - t[2];
            d[0] = (t[8] + t[10]) >> m;
            d[4 * k] = (t[8] - t[10]) >> m;
            t[8] = (short) (t[11] + t[9]) * r[10];
            d[2 * k] = t[8] + (short) ((t[9] * r[9]) >> n);
            d[6 * k] = t[8] + (short) ((t[11] * r[11]) >> n);
            t[0] = (short) (t[4] + t[7]) * r[2];
            t[1] = (short) (t[5] + t[6]) * r[0];
            t[2] = t[4] + t[6];
            t[3] = t[5] + t[7];
            t[8] = (short) (t[2] + t[3]) * r[8];
            t[2] = (short) t[2] * r[1] + t[8];
            t[3] = (short) t[3] * r[3] + t[8];
            d[7 * k] = (short) (t[4] * r[4] + t[0] + t[2]) >> n;
            d[5 * k] = (short) (t[5] * r[6] + t[1] + t[3]) >> n;
            d[3 * k] = (short) (t[6] * r[5] + t[1] + t[2]) >> n;
            d[1 * k] = (short) (t[7] * r[7] + t[0] + t[3]) >> n;
        }
    }
}

#else
void v_vec_mpy1 (short y[], const short x[], short scaler);
long int v_mac (const short *a, const short *b, long int sqr, long int *sum);
void v_fir(const short x[], const short h[], long int y[], long int n, long int order);
long int v_latsynth (short b[], const short k[], long int n, long int f);
void v_iir1 (const short *coefs, const short *input, long int *optr,
        long int *state);
void v_jpegdct (short *d, short *r);

void v_vec_mpy1(short y[], const short x[], short scaler)
{
    asm volatile (
            "        csrwi 10, 2;\n"     
            : 
            : 
            : 
            );

    size_t vl;
    vl = vsetvl_e16m2(150);

    vint16m2_t vx;
    vx = vle16_v_i16m2(&x[0], vl);
    vx = vsmul_vx_i16m2(vx, scaler, vl);

    vint16m2_t vy;
    vy = vle16_v_i16m2(&y[0], vl);
    vy = vadd_vv_i16m2_ta(vx, vy, vl);
    vse16_v_i16m2(&y[0], vy, vl);
}

long int v_mac(const short *a, const short *b, long int sqr, long int *sum)
{
    size_t vl;
    vl = vsetvl_e16m2(150);

    vint16m2_t va = vle16_v_i16m2(&a[0], vl);
    vint16m2_t vb = vle16_v_i16m2(&b[0], vl);

    vint32m4_t vm;
    vint32m1_t vr;

    vm = vwmul_vv_i32m4_ta(va, vb, vl);
    vr = vmv_s_x_i32m1_ta(*sum, vl);
    vr = vredsum_vs_i32m4_i32m1_ta(vm, vr, vl);
    *sum = vmv_x_s_i32m1_i32(vr);

    vm = vwmul_vv_i32m4_ta(vb, vb, vl);
    vr = vmv_s_x_i32m1_ta(sqr, vl);
    vr = vredsum_vs_i32m4_i32m1_ta(vm, vr, vl);
    return vmv_x_s_i32m1_i32(vr);
}

void v_fir(const short x[], const short h[], long int y[], long int n, long int order)
{
    size_t vl;
    vl = vsetvl_e16m2(order);
    vint16m2_t vh = vle16_v_i16m2(&h[0], vl);
    vint32m4_t vm;
    vint32m1_t vr;
    long int r;

    vl = vsetvl_e16m2(n);
    vint16m2_t vx = vle16_v_i16m2(&x[0], vl);
    vl = vsetvl_e32m4(n - order);
    vint32m4_t vy;
    for (long int j = 0; j < vl; ++j)
    {
        vl = vsetvl_e16m2(order);
        vm = vwmul_vv_i32m4_ta(vh, vx, vl);

        vr = vmv_s_x_i32m1_ta(0, vl);
        vr = vredsum_vs_i32m4_i32m1_ta(vm, vr, vl);
        r  = vmv_x_s_i32m1_i32(vr);

        vl = vsetvl_e16m2(n);
        vx = vslide1down_vx_i16m2_ta(vx, 0, vl);

        vl = vsetvl_e16m2(n - order);
        vy = vslide1down_vx_i32m4_ta(vy, r, vl);
    }
    vy = vsra_vx_i32m4_ta(vy, 15, vl);
    vse32_v_i32m4(&y[0], vy, vl);
}

long int v_latsynth (short b[], const short k[], long int n, long int f)
{
    asm volatile (
            "        csrwi 10, 2;\n"     
            : 
            : 
            : 
            );

    long int r;

    size_t vl;

    vl = vsetvl_e16m2(n);

    vint16m2_t vb = vle16_v_i16m2(&b[0], vl);
    vint16m2_t vk = vle16_v_i16m2(&k[0], vl);
    vint16m2_t vfn;
    vint16m2_t vt;

    vl = vsetvl_e32m4(n);

    long int x;
    vint32m4_t vx;
    vint32m4_t vf = vmv_s_x_i32m4_ta(f, vl);
    vint32m4_t vw = vwmul_vv_i32m4_ta(vb, vk, vl);
    vw = vmul_vx_i32m4_ta(vw, -1, vl);

    for (long int i = 0; i < n; ++i)
    {
        vx  = vslidedown_vx_i32m4_ta(vx, vw, vl - 1, vl);
        x   = vmv_x_s_i32m4_i32(vx) + vmv_x_s_i32m4_i32(vf);
        vf  = vslide1up_vx_i32m4_ta(vf, x, vl);
        vw  = vslide1up_vx_i32m4_ta(vw, 0, vl);
    }

    r = vmv_x_s_i32m4_i32(vf);

    vl = vsetvl_e16m2(n);

    vk  = vsra_vx_i16m2_ta(vk, 1, vl);
    vfn = vnsra_wx_i16m2_ta(vf, 16, vl);
    vt  = vsmul_vv_i16m2_ta(vk, vfn, vl);
    vb  = vadd_vv_i16m2_ta(vb, vt, vl);
    vb  = vslide1up_vx_i16m2_ta(vb, vmv_x_s_i16m2_i16(vfn), vl);

    vse16_v_i16m2(&b[0], vb, vl);

    return r;
}

void v_iir1 (const short *coefs, const short *input, long int *optr, long int *state)
{
    long int x;
    long int t;
    long int n;

    x = input[0];
    for (n = 0; n < 50; n++)
    {
        t = x + ((coefs[2] * state[0] + coefs[3] * state[1]) >> 15);
        x = t + ((coefs[0] * state[0] + coefs[1] * state[1]) >> 15);
        state[1] = state[0];
        state[0] = t;
        coefs += 4;		/* point to next filter coefs  */
        state += 2;		/* point to next filter states */
    }
    *optr++ = x;
}

void v_jpegdct (short *d, short *r)
{
    long int t[12];
    short i, j, k, m, n, p;
    for (k = 1, m = 0, n = 13, p = 8; k <= 8;
            k += 7, m += 3, n += 3, p -= 7, d -= 64)
    {
        for (i = 0; i < 8; i++, d += p)
        {
            for (j = 0; j < 4; j++)
            {
                t[j] = d[k * j] + d[k * (7 - j)];
                t[7 - j] = d[k * j] - d[k * (7 - j)];
            }
            t[8] = t[0] + t[3];
            t[9] = t[0] - t[3];
            t[10] = t[1] + t[2];
            t[11] = t[1] - t[2];
            d[0] = (t[8] + t[10]) >> m;
            d[4 * k] = (t[8] - t[10]) >> m;
            t[8] = (short) (t[11] + t[9]) * r[10];
            d[2 * k] = t[8] + (short) ((t[9] * r[9]) >> n);
            d[6 * k] = t[8] + (short) ((t[11] * r[11]) >> n);
            t[0] = (short) (t[4] + t[7]) * r[2];
            t[1] = (short) (t[5] + t[6]) * r[0];
            t[2] = t[4] + t[6];
            t[3] = t[5] + t[7];
            t[8] = (short) (t[2] + t[3]) * r[8];
            t[2] = (short) t[2] * r[1] + t[8];
            t[3] = (short) t[3] * r[3] + t[8];
            d[7 * k] = (short) (t[4] * r[4] + t[0] + t[2]) >> n;
            d[5 * k] = (short) (t[5] * r[6] + t[1] + t[3]) >> n;
            d[3 * k] = (short) (t[6] * r[5] + t[1] + t[2]) >> n;
            d[1 * k] = (short) (t[7] * r[7] + t[0] + t[3]) >> n;
        }
    }
}
#endif

static short a[200];
static short b[200];
static short c;
static long int d;
static int e;
static long int output[200];


void initialise_benchmark (void)
{
#if USE_VECTOR==0
    printf("SCALAR VERSION\n");
#else
    printf("VECTOR VERSION\n");
#endif
}


static int benchmark_body (int  rpt);

void warm_caches (int  heat)
{
    int  res = benchmark_body (heat);

    return;
}


int benchmark (void)
{
    return benchmark_body (LOCAL_SCALE_FACTOR * CPU_MHZ);
}

#include <stdio.h>
static int __attribute__ ((noinline)) benchmark_body (int rpt)
{
    int j;

    for (j = 0; j < rpt; j++)
    {
        short unsigned int in_a[200] = {
            0x0000, 0x07ff, 0x0c00, 0x0800, 0x0200, 0xf800, 0xf300, 0x0400,
            0x0000, 0x07ff, 0x0c00, 0x0800, 0x0200, 0xf800, 0xf300, 0x0400,
            0x0000, 0x07ff, 0x0c00, 0x0800, 0x0200, 0xf800, 0xf300, 0x0400,
            0x0000, 0x07ff, 0x0c00, 0x0800, 0x0200, 0xf800, 0xf300, 0x0400,
            0x0000, 0x07ff, 0x0c00, 0x0800, 0x0200, 0xf800, 0xf300, 0x0400,
            0x0000, 0x07ff, 0x0c00, 0x0800, 0x0200, 0xf800, 0xf300, 0x0400,
            0x0000, 0x07ff, 0x0c00, 0x0800, 0x0200, 0xf800, 0xf300, 0x0400,
            0x0000, 0x07ff, 0x0c00, 0x0800, 0x0200, 0xf800, 0xf300, 0x0400,
            0x0000, 0x07ff, 0x0c00, 0x0800, 0x0200, 0xf800, 0xf300, 0x0400,
            0x0000, 0x07ff, 0x0c00, 0x0800, 0x0200, 0xf800, 0xf300, 0x0400,
            0x0000, 0x07ff, 0x0c00, 0x0800, 0x0200, 0xf800, 0xf300, 0x0400,
            0x0000, 0x07ff, 0x0c00, 0x0800, 0x0200, 0xf800, 0xf300, 0x0400,
            0x0000, 0x07ff, 0x0c00, 0x0800, 0x0200, 0xf800, 0xf300, 0x0400,
            0x0000, 0x07ff, 0x0c00, 0x0800, 0x0200, 0xf800, 0xf300, 0x0400,
            0x0000, 0x07ff, 0x0c00, 0x0800, 0x0200, 0xf800, 0xf300, 0x0400,
            0x0000, 0x07ff, 0x0c00, 0x0800, 0x0200, 0xf800, 0xf300, 0x0400,
            0x0000, 0x07ff, 0x0c00, 0x0800, 0x0200, 0xf800, 0xf300, 0x0400,
            0x0000, 0x07ff, 0x0c00, 0x0800, 0x0200, 0xf800, 0xf300, 0x0400,
            0x0000, 0x07ff, 0x0c00, 0x0800, 0x0200, 0xf800, 0xf300, 0x0400,
            0x0000, 0x07ff, 0x0c00, 0x0800, 0x0200, 0xf800, 0xf300, 0x0400,
            0x0000, 0x07ff, 0x0c00, 0x0800, 0x0200, 0xf800, 0xf300, 0x0400,
            0x0000, 0x07ff, 0x0c00, 0x0800, 0x0200, 0xf800, 0xf300, 0x0400,
            0x0000, 0x07ff, 0x0c00, 0x0800, 0x0200, 0xf800, 0xf300, 0x0400,
            0x0000, 0x07ff, 0x0c00, 0x0800, 0x0200, 0xf800, 0xf300, 0x0400,
            0x0000, 0x07ff, 0x0c00, 0x0800, 0x0200, 0xf800, 0xf300, 0x0400
        };
        short unsigned int in_b[200] = {
            0x0c60, 0x0c40, 0x0c20, 0x0c00, 0xf600, 0xf400, 0xf200, 0xf000,
            0x0c60, 0x0c40, 0x0c20, 0x0c00, 0xf600, 0xf400, 0xf200, 0xf000,
            0x0c60, 0x0c40, 0x0c20, 0x0c00, 0xf600, 0xf400, 0xf200, 0xf000,
            0x0c60, 0x0c40, 0x0c20, 0x0c00, 0xf600, 0xf400, 0xf200, 0xf000,
            0x0c60, 0x0c40, 0x0c20, 0x0c00, 0xf600, 0xf400, 0xf200, 0xf000,
            0x0c60, 0x0c40, 0x0c20, 0x0c00, 0xf600, 0xf400, 0xf200, 0xf000,
            0x0c60, 0x0c40, 0x0c20, 0x0c00, 0xf600, 0xf400, 0xf200, 0xf000,
            0x0c60, 0x0c40, 0x0c20, 0x0c00, 0xf600, 0xf400, 0xf200, 0xf000,
            0x0c60, 0x0c40, 0x0c20, 0x0c00, 0xf600, 0xf400, 0xf200, 0xf000,
            0x0c60, 0x0c40, 0x0c20, 0x0c00, 0xf600, 0xf400, 0xf200, 0xf000,
            0x0c60, 0x0c40, 0x0c20, 0x0c00, 0xf600, 0xf400, 0xf200, 0xf000,
            0x0c60, 0x0c40, 0x0c20, 0x0c00, 0xf600, 0xf400, 0xf200, 0xf000,
            0x0c60, 0x0c40, 0x0c20, 0x0c00, 0xf600, 0xf400, 0xf200, 0xf000,
            0x0c60, 0x0c40, 0x0c20, 0x0c00, 0xf600, 0xf400, 0xf200, 0xf000,
            0x0c60, 0x0c40, 0x0c20, 0x0c00, 0xf600, 0xf400, 0xf200, 0xf000,
            0x0c60, 0x0c40, 0x0c20, 0x0c00, 0xf600, 0xf400, 0xf200, 0xf000,
            0x0c60, 0x0c40, 0x0c20, 0x0c00, 0xf600, 0xf400, 0xf200, 0xf000,
            0x0c60, 0x0c40, 0x0c20, 0x0c00, 0xf600, 0xf400, 0xf200, 0xf000,
            0x0c60, 0x0c40, 0x0c20, 0x0c00, 0xf600, 0xf400, 0xf200, 0xf000,
            0x0c60, 0x0c40, 0x0c20, 0x0c00, 0xf600, 0xf400, 0xf200, 0xf000,
            0x0c60, 0x0c40, 0x0c20, 0x0c00, 0xf600, 0xf400, 0xf200, 0xf000,
            0x0c60, 0x0c40, 0x0c20, 0x0c00, 0xf600, 0xf400, 0xf200, 0xf000,
            0x0c60, 0x0c40, 0x0c20, 0x0c00, 0xf600, 0xf400, 0xf200, 0xf000,
            0x0c60, 0x0c40, 0x0c20, 0x0c00, 0xf600, 0xf400, 0xf200, 0xf000,
            0x0c60, 0x0c40, 0x0c20, 0x0c00, 0xf600, 0xf400, 0xf200, 0xf000
        };
        c = 0x3;
        d = 0xAAAA;
        e = 0xEEEE;

        for (int i = 0; i < 200; i++)
        {
            a[i] = in_a[i];
            b[i] = in_b[i];
        }
        /*
         * Declared as memory variable so it doesn't get optimized out
         */


#if USE_VECTOR==0
        vec_mpy1 (a, b, c);
        c = mac (a, b, (long int) c, (long int *) output);
        fir (a, b, output);
        fir_no_red_ld (a, b, output);
        d = latsynth (a, b, N, d);
        e = d; // codebook
        iir1 (a, b, &output[100], output);
        jpegdct (a, b);
#else
        v_vec_mpy1 (a, b, c);
        c = v_mac (a, b, (long int) c, (long int *) output);
        v_fir (a, b, output, N, ORDER);
        v_fir (a, b, output, 100 + 32, 32);
        d = v_latsynth (a, b, N, d);
        e = d; // codebook
        v_iir1 (a, b, &output[100], output);
        v_jpegdct (a, b);
#endif
    }
    return 0;
}

int verify_benchmark (int unused)
{
    long int exp_output[200] =
    { 3760, 4269, 3126, 1030, 2453, -4601, 1981, -1056, 2621, 4269,
        3058, 1030, 2378, -4601, 1902, -1056, 2548, 4269, 2988, 1030,
        2300, -4601, 1822, -1056, 2474, 4269, 2917, 1030, 2220, -4601,
        1738, -1056, 2398, 4269, 2844, 1030, 2140, -4601, 1655, -1056,
        2321, 4269, 2770, 1030, 2058, -4601, 1569, -1056, 2242, 4269,
        2152, 1030, 1683, -4601, 1627, -1056, 2030, 4269, 2080, 1030,
        1611, -4601, 1555, -1056, 1958, 4269, 2008, 1030, 1539, -4601,
        1483, -1056, 1886, 4269, 1935, 1030, 1466, -4601, 1410, -1056,
        1813, 4269, 1862, 1030, 1393, -4601, 1337, -1056, 1740, 4269,
        1789, 1030, 1320, -4601, 1264, -1056, 1667, 4269, 1716, 1030,
        1968, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };

    return (0 == memcmp (output, exp_output, 200 * sizeof (output[0])))
        && (10243 == c) && (-441886230 == d) && (-441886230 == e);
}


/*
   Local Variables:
mode: C
c-file-style: "gnu"
End:
*/
