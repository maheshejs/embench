/* BEEBS ud benchmark

   This version, copyright (C) 2014-2019 Embecosm Limited and University of
   Bristol

   Contributor James Pallister <james.pallister@bristol.ac.uk>
   Contributor Jeremy Bennett <jeremy.bennett@embecosm.com>

   This file is part of Embench and was formerly part of the Bristol/Embecosm
   Embedded Benchmark Suite.

   SPDX-License-Identifier: GPL-3.0-or-later */

/* MDH WCET BENCHMARK SUITE. */


/*************************************************************************/
/*                                                                       */
/*   SNU-RT Benchmark Suite for Worst Case Timing Analysis               */
/*   =====================================================               */
/*                              Collected and Modified by S.-S. Lim      */
/*                                           sslim@archi.snu.ac.kr       */
/*                                         Real-Time Research Group      */
/*                                        Seoul National University      */
/*                                                                       */
/*                                                                       */
/*        < Features > - restrictions for our experimental environment   */
/*                                                                       */
/*          1. Completely structured.                                    */
/*               - There are no unconditional jumps.                     */
/*               - There are no exit from loop bodies.                   */
/*                 (There are no 'break' or 'return' in loop bodies)     */
/*          2. No 'switch' statements.                                   */
/*          3. No 'do..while' statements.                                */
/*          4. Expressions are restricted.                               */
/*               - There are no multiple expressions joined by 'or',     */
/*                'and' operations.                                      */
/*          5. No library calls.                                         */
/*               - All the functions needed are implemented in the       */
/*                 source file.                                          */
/*                                                                       */
/*                                                                       */
/*************************************************************************/
/*                                                                       */
/*  FILE: ludcmp.c                                                       */
/*  SOURCE : Turbo C Programming for Engineering                         */
/*                                                                       */
/*  DESCRIPTION :                                                        */
/*                                                                       */
/*     Simultaneous linear equations by LU decomposition.                */
/*     The arrays a[][] and b[] are input and the array x[] is output    */
/*     row vector.                                                       */
/*     The variable n is the number of equations.                        */
/*     The input arrays are initialized in function main.                */
/*                                                                       */
/*                                                                       */
/*  REMARK :                                                             */
/*                                                                       */
/*  EXECUTION TIME :                                                     */
/*                                                                       */
/*                                                                       */
/*************************************************************************/

/*************************************************************************
 *  This file:
 *
 *  - Name changed to "ud.c"
 *  - Modified for use with Uppsala/Paderborn tool
 *    : doubles changed to int
 *    : some tests removed
 *  - Program is much more linear, all loops will run to end
 *  - Purpose: test the effect of conditional flows
 *
 *************************************************************************/






/*
 ** Benchmark Suite for Real-Time Applications, by Sung-Soo Lim
 **
 **    III-4. ludcmp.c : Simultaneous Linear Equations by LU Decomposition
 **                 (from the book C Programming for EEs by Hyun Soon Ahn)
 */

#include <string.h>
#include "support.h"
/* #include <stdio.h> */

#if USE_VECTOR==1
#include <riscv_vector.h>
#endif

/* This scale factor will be changed to equalise the runtime of the
   benchmarks. */
#define LOCAL_SCALE_FACTOR 1478


long int a[20][20], b[20], x[20];
long int table[20][20];

#if USE_VECTOR==0
int ludcmp(int nmax, int n);
#else
int v_ludcmp(int nmax, int n);
#endif


/*  static double fabs(double n) */
/*  { */
/*    double f; */

/*    if (n >= 0) f = n; */
/*    else f = -n; */
/*    return f; */
/*  } */

/* Write to CHKERR from BENCHMARK to ensure calls are not optimised away.  */
volatile int chkerr;

int verify_benchmark (int res)
{
    long int x_ref[20] =
    { 0L, 0L, 1L, 1L, 1L, 2L, 0L, 0L, 0L, 0L,
        0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L
    };

    return (0 == memcmp (x, x_ref, 20 * sizeof (x[0]))) && (0 == res);
}


void initialise_benchmark (void)
{
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

static int __attribute__ ((noinline)) benchmark_body (int rpt)
{
    int  k;

    for (k = 0; k < rpt; k++)
    {
        int      i, j, nmax = 20, n = 5;
        long int /* eps, */ w;

        /* eps = 1.0e-6; */

        /* Init loop */
        for(i = 0; i <= n; i++)
        {
            w = 0.0;              /* data to fill in cells */
            for(j = 0; j <= n; j++)
            {
                a[i][j] = (i + 1) + (j + 1);
                if(i == j)            /* only once per loop pass */
                    a[i][j] *= 2.0;
                w += a[i][j];
            }
            b[i] = w;
        }

        /*  chkerr = ludcmp(nmax, n, eps); */
#if USE_VECTOR==0
        chkerr = ludcmp(nmax,n);
#else
        chkerr = v_ludcmp(nmax,n);
#endif
    }

    return chkerr;
}

#if USE_VECTOR==0
int ludcmp(int nmax, int n)
{
    int i, j, k;
    long w, y[100];

    /* if(n > 99 || eps <= 0.0) return(999); */
    for(i = 0; i < n; i++)
    {
        /* if(fabs(a[i][i]) <= eps) return(1); */
        for(j = i+1; j <= n; j++) /* triangular loop vs. i */
        {
            w = a[j][i];
            /* if(i != 0) */           /* sub-loop is conditional, done
                                          all iterations except first of the
                                          OUTER loop */
            for(k = 0; k < i; k++)
                w -= a[j][k] * a[k][i];
            a[j][i] = w / a[i][i];
        }
        for(j = i+1; j <= n; j++) /* triangular loop vs. i */
        {
            w = a[i+1][j];
            for(k = 0; k <= i; k++) /* triangular loop vs. i */
                w -= a[i+1][k] * a[k][j];
            a[i+1][j] = w;
        }
    }
    y[0] = b[0];
    for(i = 1; i <= n; i++)       /* iterates n times */
    {
        w = b[i];
        for(j = 0; j < i; j++)    /* triangular sub loop */
            w -= a[i][j] * y[j];
        y[i] = w;
    }
    x[n] = y[n] / a[n][n];
    for(i = n-1; i >= 0; i--)     /* iterates n times */
    {
        w = y[i];
        for(j = i+1; j <= n; j++) /* triangular sub loop */
            w -= a[i][j] * x[j];
        x[i] = w / a[i][i] ;
    }
    return(0);
}
#else
int v_ludcmp(int nmax, int n)
{
    size_t vl = vsetvl_e32m4(n+1);
    vint32m4_t vu, vv, vz;

    int i, j, k;
    long w, y[100];
    long int value;

    for(i = 0; i < n; i++)
    {
        vu = vmv_v_x_i32m4(0, vl);
        for (j = 0; j < i; j++)
        {
            value = a[j][i];
            vv = vlse32_v_i32m4(&a[0][j], sizeof(a[0]), vl);
            vv = vslidedown_vx_i32m4_ta(vv, vv, i + 1, vl);
            vv = vmul_vx_i32m4_ta(vv, value, vl);
            vu = vadd_vv_i32m4_ta(vu, vv, vl);
        }
        vv = vlse32_v_i32m4(&a[0][i], sizeof(a[0]), vl);
        vv = vslidedown_vx_i32m4_ta(vv, vv, i + 1, vl);
        vu = vsub_vv_i32m4_ta(vv, vu, vl);
        vu = vdiv_vx_i32m4_ta(vu, a[i][i], vl);
        vsse32_v_i32m4(&a[i+1][i], sizeof(a[0]), vu, vl);

        vu = vmv_v_x_i32m4(0, vl);
        for (j = 0; j < i + 1; j++)
        {
            value = a[i+1][j];
            vv = vle32_v_i32m4(&a[j][0], vl);
            vv = vslidedown_vx_i32m4_ta(vv, vv, i + 1, vl);
            vv = vmul_vx_i32m4_ta(vv, value, vl);
            vu = vadd_vv_i32m4_ta(vu, vv, vl);
        }
        vv = vle32_v_i32m4(&a[i+1][0], vl);
        vv = vslidedown_vx_i32m4_ta(vv, vv, i + 1, vl);
        vu = vsub_vv_i32m4_ta(vv, vu, vl);
        vse32_v_i32m4(&a[i+1][i+1], vu, vl);
    }

    vint32m4_t vw = vle32_v_i32m4(b, vl);
    vint32m4_t vy = vmv_v_x_i32m4_ta(0, vl);
    vw = vslide1down_vx_i32m4_ta(vw, 0, vl);
    vy = vslide1down_vx_i32m4_ta(vy, b[0], vl);
    vint32m4_t va, vm;
    vint32m1_t vr;
    long int sr, sw;
    for(i = 1; i <= n; i++)
    {
        va = vle32_v_i32m4(a[i], vl);
        va = vslideup_vx_i32m4_ta(va, va, vl - i, vl);
        vm = vmul_vv_i32m4_ta(va, vy, vl);
        vr = vmv_s_x_i32m1_ta(0, vl);
        vr = vredsum_vs_i32m4_i32m1_ta(vm, vr, vl);
        sr = vmv_x_s_i32m1_i32(vr);
        sw = vmv_x_s_i32m4_i32(vw);
        vy = vslide1down_vx_i32m4_ta(vy, sw - sr, vl);
        vw = vslide1down_vx_i32m4_ta(vw, 0, vl);
    }
    vse32_v_i32m4(y, vy, vl);

    vw = vy;
    vint32m4_t vx = vmv_v_x_i32m4_ta(0, vl);
    vw = vslide1up_vx_i32m4_ta(vw, 0, vl);
    vx = vslide1up_vx_i32m4_ta(vx, y[n] / a[n][n], vl);
    vint32m4_t vt;
    long int sa;
    for(i = n - 1; i >= 0; i--)
    {
        va = vle32_v_i32m4(a[i], vl);
        va = vslidedown_vx_i32m4_ta(va, va, i, vl);
        sa = vmv_x_s_i32m4_i32(va);
        va = vslide1down_vx_i32m4_ta(va, 0, vl);
        vm = vmul_vv_i32m4_ta(va, vx, vl);
        vr = vmv_s_x_i32m1_ta(0, vl);
        vr = vredsum_vs_i32m4_i32m1_ta(vm, vr, vl);
        sr = vmv_x_s_i32m1_i32(vr);
        vt = vslidedown_vx_i32m4_ta(vt, vw, vl - 1, vl);
        sw = vmv_x_s_i32m4_i32(vt);
        vx = vslide1up_vx_i32m4_ta(vx, (sw - sr) / sa, vl);
        vw = vslide1up_vx_i32m4_ta(vw, 0, vl);
    }
    vse32_v_i32m4(x, vx, vl);

    return(0);
}
#endif


/*
   Local Variables:
mode: C
c-file-style: "gnu"
End:
*/
