/* SAS modified this file. */
/* (C) Copyright 2003 Jens Lysgaard. All rights reserved. */
/* OSI Certified Open Source Software */
/* This software is licensed under the Common Public License Version 1.0 */

#ifndef _H_INTAP
#define _H_INTAP

typedef struct
{
  int **c;       /* Cost matrix */
  int *u;        /* Dual variables for rows */
  int *v;        /* Dual variables for columns */
  int *rall;     /* Row allocations */
  int *call;     /* Column allocations */
  char *LR;      /* Labeled rows */
  char *UC;      /* Unlabeled columns */
  int *pi;       /* Minimum reduced costs in unlabeled columns */
  int *cj;       /* pi[i] is in row cj[i] */
  int n;         /* The actual problem size  (1,...,n) */
  int Dim;       /* Allocated dimension (memory is reserved up to n=Dim) */
} INTAPRec;

typedef INTAPRec *INTAPPtr;

void INTAPInitMem(INTAPPtr *P, int n, int Dim);
void INTAPFreeMem(INTAPPtr *P);
void INTAPHungarian(INTAPPtr P);
void INTAPIterate(INTAPPtr P, int Row);
int INTAPObjValue(INTAPPtr P);

#endif
