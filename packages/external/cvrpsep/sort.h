/* SAS modified this file. */
/* (C) Copyright 2003 Jens Lysgaard. All rights reserved. */
/* OSI Certified Open Source Software */
/* This software is licensed under the Common Public License Version 1.0 */

#ifndef _H_SORT
#define _H_SORT

void SortIVInc(int *Value, int n);
void SortIVDec(int *Value, int n);
void SortDVInc(double *Value, int n);
void SortDVDec(double *Value, int n);
void SortIndexIVInc(int *Index, const int *Value, int n);
void SortIndexIVDec(int *Index, const int *Value, int n);
void SortIndexDVDec(int *Index, const double *Value, int n);
void SortIndexDVInc(int *Index, const double *Value, int n);

#endif


