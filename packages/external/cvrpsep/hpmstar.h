/* SAS modified this file. */
/* (C) Copyright 2003 Jens Lysgaard. All rights reserved. */
/* OSI Certified Open Source Software */
/* This software is licensed under the Common Public License Version 1.0 */

#ifndef _H_HPMSTAR
#define _H_HPMSTAR

#include "cnstrmgr.h"

void HPMSTAR_CreateMinVVector(int DemandSum, int CAP);

void HPMSTAR_DirectX(ReachPtr SupportPtr, /* Original support graph. */
                     ReachPtr SAdjRPtr,   /* Shrunk support graph. */
                     int NoOfCustomers,
                     const double *Demand,         /* Original demand vector. */
                     double CAP,
                     int NoOfSuperNodes,
                     double *XInSuperNode,
                     double **XMatrix,
                     double **SMatrix,
                     ReachPtr SuperNodesRPtr,
                     char SelectionRule,
                     int MaxGeneratedCuts,
                     char SearchPartialMStars,
                     CnstrMgrPointer CutsCMP,
                     double *MaxViolation);

#endif

