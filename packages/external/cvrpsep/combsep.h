/* SAS modified this file. */
/* (C) Copyright 2003 Jens Lysgaard. All rights reserved. */
/* OSI Certified Open Source Software */
/* This software is licensed under the Common Public License Version 1.0 */

#ifndef _H_COMBSEP
#define _H_COMBSEP

void COMBSEP_SeparateCombs(int NoOfCustomers,
                           int *Demand,
                           int CAP,
                           int QMin,
                           int NoOfEdges,
                           int *EdgeTail,
                           int *EdgeHead,
                           double *EdgeX,
                           int MaxNoOfCuts,
                           double *MaxViolation,
                           CnstrMgrPointer CutsCMP);

#endif

