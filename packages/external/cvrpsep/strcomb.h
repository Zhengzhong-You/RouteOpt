/* SAS modified this file. */
/* (C) Copyright 2003 Jens Lysgaard. All rights reserved. */
/* OSI Certified Open Source Software */
/* This software is licensed under the Common Public License Version 1.0 */

#ifndef _H_STRCOMB
#define _H_STRCOMB

void STRCOMB_MainStrengthenedCombs(ReachPtr SupportPtr, /* Original graph */
                                   int NoOfCustomers,
                                   int CAP,
                                   int *Demand, /* Original demand */
                                   int QMin,
                                   double **XMatrix,
                                   int MaxNoOfCuts,
                                   double *MaxViolation,
                                   CnstrMgrPointer CutsCMP);

#endif

