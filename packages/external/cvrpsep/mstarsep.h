/* SAS modified this file. */
/* (C) Copyright 2003 Jens Lysgaard. All rights reserved. */
/* OSI Certified Open Source Software */
/* This software is licensed under the Common Public License Version 1.0 */

#ifndef _H_MSTARSEP
#define _H_MSTARSEP

void MSTARSEP_SeparateMultiStarCuts(int NoOfCustomers,
                                    const double *Demand,
                                    double CAP,
                                    int NoOfEdges,
                                    const int *EdgeTail,
                                    const int *EdgeHead,
                                    const double *EdgeX,
                                    CnstrMgrPointer CMPExistingCuts,
                                    int MaxNoOfCuts,
                                    double *MaxViolation,
                                    CnstrMgrPointer CutsCMP);

#endif

