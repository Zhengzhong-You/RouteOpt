/* SAS modified this file. */
/* (C) Copyright 2003 Jens Lysgaard. All rights reserved. */
/* OSI Certified Open Source Software */
/* This software is licensed under the Common Public License Version 1.0 */

#ifndef _H_TWOMATCH
#define _H_TWOMATCH

void TWOMATCH_ExactTwoMatchings(ReachPtr SupportPtr,
                                int NoOfCustomers,
                                char *DepotEdgeBound,
                                double **XMatrix,
                                CnstrMgrPointer CutsCMP);

#endif

