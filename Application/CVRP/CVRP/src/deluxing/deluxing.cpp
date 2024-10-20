#include "deluxing.hpp"
#include <sys/time.h>
#include <iomanip>
#include "deluxing_kmeans.hpp"

#define WIDTH 10

using namespace std;


double getWallTime() {
    struct timeval time;
    return gettimeofday(&time, NULL) ? 0 : (double) time.tv_sec + (double) time.tv_usec * .000001;
}


void deLuxing(GRBmodel *lp, double UB, int NClust, int beta1, int beta2, vector<int> &Idxdel,
              double Timelimit, double Tolerance, bool Verbose) {
    double t = getWallTime();
    GRBsetintparam(GRBgetenv(lp), "Method", 2);
    GRBsetintparam(GRBgetenv(lp), "Crossover", 0);
    //    GRBsetintparam(GRBgetenv(lp), "OutputFlag", 1);
    //    GRBsetintparam(GRBgetenv(lp), "Presolve", 0);
    GRBsetintparam(GRBgetenv(lp), "Threads", 1);
    GRBoptimize(lp);
    double t1 = getWallTime() - t;
    int NRow, NCol, NN;
    GRBgetintattr(lp, "NumConstrs", &NRow);
    GRBgetintattr(lp, "NumVars", &NCol);
    NN = NCol;

    double *dj = new double[NCol];
    size_t numnzP;
    GRBXgetvars(lp, &numnzP, NULL, NULL, NULL, 0, NCol);

    int nv = max(NCol, NRow);
    double *cobj = new double[nv], *obj_real = new double[NCol];
    size_t *cmatbeg = new size_t[nv + 1], *matbeg = new size_t[NCol + 1];
    numnzP = numnzP < nv ? nv : numnzP;
    int *cmatind = new int[numnzP], *matind = new int[numnzP];
    double *cmatval = new double[numnzP], *matval = new double[numnzP];
    double *rhs_real = new double[NRow], *dual = new double[NRow];
    bool *ifdel = new bool[NCol], *delall = new bool[NCol];
    int *del_idx = new int[NCol];

    GRBgetdblattrarray(lp, "RHS", 0, NRow, rhs_real);
    GRBXgetvars(lp, &numnzP, matbeg, matind, matval, 0, NCol);
    matbeg[NCol] = numnzP;


    double LB;
    GRBgetdblattr(lp, "ObjVal", &LB);
    double gap = UB - LB + Tolerance, gap_half = (UB - LB) / 2.0 + Tolerance;
    int num_half = 0, ncol_del0 = 0;
    GRBXgetvars(lp, &numnzP, cmatbeg, cmatind, cmatval, 0, NCol);
    matbeg[NCol] = cmatbeg[NCol] = numnzP;
    GRBgetdblattrarray(lp, "RC", 0, NCol, dj);
    GRBgetdblattrarray(lp, "Obj", 0, NCol, cobj);
    vector<pair<double, int> > RCVec;
    for (int k = 0; k < NCol; ++k)
        RCVec.push_back(make_pair(dj[k], k));
    stable_sort(RCVec.begin(), RCVec.end());


    vector<int> GlobalPool;
    for (int j = 0; j < NCol; ++j)
        GlobalPool.push_back(j);
    vector<int> PoolCopy = GlobalPool;
    for (int k = 0; k < NCol; ++k) {
        GlobalPool[k] = PoolCopy[RCVec[k].second];
        if (dj[k] > gap) ++ncol_del0;
    }

    matbeg[0] = numnzP = 0;
    for (int k = 0; k < NCol - ncol_del0; ++k) {
        if (RCVec[k].first > gap_half && !num_half)
            num_half = k;
        for (auto j = cmatbeg[RCVec[k].second]; j < cmatbeg[RCVec[k].second + 1]; ++j) {
            matind[numnzP] = cmatind[j];
            matval[numnzP++] = cmatval[j];
        }
        matbeg[k + 1] = numnzP;
        obj_real[k] = cobj[RCVec[k].second];
    }
    int *idx = new int[NCol];
    vector<uint32_t> clusters(NCol, 0);
    double t2 = getWallTime();
    if (ncol_del0 >= NCol) goto Terminate;
    cout << "|R| = " << NN << ", initial del = " << ncol_del0 << ", |R_1| = " << num_half << endl;
    cout << "------------------------------------------------------------------------------------" << endl;


    iota(idx, idx + NCol, 0);
    GRBdelvars(lp, NCol, idx);
    GRBupdatemodel(lp);
    GlobalPool.erase(GlobalPool.begin() + NCol - ncol_del0, GlobalPool.end());
    NCol -= ncol_del0;
    GlobalPool.resize(NCol);
    PoolCopy.clear();
    RCVec.clear();
    GRBXaddvars(lp, num_half, matbeg[num_half], matbeg, matind, matval, obj_real, NULL, NULL, NULL, NULL);

    cmatbeg[0] = 0;
    for (int k = 0; k < NRow; ++k) {
        cmatind[k] = k;
        cmatval[k] = -1.;
        cobj[k] = UB;
        cmatbeg[k + 1] = k + 1;
    }
    GRBXaddvars(lp, NRow, NRow, cmatbeg, cmatind, cmatval, cobj, NULL, NULL, NULL, NULL);

    if (NCol <= beta1) {
        vector<vector<double> > data;
        for (int j = 0; j < NCol; ++j) {
            vector<double> col(rhs_real, rhs_real + NRow);
            for (auto k = matbeg[j]; k < matbeg[j + 1]; ++k)
                col[matind[k]] -= matval[k];
            data.emplace_back(col);
        }
        kmeans::clustering_parameters<double> parameters(NClust);
        parameters.set_random_seed(0);
        parameters.set_max_iteration(100);
        vector<vector<double> > means;
    kmeans:
        kmeans_lloyd_parallel(data, means, clusters, parameters);
    } else {
        vector<pair<double, int> > col_norm(NCol);
        for (int j = 0; j < NCol; ++j) {
            for (int k = 0; k < NRow; ++k)
                cmatval[k] = rhs_real[k];
            for (auto k = matbeg[j]; k < matbeg[j + 1]; ++k)
                if (matind[k] < NRow)
                    cmatval[matind[k]] -= matval[k];
            double sum = 0.;
            for (int k = 0; k < NRow; ++k)
                sum += cmatval[k] * cmatval[k];
            col_norm[j] = make_pair(sum, j);
        }
        stable_sort(col_norm.begin(), col_norm.end());
        double rate = 1.0 / NClust;
        for (int r = 0; r < NClust; ++r) {
            int s1 = r * NCol * rate, s2 = (r == NClust - 1 ? NCol : (r + 1) * NCol * rate);
            for (int j = s1; j < s2; ++j)
                clusters[col_norm[j].second] = r;
        }
    }
    t2 = getWallTime() - t2;


    for (int r = 0; r < NClust; ++r) {
        int ncol_del = 0;
        for (int k = 0; k < NRow; ++k) cmatval[k] = 0.;

        for (int j = 0; j < NCol; ++j)
            if (clusters[j] == r) {
                for (int k = 0; k < NRow; ++k)
                    cmatval[k] += rhs_real[k];
                for (auto k = matbeg[j]; k < matbeg[j + 1]; ++k)
                    cmatval[matind[k]] -= matval[k];
            }


        GRBsetdblattrarray(lp, "RHS", 0, NRow, cmatval);
        GRBoptimize(lp);
        GRBgetdblattrarray(lp, "PI", 0, NRow, dual);
        double threshold = UB + Tolerance;
        for (int k = 0; k < NRow; ++k)
            threshold -= dual[k] * rhs_real[k];

        double minrc = 0.;
        for (int j = num_half; j < NCol; ++j) {
            double rcj = obj_real[j];
            for (auto k = matbeg[j]; k < matbeg[j + 1]; ++k)
                rcj -= dual[matind[k]] * matval[k];
            if (rcj >= threshold) {
                ifdel[j] = true;
                ++ncol_del;
            } else {
                ifdel[j] = false;
                if (rcj < minrc) minrc = rcj;
            }
        }
        threshold -= minrc;
        GRBgetdblattrarray(lp, "RC", 0, num_half, dj);
        int dc = 0, real_half = num_half;
        vector<int> copy(real_half);
        for (int j = 0; j < num_half; ++j) {
            copy[j] = j;
            if (dj[j] >= threshold) {
                ifdel[j] = true;
                ++ncol_del;
                del_idx[dc++] = j;
            } else ifdel[j] = false;
        }


        if (dc) {
            for (int j = 0; j < num_half; ++j)
                if (ifdel[j])
                    copy.erase(std::remove(copy.begin(), copy.end(), j), copy.end());
            real_half -= dc;
            GRBdelvars(lp, dc, del_idx);
            GRBupdatemodel(lp);
        }

        for (int j = 0; j < NCol; ++j) {
            delall[j] = ifdel[j];
            ifdel[j] &= clusters[j] != r;
        }
        int extra = ncol_del, ite = 1;
        while (extra > beta2) {
            extra = 0;
            for (int k = 0; k < NRow; ++k) cmatval[k] = 0.;
            for (int j = 0; j < NCol; ++j) {
                if (ifdel[j]) {
                    for (int k = 0; k < NRow; ++k)
                        cmatval[k] += rhs_real[k];
                    for (auto k = matbeg[j]; k < matbeg[j + 1]; ++k)
                        cmatval[matind[k]] -= matval[k];
                }
            }
            GRBsetdblattrarray(lp, "RHS", 0, NRow, cmatval);
            GRBoptimize(lp);
            GRBgetdblattrarray(lp, "Pi", 0, NRow, dual);
            threshold = UB + Tolerance;
            for (int k = 0; k < NRow; ++k)
                threshold -= dual[k] * rhs_real[k];

            minrc = 0.;
            for (int j = num_half; j < NCol; ++j) {
                ifdel[j] = false;
                if (delall[j]) continue;
                double rcj = obj_real[j];
                for (auto k = matbeg[j]; k < matbeg[j + 1]; ++k)
                    rcj -= dual[matind[k]] * matval[k];
                if (rcj >= threshold) {
                    delall[j] = ifdel[j] = true;
                    ++extra;
                } else if (rcj < minrc) minrc = rcj;
            }
            threshold -= minrc;
            for (int j = 0; j < num_half; ++j) {
                ifdel[j] = false;
                if (delall[j]) continue;
                double rcj = obj_real[j];
                for (auto k = matbeg[j]; k < matbeg[j + 1]; ++k)
                    rcj -= dual[matind[k]] * matval[k];
                if (rcj >= threshold) {
                    delall[j] = ifdel[j] = true;
                    ++extra;
                }
            }

            dc = 0;
            for (int j = 0; j < num_half; ++j)
                if (ifdel[j])
                    del_idx[dc++] = std::find(copy.begin(), copy.end(), j) - copy.begin();

            if (dc) {
                for (int j = 0; j < num_half; ++j)
                    if (ifdel[j])
                        copy.erase(std::remove(copy.begin(), copy.end(), j), copy.end());
                real_half -= dc;
                GRBdelvars(lp, dc, del_idx);
                GRBupdatemodel(lp);
            }

            ncol_del += extra;
            ++ite;
            if (getWallTime() - t > Timelimit) break;
        }


        int kc = numnzP = 0;
        for (int i = 0; i < NCol; ++i) {
            if (!delall[i]) {
                for (auto k = matbeg[i]; k < matbeg[i + 1]; ++k) {
                    matind[numnzP] = matind[k];
                    matval[numnzP++] = matval[k];
                }
                obj_real[kc] = obj_real[i];
                GlobalPool[kc] = GlobalPool[i];
                clusters[kc] = clusters[i];
                matbeg[++kc] = numnzP;
            }
        }

        num_half = real_half;
        GlobalPool.resize(kc);
        NCol -= ncol_del;
        ncol_del += ncol_del0;

        if (Verbose)
            cout << "Round " << r + 1 << ", deleted " << ncol_del << " cols, ite = " << ite << ",  #cols = " <<
                    NCol << endl;
        ncol_del0 = 0;
    }

    iota(idx, idx + NCol, 0);
    GRBdelvars(lp, NCol, idx);
    GRBupdatemodel(lp);
    NCol = (int) GlobalPool.size();
    GRBXaddvars(lp, NCol, numnzP, matbeg, matind, matval, obj_real, NULL, NULL, NULL, NULL);
    GRBsetdblattrarray(lp, "RHS", 0, NRow, rhs_real);
    GRBupdatemodel(lp);


Terminate:
    delete [] cobj;
    delete [] obj_real;
    delete [] dual;
    delete [] dj;
    delete [] cmatbeg;
    delete [] cmatind;
    delete [] cmatval;
    delete [] matbeg;
    delete [] matind;
    delete [] matval;
    delete [] rhs_real;
    delete [] ifdel;
    delete [] delall;
    delete [] del_idx;
    delete [] idx;
    GRBsetintparam(GRBgetenv(lp), "Method", -1);
    //      GRBsetintparam(GRBgetenv(lp), "OutputFlag", 0);
    //	GRBsetintparam(GRBgetenv(lp), "Presolve", -1);
    stable_sort(GlobalPool.begin(), GlobalPool.end());

    if (NCol) {
        for (int j = 0; j < GlobalPool[0]; ++j)
            Idxdel.push_back(j);
        for (int k = 0; k < NCol - 1; ++k)
            for (int j = GlobalPool[k] + 1; j < GlobalPool[k + 1]; ++j)
                Idxdel.push_back(j);
        for (int j = GlobalPool[NCol - 1] + 1; j < NN; ++j)
            Idxdel.push_back(j);
    }
    cout << "------------------------------------------------------------------------------------" << endl;
    std::cout << std::setprecision(2) << std::fixed;
    cout << "Total time = " << getWallTime() - t << ", |R| = " << NN << ", NRemain = " << NCol << ", percentage = " << (
        NN - NCol) * 100.0 / NN << "%\n\n" << endl;
    //	cout << "Total time = " << getWallTime() - t << ", first LP time = " << t1 << ", clustering time  = " << t2 << ", |R| = " << NN << ", Ndel = " << NN - NCol << ", del " << (NN - NCol) *100.0 / NN << "%\n\n" << endl;
}
