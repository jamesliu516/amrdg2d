/* riemann.cpp
 *
 * Exact Riemann solver for the Euler equations in one dimension
 * Translated from the Fortran code er1pex.f and er1pex.ini
 * by Dr. E.F. Toro downloaded from
 * http://www.numeritek.com/numerica_software.html#freesample
 */

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
using namespace std;
#include"non_uniform_grid.h"

// global variables

const double gama=GAMMA;          // ratio of specific heats

static double  g1 = (gama - 1.0)/(2.0*gama);
static double  g2 = (gama + 1.0)/(2.0*gama);
static double  g3 = 2.0*gama/(gama - 1.0);
static double  g4 = 2.0/(gama - 1.0);
static double  g5 = 2.0/(gama + 1.0);
static double  g6 = (gama - 1.0)/(gama + 1.0);
static double  g7 = (gama - 1.0)/2.0;
static double  g8 = gama - 1.0;

/*
double                // density, velocity, pressure, speed of sound
    dl, ul, pl, cl,   // in left region
    dr, ur, pr, cr;   // in right region
void initialize(
    const int test,   // test number of input data set in e1rpex.ini
    double &domlen,   // domain length
    double &diaph1,   // position of diaphragm
    int &cells,       // number of cells in evaluating exact solution
    double &timeou,   // output time
    double &pscale)   // normalizing factor for pressure and energy
{
    domlen = 1.0;     // same for all tests
    cells = 1000;     // same for all tests
    pscale = 1.0;     // same for all tests

    double values[][9] = {
        // diaph1, gama, timeou, dl, ul, pl, dr, ur, pr
        { 0, 0, 0, 0, 0, 0, 0, 0, 0},  // no TEST 0
        // TEST 1 (Modified Sod)
        {0.3, 1.4, 0.20, 1.0, 0.75, 1.0, 0.125, 0.0, 0.1},
        // TEST 2 (123 problem)
        {0.5, 1.4, 0.15, 1.0, -2.0, 0.4, 1.0, 2.0, 0.4},
        // TEST 3 (Left Woodward & Colella)
        {0.5, 1.4, 0.012, 1.0, 0.0, 1000.0, 1.0, 0.0, 0.01},
        // TEST 4 (Collision of 2 shocks)
        {0.4, 1.4, 0.035, 5.99924, 19.5975, 460.894, 5.99242,
            -6.19633, 46.0950},
        // TEST 5 (Stationary contact)
        {0.8, 1.4, 0.012, 1.0, -19.59745, 1000.0, 1.0, -19.59745, 0.01}
    };

    diaph1 = values[test][0];
    gama = values[test][1];
    timeou = values[test][2];
    dl = values[test][3];
    ul = values[test][4];
    pl = values[test][5];
    dr = values[test][6];
    ur = values[test][7];
    pr = values[test][8];
}
*/

void guessp(double dl, double ul, double pl, double cl,
    double dr, double ur, double pr,double cr, double &pm)
{
    // purpose: to provide a guessed value for pressure
    //          pm in the Star Region. The choice is made
    //          according to adaptive Riemann solver using
    //          the PVRS, TRRS and TSRS approximate
    //          Riemann solvers. See Sect. 9.5 of Chapt. 9 of Ref. 1

    double cup, gel, ger, pmax, pmin, ppv, pq, ptl, ptr,
           qmax, quser, um;

    quser = 2.0;

    // compute guess pressure from PVRS Riemann solver
    cup = 0.25*(dl + dr)*(cl + cr);
    ppv = 0.5*(pl + pr) + 0.5*(ul - ur)*cup;
    ppv = max(0.0, ppv);
    pmin = min(pl,  pr);
    pmax = max(pl,  pr);
    qmax = pmax/pmin;

    if (qmax <= quser && (pmin <= ppv && ppv <= pmax))
        pm = ppv;     // select PVRS Riemann solver
    else {
        if (ppv < pmin) {
            // select Two-Rarefaction Riemann solver
            pq = pow(pl/pr, g1);
            um = (pq*ul/cl + ur/cr + g4*(pq - 1.0))/(pq/cl + 1.0/cr);
            ptl = 1.0 + g7*(ul - um)/cl;
            ptr = 1.0 + g7*(um - ur)/cr;
            pm = 0.5*(pow(pl*ptl, g3) + pow(pr*ptr, g3));
        } else {
            // select Two-Shock Riemann solver with PVRS as estimate
            gel = sqrt((g5/dl)/(g6*pl + ppv));
            ger = sqrt((g5/dr)/(g6*pr + ppv));
            pm = (gel*pl + ger*pr - (ur - ul))/(gel + ger);
        }
    }
}

void prefun(
    double &f,
    double &fd,
    double &p,
    double &dk,
    double &pk,
    double &ck)
{
    // purpose: to evaluate the pressure functions
    //          fl and fr in exact Riemann solver
    //          and their first derivatives

    double ak, bk, pratio, qrt;

    if (p <= pk) {
        // rarefaction wave
        pratio = p/pk;
        f = g4*ck*(pow(pratio, g1) - 1.0);
        fd = (1.0/(dk*ck))*pow(pratio, -g2);
    } else {
        //  shock wave
        ak = g5/dk;
        bk = g6*pk;
        qrt = sqrt(ak/(bk + p));
        f = (p - pk)*qrt;
        fd = (1.0 - 0.5*(p - pk)/(bk + p))*qrt;
    }
}

void starpu(double dl, double ul, double pl, double cl, double dr, double ur, double pr,double cr,
    double &p,
    double &u,
    const double pscale)
{
    // purpose: to compute the solution for pressure and
    //          velocity in the Star Region

    const int nriter = 20;
    const double tolpre = 1.0e-6;
    double change, fl, fld, fr, frd, pold, pstart, udiff;

    // guessed value pstart is computed
    guessp(dl,ul,pl,cl,dr,ur,pr,cr,pstart);
    pold = pstart;
    udiff = ur - ul;

//    cout << "----------------------------------------\n"
  //       << "   Iteration number     Change\n"
  //       << "----------------------------------------" << endl;

    int i = 1;
    for ( ; i <= nriter; i++) {
        prefun(fl, fld, pold, dl, pl, cl);
        prefun(fr, frd, pold, dr, pr, cr);
        p = pold - (fl + fr + udiff)/(fld + frd);
        change = 2.0*fabs((p - pold)/(p + pold));
      //  cout << '\t' << i <<  "\t\t" << change << endl;
        if (change <= tolpre)
            break;
        if (p < 0.0)
            p = tolpre;
        pold = p;
    }
    if (i > nriter)
        cout << "divergence in Newton-Raphson iteration" << endl;

    // compute velocity in star region
    u = 0.5*(ul + ur + fr - fl);
//    cout << "----------------------------------------\n"
 //        << "     Pressure           Velocity\n"
 //        << "----------------------------------------\n"
 //        << "     " << p/pscale << "\t\t" <<  u << '\n'
 //        << "----------------------------------------" << endl;
}

void sample(double dl, double ul, double pl, double cl, double dr, double ur, double pr,double cr,
    const double pm,
    const double um,
    const double s,
    double &d,
    double &u,
    double &p)
{
    // purpose: to sample the solution throughout the wave
    //          pattern. Pressure pm and velocity um in the
    //          star region are known. Sampling is performed
    //          in terms of the 'speed' s = x/t. Sampled
    //          values are d, u, p

    double c, cml, cmr, pml, pmr, shl, shr, sl, sr, stl, str;

    if (s <= um) {
        // sampling point lies to the left of the contact discontinuity
        if (pm <= pl) {
            // left rarefaction
            shl = ul - cl;
            if (s <= shl) {
                // sampled point is left data state
                d = dl;
                u = ul;
                p = pl;
            } else {
                cml = cl*pow(pm/pl, g1);
                stl = um - cml;
                if (s > stl) {
                    // sampled point is star left state
                    d = dl*pow(pm/pl, 1.0/gama);
                    u = um;
                    p = pm;
                } else {
                    // sampled point is inside left fan
                    u = g5*(cl + g7*ul + s);
                    c = g5*(cl + g7*(ul - s));
                    d = dl*pow(c/cl, g4);
                    p = pl*pow(c/cl, g3);
                }
            }
        } else {
            // left shock
            pml = pm/pl;
            sl = ul - cl*sqrt(g2*pml + g1);
            if (s <= sl) {
                // sampled point is left data state
                d = dl;
                u = ul;
                p = pl;
            } else {
                // sampled point is star left state
                d = dl*(pml + g6)/(pml*g6 + 1.0);
                u = um;
                p = pm;
            }
        }
    } else {
        // sampling point lies to the right of the contact discontinuity
        if (pm > pr) {
            // right shock
            pmr = pm/pr;
            sr  = ur + cr*sqrt(g2*pmr + g1);
            if (s >= sr) {
                // sampled point is right data state
                d = dr;
                u = ur;
                p = pr;
            } else {
                // sampled point is star right state
                d = dr*(pmr + g6)/(pmr*g6 + 1.0);
                u = um;
                p = pm;
            }
        } else {
            // right rarefaction
            shr = ur + cr;
            if (s >= shr) {
                // sampled point is right data state
                d = dr;
                u = ur;
                p = pr;
            } else {
                cmr = cr*pow(pm/pr, g1);
                str = um + cmr;
                if (s <= str) {
                    // sampled point is star right state
                    d = dr*pow(pm/pr, 1.0/gama);
                    u = um;
                    p = pm;
                } else {
                    // sampled point is inside left fan
                    u = g5*(-cr + g7*ur + s);
                    c = g5*(cr - g7*(ur - s));
                    d = dr*pow(c/cr, g4);
                    p = pr*pow(c/cr, g3);
                }
            }
        }
    }
}

int main(int argc, char *argv[])
{
    int cells;        // number of cells in evaluating exact solution
    double
        domlen,       // domain length
        diaph1,       // position of diaphragm 1
        timeou,       // output time
        pscale,       // normalizing constant
        ds, dx, pm, ps, s, um, us, xpos;

    int test = 1;
    if (argc > 1) {
        // read first command line argument as test number
        istringstream is(argv[1]);
        is >> test;
        if (test < 1 || test > 5) {
            cerr << "test number not in range 1..5" << endl;
            exit(1);
        }
    }

    initialize(test, domlen, diaph1, cells, timeou, pscale);

    // compute gamma related constants
    g1 = (gama - 1.0)/(2.0*gama);
    g2 = (gama + 1.0)/(2.0*gama);
    g3 = 2.0*gama/(gama - 1.0);
    g4 = 2.0/(gama - 1.0);
    g5 = 2.0/(gama + 1.0);
    g6 = (gama - 1.0)/(gama + 1.0);
    g7 = (gama - 1.0)/2.0;
    g8 = gama - 1.0;

    // compute sound speeds
    cl = sqrt(gama*pl/dl);
    cr = sqrt(gama*pr/dr);

    // the pressure positivity condition is tested for
    if (g4*(cl+cr) <= (ur-ul)) {

        cerr << "the initial data is such that vacuum is generated"
             << "\nstopping program" << endl;
        exit(1);
    }

    // exact solution for pressure and velocity in star region is found
    starpu(pm, um, pscale);
    dx = domlen/double(cells);

    // complete solution at time timeou is found
    ofstream file("riemann.data");
    for (int i = 0; i < cells; i++) {
        xpos = (i - 0.5)*dx;
        s  = (xpos - diaph1)/timeou;

        // solution at point (x,t) = (xpos-diaph1, timeou) is found
        sample(pm, um, s, ds, us, ps);

        // exact solution profiles are written to data file
        file << xpos << '\t' << ds << '\t' << us << '\t'
             << ps/pscale << '\t' << ps/ds/g8/pscale << '\n';
    }
}

void  solveERS(double dl,double ul, double pl,double dr,double ur, 
			  double pr, double *dil,double *dir,double *ui,double *pi)
{
    double cl,cr;
    double pscale=1.0,       // normalizing constant
           ds, dx, pm, ps, s, um, us, xpos;
    cl = sqrt(gama*pl/dl);
    cr = sqrt(gama*pr/dr);
    // the pressure positivity condition is tested for
    if (g4*(cl+cr) <= (ur-ul)) {
        *dil=dl;
        *dir=dr;
        *ui=(ul+ur);
        *pi=0.5*(pl+pr);
        return;
    }

    starpu(dl,ul,pl,cl,dr,ur,pr,cr,pm, um, pscale);

    s=0.0;
    sample(dl,ul,pl,cl,dr,ur,pr,cr,pm, um, s, ds, us, ps);
    *dir=ds;
    *ui=us;
    *pi=ps;




