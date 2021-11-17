//  ThermalEvolution1D.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
//
//  Development Plan:
//  1)  Implement a radial conduction model for a one-layer sphere (done)
//  2)  Save data into a convenient .csv format (done)
//  3)  Plot the data using a MATLAB script (done)
//  4)  Optimize current 1D model for speed and code clean up (done)
//  5)  Add ice shell (done)
//  6)  Melt/freeze the ice shell (done)
//  7)  Implement silicate hydration/dehydration
//      a) Implement latent heat tracking (done)
//      b) Incorporate specific heats (done)
//      c) Incorporate thermal conductivities
//      d) Incorporate density changes (done)
//      e) Incorporate conservation of mass (done)
//  8)  Implement mantle convection, conditionally
//  9)  Implement metallic core formation
//  10) Implement silicate melting
//

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <functional>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <array>
using namespace std;


/* linspace function as used in MATLAB/Python */
template <typename T>
vector<T> linspace(T a, T b, size_t N) {
    T h = (b - a) / static_cast<T>(N - 1);
    vector<T> xs(N);
    typename vector<T>::iterator x;
    T val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
        *x = val;
    return xs;
}

/* find minimum value in an array */
double findMin(double arr[], int N) {
    int temp = arr[0];
    for (int i = 0; i < N; i++) {
        if (temp > arr[i]) {
            temp = arr[i];
        }
    }
    return temp;
}

int main()
{
    // start recording the runtime of this simulation
    chrono::high_resolution_clock;
    chrono::seconds;
    auto t1 = chrono::high_resolution_clock::now();

    /* Define simulation parameters */
    const int n = 400;  // total number of dr layers
    const double initially_hydrated = 1;   // fraction of silicates that start hydrated

    /* Define temporal parameters */
    double yr2s = 86400 * 365;  // seconds in a year [s]
    double tstart = 3e6 * yr2s;     // formation time [s]
    double tend = 4.5e9 * yr2s; // time at present day [s]

    /* Define physical parameters */
    double initial_radius = 1609.2e3;    // thickness of each layer [m]
    double Tsurf = 273;    // surface temperature [K]
    double Tinit = 273;     // rock-metal initial temperature [K]
    double Tmelt = 273; // melting temperature of water ice [K]
    double dTmelt = 3;  // finite interval where ice melting occurs [K]
    double kc_d = 3.0; // thermal conductivity of anhydrous silicates [W/(m K)]
    double kc_h;    // temperature-dependent thermal conductivity of antigorite [W/(m K)]

    /* Initialize radial position arrays */
    double dr_layer = initial_radius / (n - 1); // radial layer step size. Adjust for multiple actual layers later. [m]
    double r_arr[n], dr_arr[n];
    r_arr[0] = 0;
    dr_arr[0] = 0;
    for (int j = 1; j < n; j++) {
        r_arr[j] = r_arr[j - 1] + dr_layer;
        dr_arr[j] = dr_layer;
    }

    /* Define thermal properties and mass arrays */
    double rho_arr[n], rho_d = 3400, rho_h = 2750, rho_w = 1000, cp_arr[n], cp_h = 1000, cp_d = 900, cp_i = 1930,
        X[n], xlhr = 3.77e5, Tdehyl = 550, Tdehyu = 900, fQl, fQs, wf = 0.13;
    for (int j = 0; j < n; j++) {
        X[j] = initially_hydrated;
        cp_arr[j] = cp_h * X[j] + cp_d * (1 - X[j]);
        rho_arr[j] = pow((X[j] / rho_h + (1 - X[j]) / rho_d), -1);
        //if (initially_hydrated) {
        //    cp_arr[j] = cp_h;   // specific heat of shell [J / (kg K)]
        //    rho_arr[j] = rho_h; // density of shell [kg]
        //    X[j] = 1;   // hydrated mass fraction
        //}
        //else {
        //    cp_arr[j] = cp_d;
        //    rho_arr[j] = rho_d;
        //    X[j] = 0;
        //}
    }
    fQl = 0.548381184555537;
    fQs = 1 - fQl;  // energy fraction going into warming each time step during dehydration

    // calculate mass [kg] and volume [m^3] of each shell
    double V[n], m_arr[n], m_init[n];
    V[0] = 0; m_arr[0] = 0;
    for (int j = 1; j < n; j++) {
        V[j] = 4 * M_PI / 3 * (pow(r_arr[j], 3) - pow(r_arr[j - 1], 3));
        m_arr[j] = V[j] * rho_arr[j];
        m_init[j] = V[j] * rho_arr[j];;
    }

    /* debug */
    double m = 0;
    for (int j = 1; j < n; j++) {
        m += V[j] * rho_arr[j];
    }

    /* initialize temperature array */
    array<double, n> T_arr, T_tmp;
    for (int j = 0; j < (n - 1); j++) {
        T_arr[j] = Tinit;
    }
    T_arr[n - 1] = Tsurf;

    /* Prepare to save a portion of temperature evolution results */
    const int tclip = 20;  // number of timesteps to save
    int ss = 0; // save counter variable

    /* Radiogenic heating parameters */
    double H0[7] = { 568.7e-6, 94.65e-6, 29.17e-6, 26.38e-6, 0.355, 7e-2, 2.7e-2 }; // specific heat production [W/kg]
    double conc[7] = { 5.4e-9, 19.9e-9, 737.9e-9, 38.7e-9, 600e-9, 10e-9, 25.7e-9 }; // concentration of radioactive isotopes [ppb]
    double lambda[7];   // radioactive time constant [1/s]
    double halflife[7] = { 7.04e8, 4.47e9, 1.28e9, 1.4e10, 7.16e5, 1.5e6, 3.7e6 }; // half-life of isotopes [yr]
    for (int i = 0; i < 7; i++) {
        H0[i] *= conc[i];
        lambda[i] = log(2) / halflife[i] / yr2s;
    }

    /* Initialize file to write to */
    ofstream out("T_evol.txt");

    /* Simulate the thermal evolution */
    int max_time_steps = 2e8;
    double t = tstart;  // time [s]
    int ifrz = 0;   // 1 = freeze, 0 = melt
    int I, i_melt;
    double cond_term[n], H_term, dTdt[n], kc_arr[n], Kappa[n], tmp, tmp2, Emelt, h1, h2, rho,
        cp, dr, drd, dru, fbr, fout, dt, qp, Qp[n], Qs, Ql, dmh, dVr, ocean_thickness, f;
    for (int i = 0; i < max_time_steps; i++) {

        // reset T_tmp and dTdt
        for (int j = 0; j < n; j++) {
            T_tmp[j] = T_arr[j];
            cond_term[j] = 0;
            dTdt[j] = 0;
        }

        // update material thermal properties
        for (int j = 0; j < n; j++) {
            kc_h = 1 / (0.404 + 0.000246 * T_arr[j]);   // thermal conductivity of antigorite [W / (m K)]
            f = X[j] * rho_arr[j] / rho_h;  // volume fraction of hydrated silicates
            kc_arr[j] = f * kc_h + (1 - f) * kc_d;  // thermal conductivity of shell [W / (m K)]
            cp_arr[j] = cp_h * X[j] + cp_d * (1 - X[j]);    // specific heat of shell [J / (kg K)]
            Kappa[j] = kc_arr[j] / rho_arr[j] / cp_arr[j];  // thermal diffusivity of shell [m^2 / s]
        }

        // determine maximum timestep, dt [s]
        dt = 0.3 * dr_arr[1] * dr_arr[1] / Kappa[1];
        for (int j = 1; j < n; j++) {
            tmp = 0.3 * dr_arr[j] * dr_arr[j] / Kappa[j];
            if (dt > tmp) dt = tmp;
        }

        // conductive and radiogenic heating solution
        for (int j = 1; j < (n - 1); j++) {

            // determine appropriate density and specific heat
            rho = rho_arr[j];
            cp = cp_arr[j];

            // determine current and neighboring dr step sizes
            dr = dr_arr[j];
            drd = dr_arr[j - 1];
            dru = dr_arr[j + 1];

            // solve for conductive solution
            h1 = pow(r_arr[j] - dr / 2, 2) * (T_arr[j] - T_arr[j - 1]) / (dr / kc_arr[j] + drd / kc_arr[j - 1]);
            h2 = pow(r_arr[j] + dr / 2, 2) * (T_arr[j + 1] - T_arr[j]) / (dru / kc_arr[j + 1] + dr / kc_arr[j]);
            tmp = (2 / rho / cp / dr / r_arr[j] / r_arr[j]) * (h2 - h1);
            cond_term[j] = tmp;
        }

        // solve for radiogenic specific heating (W / kg)
        H_term = 0;
        for (int h = 0; h < 7; h++) {
            H_term += H0[h] * exp(-lambda[h] * t);
        }

        // update temperature profile
        cond_term[0] = cond_term[1];
        for (int j = 1; j < (n - 1); j++) {
            dTdt[j] = cond_term[j] + H_term / cp_arr[j];
            T_tmp[j] += dTdt[j] * dt;
        }
        T_tmp[0] = T_tmp[1];
        T_tmp[n - 1] = Tsurf;

        // silicate dehydration
        if (initially_hydrated) {
            for (int j = 1; j < n; j++) {
                if (T_tmp[j] >= Tdehyl && T_tmp[j] <= Tdehyu && dTdt[j] > 0) {

                    qp = m_arr[j] * cp_arr[j] * dTdt[j] * dt;
                    T_tmp[j] = T_arr[j] + (fQs * qp) / (m_arr[j] * cp_arr[j]);
                    dmh = fQl * qp / xlhr;
                    X[j] = (X[j] * m_arr[j] - dmh) / (m_arr[j] - wf * dmh);
                    m_arr[j] -= wf * dmh;

                    // handle excess energies (for some reason, shells never fully dehydrate)
                    if (T_tmp[j] >= Tdehyu || X[j] <= 0 || m_arr[j] <= (m_init[j] * (1 - wf))) {
                        T_tmp[j] = Tdehyu;
                        X[j] = 0;
                        m_arr[j] = m_init[j] * (1 - wf);
                    }
                }
            }
            X[0] = X[1];
            T_tmp[0] = T_tmp[1];

            // restructure radial profile
            for (int j = 0; j < (n - 1); j++) {
                r_arr[j + 1] = pow(pow(r_arr[j], 3) + (3. / 4. / M_PI) * (X[j + 1] *
                    m_arr[j + 1] / rho_h + (1 - X[j + 1]) * m_arr[j + 1] / rho_d), 1.0 / 3.0);
            }

            // update volume and density
            for (int j = 1; j < n; j++) {
                V[j] = 4 * M_PI / 3 * (pow(r_arr[j], 3) - pow(r_arr[j - 1], 3));
                rho_arr[j] = m_arr[j] / V[j];
                dr_arr[j] = r_arr[j] - r_arr[j - 1];
            }
            rho_arr[0] = rho_arr[1];
        }

        // update temperature array
        for (int j = 0; j < n; j++) T_arr[j] = T_tmp[j];

        /* Save results occasionally */
        if (i % tclip == 0) {

            // write results
            for (int j = 0; j < n; j++) {
                out << log10(t) << "," << r_arr[j] << "," << T_arr[j] << "," << ocean_thickness <<
                    "," << fbr << "," << fout << "," << i_melt << "," << kc_arr[j] << "," << ifrz <<
                    "," << m_arr[j] << "," << i << "," << X[j] << "," << rho_arr[j] << "," << dTdt[j] << ",\n";
            }
            ss++;
        }


        // update time for next iteration
        t += dt;
        if (t > tend) {
            cout << "Break at iteration i = " << i << endl;
            break;
        }
    }

    out << ss << "," << n << "," << n << "," << initially_hydrated << ",0,0,0,0,\n";
    out.close();

    // end time of simulation
    auto t2 = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(t2 - t1);
    cout << "Total run time = " << duration.count() << " s" << endl;
}
