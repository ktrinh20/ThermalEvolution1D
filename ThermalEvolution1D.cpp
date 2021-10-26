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
//      b) Incorporate specific heats (done, but revisit later for more comprehensive treatment)
//      c) Incorporate thermal conductivities
//      d) Incorporate density changes
//      e) Incorporate conservation of mass
//  8)  Implement mantle convection, conditionally
//  9)  Implement metallic core formation
//  10)  Implement silicate melting
//
//  Description of current state:
//      I'm currently adapting cp_arr[3] to cp_arr[n]

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
    const int n1 = 287 + 1;  // # of shells in deep interior, plus center point
    const int n2 = 106 * 2;  // # of shells in hydrosphere
    const int n = 500;  // total number of dr layers

    /* Define physical parameters */
    double thicknesses[2] = { 1455e3, 106e3 };    // thickness of each layer [m]
    double Tsurf = 100;    // surface temperature [K]
    double Tinit = 273;     // rock-metal initial temperature [K]
    double Tinit2 = 100;    // hydrosphere initial temperature [K]
    double Tmelt = 273; // melting temperature of water ice [K]
    double dTmelt = 3;  // finite interval where ice melting occurs [K]
    double kc_d = 3.0; // thermal conductivity of anhydrous silicates [W/(m K)]

    /* Initialize radial position arrays */
    double dr_layer[2] = { thicknesses[0] / n1, thicknesses[1] / n2 }; // radial layer step size. Adjust for multiple actual layers later. [m]
    double r_arr[n], dr_arr[n];
    r_arr[0] = 0;
    dr_arr[0] = dr_layer[0];
    for (int j = 1; j < n; j++) {
        if (j <= n1) {
            r_arr[j] = r_arr[j - 1] + dr_layer[0];
            dr_arr[j] = dr_layer[0];
        }
        else {
            r_arr[j] = r_arr[j - 1] + dr_layer[1];
            dr_arr[j] = dr_layer[1];
        }
    }

    /* Define temporal parameters */
    double yr2s = 86400 * 365;  // seconds in a year [s]
    double tstart = 5e6 * yr2s;     // formation time [s]
    double tend = 4.5e9 * yr2s; // time at present day [s]

    /* Define thermal properties and mass arrays */
    double rho_arr[n], rho_d = 3500, rho_h = 3500, rho_w = 1000, cp_arr[n], cp_h = 1000, cp_d = 900, cp_i = 1930,
        X[n], xlhr = 3.77e5, Tdehyl = 550, Tdehyu = 900, fQl, fQs;
    for (int j = 0; j < n; j++) {
        if (j <= n1) {
            cp_arr[j] = cp_h;   // specific heat of shell [J / (kg K)]
            rho_arr[j] = rho_h; // density of shell [kg]
            X[j] = 1;   // hydrated mass fraction
        }
        else {
            cp_arr[j] = cp_i;
            rho_arr[j] = rho_w;
            X[j] = 1;
        }
    }
    fQl = xlhr / (xlhr + cp_h * (Tdehyu - Tdehyl)); // energy fraction going into phase change each time step during dehydration
    fQs = 1 - fQl;  // energy fraction going into warming each time step during dehydration

    // recalculate mass [kg] and volume [m^3] of each shell
    double V[n], m_arr[n];
    V[0] = 0; m_arr[0] = 0;
    for (int j = 1; j < n; j++) {
        V[j] = 4 * M_PI / 3 * (pow(r_arr[j], 3) - pow(r_arr[j - 1], 3));
        m_arr[j] = V[j] * rho_arr[j];
    }

    /* initialize temperature array */
    array<double, n> T_arr, T_tmp;
    for (int j = 0; j < (n - 1); j++) {
        if (j <= n1) {
            T_arr[j] = Tinit;
        }
        else {
            T_arr[j] = Tinit2;
        }
    }
    T_arr[n - 1] = Tsurf;

    /* Set energetic constants */
    double xlhi = 3.33e5;    // latent heat density of water ice [J / kg]
    //double xlhr = 3.77e5;   // latent heat density of hydrated silicates [J / kg]
    //double Tdehyl = 550;    // lower bound temp. of silicate dehydration [K]
    //double Tdehyu = 900;    // upper bound temp. of silicate dehydration [K]
    double iceheat0 = rho_w * (cp_i * dTmelt + xlhi); // heat needed for each layer to fully melt [J / m^3]
    //double rockheat0 = rho_h * xlhr; // latent heat needed for each layer to fully dehydration [J / m^3]
    double LI_arr[n]; // heat needed to melt ice layer [J / m^3]
    LI_arr[0] = 0;
    for (int j = 1; j < n; j++) {
        if (j <= n1) {
            //LR_arr[j] = rockheat0;
        }
        else {
            LI_arr[j] = iceheat0;
        }
    }

    /* Prepare to save a portion of temperature evolution results */
    const int tclip = 3000;  // number of timesteps to save
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

    /* Initialize ocean thickness */
    double frac_melt[n]; // percent melted
    for (int j = 0; j < n - 2; j++) {
        frac_melt[j] = 0; // initially liquid hydrosphere
    }
    frac_melt[n - 1] = 0; // frozen surface
    double ocean_thickness; // track ocean thickness over time

    /* Simulate the thermal evolution */
    int max_time_steps = 2e8;
    double t = tstart;  // time [s]
    int ifrz = 0;   // 1 = freeze, 0 = melt
    int I, i_melt;
    double cond_term[n], H_term, dTdt[n], kc_arr[n], Kappa[n], tmp, Emelt, h1, h2, rho,
        cp, dr, drd, dru, isOcean, fbr, fout, dt, qp, Qp[n], Qs, Ql, dmh, dVr;
    for (int i = 0; i < max_time_steps; i++) {

        // reset T_tmp and dTdt
        for (int j = 0; j < n; j++) {
            T_tmp[j] = T_arr[j];
            cond_term[j] = 0;
            dTdt[j] = 0;
        }

        // find index where melting occurs if it exists
        i_melt = 0;
        for (int j = n1 + 1; j < n; j++) {
            if (frac_melt[j] >= 1 && frac_melt[j + 1] <= 0 && ifrz == 1) {
                // completely liquid and solid interface with intent to freeze
                i_melt = j;
                break;
            }
            else if (frac_melt[j] < 1) {
                // above layer is partially melted or completely frozen
                i_melt = j;
                break;
            }
        }

        // recalculate volume of each shell
        V[0] = 0;
        for (int j = 1; j < n; j++) {
            V[j] = 4 * M_PI / 3 * (pow(r_arr[j], 3) - pow(r_arr[j - 1], 3));
            if (j <= n1) {
                V[j] = X[j] * m_arr[j] / rho_h + (1 - X[j]) * m_arr[j] / rho_d;
            }
            else {
                V[j] = m_arr[j] / rho_w;
            }
        }


        // update material thermal properties
        for (int j = 0; j < n; j++) {
            rho_arr[j] = 1 / (X[j] / rho_h + (1 - X[j]) / rho_d);

            if (j <= n1) {
                // specific heat of silicates
                if (T_arr[j] <= Tdehyl) {
                    cp_arr[j] = cp_h;
                    kc_arr[j] = kc_d;
                }
                else if (T_arr[j] < Tdehyu) {
                    cp_arr[j] = cp_h * X[j] + cp_d * (1 - X[j]);
                    kc_arr[j] = kc_d;
                    //kc_arr[j] = kc_d * (1 - heat_arr[j] / rockheat0) + (1 / (0.404 + 0.000246 * T_arr[j])) * (heat_arr[j] / rockheat0);
                }
                else {
                    cp_arr[j] = cp_d;
                    kc_arr[j] = 1 / (0.404 + 0.000246 * T_arr[j]);
                    //kc_arr[j] = kc_d;
                }
            }
            else {
                kc_arr[j] = 0.4685 + 488.12 / T_arr[j];
            }
            Kappa[j] = kc_arr[j] / rho_arr[j] / cp_arr[j];
        }

        // determine maximum timestep, dt [s]
        dt = 0.3 * dr_arr[0] * dr_arr[0] / Kappa[0];
        for (int j = 0; j < n; j++) {
            tmp = 0.3 * dr_arr[j] * dr_arr[j] / Kappa[j];
            if (dt > tmp) dt = tmp;
        }

        // conductive and radiogenic heating solution
        for (int j = 1; j < (n - 1); j++) {


            // skip conduction for ocean and ocean boundaries
            if (j >= n1 && j <= i_melt + 1) continue;

            // determine appropriate density and specific heat
            rho = rho_arr[j];
            if (j <= n1) {
                cp = cp_arr[j];
            }
            else {
                cp = cp_i;
            }

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

        // check to see if there is an ocean
        isOcean = 0;
        for (int j = n1 + 1; j < n; j++) {
            if (frac_melt[j] > 0) {
                isOcean = 1;
                break;
            }
        }

        // conduction at outermost silicates
        I = n1;
        dru = dr_arr[I + 1];
        dr = dr_arr[I];
        drd = dr_arr[I - 1];
        h1 = pow((r_arr[I] - dr / 2), 2) * (T_arr[I] - T_arr[I - 1]) / (dr / kc_arr[I] + drd / kc_arr[I - 1]);
        if (isOcean) {
            // there is an ocean
            h2 = pow((r_arr[I] + (dr / 2)), 2) * (T_arr[I + 1] - T_arr[I]) * kc_arr[I] / (dr / 2);
        }
        else {
            // no ocean
            h2 = pow((r_arr[I] + (dr / 2)), 2) * (T_arr[I + 1] - T_arr[I]) / (dr / kc_arr[I] + dru / kc_arr[I + 1]);
        }
        cond_term[I] = 2 / (rho_arr[I] * cp_arr[I] * dr * r_arr[I] * r_arr[I]) * (h2 - h1);

        // heat flux out of base of the ocean
        fbr = 2 * h2 / pow(r_arr[i_melt], 2);

        // conduction at base of ice shell
        I = i_melt + 1; // last completely melted layer
        dru = dr_arr[I + 1];
        dr = dr_arr[I];
        drd = dr_arr[I - 1];
        if (isOcean) {
            // there is an ocean
            h1 = pow((r_arr[I] - (dr / 2)), 2) * (T_arr[I] - T_arr[I - 1]) * kc_arr[I] / (dr / 2);
        }
        else {
            h1 = pow((r_arr[I] - (dr / 2)), 2) * (T_arr[I] - T_arr[I - 1]) / (dr / kc_arr[I] + drd / kc_arr[I - 1]);
        }
        h2 = pow(r_arr[I] + (dr / 2), 2) * (T_arr[I + 1] - T_arr[I]) / (dru / kc_arr[I + 1] + dr / kc_arr[I]);
        cond_term[I] = 2 / (rho_arr[I] * cp_arr[I] * dr * r_arr[i_melt] * r_arr[i_melt]) * (h2 - h1);

        // heat flux out of the ocean
        fout = 2 * h1 / pow(r_arr[i_melt], 2);

        // power and temperature of partially melted layer
        cond_term[i_melt] = (fout - fbr) / (dr_arr[i_melt] * rho_w * cp_i);

        // update temperature profile
        cond_term[0] = cond_term[1];
        for (int j = 0; j < n; j++) {
            dTdt[j] = cond_term[j];
            if (j <= n1) dTdt[j] += H_term / cp_arr[j];
            T_tmp[j] += dTdt[j] * dt;
        }
        //T_tmp[0] = T_tmp[1];
        T_tmp[n - 1] = Tsurf;

        // simple dehydration for T only (wait, are hbr and hout affected?!)
        for (int j = 1; j <= n1; j++) {
            if (T_tmp[j] >= Tdehyl && T_tmp[j] <= Tdehyu) {
                
                qp = m_arr[j] * cp_arr[j] * dTdt[j] * dt;
                T_tmp[j] = T_arr[j] + (fQs * qp) / (m_arr[j] * cp_arr[j]); // including m_arr makes the code take forever
                dmh = fQl * qp / xlhr;
                X[j] = (X[j] * m_arr[j] - dmh) / (m_arr[j] - 0.13 * dmh);
                m_arr[j] -= 0.13 * dmh;

                // handle excess energies
                if (T_tmp[j] >= Tdehyu || X[j] <= 0) {
                    T_tmp[j] = Tdehyu;
                    X[j] = 0;
                    //m_arr[j] = V[j] * rho_d;
                }
                else if (T_tmp[j] <= Tdehyl || X[j] >= 1) {
                    T_tmp[j] = Tdehyl;
                    X[j] = 1;
                    //m_arr[j] = V[j] * rho_h;
                }
            }
        }
        X[0] = X[1];
        T_tmp[0] = T_tmp[1];

        // restructure radial profile


        // melting or freezing occurs
        if ((T_tmp[i_melt] >= (Tmelt - dTmelt)) || frac_melt[i_melt] > 0) {

            // if ice-layer is at least partially melted...
            if (frac_melt[i_melt] > 0) {
                Emelt = -(fbr - fout) * dt / dr;
            }
            else if (T_tmp[i_melt] > (Tmelt - dTmelt)) {
                // i_melt should completely melt
                Emelt = (T_tmp[i_melt] - (Tmelt - dTmelt)) * rho_w * cp_i;
            }
            else {
                cout << "Error occured when calculating energy for ice-water interface!" << endl;
            }

            if (Emelt < 0) ifrz = 1;    // freeze next time step
            else if (Emelt > 0) ifrz = 0;   // melt next time step

            // melt/freeze stuff!
            LI_arr[i_melt] = LI_arr[i_melt] - Emelt;
            frac_melt[i_melt] = 1 - LI_arr[i_melt] / iceheat0;

            // handle excess/deficit heat
            if (frac_melt[i_melt] > 1) {
                LI_arr[i_melt] = 0;
                frac_melt[i_melt] = 1;
            }
            else if (frac_melt[i_melt] < 0) {
                LI_arr[i_melt] = iceheat0;
                frac_melt[i_melt] = 0;
            }
            T_tmp[i_melt] = Tmelt - dTmelt + frac_melt[i_melt] * dTmelt;
        }

        // update temperature array
        for (int j = 0; j < n; j++) T_arr[j] = T_tmp[j];

        /* Save results occasionally */
        if (i % tclip == 0) {

            // calculate ocean thickness
            ocean_thickness = 0;
            for (int j = n1; j < n; j++) {
                ocean_thickness = ocean_thickness + frac_melt[j] * dr_layer[1];
            }

            // write results
            for (int j = 0; j < n; j++) {

                // scale heat fluxes to surfaces
                fbr = fbr * r_arr[i_melt] * r_arr[i_melt] / r_arr[n - 1] / r_arr[n - 1];
                fout = fbr * r_arr[i_melt] * r_arr[i_melt] / r_arr[n - 1] / r_arr[n - 1];

                // write to file
                out << log10(t) << "," << r_arr[j] << "," << T_arr[j] << "," << ocean_thickness <<
                    "," << fbr << "," << fout << "," << i_melt << "," << kc_arr[j] << "," << ifrz <<
                    "," << isOcean << "," << i << "," << X[j] << "," << rho_arr[j] << "," << dTdt[j] << ",\n";
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

    out << ss << "," << n << "," << n1 << ",0,0,0,0,0,\n";
    out.close();

    // end time of simulation
    auto t2 = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(t2 - t1);
    cout << "Total run time = " << duration.count() << " s" << endl;
}
