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

    /* Define physical parameters */
    double thicknesses[2] = { 1455e3, 106e3 };    // thickness of each layer [m]
    double Tsurf = 100;    // surface temperature [K]
    double Tinit = 273;     // rock-metal initial temperature [K]
    double Tinit2 = 100;    // hydrosphere initial temperature [K]
    double Tmelt = 273; // melting temperature of water ice [K]
    double dTmelt = 3;  // finite interval where ice melting occurs [K]
    double rho_arr[2] = { 3500, 1000 };  // density of layers [kg/m^3]
    //double cp_arr[3] = { 1000, 840, 1930 }; // specificy heat of layers [J/kg/K]
    double kc_s = 3.0; // thermal conductivity of silicates [W/(m K)]

    /* Define temporal parameters */
    double yr2s = 86400 * 365;  // seconds in a year [s]
    double tstart = 5e6 * yr2s;     // formation time [s]
    double tend = 4.5e9 * yr2s; // time at present day [s]

    /* Define simulation parameters */
    const int n1 = 1455 + 1;  // # of layers in deep interior
    const int n2 = 106 * 2;  // # of layers in hydrosphere
    const int n = 1668;  // total number of dr layers

    /* Define specific heat stuff */
    double cp_h = 1000; // specific heat of hydrated silicates [J/kg/K]
    double cp_d = 840;  // specific heat of dehydrated silicates [J/kg/K]
    double cp_i = 1930; // specific heat of water ice [J/kg/K]
    double cp_arr[n];
    for (int j = 0; j < n; j++) {
        if (j <= n1) {
            cp_arr[j] = cp_h;
        }
        else {
            cp_arr[j] = cp_i;
        }
    }

    /* initialize temperature array*/
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

    /* Set energetic constants */
    double xlhi = 3.33e5;    // latent heat density of water ice [J / kg]
    double xlhr = 3.77e5;   // latent heat density of hydrated silicates [J / kg]
    double Tdehyl = 550;    // lower bound temp. of silicate dehydration [K]
    double Tdehyu = 900;    // upper bound temp. of silicate dehydration [K]
    double iceheat0 = rho_arr[1] * (cp_i * dTmelt + xlhi); // heat needed for each layer to fully melt [J]
    double rockheat0 = rho_arr[0] * (cp_h * (Tdehyu - Tdehyl) + xlhr); // heat needed for each layer to fully dehydration [J]
    double heat_arr[n]; // heat needed to complete a phase change [J]
    for (int j = n1 + 1; j < n; j++) {
        if (j <= n1) {
            heat_arr[j] = rockheat0;
        }
        else {
            heat_arr[j] = iceheat0;
        }
    }

    /* Partition heat production for silicate dehydration/hydration */
    double dehys = rho_arr[0] * cp_h * (Tdehyu - Tdehyl); // total heat needed to warm throughout dehydration [J]
    double dehyl = rho_arr[0] * xlhr; // total latent heat needed to dehydrate a layer [J]
    double fracs = dehys / rockheat0; // fraction of heat dedicated to secular warming during dehydration [unitless]
    double fracl = dehyl / rockheat0; // fraction of heat dedicated to latent heat of dehydration [unitless]

    /* Prepare to save a portion of temperature evolution results */
    const int tclip = 4000;  // number of timesteps to save
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
        cp, dr, drd, dru, isOcean, fbr, fout, dt, Qprod, Qsec, Ql;
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

        // update rock thermal properties
        for (int j = 0; j < n; j++) {
            if (j <= n1) {
                kc_arr[j] = kc_s;
                Kappa[j] = kc_s / rho_arr[0] / cp_arr[j];

                // specific heat of silicates
                if (T_arr[j] < Tdehyl) {
                    cp_arr[j] = cp_h;
                }
                else if (T_arr[j] <= Tdehyu) {
                    cp_arr[j] = cp_h - (cp_h - cp_d) * (Tdehyu - T_arr[j]) / (Tdehyu - Tdehyl);
                }
                else {
                    cp_arr[j] = cp_d;
                }
            }
            else {
                kc_arr[j] = 0.4685 + 488.12 / T_arr[j];
                Kappa[j] = kc_arr[j] / rho_arr[1] / cp_i;
            }
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
            if (j <= n1) {
                rho = rho_arr[0];
                cp = cp_arr[j];
            }
            else {
                rho = rho_arr[1];
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
        cond_term[I] = 2 / (rho_arr[0] * cp_arr[I] * dr * r_arr[I] * r_arr[I]) * (h2 - h1);

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
        cond_term[I] = 2 / (rho_arr[1] * cp_arr[I] * dr * r_arr[i_melt] * r_arr[i_melt]) * (h2 - h1);

        // heat flux out of the ocean
        fout = 2 * h1 / pow(r_arr[i_melt], 2);

        // power and temperature of partially melted layer
        cond_term[i_melt] = (fout - fbr) / (dr_arr[i_melt] * rho_arr[1] * cp_i);

        // update temperature profile
        for (int j = 0; j < n; j++) {
            dTdt[j] = cond_term[j];
            if (j <= n1) dTdt[j] += H_term / cp_arr[j];
            T_tmp[j] += dTdt[j] * dt;
        }
        T_tmp[0] = T_tmp[1];
        T_tmp[n - 1] = Tsurf;

        // hydration or dehydration occurs
        for (int j = 0; j <= n1; j++) {

            // if dehydration should occur (Tdehyl <= T <= Tdehyu)
            if (T_tmp[j] >= Tdehyl && T_tmp[j] <= Tdehyu) {

                // partition the heat produced to secular warming and latent heat
                Qprod = rho_arr[0] * cp_arr[j] * dTdt[j] * dt; // heat produced this timestep [J]
                Qsec = Qprod * fracs; // heat dedicated to secular warming [J]
                Ql = Qprod * fracl; // heat dedicated to latent heat of dehydration [J]

                // add/subtract from heat_arr[j]
                heat_arr[j] += Ql;

                // adjust temperature
                T_tmp[j] = T_arr[j] + Qsec / rho_arr[0] / cp_arr[j]; // update this when we have hydration dependent rho, cp

                // account for excess energy
                if (heat_arr[j] < 0 && Qprod > 0) {
                    heat_arr[j] = 0;
                    T_tmp[j] = Tdehyu;
                }
                else if (heat_arr[j] > rockheat0 && Qprod < 0) {
                    heat_arr[j] = rockheat0;
                    T_tmp[j] = Tdehyl;
                }

            }
        }

        // melting or freezing occurs
        if ((T_tmp[i_melt] >= (Tmelt - dTmelt)) || frac_melt[i_melt] > 0) {

            // if ice-layer is at least partially melted...
            if (frac_melt[i_melt] > 0) {
                Emelt = -(fbr - fout) * dt / dr;
            }
            else if (T_tmp[i_melt] > (Tmelt - dTmelt)) {
                // i_melt should completely melt
                Emelt = (T_tmp[i_melt] - (Tmelt - dTmelt)) * rho_arr[1] * cp_i;
            }
            else {
                cout << "Error occured when calculating energy for ice-water interface!" << endl;
            }

            if (Emelt < 0) ifrz = 1;    // freeze next time step
            else if (Emelt > 0) ifrz = 0;   // melt next time step

            // melt/freeze stuff!
            heat_arr[i_melt] = heat_arr[i_melt] - Emelt;
            frac_melt[i_melt] = 1 - heat_arr[i_melt] / iceheat0;

            // handle excess/deficit heat
            if (frac_melt[i_melt] > 1) {
                heat_arr[i_melt] = 0;
                frac_melt[i_melt] = 1;
            }
            else if (frac_melt[i_melt] < 0) {
                heat_arr[i_melt] = iceheat0;
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
                    "," << isOcean << "," << i << ",\n";
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

