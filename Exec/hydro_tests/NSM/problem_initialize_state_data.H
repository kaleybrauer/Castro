#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H
#define _USE_MATH_DEFINES

#include <prob_parameters.H>
#include <eos.H>
#include <math.h>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k, Array4<Real> const& state, const GeometryData& geomdata)
{
    int coord_type = geomdata.Coord();

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    // current coordinate location of the middle of this grid cell
    Real x = problem::center[0] - (problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt));
    Real y = 0.0;
    Real z = 0.0;

#if AMREX_SPACEDIM >= 2
    y = problem::center[1] - (problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt));
#endif

#if AMREX_SPACEDIM == 3
    z = problem::center[2] - (problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt));
#endif

    Real r = sqrt(x * x + y * y + z * z);
    Real v = r / problem::t_init; // homologous (free) expansion

    Real zetaV = 1.69;
    Real zetaRho = 0.12;
    Real M_chand = 1.4; // chandrasekhar mass in Msun
    Real v_t = 6e8 * zetaV * sqrt(pow(10.0, problem::exp_energy - 51.0) / (problem::M_ej / M_chand)); // approx 0.13c

    Real dens_exp;
    if (v < v_t)
        dens_exp = 1;
    else
        dens_exp = 10;

    Real Msun = 1.99e33; // solar mass in grams
    Real dens_ej = zetaRho * (problem::M_ej * Msun / pow((v_t * problem::t_init), 3.0)) * pow(r / (v_t * problem::t_init), -dens_exp);

    Real dens_zone;
    Real vx = 0.0;
    Real vy = 0.0;
    Real vz = 0.0;

    // to identify if the cell is ejecta or not
    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = 0.0_rt;
    }

    if (dens_ej > problem::dens_ISM)
    {
        dens_zone = dens_ej;
        vx = x / problem::t_init;
        vy = y / problem::t_init;
        vz = z / problem::t_init;

        state(i,j,k,UFS) = dens_zone;
        // state(i,j,k,UFS+1) = 0.0_rt;
    }
    else {
        dens_zone = problem::dens_ISM; // only extend broken power law until it equals rho_ISM
        // state(i,j,k,UFS) = 0.0_rt;
        state(i,j,k,UFS+1) = dens_zone;
    }

    if ((i == 21) and (j == 0)) {
        printf("i, j, k: %d, %d, %d; r:%E\n", i,j,k,r);
        printf("UFS0:%E, UFS1:%E\n", state(i,j,k,UFS), state(i,j,k,UFS+1));
    }
    // printf("vx:%E\n", vx);
    // printf("vy:%E\n", vy);

    if (coord_type != 0) {
        vx = sqrt(vx * vx);
        vy = sqrt(vy * vy);
        vz = sqrt(vz * vz);
    }

    Real v_min = 1e7;
    Real eint_min = 0.5e0_rt * problem::dens_ISM * v_min * v_min;
    Real ekin = 0.5e0_rt * (vx * vx + vy * vy + vz * vz) * dens_zone;
    Real eint = max(0.05 * ekin, eint_min); // setting E_int to 5% of the kinetic energy

    state(i,j,k,URHO) = dens_zone;
    state(i,j,k,UFS) = state(i,j,k,URHO) - state(i,j,k,UFS+1);
    state(i,j,k,UMX) = dens_zone * vx;
    state(i,j,k,UMY) = dens_zone * vy;
    state(i,j,k,UMZ) = dens_zone * vz;

    // state(i,j,k,UEDEN) = eint +
    //                      0.5e0_rt * (state(i,j,k,UMX) * state(i,j,k,UMX) / state(i,j,k,URHO) +
    //                                  state(i,j,k,UMY) * state(i,j,k,UMY) / state(i,j,k,URHO) +
    //                                  state(i,j,k,UMZ) * state(i,j,k,UMZ) / state(i,j,k,URHO));

    state(i,j,k,UEDEN) = ekin + eint;
    state(i,j,k,UEINT) = eint;


    // The temperature initialization will be done later
    // in the call to clean_state. We want to avoid
    // EOS calls in ca_initdata. So that there is valid
    // data for an initial guess for the EOS, we initialize
    // to the ISM temperature everywhere.

    state(i,j,k,UTEMP) = problem::temp_ISM;

    state(i,j,k,UFS) = state(i,j,k,URHO);
}

#endif
