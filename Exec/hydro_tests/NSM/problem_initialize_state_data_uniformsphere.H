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
    
#if AMREX_SPACEDIM == 1
    int nsubx = problem::nsub;
    int nsuby = 1;
    int nsubz = 1;
#elif AMREX_SPACEDIM == 2
    int nsubx = problem::nsub;
    int nsuby = problem::nsub;
    int nsubz = 1;
#else
    int nsubx = problem::nsub;
    int nsuby = problem::nsub;
    int nsubz = problem::nsub;
#endif

    Real ds[3] = {0.0, 0.0, 0.0};
    Real xmin = 0.0;
    Real ymin = 0.0;
    Real zmin = 0.0;

    ds[0] = dx[0] / nsubx;
    xmin = problo[0] + dx[0] * static_cast<Real>(i);

#if AMREX_SPACEDIM >= 2
    ds[1] = dx[1] / nsuby;
    ymin = problo[1] + dx[1] * static_cast<Real>(j);
#endif

#if AMREX_SPACEDIM == 3
    ds[2] = dx[2] / nsubz;
    zmin = problo[2] + dx[2] * static_cast<Real>(k);
#endif

    Real erg_per_Msol_cgs = 5.03e-34;
    Real v_max = sqrt( pow (10.0, problem::exp_energy) * erg_per_Msol_cgs / (3.0 * problem::M_ej / 10.0));
    printf("vmax: %.2e\n", v_max);
    Real r_init = v_max * problem::t_init;
    Real dens_ej = problem::M_ej * 1.988435e33 / (4.0 * M_PI * pow(r_init, 3.0) / 3.0); // uniform sphere

    Real vol_ej = 0.0_rt;
    Real vol_ISM = 0.0_rt;

    for (int kk = 0; kk <= nsubz-1; ++kk) {
        Real zz = zmin + ds[2] * (static_cast<Real>(kk) + 0.5_rt);

        for (int jj = 0; jj <= nsuby-1; ++jj) {
            Real yy = ymin + ds[1] * (static_cast<Real>(jj) + 0.5_rt);

            for (int ii = 0; ii <= nsubx-1; ++ii) {
                Real xx = xmin + ds[0] * (static_cast<Real>(ii) + 0.5_rt);

                Real dist = (problem::center[0] - xx) * (problem::center[0] - xx) +
                            (problem::center[1] - yy) * (problem::center[1] - yy) +
                            (problem::center[2] - zz) * (problem::center[2] - zz);

                bool ej = false;

                if (dist <= r_init * r_init) {
                    ej = true;
                }

                if (coord_type == 1) {

                    // The volume of a cell is a annular cylindrical region.
                    // The main thing that matters is the distance from the
                    // symmetry axis.
                    //   V = pi*dy*(x_r**2 - x_l**2) = pi*dy*dx*0.5*xx
                    // (where x_r is the coordinate of the x right edge,
                    //        x_l is the coordinate of the x left edge,
                    //    and xx  is the coordinate of the x center of the cell)
                    //
                    // since dx and dy are constant, they cancel out

                    if (ej) {
                        vol_ej = vol_ej + xx;
                    }
                    else {
                        vol_ISM = vol_ISM + xx;
                    }

                }
                else if (coord_type == 2) {

                    Real xl = xx - 0.5_rt * ds[0];
                    Real xr = xx + 0.5_rt * ds[0];

                    // the volume of a subzone is (4/3) pi (xr^3 - xl^3).
                    // we can factor this as: (4/3) pi dr (xr^2 + xl*xr + xl^2)
                    // The (4/3) pi dr factor is common, so we can neglect it.

                    if (ej) {
                        vol_ej = vol_ej + (xr * xr + xl * xr + xl * xl);
                    }
                    else {
                        vol_ISM = vol_ISM + (xr * xr + xl * xr + xl * xl);
                    }

                }
                else {

                    // Cartesian -- equal volume of dx * dy * dz is factored out.

                    if (ej) {
                        vol_ej = vol_ej + 1.0_rt;
                    }
                    else {
                        vol_ISM = vol_ISM + 1.0_rt;
                    }

                }

            }
        }
    }

    Real dens_zone = (vol_ej * dens_ej + vol_ISM * problem::dens_ISM) / (vol_ej + vol_ISM);

    Real curr_x = problem::center[0] - (problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt));
    Real curr_y = problem::center[1] - (problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt));
    Real curr_z = 0;
    Real curr_r = sqrt(curr_x * curr_x + curr_y * curr_y + curr_z * curr_z);

    Real curr_v = v_max * (curr_r / r_init) * vol_ej / (vol_ej + vol_ISM); // zero velocity if vol_ej = 0
    Real vx = curr_v * curr_x / curr_r;
    Real vy = curr_v * curr_y / curr_r;

    Real v_min = 1e7;
    Real eint_min = 0.5e0_rt * problem::dens_ISM * v_min * v_min;
    Real ekin = 0.5e0_rt * curr_v * curr_v * dens_zone;
    Real eint = max(0.05 * ekin, eint_min); // setting E_int to 5% of the kinetic energy

    state(i,j,k,URHO) = dens_zone;
    printf ("r_init, r_curr: %.2e, %.2e\n", r_init,curr_r);
    printf ("rho: %.2e\n", dens_zone);
    state(i,j,k,UMX) = dens_zone * vx;
    state(i,j,k,UMY) = dens_zone * vy;
    state(i,j,k,UMZ) = 0.e0_rt;

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
