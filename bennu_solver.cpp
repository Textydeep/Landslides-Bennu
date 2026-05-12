// =============================================================================
// bennu_solver.cpp
//
// Finite-volume solver for shallow granular flow on a rotating top-shaped
// asteroid (modelled as a double cone). Solves a hyperbolic system in
// downslope coordinates using a non-oscillatory central (Kurganov-Tadmor)
// scheme with a minmod slope limiter and a two-stage Runge-Kutta integrator.
//
// Conserved variables (per unit area in the meridional plane):
//   p = h * tau          (mass)
//   q = h * u * tau      (downslope momentum)
//   r = h * v * tau      (azimuthal momentum)
// where h = flow height, u = downslope velocity, v = azimuthal velocity,
// and tau = x * sin(zeta) is the radial coordinate.
//
// Runtime parameters are read from bennu_input.txt.
// The gravity field (grav_file, eta_file) must be generated beforehand
// using gravity_field.ipynb.
// Compile-time grid size RES is set in config.h.
// =============================================================================

#include <iostream>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include "config.h"

using namespace std;

// =============================================================================
// Global arrays  (static size: RES active cells + 4 ghost cells)
// =============================================================================

double grav[RES + 4];   // normalised gravity magnitude at each surface point
double eta[RES + 4];    // angle of gravity vector from the axial direction (rad)

// Primitive variables
double h[RES + 4], htemp[RES + 4];   // flow height
double u[RES + 4], utemp[RES + 4];   // downslope velocity
double v[RES + 4], vtemp[RES + 4];   // azimuthal velocity

// Conserved variables and their RK stage copies
double p[RES + 4], ptemp[RES + 4], ptemphat[RES + 4];
double q[RES + 4], qtemp[RES + 4], qtemphat[RES + 4];
double r[RES + 4], rtemp[RES + 4], rtemphat[RES + 4];

// Maximum wave speed at each cell face (for CFL)
double face_speed[RES + 4];

// =============================================================================
// Global scalars — set by load_parameters() then used throughout
// =============================================================================

// Numerical
int    output_interval;
double t_final;
double limiter_theta;
double dx, dt;
const  double rk_weight = 0.5;    // weight for the RK corrector (0.5 = standard TVD-RK2)

// Geometry
double apex_half_angle;    // xi:  cone half-apex angle (rad)
double semi_apex_angle;    // zeta = pi/2 - xi
double domain_length;      // normalised downslope extent (xlim)
double domain_start;       // = 1 - domain_length (xstart)

// Physical
double friction_angle;     // Coulomb friction angle delta (rad)
double omega;              // current non-dimensional spin rate
double aspect_ratio;       // shallowness parameter epsilon
double regolith_depth;     // initial regolith layer thickness (non-dim)
double density_ratio;      // core-to-regolith density ratio rho_bar
double gravity_scale;      // scale factor for the loaded gravity field
double fillet_length;      // cone-junction fillet radius Dc
double cone_slant_height;  // current cone slant height (non-dim)

// Derived
double earth_pressure_k;   // lateral earth pressure coefficient Kx
int    n_avalanches;
string grav_file, eta_file, output_prefix;

// Precomputed trig constants (set once in main after parameter loading)
double SIN_ZETA, COS_ZETA;    // sin/cos of semi_apex_angle
double SIN_XI;                 // sin of apex_half_angle
double TAN_DELTA;              // tan of friction_angle

// Per-cell trig arrays (set in load_gravity_field after eta[] is loaded)
double cos_xi_eta[RES + 4];   // cos(apex_half_angle - eta[j])
double sin_xi_eta[RES + 4];   // sin(apex_half_angle - eta[j])

// Time / step tracking
int    step_count    = 0;
double t_cumulative  = 0.0;
double d_omega       = 0.0;   // change in omega over one time step
double ang_accel     = 0.0;   // dOmega/dt
double ang_accel_moi = 0.0;   // contribution from MOI change
double ang_accel_mom = 0.0;   // contribution from angular momentum flux

// Diagnostics
double ke_current    = 0.0;
double total_mass    = 0.0;
double total_momentum = 0.0;

// Angular momentum integrals
double moi_regolith         = 0.0, moi_regolith_prev    = 0.0;
double ang_mom_regolith     = 0.0, ang_mom_regolith_prev = 0.0;
double total_mass_prev      = 0.0;
double moi_core;    // moment of inertia of solid cone core

// Local state — set inside loops before calling flux/source functions.
// These globals are side-effects of compute_stress_state() and are read by
// flux(), S2(), S3() in the same sequential evaluation context.
double x_loc;           // current downslope coordinate
double tau_loc;         // current radial coordinate tau = x * sin(zeta)
double v_azim;          // local azimuthal velocity at evaluation point
double basal_stress;    // effective normal stress at flow base (phi)
double normal_gravity;  // gravity component normal to slope surface (Gn)
double tangential_gravity; // gravity component along slope (Gt)

// Interface reconstruction values — set by reconstruct_face_values()
double Qp_p, Qm_p;   // p reconstructed to face: plus (right) and minus (left)
double Qp_q, Qm_q;
double Qp_r, Qm_r;

// Flux temporaries
double f1, f2, f3, fp, fm;
double max_wave_speed;

// =============================================================================
// Function declarations
// =============================================================================

void   load_parameters(const string& filename);
void   load_gravity_field();

double interp_p(int jp);
double interp_q(int jp);
double interp_r(int jp);

void   compute_stress_state(int side, int jj);
void   reconstruct_face_values(int jj);

double max_face_speed(int jj);
double max_eigenvalue(int jj, int side);

double flux(double pp, double qq, double rr, int jj, int component, int side);
double central_flux(int jj, int component);
double limited_slope(int jp, int component);
double minmod(double a, double b, double c);

double S1(int jp);
double S2(int jp);
double S3(int jp);

// =============================================================================
// main
// =============================================================================

int main()
{
    load_parameters("bennu_input.txt");

    // Derived constants (computed once after loading)
    domain_start     = 1.0 - domain_length;
    dx               = domain_length / RES;
    dt               = dx / 4.0;
    semi_apex_angle  = PI / 2.0 - apex_half_angle;
    earth_pressure_k = 2.0 / (cos(friction_angle) * cos(friction_angle)) - 1.0;

    // Cache trig values — these appear in every inner-loop call
    SIN_ZETA  = sin(semi_apex_angle);
    COS_ZETA  = cos(semi_apex_angle);
    SIN_XI    = sin(apex_half_angle);
    TAN_DELTA = tan(friction_angle);

    load_gravity_field();

    // -------------------------------------------------------------------------
    // Initial conditions: uniform regolith layer at rest
    // -------------------------------------------------------------------------
    total_mass = 0.0;
    for (int j = 2; j <= RES + 1; j++)
    {
        x_loc   = domain_start + domain_length * (j - 1.5) / RES;
        tau_loc = x_loc * SIN_ZETA;
        h[j] = regolith_depth;
        u[j] = 0.0;
        v[j] = 0.0;
        p[j] = h[j] * tau_loc;
        q[j] = h[j] * u[j] * tau_loc;
        r[j] = h[j] * v[j] * tau_loc;
        total_mass += 2.0 * PI * p[j] * dx;
    }

    // -------------------------------------------------------------------------
    // Open output files
    // -------------------------------------------------------------------------
    ofstream file_h, file_u, file_v, file_omega, file_time, file_mass, file_ke;
    file_h    .open(output_prefix + "_h.txt");
    file_u    .open(output_prefix + "_u.txt");
    file_v    .open(output_prefix + "_v.txt");
    file_omega.open(output_prefix + "_omega.txt");
    file_time .open(output_prefix + "_time.txt");
    file_mass .open(output_prefix + "_mass.txt");
    file_ke   .open(output_prefix + "_ke.txt");

    // =========================================================================
    // Outer loop: successive avalanche / landslide events
    // =========================================================================
    for (int i_layer = 0; i_layer < n_avalanches; i_layer++)
    {
        // Reset local time; omega and t_cumulative persist across events
        double t = 0.0;

        // Add a fresh regolith layer on top of any residual from the previous event
        if (i_layer != 0)
        {
            for (int j = 2; j <= RES + 1; j++)
            {
                x_loc   = domain_start + domain_length * (j - 1.5) / RES;
                tau_loc = x_loc * SIN_ZETA;
                h[j] = regolith_depth + h[j];
                u[j] = 0.0;
                v[j] = 0.0;
                p[j] = h[j] * tau_loc;
                q[j] = h[j] * u[j] * tau_loc;
                r[j] = h[j] * v[j] * tau_loc;
            }
        }

        // Update core MOI for the eroded cone (slant height shrinks each event)
        double base_radius  = cone_slant_height * cos(apex_half_angle);
        double br2 = base_radius * base_radius, br4 = br2 * br2;
        static const double inner_r5 = []{
            double r = 241.8 / 376.0 * cos(PI / 4.0); double r2 = r*r; return r2*r2*r;
        }();
        moi_core = (density_ratio - 1.0) * PI * inner_r5 / 10.0
                   + PI * br4 * base_radius / 10.0;
        cone_slant_height -= 2.0 * aspect_ratio * regolith_depth;

        d_omega = 0.0;   // no spin-rate change on the first step of each event

        // =====================================================================
        // Time-stepping loop
        // =====================================================================
        while (t < t_final)
        {
            omega += d_omega;

            // --- Periodic output -------------------------------------------
            if (step_count % output_interval == 0)
            {
                for (int j = 2; j <= RES + 1; j++)
                {
                    file_h << h[j] << "\n";
                    file_u << u[j] << "\n";
                    file_v << v[j] << "\n";
                }
                file_mass  << total_mass    << "\n";
                file_ke    << ke_current    << "\n";
                file_omega << omega         << "\n";
                file_time  << t_cumulative  << "\n";
            }

            // --- Build temporary conserved / primitive variables -----------
            for (int j = 2; j <= RES + 1; j++)
            {
                x_loc   = domain_start + domain_length * (j - 1.5) / RES;
                tau_loc = x_loc * SIN_ZETA;
                utemp[j]     = q[j] / p[j];
                vtemp[j]     = r[j] / p[j];
                ptemp[j]     = max(p[j], 1.0e-16 * tau_loc);
                htemp[j]     = ptemp[j] / tau_loc;
                qtemp[j]     = ptemp[j] * utemp[j];
                rtemp[j]     = ptemp[j] * vtemp[j];
                ptemphat[j]  = ptemp[j];
                qtemphat[j]  = qtemp[j];
                rtemphat[j]  = rtemp[j];
            }

            // --- Ghost cells for predictor (wall BC at pole and equator) ---
            // Pole (j = 0, 1): reflect u (normal), extrapolate v (tangential)
            ptemp[0] = ptemp[3]; qtemp[0] = -qtemp[3]; rtemp[0] = rtemp[3];
            htemp[0] = htemp[3]; utemp[0] = -utemp[3];
            vtemp[0] = (1.0 + 2.0 * limiter_theta) * vtemp[2] - 2.0 * limiter_theta * vtemp[3];

            ptemp[1] = ptemp[2]; qtemp[1] = -qtemp[2]; rtemp[1] = rtemp[2];
            htemp[1] = htemp[2]; utemp[1] = -utemp[2];
            vtemp[1] = (1.0 + limiter_theta) * vtemp[2] - limiter_theta * vtemp[3];

            // Equator (j = RES+2, RES+3): reflect u, extrapolate v
            ptemp[RES + 2] = ptemp[RES + 1]; qtemp[RES + 2] = -qtemp[RES + 1]; rtemp[RES + 2] = rtemp[RES + 1];
            htemp[RES + 2] = htemp[RES + 1]; utemp[RES + 2] = -utemp[RES + 1];
            vtemp[RES + 2] = (1.0 + limiter_theta) * vtemp[RES + 1] - limiter_theta * vtemp[RES];

            ptemp[RES + 3] = ptemp[RES]; qtemp[RES + 3] = -qtemp[RES]; rtemp[RES + 3] = rtemp[RES];
            htemp[RES + 3] = htemp[RES]; utemp[RES + 3] = -utemp[RES];
            vtemp[RES + 3] = (1.0 + 2.0 * limiter_theta) * vtemp[RES + 1] - 2.0 * limiter_theta * vtemp[RES];

            // Store previous-step angular momentum integrals
            moi_regolith_prev = ang_mom_regolith_prev = total_mass_prev = 0.0;
            for (int j = 2; j <= RES + 1; j++)
            {
                x_loc   = domain_start + domain_length * (j - 1.5) / RES;
                tau_loc = x_loc * SIN_ZETA;
                double xs   = x_loc * SIN_XI;
                double xs2  = xs * xs, xs3 = xs2 * xs;
                total_mass_prev       += 2.0 * PI * p[j] * dx;
                moi_regolith_prev     += 2.0 * PI * xs3 * aspect_ratio * p[j] / tau_loc * dx;
                ang_mom_regolith_prev += 2.0 * PI * xs2 * aspect_ratio * r[j] / tau_loc * dx;
            }

            // --- Predictor step (first RK stage) ---------------------------
            for (int j = 2; j <= RES + 1; j++)
            {
                x_loc   = domain_start + domain_length * (j - 1.5) / RES;
                tau_loc = x_loc * SIN_ZETA;
                p[j] = ptemp[j] - ((central_flux(j + 1, 1) - central_flux(j, 1)) / dx - S1(j)) * dt;
                q[j] = qtemp[j] - ((central_flux(j + 1, 2) - central_flux(j, 2)) / dx - S2(j)) * dt;
                r[j] = rtemp[j] - ((central_flux(j + 1, 3) - central_flux(j, 3)) / dx - S3(j)) * dt;
            }

            // Update temporary primitives after predictor
            for (int j = 2; j <= RES + 1; j++)
            {
                x_loc   = domain_start + domain_length * (j - 1.5) / RES;
                tau_loc = x_loc * SIN_ZETA;
                utemp[j] = q[j] / p[j];
                vtemp[j] = r[j] / p[j];
                ptemp[j] = max(p[j], 1.0e-16 * tau_loc);
                htemp[j] = ptemp[j] / tau_loc;
                qtemp[j] = ptemp[j] * utemp[j];
                rtemp[j] = ptemp[j] * vtemp[j];
            }

            // Ghost cells after predictor (wall BCs, no v extrapolation needed)
            ptemp[0] = ptemp[3]; qtemp[0] = -qtemp[3]; rtemp[0] = rtemp[3];
            ptemp[1] = ptemp[2]; qtemp[1] = -qtemp[2]; rtemp[1] = rtemp[2];
            ptemp[RES + 2] = ptemp[RES + 1]; qtemp[RES + 2] = -qtemp[RES + 1]; rtemp[RES + 2] = rtemp[RES + 1];
            ptemp[RES + 3] = ptemp[RES];     qtemp[RES + 3] = -qtemp[RES];     rtemp[RES + 3] = rtemp[RES];

            // --- Corrector step (second RK stage) --------------------------
            for (int j = 2; j <= RES + 1; j++)
            {
                x_loc   = domain_start + domain_length * (j - 1.5) / RES;
                tau_loc = x_loc * SIN_ZETA;
                p[j] = rk_weight * ptemphat[j]
                     + (1.0 - rk_weight) * (ptemp[j] - ((central_flux(j + 1, 1) - central_flux(j, 1)) / dx - S1(j)) * dt);
                q[j] = rk_weight * qtemphat[j]
                     + (1.0 - rk_weight) * (qtemp[j] - ((central_flux(j + 1, 2) - central_flux(j, 2)) / dx - S2(j)) * dt);
                r[j] = rk_weight * rtemphat[j]
                     + (1.0 - rk_weight) * (rtemp[j] - ((central_flux(j + 1, 3) - central_flux(j, 3)) / dx - S3(j)) * dt);

                // Damp velocity sign-reversals that exceed machine precision
                if (q[j] * qtemphat[j] < 0.0)
                    q[j] = (q[j] > 0.0 ? 1.0 : -1.0) * p[j] * 1.0e-8;
                if (r[j] * rtemphat[j] < 0.0)
                    r[j] = (r[j] > 0.0 ? 1.0 : -1.0) * p[j] * 1.0e-8;
            }

            // --- Angular momentum conservation: compute d_omega ------------
            moi_regolith = ang_mom_regolith = total_mass = 0.0;
            for (int j = 2; j <= RES + 1; j++)
            {
                x_loc   = domain_start + domain_length * (j - 1.5) / RES;
                tau_loc = x_loc * SIN_ZETA;
                double xs   = x_loc * SIN_XI;
                double xs2  = xs * xs, xs3 = xs2 * xs;
                total_mass       += 2.0 * PI * p[j] * dx;
                moi_regolith     += 2.0 * PI * xs3 * aspect_ratio * p[j] / tau_loc * dx;
                ang_mom_regolith += 2.0 * PI * xs2 * aspect_ratio * r[j] / tau_loc * dx;
            }
            d_omega = -(omega * (moi_regolith - moi_regolith_prev)
                        + (ang_mom_regolith - ang_mom_regolith_prev))
                       / (moi_core + moi_regolith);
            ang_accel     = d_omega / dt;
            ang_accel_moi = -(omega * (moi_regolith - moi_regolith_prev))
                             / (moi_core + moi_regolith) / dt;
            ang_accel_mom = -(ang_mom_regolith - ang_mom_regolith_prev)
                             / (moi_core + moi_regolith) / dt;

            // --- Mass shedding: drop cells where basal stress > 0 ---------
            // Interior cells
            for (int j = 2; j <= RES; j++)
            {
                x_loc   = domain_start + domain_length * (j - 1.5) / RES;
                tau_loc = x_loc * SIN_ZETA;
                compute_stress_state(3, j);
                h[j] = p[j] / tau_loc;
                u[j] = q[j] / p[j];
                v[j] = r[j] / p[j];
                if (basal_stress > 0.0)
                {
                    p[j] = 1.0e-8 * tau_loc;
                    q[j] = p[j] * u[j];
                    r[j] = p[j] * v[j];
                }
            }
            // Equatorial cell: additional longitudinal curvature contribution
            {
                int j  = RES + 1;
                x_loc   = domain_start + domain_length * (j - 1.5) / RES;
                tau_loc = x_loc * SIN_ZETA;
                compute_stress_state(3, j);
                double u_eq = q[j] / p[j];
                double v_eq = r[j] / p[j];
                if (basal_stress + 2.0 * u_eq * u_eq * SIN_ZETA / fillet_length > 0.0)
                {
                    p[j] = 1.0e-8 * tau_loc;
                    q[j] = p[j] * u_eq;
                    r[j] = p[j] * v_eq;
                }
            }

            // --- Reconstruct primitives from conserved variables -----------
            ke_current = total_momentum = 0.0;
            for (int j = 2; j <= RES + 1; j++)
            {
                x_loc   = domain_start + domain_length * (j - 1.5) / RES;
                tau_loc = x_loc * SIN_ZETA;
                h[j] = p[j] / tau_loc;
                u[j] = q[j] / p[j];
                v[j] = r[j] / p[j];
                ke_current    += 0.5 * (u[j] * u[j] + v[j] * v[j]);
                total_momentum += abs(h[j] * u[j] * tau_loc)
                                + abs(h[j] * (abs(v[j]) > 1.0e-5 ? v[j] : 0.0) * tau_loc);
            }

            // --- CFL condition: update dt from maximum wave speed ----------
            max_wave_speed = 1.0e-10;
            for (int j = 1; j <= RES + 1; j++)
            {
                face_speed[j] = max_face_speed(j);
                if (abs(face_speed[j]) > max_wave_speed)
                    max_wave_speed = abs(face_speed[j]);
            }
            dt = min(dx / 4.0, dx / 4.0 / max_wave_speed);

            // --- Advance time counters ------------------------------------
            t            += dt;
            t_cumulative += dt;
            step_count++;

            cout << fixed << setprecision(10)
                 << "t=" << t_cumulative
                 << "  step=" << step_count
                 << "  alpha=" << ang_accel
                 << "  omega=" << omega
                 << "  dt=" << dt << "\n";

        } // end time-stepping loop
    } // end avalanche loop

    file_h.close(); file_u.close(); file_v.close();
    file_omega.close(); file_time.close(); file_mass.close(); file_ke.close();
    return 0;
}

// =============================================================================
// Gravity field loader
// =============================================================================

void load_gravity_field()
{
    ifstream eta_in(eta_file);
    ifstream grav_in(grav_file);
    string line;

    grav[0] = eta[0] = grav[RES + 3] = eta[RES + 3] = 0.0;

    if (!eta_in.is_open() || !grav_in.is_open())
    {
        cerr << "Error: cannot open gravity field files ("
             << grav_file << ", " << eta_file << ")\n";
        return;
    }

    int count = 2;
    while (getline(eta_in, line))
        eta[count++] = stof(line);

    count = 2;
    while (getline(grav_in, line))
        grav[count++] = gravity_scale * stof(line);

    // Fill ghost cells by copying nearest interior values
    grav[1] = grav[2];    grav[RES + 2] = grav[RES + 1];
    eta[1]  = eta[2];     eta[RES + 2]  = eta[RES + 1];

    // Precompute per-cell trig combinations (called O(steps*cells) times)
    for (int j = 0; j <= RES + 3; j++)
    {
        cos_xi_eta[j] = cos(apex_half_angle - eta[j]);
        sin_xi_eta[j] = sin(apex_half_angle - eta[j]);
    }

    eta_in.close();
    grav_in.close();
}

// =============================================================================
// Input file parser
// =============================================================================

void load_parameters(const string& filename)
{
    ifstream fin(filename);
    if (!fin.is_open())
    {
        cerr << "Warning: cannot open " << filename << "; using built-in defaults.\n";
        // Set defaults matching original code
        output_interval  = 50;
        t_final          = 10.0;
        limiter_theta    = 1.0;
        apex_half_angle  = 45.0 * PI / 180.0;
        domain_length    = 0.99;
        fillet_length    = 0.005;
        friction_angle   = 20.0 * PI / 180.0;
        omega            = 1.5;
        aspect_ratio     = 0.005;
        regolith_depth   = 1.0;
        density_ratio    = 1.2 / 0.75;
        gravity_scale    = 1.0;
        cone_slant_height = 470.0 / 376.0;
        n_avalanches     = 25;
        grav_file        = "grav45.txt";
        eta_file         = "eta45.txt";
        output_prefix    = "run1";
        return;
    }

    string line;
    while (getline(fin, line))
    {
        // Strip comments
        auto hash = line.find('#');
        if (hash != string::npos) line = line.substr(0, hash);

        // Split on '='
        auto eq = line.find('=');
        if (eq == string::npos) continue;

        string key = line.substr(0, eq);
        string val = line.substr(eq + 1);

        // Trim whitespace
        auto trim = [](string& s) {
            size_t start = s.find_first_not_of(" \t\r\n");
            size_t end   = s.find_last_not_of(" \t\r\n");
            s = (start == string::npos) ? "" : s.substr(start, end - start + 1);
        };
        trim(key); trim(val);
        if (key.empty() || val.empty()) continue;

        if      (key == "output_interval")    output_interval   = stoi(val);
        else if (key == "t_final")            t_final           = stod(val);
        else if (key == "limiter_theta")      limiter_theta     = stod(val);
        else if (key == "apex_half_angle_deg")apex_half_angle   = stod(val) * PI / 180.0;
        else if (key == "domain_length")      domain_length     = stod(val);
        else if (key == "fillet_length")      fillet_length     = stod(val);
        else if (key == "friction_angle_deg") friction_angle    = stod(val) * PI / 180.0;
        else if (key == "omega_init")         omega             = stod(val);
        else if (key == "aspect_ratio")       aspect_ratio      = stod(val);
        else if (key == "regolith_depth")     regolith_depth    = stod(val);
        else if (key == "density_ratio")      density_ratio     = stod(val);
        else if (key == "gravity_scale")      gravity_scale     = stod(val);
        else if (key == "cone_slant_height")  cone_slant_height = stod(val);
        else if (key == "n_avalanches")       n_avalanches      = stoi(val);
        else if (key == "grav_file")          grav_file         = val;
        else if (key == "eta_file")           eta_file          = val;
        else if (key == "output_prefix")      output_prefix     = val;
    }
    fin.close();
}

// =============================================================================
// Polynomial interpolants (piecewise-constant; hooks for future extension)
// =============================================================================

double interp_p(int jp) { return ptemp[jp]; }
double interp_q(int jp) { return qtemp[jp]; }
double interp_r(int jp) { return rtemp[jp]; }

// =============================================================================
// compute_stress_state
//
// Sets globals: normal_gravity (Gn), v_azim (local v), basal_stress (phi).
// side = 1 : plus (right) reconstruction at face jj
// side = 2 : minus (left) reconstruction at face jj
// side = 3 : cell-centre evaluation for cell jj
// =============================================================================

void compute_stress_state(int side, int jj)
{
    double xx  = domain_start + domain_length * (jj - 2) / RES;
    double tau = xx * SIN_ZETA;   // face coordinate (local)

    if (side == 1)
    {
        normal_gravity = omega * omega * COS_ZETA * tau
                       - grav[jj] * cos_xi_eta[jj];
        v_azim         = Qp_r / Qp_p;
        basal_stress   = normal_gravity
                       + 2.0 * omega * v_azim * COS_ZETA
                       + COS_ZETA * v_azim * v_azim / tau;
    }
    else if (side == 2)
    {
        normal_gravity = omega * omega * COS_ZETA * tau
                       - grav[jj - 1] * cos_xi_eta[jj - 1];
        v_azim         = Qm_r / Qm_p;
        basal_stress   = normal_gravity
                       + 2.0 * omega * v_azim * COS_ZETA
                       + COS_ZETA * v_azim * v_azim / tau;
    }
    else   // side == 3: cell centre
    {
        tau            = (xx + domain_length / 2.0 / RES) * SIN_ZETA;
        normal_gravity = omega * omega * COS_ZETA * tau
                       - grav[jj] * cos_xi_eta[jj];
        v_azim         = r[jj] / p[jj];
        basal_stress   = normal_gravity
                       + 2.0 * omega * v_azim * COS_ZETA
                       + COS_ZETA * v_azim * v_azim / tau;
    }
}

// =============================================================================
// reconstruct_face_values
//
// MUSCL reconstruction: sets Qp_* (right state) and Qm_* (left state)
// for the interface between cells jj-1 and jj.
// =============================================================================

void reconstruct_face_values(int jj)
{
    Qp_p = interp_p(jj)     - dx * limited_slope(jj,     1) / 2.0;
    Qm_p = interp_p(jj - 1) + dx * limited_slope(jj - 1, 1) / 2.0;
    Qp_q = interp_q(jj)     - dx * limited_slope(jj,     2) / 2.0;
    Qm_q = interp_q(jj - 1) + dx * limited_slope(jj - 1, 2) / 2.0;
    Qp_r = interp_r(jj)     - dx * limited_slope(jj,     3) / 2.0;
    Qm_r = interp_r(jj - 1) + dx * limited_slope(jj - 1, 3) / 2.0;
}

// =============================================================================
// max_face_speed  —  max wave speed across both reconstructed states at face jj
// =============================================================================

double max_face_speed(int jj)
{
    return max(max_eigenvalue(jj, 1), max_eigenvalue(jj, 2));
}

// =============================================================================
// max_eigenvalue  —  largest eigenvalue magnitude for side 1 (plus) or 2 (minus)
// =============================================================================

double max_eigenvalue(int jj, int side)
{
    double e1, e2, e3, root;
    double xx  = domain_start + domain_length * (jj - 1) / RES;
    reconstruct_face_values(jj + 1);
    compute_stress_state(side, jj + 1);

    if (side == 1)
    {
        double tau  = xx * SIN_ZETA;
        double tau2 = tau * tau,  tau3 = tau2 * tau, tau4 = tau2 * tau2;
        double Pp2  = Qp_p * Qp_p, Pp3 = Pp2 * Qp_p, Pp4 = Pp2 * Pp2, Pp5 = Pp4 * Qp_p;
        double Pr2  = Qp_r * Qp_r;
        root = aspect_ratio * earth_pressure_k * tau
             * (-normal_gravity * Pp5 * tau4
                - Pp3 * Pr2 * tau3 * COS_ZETA
                - 2.0 * Pp4 * Qp_r * tau4 * omega * COS_ZETA);
        root  = (root < 0.0) ? 0.0 : sqrt(root) / (Pp2 * tau3);
        e1    = abs(Qp_q / Qp_p);
        e2    = abs(e1 - root);
        e3    = abs(e1 + root);
    }
    else
    {
        double tau  = xx * SIN_ZETA;
        double tau2 = tau * tau,  tau3 = tau2 * tau, tau4 = tau2 * tau2;
        double Pm2  = Qm_p * Qm_p, Pm3 = Pm2 * Qm_p, Pm4 = Pm2 * Pm2, Pm5 = Pm4 * Qm_p;
        double Mr2  = Qm_r * Qm_r;
        root = aspect_ratio * earth_pressure_k * tau
             * (-normal_gravity * Pm5 * tau4
                - Pm3 * Mr2 * tau3 * COS_ZETA
                - 2.0 * Pm4 * Qm_r * tau4 * omega * COS_ZETA);
        root  = (root < 0.0) ? 0.0 : sqrt(root) / (Pm2 * tau3);
        e1    = abs(Qm_q / Qm_p);
        e2    = abs(e1 - root);
        e3    = abs(e1 + root);
    }
    return max(e1, max(e2, e3));
}

// =============================================================================
// flux  —  physical flux vector component for the given reconstructed state
// component: 1 = f1 (mass), 2 = f2 (downslope momentum), 3 = f3 (azimuthal momentum)
// =============================================================================

double flux(double pp, double qq, double rr, int jj, int component, int side)
{
    compute_stress_state(side, jj);
    f1 = qq;
    double h_sq = (pp / tau_loc) * (pp / tau_loc);
    f2 = qq * qq / pp
         - aspect_ratio * earth_pressure_k * h_sq * basal_stress * tau_loc / 2.0;
    f3 = qq * rr / pp;

    if (component == 1) return f1;
    if (component == 2) return f2;
    return f3;
}

// =============================================================================
// central_flux  —  Kurganov-Tadmor central numerical flux at face jj
// =============================================================================

double central_flux(int jj, int component)
{
    reconstruct_face_values(jj);
    fp = flux(Qp_p, Qp_q, Qp_r, jj, component, 1);
    fm = flux(Qm_p, Qm_q, Qm_r, jj, component, 2);

    double a = max_face_speed(jj - 1);   // local dissipation coefficient
    if (component == 1) return (fp + fm) / 2.0 - a * (Qp_p - Qm_p) / 2.0;
    if (component == 2) return (fp + fm) / 2.0 - a * (Qp_q - Qm_q) / 2.0;
    return                     (fp + fm) / 2.0 - a * (Qp_r - Qm_r) / 2.0;
}

// =============================================================================
// limited_slope  —  minmod-limited slope estimate for cell jp
// component: 1 = p, 2 = q, 3 = r
// =============================================================================

double limited_slope(int jp, int component)
{
    if (jp == RES + 3)
        cerr << "Error: limited_slope called with out-of-bounds index " << jp << "\n";

    if (component == 1)
        return minmod(limiter_theta * (ptemp[jp] - ptemp[jp - 1]) / dx,
                      (ptemp[jp + 1] - ptemp[jp - 1]) / (2.0 * dx),
                      limiter_theta * (ptemp[jp + 1] - ptemp[jp]) / dx);
    if (component == 2)
        return minmod(limiter_theta * (qtemp[jp] - qtemp[jp - 1]) / dx,
                      (qtemp[jp + 1] - qtemp[jp - 1]) / (2.0 * dx),
                      limiter_theta * (qtemp[jp + 1] - qtemp[jp]) / dx);
    return     minmod(limiter_theta * (rtemp[jp] - rtemp[jp - 1]) / dx,
                      (rtemp[jp + 1] - rtemp[jp - 1]) / (2.0 * dx),
                      limiter_theta * (rtemp[jp + 1] - rtemp[jp]) / dx);
}

// =============================================================================
// minmod limiter
// =============================================================================

double minmod(double a, double b, double c)
{
    if      (a < 0.0 && b < 0.0 && c < 0.0) return max(a, max(b, c));
    else if (a > 0.0 && b > 0.0 && c > 0.0) return min(a, min(b, c));
    else                                      return 0.0;
}

// =============================================================================
// Source terms  (uses globals x_loc, tau_loc, basal_stress set before calling)
// =============================================================================

// S1: source for p (mass) — zero by axisymmetry
double S1(int /*jp*/) { return 0.0; }

// S2: source for q (downslope momentum)
double S2(int jp)
{
    double ht = htemp[jp];
    double ut = utemp[jp];
    double vt = vtemp[jp];
    int    sgn_u = (ut > 0.0) ? 1 : ((ut < 0.0) ? -1 : 0);
    compute_stress_state(3, jp);
    tangential_gravity = omega * omega * SIN_ZETA * tau_loc
                       + grav[jp] * sin_xi_eta[jp];
    return SIN_ZETA * vt * vt * ht
         + ht * tau_loc * basal_stress * TAN_DELTA * sgn_u
         + tau_loc * ht * (tangential_gravity + 2.0 * vt * SIN_ZETA * omega);
}

// S3: source for r (azimuthal momentum)
// Note: basal_stress and tau_loc must already be set by the preceding S2 call.
double S3(int jp)
{
    double ht = htemp[jp];
    double ut = utemp[jp];
    double vt = vtemp[jp];
    int    sgn_v = (vt > 0.0) ? 1 : ((vt < 0.0) ? -1 : 0);
    return -SIN_ZETA * ut * vt * ht
         + ht * tau_loc * basal_stress * TAN_DELTA * sgn_v
         - tau_loc * ht * (tau_loc * ang_accel + 2.0 * ut * SIN_ZETA * omega);
}
