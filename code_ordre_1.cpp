#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>

//  -----------------------------------------------------
//  Petit utilitaire pour indexer (i,j) dans un tableau 2D
//  en 1D: idx(i,j) = i + j*Nx
inline int idx(int i, int j, int Nx) {
    return i + j * Nx;
}

// ==================================================
// === Fonctions de calcul de flux de Rusanov (x,y)
// ==================================================

void compute_alpha_beta(
    double u1,
    double u2,
    double rho,
    double p,
    const std::vector<double>& B, // B = {B1, B2, B3}
    double gamma_,
    double& alpha_out,
    double& beta_out
) {
    // a^2 = gamma * p / rho
    double a2 = gamma_ * p / rho;
    // |B|^2
    double b2 = B[0] * B[0] + B[1] * B[1] + B[2] * B[2];

    // c1, c2
    double disc1 = (a2 + b2) * (a2 + b2) - 4.0 * a2 * (B[0] * B[0]);
    double disc2 = (a2 + b2) * (a2 + b2) - 4.0 * a2 * (B[1] * B[1]);
    if (disc1 < 0.0) disc1 = 0.0;
    if (disc2 < 0.0) disc2 = 0.0;

    double c1 = 0.5 * ((a2 + b2) + std::sqrt(disc1));
    double c2 = 0.5 * ((a2 + b2) + std::sqrt(disc2));

    alpha_out = std::fabs(u1) + std::sqrt(std::fabs(c1));
    beta_out = std::fabs(u2) + std::sqrt(std::fabs(c2));
}

// -----------------------------------------------------
// Flux de Rusanov en X
std::vector<double> rusanov_flux_x(
    const std::vector<double>& UL,
    const std::vector<double>& UR,
    double gamma_
) {
    //  U = (rho, rho_u1, rho_u2, rho_u3, B1, B2, B3, E)
    //   UL[0] = rho_L; etc.
    double rhoL = UL[0], rhoU1L = UL[1], rhoU2L = UL[2], rhoU3L = UL[3];
    double B1L = UL[4], B2L = UL[5], B3L = UL[6], EL = UL[7];

    double rhoR = UR[0], rhoU1R = UR[1], rhoU2R = UR[2], rhoU3R = UR[3];
    double B1R = UR[4], B2R = UR[5], B3R = UR[6], ER = UR[7];

    // vitesses
    double u1L = rhoU1L / rhoL, u2L = rhoU2L / rhoL, u3L = rhoU3L / rhoL;
    double u1R = rhoU1R / rhoR, u2R = rhoU2R / rhoR, u3R = rhoU3R / rhoR;

    // pressions
    double pL = (gamma_ - 1.0) * (EL - 0.5 * rhoL * (u1L * u1L + u2L * u2L + u3L * u3L)
        - 0.5 * (B1L * B1L + B2L * B2L + B3L * B3L));
    double pR = (gamma_ - 1.0) * (ER - 0.5 * rhoR * (u1R * u1R + u2R * u2R + u3R * u3R)
        - 0.5 * (B1R * B1R + B2R * B2R + B3R * B3R));

    // flux fL, fR
    std::vector<double> fL(8), fR(8);

    fL[0] = rhoU1L;
    fL[1] = rhoU1L * u1L + pL - 0.5 * B1L * B1L;
    fL[2] = rhoU1L * u2L - B1L * B2L;
    fL[3] = rhoU1L * u3L - B1L * B3L;
    fL[4] = 0.0;
    fL[5] = u1L * B2L - u2L * B1L;
    fL[6] = u1L * B3L - u3L * B1L;
    fL[7] = (EL + pL) * u1L - (u1L * B1L + u2L * B2L + u3L * B3L) * B1L;

    fR[0] = rhoU1R;
    fR[1] = rhoU1R * u1R + pR - 0.5 * B1R * B1R;
    fR[2] = rhoU1R * u2R - B1R * B2R;
    fR[3] = rhoU1R * u3R - B1R * B3R;
    fR[4] = 0.0;
    fR[5] = u1R * B2R - u2R * B1R;
    fR[6] = u1R * B3R - u3R * B1R;
    fR[7] = (ER + pR) * u1R - (u1R * B1R + u2R * B2R + u3R * B3R) * B1R;

    // alpha
    double alphaL, betaL, alphaR, betaR;
    {
        std::vector<double> BL = { B1L, B2L, B3L };
        std::vector<double> BR = { B1R, B2R, B3R };
        compute_alpha_beta(u1L, u2L, rhoL, pL, BL, gamma_, alphaL, betaL);
        compute_alpha_beta(u1R, u2R, rhoR, pR, BR, gamma_, alphaR, betaR);
    }
    double alpha = std::max(alphaL, alphaR);

    // difference (UR-UL)
    std::vector<double> dU(8);
    for (int k = 0; k < 8; k++) {
        dU[k] = UR[k] - UL[k];
    }

    // flux
    std::vector<double> flux(8);
    for (int k = 0; k < 8; k++) {
        flux[k] = 0.5 * (fL[k] + fR[k]) - 0.5 * alpha * dU[k];
    }
    return flux;
}

// -----------------------------------------------------
// Flux de Rusanov en Y
std::vector<double> rusanov_flux_y(
    const std::vector<double>& UL,
    const std::vector<double>& UR,
    double gamma_
) {
    double rhoL = UL[0], rhoU1L = UL[1], rhoU2L = UL[2], rhoU3L = UL[3];
    double B1L = UL[4], B2L = UL[5], B3L = UL[6], EL = UL[7];

    double rhoR = UR[0], rhoU1R = UR[1], rhoU2R = UR[2], rhoU3R = UR[3];
    double B1R = UR[4], B2R = UR[5], B3R = UR[6], ER = UR[7];

    double u1L = rhoU1L / rhoL, u2L = rhoU2L / rhoL, u3L = rhoU3L / rhoL;
    double u1R = rhoU1R / rhoR, u2R = rhoU2R / rhoR, u3R = rhoU3R / rhoR;

    double pL = (gamma_ - 1.0) * (EL - 0.5 * rhoL * (u1L * u1L + u2L * u2L + u3L * u3L)
        - 0.5 * (B1L * B1L + B2L * B2L + B3L * B3L));
    double pR = (gamma_ - 1.0) * (ER - 0.5 * rhoR * (u1R * u1R + u2R * u2R + u3R * u3R)
        - 0.5 * (B1R * B1R + B2R * B2R + B3R * B3R));

    std::vector<double> gL(8), gR(8);

    gL[0] = rhoU2L;
    gL[1] = rhoU2L * u1L - B1L * B2L;
    gL[2] = rhoU2L * u2L + pL - 0.5 * B2L * B2L;
    gL[3] = rhoU2L * u3L - B2L * B3L;
    gL[4] = u2L * B1L - u1L * B2L;
    gL[5] = 0.0;
    gL[6] = u2L * B3L - u3L * B2L;
    gL[7] = (EL + pL) * u2L - (u1L * B1L + u2L * B2L + u3L * B3L) * B2L;

    gR[0] = rhoU2R;
    gR[1] = rhoU2R * u1R - B1R * B2R;
    gR[2] = rhoU2R * u2R + pR - 0.5 * B2R * B2R;
    gR[3] = rhoU2R * u3R - B2R * B3R;
    gR[4] = u2R * B1R - u1R * B2R;
    gR[5] = 0.0;
    gR[6] = u2R * B3R - u3R * B2R;
    gR[7] = (ER + pR) * u2R - (u1R * B1R + u2R * B2R + u3R * B3R) * B2R;

    double alphaL, betaL, alphaR, betaR;
    {
        std::vector<double> BL = { B1L, B2L, B3L };
        std::vector<double> BR = { B1R, B2R, B3R };
        compute_alpha_beta(u1L, u2L, rhoL, pL, BL, gamma_, alphaL, betaL);
        compute_alpha_beta(u1R, u2R, rhoR, pR, BR, gamma_, alphaR, betaR);
    }
    double beta = std::max(betaL, betaR);

    std::vector<double> dU(8);
    for (int k = 0; k < 8; k++) {
        dU[k] = UR[k] - UL[k];
    }
    std::vector<double> flux(8);
    for (int k = 0; k < 8; k++) {
        flux[k] = 0.5 * (gL[k] + gR[k]) - 0.5 * beta * dU[k];
    }
    return flux;
}

// -----------------------------------------------------
//   Calcul de la vitesse d'onde max sur le domaine
//   => pour adapter dt = CFL * dx / alpha_max
// -----------------------------------------------------
double compute_max_wave_speed(
    const std::vector<double>& U,
    int Nx, int Ny,
    double gamma_
) {
    double alpha_max = 0.0;
    // On parcourt toutes les cellules
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            int k = 8 * idx(i, j, Nx);
            double rho_ = U[k + 0];
            double rhoU1_ = U[k + 1];
            double rhoU2_ = U[k + 2];
            double rhoU3_ = U[k + 3];
            double B1_ = U[k + 4];
            double B2_ = U[k + 5];
            double B3_ = U[k + 6];
            double E_ = U[k + 7];

            double u1_ = rhoU1_ / rho_;
            double u2_ = rhoU2_ / rho_;
            double u3_ = rhoU3_ / rho_;

            double p_ = (gamma_ - 1.0) * (E_
                - 0.5 * rho_ * (u1_ * u1_ + u2_ * u2_ + u3_ * u3_)
                - 0.5 * (B1_ * B1_ + B2_ * B2_ + B3_ * B3_));

            // On approxime la vitesse rapide (fast magnetosonic).
            // compute_alpha_beta donne alpha = |u1| + sqrt(c1)
            //            et beta  = |u2| + sqrt(c2)
            // On va juste prendre un max(|u| + ???)
            std::vector<double> B = { B1_, B2_, B3_ };
            double alpha_, beta_;
            compute_alpha_beta(u1_, u2_, rho_, p_, B, gamma_, alpha_, beta_);

            // On prend le max( alpha_, beta_, |u3| + c_approx ) 
            // mais pour simplifier on prend le plus grand entre alpha_ et beta_:
            double wave_speed_cell = std::max(alpha_, beta_);
            // On compare avec la composante u3 : on pourrait affiner si on voulait
            if (std::fabs(u3_) > wave_speed_cell)
                wave_speed_cell = std::fabs(u3_);

            alpha_max = std::max(alpha_max, wave_speed_cell);
        }
    }
    return alpha_max;
}


























// -----------------------------------------------------
//   Mise à jour "Forward Euler" sur un pas dt,
//   avec splitting en x puis y, en reprenant
//   votre logique de flux_plus, flux_minus, ...
// -----------------------------------------------------
std::vector<double> forward_euler_update(
    const std::vector<double>& U_in,
    double dt,
    int Nx, int Ny,
    double dx, double dy,
    double gamma_
) {
    // Copie initiale
    std::vector<double> U_new = U_in;

    // ========== Flux en X ==========
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            int ip = (i + 1) % Nx;
            int im = (i - 1 + Nx) % Nx;
            int jp = (j + 1) % Ny;
            int jm = (j - 1 + Ny) % Ny;

            int k0 = 8 * idx(i, j, Nx);
            int k_ip = 8 * idx(ip, j, Nx);
            int k_im = 8 * idx(im, j, Nx);
            int k_jp = 8 * idx(i, jp, Nx);
            int k_jm = 8 * idx(i, jm, Nx);
            int k_jp_ip = 8 * idx(ip, jp, Nx);
            int k_jm_im = 8 * idx(im, jm, Nx);
            int k_jm_ip = 8 * idx(ip, jm, Nx);
            int k_jp_im = 8 * idx(im, jp, Nx);

            // Construction de U[i,j], U[ip,j], etc.
            std::vector<double> Uij(8), Uipj(8), Uimj(8), Ujpj(8), Ujpj_ip(8), Ujmj_im(8), Ujmj_ip(8), Ujpj_im(8);
            for (int c = 0; c < 8; c++) {
                Uij[c] = U_in[k0 + c];
                Uipj[c] = U_in[k_ip + c];
                Uimj[c] = U_in[k_im + c];
                Ujpj[c] = U_in[k_jp + c];
                Ujpj_ip[c] = U_in[k_jp_ip + c];
                Ujmj_im[c] = U_in[k_jm_im + c];
                Ujmj_ip[c] = U_in[k_jm_ip + c];
                Ujpj_im[c] = U_in[k_jp_im + c];
            }

            // flux_plus  = rusanov_flux_x(Uij, Ujpj_ip) - rusanov_flux_x(Ujmj_im, Uij)
            auto flux_x_plus = rusanov_flux_x(Uij, Ujpj_ip, gamma_);
            auto flux_x_plus2 = rusanov_flux_x(Ujmj_im, Uij, gamma_);
            for (int c = 0; c < 8; c++) {
                flux_x_plus[c] -= flux_x_plus2[c];
            }

            // flux_minus = rusanov_flux_x(Uij, Ujmj_ip) - rusanov_flux_x(Ujpj_im, Uij)
            auto flux_x_minus = rusanov_flux_x(Uij, Ujmj_ip, gamma_);
            auto flux_x_minus2 = rusanov_flux_x(Ujpj_im, Uij, gamma_);
            for (int c = 0; c < 8; c++) {
                flux_x_minus[c] -= flux_x_minus2[c];
            }

            // flux_standard_x = rusanov_flux_x(Uij, Uipj) - rusanov_flux_x(Uimj, Uij)
            auto flux_standard_x = rusanov_flux_x(Uij, Uipj, gamma_);
            auto flux_standard_x2 = rusanov_flux_x(Uimj, Uij, gamma_);
            for (int c = 0; c < 8; c++) {
                flux_standard_x[c] -= flux_standard_x2[c];
            }

            for (int c = 0; c < 8; c++) {
                U_new[k0 + c] -= (dt / (4.0 * dx)) * (
                    flux_x_plus[c]
                    + 2.0 * flux_standard_x[c]
                    + flux_x_minus[c]);
            }
           /* for (int c = 0; c < 8; c++) {
                U_new[k0 + c] -= (dt / dx) * (flux_standard_x[c]);
            }*/
        }
    }

    // ========== Flux en Y ==========
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            int ip = (i + 1) % Nx;
            int im = (i - 1 + Nx) % Nx;
            int jp = (j + 1) % Ny;
            int jm = (j - 1 + Ny) % Ny;

            int k0 = 8 * idx(i, j, Nx);
            int k_ip = 8 * idx(ip, j, Nx);
            int k_im = 8 * idx(im, j, Nx);
            int k_jp = 8 * idx(i, jp, Nx);
            int k_jm = 8 * idx(i, jm, Nx);
            int k_jp_ip = 8 * idx(ip, jp, Nx);
            int k_jm_im = 8 * idx(im, jm, Nx);
            int k_jm_ip = 8 * idx(ip, jm, Nx);
            int k_jp_im = 8 * idx(im, jp, Nx);

            std::vector<double> Uij(8), Ujpip(8), Ujmim(8), Ujmip(8), Ujpim(8);
            for (int c = 0; c < 8; c++) {
                Uij[c] = U_new[k0 + c];  // On prend dans U_new (car on a déjà mis à jour en x)
            }
            // Idem pour accès diag
            std::vector<double> Ujpj_ip(8), Ujmj_im(8), Ujmj_ip(8), Ujpj_im(8);
            for (int c = 0; c < 8; c++) {
                Ujpj_ip[c] = U_new[k_jp_ip + c];
                Ujmj_im[c] = U_new[k_jm_im + c];
                Ujmj_ip[c] = U_new[k_jm_ip + c];
                Ujpj_im[c] = U_new[k_jp_im + c];
            }
            std::vector<double> Ujpj_i(8), Ujmj_i(8), Uij_ip(8), Uij_im(8);
            for (int c = 0; c < 8; c++) {
                Ujpj_i[c] = U_new[k_jp + c];
                Ujmj_i[c] = U_new[k_jm + c];
                Uij_ip[c] = U_new[k_ip + c];
                Uij_im[c] = U_new[k_im + c];
            }

            // flux_y_plus  = rusanov_flux_y(Uij, Ujpj_ip) - rusanov_flux_y(Ujmj_im, Uij)
            auto flux_y_plus = rusanov_flux_y(Uij, Ujpj_ip, gamma_);
            auto flux_y_plus2 = rusanov_flux_y(Ujmj_im, Uij, gamma_);
            for (int c = 0; c < 8; c++) {
                flux_y_plus[c] -= flux_y_plus2[c];
            }

            // flux_y_minus = rusanov_flux_y(Ujmj_ip, Uij) - rusanov_flux_y(Uij, Ujpj_im)
            auto flux_y_minus = rusanov_flux_y(Ujmj_ip, Uij, gamma_);
            auto flux_y_minus2 = rusanov_flux_y(Uij, Ujpj_im, gamma_);
            for (int c = 0; c < 8; c++) {
                flux_y_minus[c] -= flux_y_minus2[c];
            }

            // flux_standard_y = rusanov_flux_y(Uij, Ujpj_i) - rusanov_flux_y(Ujmj_i, Uij)
            auto flux_standard_y = rusanov_flux_y(Uij, Ujpj_i, gamma_);
            auto flux_standard_y2 = rusanov_flux_y(Ujmj_i, Uij, gamma_);
            for (int c = 0; c < 8; c++) {
                flux_standard_y[c] -= flux_standard_y2[c];
            }

            for (int c = 0; c < 8; c++) {
                U_new[k0 + c] -= (dt / (4.0 * dy)) * (flux_y_plus[c]
                    + 2.0 * flux_standard_y[c]
                    + flux_y_minus[c]);
            }
            /*for (int c = 0; c < 8; c++) {
                U_new[k0 + c] -= (dt / dy) * ( flux_standard_y[c]);
            }*/
        }
    }

    return U_new;
}








// =============================
// ==========   MAIN   =========
// =============================
int main() {
    // --- Paramètres de la grille
    int Nx = 400, Ny = 400;
    double Lx = 4.0 * M_PI, Ly = 4.0 * M_PI;
    double dx = Lx / Nx;
    double dy = Ly / Ny;

    // Coordonnées x,y
    std::vector<double> x(Nx), y(Ny);
    for (int i = 0; i < Nx; i++) x[i] = i * dx;
    for (int j = 0; j < Ny; j++) y[j] = j * dy;

    double gamma_ = 5.0 / 3.0;

    // On stocke 8 variables par maille: (rho, rho*u1, rho*u2, rho*u3, B1, B2, B3, E)
    std::vector<double> U(8 * Nx * Ny, 0.0);

    // ---------------------
    // Init
    // ---------------------
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            double xx = x[i];
            double yy = y[j];
            double rho0 = gamma_ * gamma_;   // gamma^2
            double u1_0 = -std::sin(yy);
            double u2_0 = std::sin(xx);
            double u3_0 = 0.0;
            double B1_0 = -std::sin(yy);
            double B2_0 = std::sin(2.0 * xx);
            double B3_0 = 0.0;
            double p0 = gamma_;

            double E0 = p0 / (gamma_ - 1.0)
                + 0.5 * rho0 * (u1_0 * u1_0 + u2_0 * u2_0 + u3_0 * u3_0)
                + 0.5 * (B1_0 * B1_0 + B2_0 * B2_0 + B3_0 * B3_0);

            int k = 8 * idx(i, j, Nx);
            U[k + 0] = rho0;
            U[k + 1] = rho0 * u1_0;
            U[k + 2] = rho0 * u2_0;
            U[k + 3] = rho0 * u3_0;
            U[k + 4] = B1_0;
            U[k + 5] = B2_0;
            U[k + 6] = B3_0;
            U[k + 7] = E0;
        }
    }

    // ---------------------
    // Paramètres temporels
    // ---------------------
    double CFL = 0.45;
    double t = 0.0;
    double final_time = M_PI;

    // On fait une boucle "while" tant que t < final_time
    // et on recalcule dt a chaque fois
    int nstep = 0;
    while (t < final_time) {
        // 1) Calculer la vitesse d'onde max
        double alpha_max = compute_max_wave_speed(U, Nx, Ny, gamma_);
        // 2) dt = CFL * min(dx,dy) / alpha_max
        double dt = CFL * std::min(dx, dy) / 10.0;
        // on s'assure de ne pas dépasser final_time
        if (t + dt > final_time) {
            dt = final_time - t;
        }

        // -- RK2 STAGE 1
        std::vector<double> U_star = forward_euler_update(U, dt, Nx, Ny, dx, dy, gamma_);

        // -- RK2 STAGE 2
        std::vector<double> U_star2 = forward_euler_update(U_star, dt, Nx, Ny, dx, dy, gamma_);

        // -- Combinaison
        //    U^{n+1} = 0.5*( U^n + U_star2 )
        // (forme simplifiée du TVD RK2)
        for (size_t c = 0; c < U.size(); c++) {
            U[c] = 0.5 * (U[c] + U_star2[c]);
        }

        // Incrément
        t += dt;
        nstep++;
        std::cout << "Step " << nstep << ", time=" << t << "/" << final_time
            << ", dt=" << dt << ", alpha_max=" << alpha_max << std::endl;
    }

    // ==================================
    // Ecriture du résultat final
    // ==================================
    std::ofstream ofs("solution.dat");
    if (!ofs) {
        std::cerr << "Impossible d'ouvrir solution.dat en écriture.\n";
        return 1;
    }

    // On calcule la pression finale p = (gamma-1)(E - 0.5*rho u^2 - 0.5 B^2)
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            int k = 8 * idx(i, j, Nx);
            double rho_ = U[k + 0];
            double rhoU1_ = U[k + 1];
            double rhoU2_ = U[k + 2];
            double rhoU3_ = U[k + 3];
            double B1_ = U[k + 4];
            double B2_ = U[k + 5];
            double B3_ = U[k + 6];
            double E_ = U[k + 7];

            double u1_ = rhoU1_ / rho_;
            double u2_ = rhoU2_ / rho_;
            double u3_ = rhoU3_ / rho_;
            double p_ = (gamma_ - 1.0) * (E_
                - 0.5 * rho_ * (u1_ * u1_ + u2_ * u2_ + u3_ * u3_)
                - 0.5 * (B1_ * B1_ + B2_ * B2_ + B3_ * B3_));

            ofs << x[i] << " " << y[j] << " " << p_ << "\n";
        }
        ofs << "\n";
    }
    ofs.close();
    std::cout << "Simulation terminée, resultat dans solution.dat\n";

    return 0;
}
