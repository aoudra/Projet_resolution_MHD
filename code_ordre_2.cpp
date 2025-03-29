#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm> // pour std::max

// -----------------------------------------------------------------------------
// 1) Indice linéaire et périodicité
// -----------------------------------------------------------------------------
inline int idx(int i, int j, int Nx) {
    return i + j * Nx;
}
inline int periodic(int i, int N) {
    return (i % N + N) % N;
}

// -----------------------------------------------------------------------------
// 2) minmod pour la reconstruction MUSCL
// -----------------------------------------------------------------------------
inline double minmod(double a, double b, double c) {
    auto sgn = [](double x) { return (x > 0) - (x < 0); };
    int sa = sgn(a), sb = sgn(b), sc = sgn(c);
    if (sa != 0 && sa == sb && sb == sc) {
        double m = std::min({ std::fabs(a), std::fabs(b), std::fabs(c) });
        return sa * m;
    }
    return 0.0;
}

// -----------------------------------------------------------------------------
// 3) Calcul de alpha, beta (vitesse max en x et y) et flux de Rusanov
// -----------------------------------------------------------------------------
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
    // c_s^2 = gamma * p / rho
    double a2 = gamma_ * p / rho;
    // b^2 = B1^2 + B2^2 + B3^2
    double b2 = B[0] * B[0] + B[1] * B[1] + B[2] * B[2];
    // discriminants
    double disc1 = (a2 + b2) * (a2 + b2) - 4.0 * a2 * (B[0] * B[0]);
    double disc2 = (a2 + b2) * (a2 + b2) - 4.0 * a2 * (B[1] * B[1]);
    if (disc1 < 0.0) disc1 = 0.0;
    if (disc2 < 0.0) disc2 = 0.0;
    double c1 = 0.5 * ((a2 + b2) + std::sqrt(disc1));
    double c2 = 0.5 * ((a2 + b2) + std::sqrt(disc2));
    alpha_out = std::fabs(u1) + std::sqrt(c1); // max vitesse en x
    beta_out = std::fabs(u2) + std::sqrt(c2); // max vitesse en y
}

// Flux de Rusanov en direction x
std::vector<double> rusanov_flux_x(
    const std::vector<double>& U_L,
    const std::vector<double>& U_R,
    double gamma_
) {
    // Gauche
    double rhoL = U_L[0];
    double rho_u1L = U_L[1];
    double rho_u2L = U_L[2];
    double rho_u3L = U_L[3];
    double B1L = U_L[4];
    double B2L = U_L[5];
    double B3L = U_L[6];
    double EL = U_L[7];

    double u1L = rho_u1L / rhoL;
    double u2L = rho_u2L / rhoL;
    double u3L = rho_u3L / rhoL;
    double pL = (gamma_ - 1.0) * (
        EL
        - 0.5 * rhoL * (u1L * u1L + u2L * u2L + u3L * u3L)
        - 0.5 * (B1L * B1L + B2L * B2L + B3L * B3L)
        );

    // Droite
    double rhoR = U_R[0];
    double rho_u1R = U_R[1];
    double rho_u2R = U_R[2];
    double rho_u3R = U_R[3];
    double B1R = U_R[4];
    double B2R = U_R[5];
    double B3R = U_R[6];
    double ER = U_R[7];

    double u1R = rho_u1R / rhoR;
    double u2R = rho_u2R / rhoR;
    double u3R = rho_u3R / rhoR;
    double pR = (gamma_ - 1.0) * (
        ER
        - 0.5 * rhoR * (u1R * u1R + u2R * u2R + u3R * u3R)
        - 0.5 * (B1R * B1R + B2R * B2R + B3R * B3R)
        );

    // flux en x
    std::vector<double> fL(8), fR(8);
    fL[0] = rho_u1L;
    fL[1] = rho_u1L * u1L + pL - 0.5 * B1L * B1L;
    fL[2] = rho_u1L * u2L - B1L * B2L;
    fL[3] = rho_u1L * u3L - B1L * B3L;
    fL[4] = 0.0;
    fL[5] = u1L * B2L - u2L * B1L;
    fL[6] = u1L * B3L - u3L * B1L;
    fL[7] = (EL + pL) * u1L - (u1L * B1L + u2L * B2L + u3L * B3L) * B1L;

    fR[0] = rho_u1R;
    fR[1] = rho_u1R * u1R + pR - 0.5 * B1R * B1R;
    fR[2] = rho_u1R * u2R - B1R * B2R;
    fR[3] = rho_u1R * u3R - B1R * B3R;
    fR[4] = 0.0;
    fR[5] = u1R * B2R - u2R * B1R;
    fR[6] = u1R * B3R - u3R * B1R;
    fR[7] = (ER + pR) * u1R - (u1R * B1R + u2R * B2R + u3R * B3R) * B1R;

    // alpha max
    double alphaL, betaL;
    double alphaR, betaR;
    {
        std::vector<double> BL = { B1L,B2L,B3L };
        std::vector<double> BR = { B1R,B2R,B3R };
        compute_alpha_beta(u1L, u2L, rhoL, pL, BL, gamma_, alphaL, betaL);
        compute_alpha_beta(u1R, u2R, rhoR, pR, BR, gamma_, alphaR, betaR);
    }
    double alpha = std::max(alphaL, alphaR);

    // différence d'état
    std::vector<double> dU(8);
    for (int c = 0; c < 8; c++) {
        dU[c] = U_R[c] - U_L[c];
    }

    // flux Rusanov
    std::vector<double> flux(8);
    for (int c = 0; c < 8; c++) {
        flux[c] = 0.5 * (fL[c] + fR[c]) - 0.5 * alpha * dU[c];
    }
    return flux;
}

// Flux de Rusanov en direction y (similaire)
std::vector<double> rusanov_flux_y(
    const std::vector<double>& U_D, // "down"
    const std::vector<double>& U_U, // "up"
    double gamma_
) {
    double rhoD = U_D[0];
    double rho_u1D = U_D[1];
    double rho_u2D = U_D[2];
    double rho_u3D = U_D[3];
    double B1D = U_D[4];
    double B2D = U_D[5];
    double B3D = U_D[6];
    double ED = U_D[7];

    double u1D = rho_u1D / rhoD;
    double u2D = rho_u2D / rhoD;
    double u3D = rho_u3D / rhoD;
    double pD = (gamma_ - 1.0) * (
        ED
        - 0.5 * rhoD * (u1D * u1D + u2D * u2D + u3D * u3D)
        - 0.5 * (B1D * B1D + B2D * B2D + B3D * B3D)
        );

    double rhoU = U_U[0];
    double rho_u1U = U_U[1];
    double rho_u2U = U_U[2];
    double rho_u3U = U_U[3];
    double B1U = U_U[4];
    double B2U = U_U[5];
    double B3U = U_U[6];
    double EU = U_U[7];

    double u1U = rho_u1U / rhoU;
    double u2U = rho_u2U / rhoU;
    double u3U = rho_u3U / rhoU;
    double pU = (gamma_ - 1.0) * (
        EU
        - 0.5 * rhoU * (u1U * u1U + u2U * u2U + u3U * u3U)
        - 0.5 * (B1U * B1U + B2U * B2U + B3U * B3U)
        );

    // flux en y
    std::vector<double> gD(8), gU(8);
    gD[0] = rho_u2D;
    gD[1] = rho_u2D * u1D - B2D * B1D;
    gD[2] = rho_u2D * u2D + pD - 0.5 * B2D * B2D;
    gD[3] = rho_u2D * u3D - B2D * B3D;
    gD[4] = u2D * B1D - u1D * B2D;
    gD[5] = 0.0;
    gD[6] = u2D * B3D - u3D * B2D;
    gD[7] = (ED + pD) * u2D - (u1D * B1D + u2D * B2D + u3D * B3D) * B2D;

    gU[0] = rho_u2U;
    gU[1] = rho_u2U * u1U - B2U * B1U;
    gU[2] = rho_u2U * u2U + pU - 0.5 * B2U * B2U;
    gU[3] = rho_u2U * u3U - B2U * B3U;
    gU[4] = u2U * B1U - u1U * B2U;
    gU[5] = 0.0;
    gU[6] = u2U * B3U - u3U * B2U;
    gU[7] = (EU + pU) * u2U - (u1U * B1U + u2U * B2U + u3U * B3U) * B2U;

    double alphaD, betaD;
    double alphaU, betaU;
    {
        std::vector<double> BD = { B1D,B2D,B3D }, BU = { B1U,B2U,B3U };
        compute_alpha_beta(u1D, u2D, rhoD, pD, BD, gamma_, alphaD, betaD);
        compute_alpha_beta(u1U, u2U, rhoU, pU, BU, gamma_, alphaU, betaU);
    }
    double bet = std::max(betaD, betaU);

    std::vector<double> dU(8);
    for (int c = 0; c < 8; c++) {
        dU[c] = U_U[c] - U_D[c];
    }

    std::vector<double> flux(8);
    for (int c = 0; c < 8; c++) {
        flux[c] = 0.5 * (gD[c] + gU[c]) - 0.5 * bet * dU[c];
    }
    return flux;
}

// -----------------------------------------------------------------------------
// 4) Reconstruction bilinéaire: calcul de slope_x, slope_y et "coin" SW, NE, ...
// -----------------------------------------------------------------------------

// Calcul d'une "slope" MUSCL en x et y, composante par composante.
void compute_slopes(
    const std::vector<double>& U,
    int Nx, int Ny,
    std::vector<double>& slope_x,
    std::vector<double>& slope_y
) {
    for (int j = 0; j < Ny; j++) {
        int jm = periodic(j - 1, Ny);
        int jp = periodic(j + 1, Ny);
        for (int i = 0; i < Nx; i++) {
            int im = periodic(i - 1, Nx);
            int ip = periodic(i + 1, Nx);

            for (int c = 0; c < 8; c++) {
                double Ui = U[8 * idx(i, j, Nx) + c];
                double Uim = U[8 * idx(im, j, Nx) + c];
                double Uip = U[8 * idx(ip, j, Nx) + c];
                double Ujm = U[8 * idx(i, jm, Nx) + c];
                double Ujp = U[8 * idx(i, jp, Nx) + c];

                // slope x
                double dR = (Uip - Ui);
                double dL = (Ui - Uim);
                double dM = 0.5 * (Uip - Uim);
                slope_x[8 * idx(i, j, Nx) + c] = minmod(dR, dM, dL);

                // slope y
                double dU = (Ujp - Ui);
                double dD = (Ui - Ujm);
                double dMy = 0.5 * (Ujp - Ujm);
                slope_y[8 * idx(i, j, Nx) + c] = minmod(dU, dMy, dD);
            }
        }
    }
}

// Bili. reconstruction "coin" : p_{i,j}( x_i + sign_x*dx/2, y_j + sign_y*dy/2 ).
inline void bilinear_corner(
    const std::vector<double>& U,
    const std::vector<double>& sx,
    const std::vector<double>& sy,
    int i, int j,
    int Nx, int Ny,
    double sign_x,
    double sign_y,
    std::vector<double>& out // taille 8
) {
    int k = 8 * idx(i, j, Nx);
    for (int c = 0; c < 8; c++) {
        double U0 = U[k + c];
        double sx0 = sx[k + c];
        double sy0 = sy[k + c];
        // p_{i,j}(...,...) = U0 + 0.5*sign_x*sx0 + 0.5*sign_y*sy0
        out[c] = U0 + 0.5 * sign_x * sx0 + 0.5 * sign_y * sy0;
    }
}


// -----------------------------------------------------------------------------
// 5) Programme principal
// -----------------------------------------------------------------------------
int main() {
    // -------------------------
    // 5.1) Paramètres
    // -------------------------
    int Nx = 400;
    int Ny = 400;
    double Lx = 4.0 * M_PI;
    double Ly = 4.0 * M_PI;
    double dx = Lx / Nx;
    double dy = Ly / Ny;

    double gamma_ = 5.0 / 3.0;

    // Coordonnées (optionnel, pour sortie)
    std::vector<double> xC(Nx), yC(Ny);
    for (int i = 0; i < Nx; i++) {
        xC[i] = (i + 0.5) * dx; // centre
    }
    for (int j = 0; j < Ny; j++) {
        yC[j] = (j + 0.5) * dy; // centre
    }

    // -------------------------
    // 5.2) Stockage U : 8*Nx*Ny
    // -------------------------
    std::vector<double> U(8 * Nx * Ny, 0.0);

    // Initialisation
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            double xx = xC[i];
            double yy = yC[j];

            double rho0 = std::pow(gamma_, 2.0);
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

    // Temps
    double max_wave_speed = 6.0;
    double CFL = 0.4;
    double dt = CFL * std::min(dx, dy) / max_wave_speed;
    double final_time = M_PI;
    int nt = (int)(final_time / dt);

    // -------------------------
    // 5.3) Boucle en temps
    // -------------------------
    for (int n = 0; n < nt; n++) {
        std::vector<double> U_new = U;

        // 1) calcul slopes
        std::vector<double> slope_x(8 * Nx * Ny, 0.0), slope_y(8 * Nx * Ny, 0.0);
        compute_slopes(U, Nx, Ny, slope_x, slope_y);

        // ---------------------------------------------------------------------
        // 2) Schéma en X
        // ---------------------------------------------------------------------
        for (int j = 0; j < Ny; j++) {
            int jm = periodic(j - 1, Ny);
            int jp = periodic(j + 1, Ny);
            for (int i = 0; i < Nx; i++) {
                int im = periodic(i - 1, Nx);
                int ip = periodic(i + 1, Nx);

                int k0 = 8 * idx(i, j, Nx);
                int kJpIp = 8 * idx(ip, jp, Nx);
                int kJmIm = 8 * idx(im, jm, Nx);
                int kJmIp = 8 * idx(ip, jm, Nx);
                int kJpIm = 8 * idx(im, jp, Nx);

                // On va définir Uji, Ujp_ip, etc. comme "coins" bilinéaires :
                // Uji      = coin SW de la maille (i,j)
                // Ujp_ip   = coin NE de la maille (i+1, j+1)
                // Ujm_im   = coin SW de la maille (i-1, j-1)
                // Ujm_ip   = coin SW de (i+1, j-1)
                // Ujp_im   = coin SW de (i-1, j+1)
                // etc.
                std::vector<double> Uji(8), Ujp_ip(8), Ujm_im(8), Ujm_ip(8), Ujp_im(8);

                // SW corner (i,j)
                bilinear_corner(U, slope_x, slope_y, i, j, Nx, Ny, -1.0, -1.0, Uji);
                // NE corner (i+1, j+1)
                {
                    int i2 = ip;
                    int j2 = jp;
                    bilinear_corner(U, slope_x, slope_y, i2, j2, Nx, Ny, +1.0, +1.0, Ujp_ip);
                }
                // SW corner (i-1, j-1)
                {
                    int i2 = im;
                    int j2 = jm;
                    bilinear_corner(U, slope_x, slope_y, i2, j2, Nx, Ny, -1.0, -1.0, Ujm_im);
                }
                // SW corner (i+1, j-1)
                {
                    int i2 = ip;
                    int j2 = jm;
                    bilinear_corner(U, slope_x, slope_y, i2, j2, Nx, Ny, -1.0, -1.0, Ujm_ip);
                }
                // SW corner (i-1, j+1)
                {
                    int i2 = im;
                    int j2 = jp;
                    bilinear_corner(U, slope_x, slope_y, i2, j2, Nx, Ny, -1.0, -1.0, Ujp_im);
                }

                // flux_x_plus   = rusanov_flux_x(Uji, Ujp_ip) - rusanov_flux_x(Ujm_im, Uji)
                auto flux_x_plus = rusanov_flux_x(Uji, Ujp_ip, gamma_);
                auto flux_x_plus2 = rusanov_flux_x(Ujm_im, Uji, gamma_);
                for (int c = 0; c < 8; c++) {
                    flux_x_plus[c] -= flux_x_plus2[c];
                }

                // flux_x_minus  = rusanov_flux_x(Uji, Ujm_ip) - rusanov_flux_x(Ujp_im, Uji)
                auto flux_x_minus = rusanov_flux_x(Uji, Ujm_ip, gamma_);
                auto flux_x_minus2 = rusanov_flux_x(Ujp_im, Uji, gamma_);
                for (int c = 0; c < 8; c++) {
                    flux_x_minus[c] -= flux_x_minus2[c];
                }

                // flux_standard_x = rusanov_flux_x(Uji, Uj_ip) - rusanov_flux_x(Uj_im, Uji)
                // => on reconstruit "Uj_ip" (maille (i, j), coin orienté ?),
                //    et "Uj_im", etc. 
                // Pour coller à ton code, on considère :
                // Uj_ip = coin SW de (i+1, j)
                // Uj_im = coin SW de (i-1, j)
                std::vector<double> Uj_ip(8), Uj_im(8);
                {
                    int i2 = ip;
                    int j2 = j;
                    bilinear_corner(U, slope_x, slope_y, i2, j2, Nx, Ny, -1.0, -1.0, Uj_ip);
                }
                {
                    int i2 = im;
                    int j2 = j;
                    bilinear_corner(U, slope_x, slope_y, i2, j2, Nx, Ny, -1.0, -1.0, Uj_im);
                }
                auto flux_standard_x = rusanov_flux_x(Uji, Uj_ip, gamma_);
                auto flux_standard_x2 = rusanov_flux_x(Uj_im, Uji, gamma_);
                for (int c = 0; c < 8; c++) {
                    flux_standard_x[c] -= flux_standard_x2[c];
                }

                // Mise à jour
                for (int c = 0; c < 8; c++) {
                    U_new[k0 + c] -= (dt / (4.0 * dx)) * (
                        flux_x_plus[c]
                        + 2.0 * flux_standard_x[c]
                        + flux_x_minus[c]);
                }
            }
        }

        // ---------------------------------------------------------------------
        // 3) Schéma en Y
        //    => flux_y_plus, flux_y_minus, flux_standard_y
        // ---------------------------------------------------------------------
        for (int j = 0; j < Ny; j++) {
            int jm = periodic(j - 1, Ny);
            int jp = periodic(j + 1, Ny);
            for (int i = 0; i < Nx; i++) {
                int im = periodic(i - 1, Nx);
                int ip = periodic(i + 1, Nx);

                int k0 = 8 * idx(i, j, Nx);

                // Reconstruire "coins" pour la direction Y
                // (même idée que ci-dessus, en t'assurant de bien coller aux
                //  indices U[:,j,i] => Uji, U[:,jp,ip] => Ujp_ip, etc.)
                std::vector<double> Uji(8), Ujp_ip(8), Ujm_im(8), Ujm_ip(8), Ujp_im(8);

                // On associe Uji = coin SW(i,j)
                bilinear_corner(U, slope_x, slope_y, i, j, Nx, Ny, -1.0, -1.0, Uji);
                // coin NE(i+1, j+1)
                {
                    int i2 = ip; int j2 = jp;
                    bilinear_corner(U, slope_x, slope_y, i2, j2, Nx, Ny, +1.0, +1.0, Ujp_ip);
                }
                // coin SW(i-1, j-1)
                {
                    int i2 = im; int j2 = jm;
                    bilinear_corner(U, slope_x, slope_y, i2, j2, Nx, Ny, -1.0, -1.0, Ujm_im);
                }
                // coin SW(i+1, j-1)
                {
                    int i2 = ip; int j2 = jm;
                    bilinear_corner(U, slope_x, slope_y, i2, j2, Nx, Ny, -1.0, -1.0, Ujm_ip);
                }
                // coin SW(i-1, j+1)
                {
                    int i2 = im; int j2 = jp;
                    bilinear_corner(U, slope_x, slope_y, i2, j2, Nx, Ny, -1.0, -1.0, Ujp_im);
                }

                // flux_y_plus  = rusanov_flux_y(Uji, Ujp_ip) - rusanov_flux_y(Ujm_im, Uji)
                auto flux_y_plus = rusanov_flux_y(Uji, Ujp_ip, gamma_);
                auto flux_y_plus2 = rusanov_flux_y(Ujm_im, Uji, gamma_);
                for (int c = 0; c < 8; c++) {
                    flux_y_plus[c] -= flux_y_plus2[c];
                }

                // flux_y_minus = rusanov_flux_y(Uji, Ujm_ip) - rusanov_flux_y(Ujp_im, Uji)
                auto flux_y_minus = rusanov_flux_y(Uji, Ujm_ip, gamma_);
                auto flux_y_minus2 = rusanov_flux_y(Ujp_im, Uji, gamma_);
                for (int c = 0; c < 8; c++) {
                    flux_y_minus[c] -= flux_y_minus2[c];
                }

                // flux_standard_y = rusanov_flux_y(Uji, ???) - rusanov_flux_y(???, Uji)
                // => idem X
                std::vector<double> Uj_jp(8), Uj_jm(8);
                {
                    int i2 = i; int j2 = jp;
                    bilinear_corner(U, slope_x, slope_y, i2, j2, Nx, Ny, -1.0, -1.0, Uj_jp);
                }
                {
                    int i2 = i; int j2 = jm;
                    bilinear_corner(U, slope_x, slope_y, i2, j2, Nx, Ny, -1.0, -1.0, Uj_jm);
                }
                auto flux_standard_y = rusanov_flux_y(Uji, Uj_jp, gamma_);
                auto flux_standard_y2 = rusanov_flux_y(Uj_jm, Uji, gamma_);
                for (int c = 0; c < 8; c++) {
                    flux_standard_y[c] -= flux_standard_y2[c];
                }

                // Mise à jour
                for (int c = 0; c < 8; c++) {
                    U_new[k0 + c] -= (dt / (4.0 * dy)) * (
                        flux_y_plus[c] +
                        2.0 * flux_standard_y[c]
                        + flux_y_minus[c]);
                }
            }
        }

        // validation
        U = U_new;
        std::cout << "Iteration " << n + 1 << " / " << nt << std::endl;
    }

    // -------------------------------------------------------------------------
    // 5.4) Sauvegarde : on calcule la pression finale et on écrit solution.dat
    // -------------------------------------------------------------------------
    std::ofstream ofs("solution.dat");
    if (!ofs) {
        std::cerr << "Impossible d'ouvrir solution.dat\n";
        return -1;
    }
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            int k = 8 * idx(i, j, Nx);
            double rho_ = U[k + 0];
            double rho_u1_ = U[k + 1];
            double rho_u2_ = U[k + 2];
            double rho_u3_ = U[k + 3];
            double B1_ = U[k + 4];
            double B2_ = U[k + 5];
            double B3_ = U[k + 6];
            double E_ = U[k + 7];

            double u1_ = rho_u1_ / rho_;
            double u2_ = rho_u2_ / rho_;
            double u3_ = rho_u3_ / rho_;
            double p_ = (gamma_ - 1.0) * (
                E_
                - 0.5 * rho_ * (u1_ * u1_ + u2_ * u2_ + u3_ * u3_)
                - 0.5 * (B1_ * B1_ + B2_ * B2_ + B3_ * B3_)
                );

            ofs << xC[i] << " " << yC[j] << " " << p_ << "\n";
        }
        ofs << "\n";
    }
    ofs.close();

    std::cout << "Simulation terminee. Donnees dans solution.dat\n";
    return 0;
}
