#include "Fluid.h"
#include <algorithm>

Fluid::Fluid(int N, float dt, float diffusion, float viscosity)
{
    this->size = N;
    this->dt = dt;
    this->diff = diffusion;
    this->visc = viscosity;

    this->pre_density = vector<float>(N * N, 0);
    this->density = vector<float>(N * N, 0);

    this->Vx = vector<float>(N * N, 0);
    this->Vy = vector<float>(N * N, 0);

    this->pre_Vx = vector<float>(N * N, 0);
    this->pre_Vy = vector<float>(N * N, 0);
}

void Fluid::addDensity(int x, int y, float amount)
{
    int index = IX(x, y);
    this->density[index] += amount;
}

void Fluid::addVelocity(int x, int y, float amountX, float amountY)
{
    int index = IX(x, y);
    this->Vx[index] += amountX;
    this->Vy[index] += amountY;
}

int Fluid::IX(int x, int y)
{   
    return x + y * this->size;
}

void Fluid::diffuse(int b, vector<float>& x, vector<float>& x0, float diff, int iter)
{
    int N = this->size;
    int dt = this->dt;
    float a = dt * diff * (N - 2) * (N - 2);
    lin_solve(b, x, x0, a, 1 + 6 * a, iter);
}

void Fluid::lin_solve(int b, vector<float>& x, vector<float>& x0, float a, float c, int iter)
{
    int N = this->size;
    float cRecip = 1.0 / c;
    for (int k = 0; k < iter; k++) 
    {
        for (int j = 1; j < N - 1; j++) 
        {
            for (int i = 1; i < N - 1; i++) 
            {
                x[IX(i, j)] = 
                    (x0[IX(i, j)]
                    + a * (x[IX(i + 1, j)]
                         + x[IX(i - 1, j)]
                         + x[IX(i, j + 1)]
                         + x[IX(i, j - 1)]
                         )) * cRecip;
            }
        }
        set_bnd(b, x);
    }
}

void Fluid::project(vector<float>& velocX, vector<float>& velocY, vector<float>& p, vector<float>& div, int iter)
{
    int N = this->size;
    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            div[IX(i, j)] = -0.5f * (
                velocX[IX(i + 1, j)]
                - velocX[IX(i - 1, j)]
                + velocY[IX(i, j + 1)]
                - velocY[IX(i, j - 1)]
                ) / N;
            p[IX(i, j)] = 0;
        }
    }
    set_bnd(0, div);
    set_bnd(0, p);
    lin_solve(0, p, div, 1, 6, iter);

    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            velocX[IX(i, j)] -= 0.5f * (p[IX(i + 1, j)]
                - p[IX(i - 1, j)]) * N;
            velocY[IX(i, j)] -= 0.5f * (p[IX(i, j + 1)]
                - p[IX(i, j - 1)]) * N;
        }
    }
    set_bnd(1, velocX);
    set_bnd(2, velocY);
}

void Fluid::advect(int b, vector<float>& d, vector<float>& d0, vector<float>& velocX, vector<float>& velocY)
{
    int N = this->size;
    int i, j, i0, j0, i1, j1;
    float x, y, s0, t0, s1, t1, dt0;
    dt0 = this->dt * N;
    for (i = 1; i < N - 1; i++) 
    {
        for (j = 1; j < N - 1; j++) 
        {
            x = i - dt0 * velocX[IX(i, j)];
            y = j - dt0 * velocY[IX(i, j)];

            if (x < 0.5f) x = 0.5f;
            if (x > N - 1.5f) x = N - 1.5f;
            i0 = (int)x; 
            i1 = i0 + 1;
            if (y < 0.5f) y = 0.5f; 
            if (y > N - 1.5) y = N - 1.5f;
            j0 = (int)y; 
            j1 = j0 + 1;

            s1 = x - i0; 
            s0 = 1 - s1; 
            t1 = y - j0; 
            t0 = 1 - t1;

            d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) +
                s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
        }
    }
    set_bnd(b, d);
}


void Fluid::set_bnd(int b, vector<float>& x)
{
    int N = this->size;

    for (int k = 1; k < N - 1; k++) {
        for (int i = 1; i < N - 1; i++) {
            x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
            x[IX(i, N - 1)] = b == 2 ? -x[IX(i, N - 2)] : x[IX(i, N - 2)];
        }
    }
    for (int k = 1; k < N - 1; k++) {
        for (int j = 1; j < N - 1; j++) {
            x[IX(0, j)] = b == 1 ? -x[IX(1, j)] : x[IX(1, j)];
            x[IX(N - 1, j)] = b == 1 ? -x[IX(N - 2, j)] : x[IX(N - 2, j)];
        }
    }

    x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(0, N - 1)] = 0.5f * (x[IX(1, N - 1)] + x[IX(0, N - 2)]);
    x[IX(N - 1, 0)] = 0.5f * (x[IX(N - 2, 0)] + x[IX(N - 1, 1)]);
    x[IX(N - 1, N - 1)] = 0.5f * (x[IX(N - 2, N - 1)] + x[IX(N - 1, N - 2)]);
}


void Fluid::fluidCubeStep()
{
    int N = this->size;
    float visc = this->visc;
    float diff = this->diff;
    float dt = this->dt;
    int iter = 4;
    vector<float>& Vx = this->Vx;
    vector<float>& Vy = this->Vy;
    vector<float>& Vx0 = this->pre_Vx;
    vector<float>& Vy0 = this->pre_Vy;
    vector<float>& s = this->pre_density;
    vector<float>& density = this->density;

    diffuse(1, Vx0, Vx, visc, iter);
    diffuse(2, Vy0, Vy, visc, iter);

    project(Vx0, Vy0, Vx, Vy, iter);

    advect(1, Vx, Vx0, Vx0, Vy0);
    advect(2, Vy, Vy0, Vx0, Vy0);


    project(Vx, Vy, Vx0, Vy0, iter);

    diffuse(0, s, density, diff, iter);
    advect(0, density, s, Vx, Vy);

    this->pre_density = this->density;
    this->pre_Vx = this->Vx;
    this->pre_Vy = this->Vy;
}

void Fluid::fadeDensity()
{
    for (int i = 0; i < this->density.size(); ++i)
    {
        float d = this->density[i];
        this->density[i] = d > 0.1 ? d - 0.1 : 0;
    }
}