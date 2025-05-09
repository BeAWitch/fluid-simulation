#pragma once

#include <vector>
using namespace std;

class Fluid
{
public:
	Fluid() = default;
    Fluid(int N, float dt, float diffusion, float viscosity);

    void addDensity(int x, int y, float amount);
    void addVelocity(int x, int y, float amountX, float amountY);
    void diffuse(int b, vector<float>& x, vector<float>& x0, float diff, int iter);
    void lin_solve(int b, vector<float>& x, vector<float>& x0, float a, float c, int iter); // 求解线性方程
    void project(vector<float>& velocX, vector<float>& velocY, vector<float>& p, vector<float>& div, int iter);
    void advect(int b, vector<float>& d, vector<float>& d0, vector<float>& velocX, vector<float>& velocY);
    void set_bnd(int b, vector<float>& x);

    void fluidCubeStep();

    void fadeDensity();

    int IX(int x, int y); // 计算索引
    vector<float> getDensity() const { return this->density; }

private:
    int size;
    float dt; // 时间步长
    float diff; // 扩散
    float visc; // 粘性

    vector<float> pre_density; // 先前的密度
    vector<float> density; // 密度

    vector<float> Vx; // 速度分量 x
    vector<float> Vy; // 速度分量 y

    vector<float> pre_Vx; // 先前的速度分量 x
    vector<float> pre_Vy; // 先前的速度分量 y
};
