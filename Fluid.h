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
    void lin_solve(int b, vector<float>& x, vector<float>& x0, float a, float c, int iter); // ������Է���
    void project(vector<float>& velocX, vector<float>& velocY, vector<float>& p, vector<float>& div, int iter);
    void advect(int b, vector<float>& d, vector<float>& d0, vector<float>& velocX, vector<float>& velocY);
    void set_bnd(int b, vector<float>& x);

    void fluidCubeStep();

    void fadeDensity();

    int IX(int x, int y); // ��������
    vector<float> getDensity() const { return this->density; }

private:
    int size;
    float dt; // ʱ�䲽��
    float diff; // ��ɢ
    float visc; // ճ��

    vector<float> pre_density; // ��ǰ���ܶ�
    vector<float> density; // �ܶ�

    vector<float> Vx; // �ٶȷ��� x
    vector<float> Vy; // �ٶȷ��� y

    vector<float> pre_Vx; // ��ǰ���ٶȷ��� x
    vector<float> pre_Vy; // ��ǰ���ٶȷ��� y
};
