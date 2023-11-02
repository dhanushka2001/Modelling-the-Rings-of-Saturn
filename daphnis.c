#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// declare functions
double min(double x, double y);
void swapArray(double *a, double *b, int n);
void add_boundary_velocity(int N, int M, double *vx, double *vy, double *S_vx, double *S_vy);
void add_boundary_source(int N, int M, double *p, int S);
void add_source(int N, int M, double *v, double *a, double dt);
void add_non_inertial_source(int N, int M, double *vx, double *vy, double dx, double dt);
void diffuse(int N, int M, double *p, double *p0, double D, double dx, double dt);
void advect(int N, int M, double *p, double *p0, double *vx, double *vy, double dx, double dt);
void new_advect(int N, int M, double *p, double *p0, double *vx, double *vy, double dx, double dt);
void new_advect_x_y(int N, int M, double *vx, double *vy, double *vx0, double *vy0, double *p, double *T_p, double dx, double dt);
void vel_step(int N, int M, double *vx, double *vy, double *vx0, double *vy0, double *ax, double *ay, double *S_vx, double *S_vy, double *p, double *T_p, double visc, double dx, double dt);
void dens_step(int N, int M, double *p, double *p0, double *vx, double *vy, double D, double dx, double dt);
void write_out(int N, int M, double *p);
void write_out_u(int N, int M, double *p);
void write_out_v(int N, int M, double *p);

double min(double x, double y)
{
    if (x <= y)
    {
        return x;
    }
    else
    {
        return y;
    }
}

void swapArray(double *a, double *b, int n)
{
    for (int i = 0; i < n; i++)
    {
        double tmp = a[i];
        a[i] = b[i];
        b[i] = tmp;
    }
    // double *tmp;
    // tmp = a;
    // a = b;
    // b = tmp;
}

int IX(int i, int j, int N)
{
    return i + N * j;
}

void add_boundary_velocity(int N, int M, double *vx, double *vy, double *S_vx, double *S_vy)
{
    int i0 = N / 2;
    int j0 = M / 2;
    for (int j = 0; j < M; j++)
    {
        if (j < j0)
        {
            vx[0 + N * j] = S_vx[j];
            vy[0 + N * j] = S_vy[j];
            //vx[1 + N*j] = S_vx[j];
        }
        else
        {
            vx[N - 1 + N * j] = S_vx[j];
            vy[N - 1 + N * j] = S_vx[j];
            //vx[N-2 + N*j] = S_vx[j];
        }
    }
    // v=0 at Daphnis
    vx[IX(i0, j0, N)] = 0;
    vx[IX(i0 - 1, j0, N)] = 0;
    vx[IX(i0, j0 - 1, N)] = 0;
    vx[IX(i0 - 1, j0 - 1, N)] = 0;

    vy[IX(i0, j0, N)] = 0;
    vy[IX(i0 - 1, j0, N)] = 0;
    vy[IX(i0, j0 - 1, N)] = 0;
    vy[IX(i0 - 1, j0 - 1, N)] = 0;
}

void add_boundary_source(int N, int M, double *p, int S)
{
    int i0 = N / 2;
    int j0 = M / 2;
    for (int j = 0; j < M; j++)
    {
        // if ((j0 - 1 - 4 <= j) && (j <= j0 + 4))
        // {
        //     continue;
        // }
        // else
        // {
        //     if (j < j0 - 1 - 4)
        //     {
        //         p[0 + N * j] = (double)S;
        //         //p[1 + N*j] = (double) S;
        //         //p[2 + N*j] = (double) S;
        //         //p[3 + N*j] = (double) S;
        //     }
        //     else
        //     {
        //         //p[N-4 + N*j] = (double) S;
        //         //p[N-3 + N*j] = (double) S;
        //         //p[N-2 + N*j] = (double) S;
        //         p[N - 1 + N * j] = (double)S;
        //     }
        // }
        if (j <= j0 - 1)
        {
            p[0 + N * j] = (double)S;
        }
        else
        {
            p[N - 1 + N * j] = (double)S;
        }

        p[IX(i0, j0, N)] = 0;
        p[IX(i0 - 1, j0, N)] = 0;
        p[IX(i0, j0 - 1, N)] = 0;
        p[IX(i0 - 1, j0 - 1, N)] = 0;

        // if (j < j0)
        // {
        //     p[0 + N*j] = (double) S;
        //     //p[1 + N*j] = (double) S;
        //     //p[2 + N*j] = (double) S;
        //     //p[3 + N*j] = (double) S;
        // }
        // else
        // {
        //     //p[N-4 + N*j] = (double) S;
        //     //p[N-3 + N*j] = (double) S;
        //     //p[N-2 + N*j] = (double) S;
        //     p[N-1 + N*j] = (double) S;
        // }
    }
}

void add_source(int N, int M, double *v, double *a, double dt)
{
    // int i = (int) N/2;
    // int j = (int) M/2;
    // p[IX(i,j,N)] -= dt*S;
    for (int i = 0; i < N * M; i++)
    {
        v[i] += dt * a[i];
        //p[i] -= p[i]*S[i];
    }
}

void add_non_inertial_source(int N, int M, double *vx, double *vy, double dx, double dt)
{
    double G = 6.67E-20;           //e-20
    double saturn_mass = 5.683e26; //e26
    double daphnis_radius = 136505;
    double x, y;
    double N0 = (double)(N - 1) / 2;
    double M0 = (double)(M - 1) / 2;
    double w = sqrt(G * saturn_mass / (daphnis_radius * daphnis_radius * daphnis_radius));
    double c1 = 2 * w; // Coriolis
    double c2 = w * w; // Centrifugal
    for (int i = 0; i < N * M; i++)
    {
        x = dx * (i % N - N0);
        y = dx * (i / N - M0);
        vx[i] += dt * (-1 * c1 * vy[i] + c2 * x);
        vy[i] += dt * (c1 * vx[i] + c2 * y);
    }
}

void diffuse(int N, int M, double *p, double *p0, double D, double dx, double dt)
{
    double a = D * dt / (dx * dx);
    int n = 100;

    for (int k = 0; k < n; k++)
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < M; j++)
            {
                // boundary conditions
                if (i == 0)
                {
                    if (j == 0)
                    {
                        //p[IX(i,j,N)] = (p0[IX(i,j,N)] + a*(p[IX(N-1,j,N)]+p[IX(i+1,j,N)]+dx*p[IX(i,j+1,N)]))/(1+2*a+a*dx);
                        p[IX(i, j, N)] = (p0[IX(i, j, N)] + a * dx * (p[IX(i + 1, j, N)] + p[IX(i, j + 1, N)])) / (1 + 2 * a * dx);
                    }
                    else if (j == M - 1)
                    {
                        //p[IX(i,j,N)] = (p0[IX(i,j,N)] + a*(p[IX(N-1,j,N)]+p[IX(i+1,j,N)]+dx*p[IX(i,j-1,N)]))/(1+2*a+a*dx);
                        p[IX(i, j, N)] = (p0[IX(i, j, N)] + a * dx * (p[IX(i + 1, j, N)] + p[IX(i, j - 1, N)])) / (1 + 2 * a * dx);
                    }
                    else
                    {
                        //p[IX(i,j,N)] = (p0[IX(i,j,N)] + a*(p[IX(N-1,j,N)]+p[IX(i+1,j,N)]+p[IX(i,j-1,N)]+p[IX(i,j+1,N)]))/(1+4*a);
                        p[IX(i, j, N)] = (p0[IX(i, j, N)] + a * (dx * p[IX(i + 1, j, N)] + p[IX(i, j - 1, N)] + p[IX(i, j + 1, N)])) / (1 + 2 * a + a * dx);
                    }
                }
                else if (i == N - 1)
                {
                    if (j == 0)
                    {
                        //p[IX(i,j,N)] = (p0[IX(i,j,N)] + a*(p[IX(i-1,j,N)]+p[IX(0,j,N)]+dx*p[IX(i,j+1,N)]))/(1+2*a+a*dx);
                        p[IX(i, j, N)] = (p0[IX(i, j, N)] + a * dx * (p[IX(i - 1, j, N)] + dx * p[IX(i, j + 1, N)])) / (1 + 2 * a * dx);
                    }
                    else if (j == M - 1)
                    {
                        //p[IX(i,j,N)] = (p0[IX(i,j,N)] + a*(p[IX(i-1,j,N)]+p[IX(0,j,N)]+dx*p[IX(i,j-1,N)]))/(1+2*a+a*dx);
                        p[IX(i, j, N)] = (p0[IX(i, j, N)] + a * dx * (p[IX(i - 1, j, N)] + dx * p[IX(i, j - 1, N)])) / (1 + 2 * a * dx);
                    }
                    else
                    {
                        //p[IX(i,j,N)] = (p0[IX(i,j,N)] + a*(p[IX(i-1,j,N)]+p[IX(0,j,N)]+p[IX(i,j-1,N)]+p[IX(i,j+1,N)]))/(1+4*a);
                        p[IX(i, j, N)] = (p0[IX(i, j, N)] + a * (dx * p[IX(i - 1, j, N)] + p[IX(i, j - 1, N)] + p[IX(i, j + 1, N)])) / (1 + 2 * a + a * dx);
                    }
                }
                else
                {
                    if (j == 0)
                    {
                        p[IX(i, j, N)] = (p0[IX(i, j, N)] + a * (p[IX(i - 1, j, N)] + p[IX(i + 1, j, N)] + dx * p[IX(i, j + 1, N)])) / (1 + 2 * a + a * dx);
                    }
                    else if (j == M - 1)
                    {
                        p[IX(i, j, N)] = (p0[IX(i, j, N)] + a * (p[IX(i - 1, j, N)] + p[IX(i + 1, j, N)] + dx * p[IX(i, j - 1, N)])) / (1 + 2 * a + a * dx);
                    }
                    else
                    {
                        p[IX(i, j, N)] = (p0[IX(i, j, N)] + a * (p[IX(i - 1, j, N)] + p[IX(i + 1, j, N)] + p[IX(i, j - 1, N)] + p[IX(i, j + 1, N)])) / (1 + 4 * a);
                    }
                }
            }
        }
    }
}

void advect(int N, int M, double *p, double *p0, double *vx, double *vy, double dx, double dt)
{
    double x, y, s0, s1, t1, t0;
    int i0, j0, i1, j1;
    for (int j = 0; j < M; j++)
    {
        for (int i = 0; i < N; i++)
        {
            x = i - dt * vx[IX(i, j, N)] / dx; // needs to be old v, currently v is constant though
            y = j - dt * vy[IX(i, j, N)] / dx; // for small dx this should be okay
            //printf("%i %lf %lf\n", j, vy[IX(i,j,N)], y);

            // if x was initially < 0
            while (x < 0)
            {
                // after this x will be somewhere between [0,N]
                x += N;
            }
            // if x was initially > N
            while (x > N)
            {
                // after this x will be somewhere between [0,N]
                x -= N;
            }
            if (y < 0)
            {
                y = 0;
            }
            if (y > M - 1)
            {
                y = M - 2;
            }
            i0 = (int)x;
            i1 = (int)i0 + 1;
            if (i0 == N - 1)
            {
                i1 = 0;
            }
            j0 = (int)y;
            j1 = (int)j0 + 1;
            s1 = x - i0;
            s0 = 1 - s1; // = i1-x
            t1 = y - j0;
            t0 = 1 - t1; // = t1-x
            // if (i<5 && j==49)
            // {
            //     printf("%i %i %i %i\n",i0,i1,j0,j1);
            // }
            p[IX(i, j, N)] = s0 * (t0 * p0[IX(i0, j0, N)] + t1 * p0[IX(i0, j1, N)]) + s1 * (t0 * p0[IX(i1, j0, N)] + t1 * p0[IX(i1, j1, N)]);
        }
    }
}

void new_advect(int N, int M, double *p, double *p0, double *vx, double *vy, double dx, double dt)
{
    for (int j = 0; j < M; j++)
    {
        for (int i = 0; i < N; i++)
        {
            p[IX(i, j, N)] = 0;
        }
    }

    double x, y, s0, s1, t1, t0;
    int i0, j0, i1, j1;
    for (int j = 0; j < M; j++)
    {
        for (int i = 0; i < N; i++)
        {
            //printf("%i/%i\n",(i+1)+(j)*N,M*N);
            // x and y are the (double) gridpoints that (i,j) points to.
            x = i + dt * vx[IX(i, j, N)] / dx; // needs to be old v, currently v is constant though
            y = j + dt * vy[IX(i, j, N)] / dx; // for small dx this should be okay
            //printf("%i %lf %lf\n", j, vy[IX(i,j,N)], y);
            //printf("%lf %lf %lf %lf\n",x,y,dt,dx);

            // if x was initially < 0
            // while (x < 0)
            // {
            //     // after this x will be somewhere between [0,N]
            //     x += N;
            // }
            // // if x was initially > N
            // while (x > N)
            // {
            //     // after this x will be somewhere between [0,N]
            //     x -= N;
            // }
            if (x < 0 || x > N - 1)
            {
                //printf("%i %i %lf %lf\n",i, j, x, y);
                continue;
            }
            if (y < 0)
            {
                //y = 0;
                continue;
            }
            if (y > M - 1)
            {
                //y = M-2;
                continue;
            }
            i0 = (int)x;
            i1 = (int)i0 + 1;
            j0 = (int)y;
            j1 = (int)j0 + 1;
            if (i0 == N - 1)
            {
                i0 = N - 2;
                i1 = N - 1;
            }
            if (j0 == M - 1)
            {
                j0 = M - 2;
                j1 = M - 1;
            }
            s1 = x - i0;
            s0 = 1 - s1; // = i1-x
            t1 = y - j0;
            t0 = 1 - t1; // = t1-x

            // if (i == 0 || i == N-1)
            // {
            //     printf("DEN %i %i %lf %0.15lf\n",i, j, x-(double)i, y-(double)j);
            // }
            double k = p0[IX(i, j, N)];
            // if (j > 34)
            // {
            //     printf("DEN: i=%i j=%i dx=%lf dy=%0.15lf vx=%0.15lf vy=%0.15lf\n",i, j, x-(double)i, y-(double)j, vx[IX(i,j,N)], vy[IX(i,j,N)]);
            // }
            p[IX(i0, j0, N)] += s0 * t0 * k;
            p[IX(i0, j1, N)] += s0 * t1 * k;
            p[IX(i1, j0, N)] += s1 * t0 * k;
            p[IX(i1, j1, N)] += s1 * t1 * k;
        }
    }
}

void new_advect_x_y(int N, int M, double *vx, double *vy, double *vx0, double *vy0, double *p, double *T_p, double dx, double dt)
{
    for (int j = 0; j < M; j++)
    {
        for (int i = 0; i < N; i++)
        {
            vx[IX(i, j, N)] = 0;
            vy[IX(i, j, N)] = 0;
            T_p[IX(i, j, N)] = 0;
        }
    }

    double x, y, s0, s1, t1, t0;
    int i0, j0, i1, j1;
    for (int j = 0; j < M; j++)
    {
        for (int i = 0; i < N; i++)
        {
            // x and y are the (double) gridpoints that (i,j) points to.
            x = i + dt * vx0[IX(i, j, N)] / dx; // needs to be old v, currently v is constant though
            y = j + dt * vy0[IX(i, j, N)] / dx; // for small dx this should be okay

            // if x was initially < 0
            // while (x < 0)
            // {
            //     // after this x will be somewhere between [0,N]
            //     x += N;
            // }
            // // if x was initially > N
            // while (x > N)
            // {
            //     // after this x will be somewhere between [0,N]
            //     x -= N;
            // }
            if (x < 0 || x > N - 1)
            {
                //printf("%i %i %lf %lf\n",i, j, x, y);
                continue;
            }
            if (y < 0)
            {
                //y = 0;
                continue;
            }
            if (y > M - 1)
            {
                //y = M-2;
                continue;
            }
            i0 = (int)x;
            i1 = (int)i0 + 1;
            j0 = (int)y;
            j1 = (int)j0 + 1;
            if (i0 == N - 1)
            {
                i0 = N - 2;
                i1 = N - 1;
            }
            if (j0 == M - 1)
            {
                j0 = M - 2;
                j1 = M - 1;
            }
            s1 = x - i0;
            s0 = 1 - s1; // = i1-x
            t1 = y - j0;
            t0 = 1 - t1; // = t1-x

            //printf("VEL %i %i %lf %0.15lf\n",i, j, x-(double)i, y-(double)j);

            double k1 = vx0[IX(i, j, N)];
            double k2 = vy0[IX(i, j, N)];
            double p1 = p[IX(i, j, N)];

            /* total sums of m_i*u_i */
            vx[IX(i0, j0, N)] += s0 * t0 * k1 * p1;
            vx[IX(i0, j1, N)] += s0 * t1 * k1 * p1;
            vx[IX(i1, j0, N)] += s1 * t0 * k1 * p1;
            vx[IX(i1, j1, N)] += s1 * t1 * k1 * p1;

            vy[IX(i0, j0, N)] += s0 * t0 * k2 * p1;
            vy[IX(i0, j1, N)] += s0 * t1 * k2 * p1;
            vy[IX(i1, j0, N)] += s1 * t0 * k2 * p1;
            vy[IX(i1, j1, N)] += s1 * t1 * k2 * p1;

            /* total sums of "masses", m_i */
            T_p[IX(i0, j0, N)] += s0 * t0 * p1;
            T_p[IX(i0, j1, N)] += s0 * t1 * p1;
            T_p[IX(i1, j0, N)] += s1 * t0 * p1;
            T_p[IX(i1, j1, N)] += s1 * t1 * p1;
        }
    }

    for (int j = 0; j < M; j++)
    {
        for (int i = 0; i < N; i++)
        {
            double total_p = T_p[IX(i, j, N)];
            if (total_p < 1)
            {
                vx[IX(i, j, N)] = vx0[IX(i, j, N)];
                vy[IX(i, j, N)] = vy0[IX(i, j, N)];
            }
            else
            {
                vx[IX(i, j, N)] /= total_p;
                vy[IX(i, j, N)] /= total_p;
            }
        }
    }
}

void vel_step(int N, int M, double *vx, double *vy, double *vx0, double *vy0, double *ax, double *ay, double *S_vx, double *S_vy, double *p, double *T_p, double visc, double dx, double dt)
{
    //write_out(N,M,vy);
    add_source(N, M, vx, ax, dt);
    add_source(N, M, vy, ay, dt);
    //add_non_inertial_source(N, M, vx, vy, dx, dt);
    swapArray(vx0, vx, N * M);
    swapArray(vy0, vy, N * M);
    diffuse(N, M, vx, vx0, visc, dx, dt);
    diffuse(N, M, vy, vy0, visc, dx, dt);
    swapArray(vx0, vx, N * M);
    swapArray(vy0, vy, N * M);
    new_advect_x_y(N, M, vx, vy, vx0, vy0, p, T_p, dx, dt);
    add_boundary_velocity(N, M, vx, vy, S_vx, S_vy);
    //printf("vel %0.15lf %0.15lf %0.15lf %0.15lf\n",vx[IX(49,24,N)],vy[IX(49,24,N)],vx0[IX(49,24,N)],vy0[IX(49,24,N)]);
}

void dens_step(int N, int M, double *p, double *p0, double *vx, double *vy, double D, double dx, double dt)
{
    // p ALREADY has boundary source initially.
    //printf("den %0.15lf %0.15lf\n",vx[IX(49,24,N)],vy[IX(49,24,N)]);
    swapArray(p0, p, N * M);                 // swap p_old and p_new
    diffuse(N, M, p, p0, D, dx, dt);         // diffuse p_old iteratively and store in p_new
    swapArray(p0, p, N * M);                 // swap p_old and p_new
    new_advect(N, M, p, p0, vx, vy, dx, dt); // advect p_old and store in p_new
    add_boundary_source(N, M, p, 100);       // add source to p_new
}

void write_out(int N, int M, double *p)
{
    FILE *output_file = fopen("E:/Warwick/saturn/daphnis_e10.txt", "a");
    for (int j = 0; j < M; j++)
    {
        for (int i = 0; i < N; i++)
        {
            fprintf(output_file, "%0.10lf ", p[IX(i, j, N)]);
        }
        fprintf(output_file, "\n");
    }
    fclose(output_file);
}

void write_out_u(int N, int M, double *p)
{
    FILE *output_file = fopen("C:/C Projects/daphnis_u.txt", "a");
    for (int j = 0; j < M; j++)
    {
        for (int i = 0; i < N; i++)
        {
            fprintf(output_file, "%0.10lf ", p[IX(i, j, N)]);
        }
        fprintf(output_file, "\n");
    }
    fclose(output_file);
}

void write_out_v(int N, int M, double *p)
{
    FILE *output_file = fopen("C:/C Projects/daphnis_v.txt", "a");
    for (int j = 0; j < M; j++)
    {
        for (int i = 0; i < N; i++)
        {
            fprintf(output_file, "%0.10lf ", p[IX(i, j, N)]);
        }
        fprintf(output_file, "\n");
    }
    fclose(output_file);
}

int main()
{
    // we are using KM and KG
    double G = 6.67E-20;           //e-20
    double saturn_mass = 5.683e26; //e26
    double daphnis_mass = 7.7e10;  //e13
    double daphnis_radius = 136505;
    double C = sqrt(G * saturn_mass);

    int N = 800;      // N and M are # of gridpoints, L and W (calculated) are real lengths.
    int M = 100;        // width, Keepler Gap (Daphnis) is 35km. N=800,M=100,L=3200
    double L = 3200;   // width = M * L/N
    double dx = L / N; // 4km is ideal
    double dt = 100;
    double D = 1e-12;
    double visc = 1e-12;
    int i0 = N / 2;
    int j0 = M / 2;
    double N0 = (double)(N - 1) / 2;
    double M0 = (double)(M - 1) / 2;

    double *p0 = (double *)malloc(sizeof(double) * N * M);
    double *p = (double *)malloc(sizeof(double) * N * M);
    double *S = (double *)malloc(sizeof(double) * N * M);
    double *S_vx = (double *)malloc(sizeof(double) * M);
    double *S_vy = (double *)malloc(sizeof(double) * M);
    double *vx = (double *)malloc(sizeof(double) * N * M);
    double *vy = (double *)malloc(sizeof(double) * N * M);
    double *vx0 = (double *)malloc(sizeof(double) * N * M);
    double *vy0 = (double *)malloc(sizeof(double) * N * M);
    double *ax = (double *)malloc(sizeof(double) * N * M);
    double *ay = (double *)malloc(sizeof(double) * N * M);
    double *T_p = (double *)malloc(sizeof(double) * N * M);

    double v0 = sqrt(G * saturn_mass / daphnis_radius);

    // INITIAL CONDITIONS
    for (int j = 0; j < M; j++)
    {
        for (int i = 0; i < N; i++)
        {
            int index = IX(i, j, N);
            double i1 = i - N0;
            double j1 = j - M0;
            double x = i1 * dx; // x-distance from Daphnis (left is neg, right is pos)
            double y = j1 * dx; // y-distance from Daphnis (up is neg, down is pos)
            double y1 = daphnis_radius + j1 * dx; // y-distance from Saturn (always positive, force is pointing up so need to convert to negative)
            double radius = sqrt(x * x + y1 * y1); // distance to Saturn
            double v = C / sqrt(radius); // velocity
            double vx1 = v * (y1 / radius) - v0; // x-velocity (in frame of Daphnis, subtracting v0)
            //double vy1 = v * (-1*x / radius); // y-velocity
            // R3 is |r|^3
            double R3 = sqrt(x * x + y * y) * sqrt(x * x + y * y) * sqrt(x * x + y * y);
            //printf("%i %lf\n",j,v);
            if (i == 0)
            {
                S_vx[j] = vx1;
                S_vy[j] = 0;//vy1;
            }
            vx[index] = vx1;
            vy[index] = 0;//vy1;
            vx0[index] = 0;
            vy0[index] = 0;
            p0[index] = 0;
            p[index] = 100;
            // if ((j0 - 1 - 4 <= j) && (j <= j0 + 4))
            // {
            //     p[index] = 0;
            // }
            ax[index] = -1 * G * daphnis_mass * (x) / (R3);                                         // F_g(Daphnis) in x
            ay[index] = -1 * G * daphnis_mass * (y) / (R3);                                         // F_g(Daphnis) in y
            //ax[index] += -1 * G * saturn_mass * (x) / (radius*radius*radius);                       // F_g(Saturn) in x
            //ay[index] += -1 * G * saturn_mass * (y1) / (radius*radius*radius);                      // F_g(Saturn) in y
            //ay[index] += 1 * G * saturn_mass / (daphnis_radius*daphnis_radius);                     // F_fictitious (orbitting)
            //double w = sqrt(G * saturn_mass / (daphnis_radius * daphnis_radius * daphnis_radius));  // = v/r
            //double c1 = 2 * w; // Coriolis
            //double c2 = w * w; // Centrifugal
            //ax[index] += (-1 * c1 * vy[index] + c2 * x);                                                // Coriolis + Centrifugal in x
            //ay[index] += (c1 * vx[index] + c2 * y);                                                     // Coriolis + Centrifugal in y

            

            //printf("%i %lf %lf %0.15lf %0.15lf\n", index, vx[index], vy[index], ax[index], ay[index]);
            // if ((N % 2 == 0) && (M % 2 == 0))
            // {
            //     if ((i == i0 && j == j0) || (i == i0 - 1 && j == j0) || (i == i0 && j == j0 - 1) || (i == i0 - 1 && j == j0 - 1))
            //     {
            //         ax[index] = 0;
            //         ay[index] = 0;
            //     }
            // }
            // else if ((N % 2 == 0) && (M % 2 != 0))
            // {
            //     if ((i == i0 && j == j0) || (i == i0 - 1 && j == j0))
            //     {
            //         ax[index] = 0;
            //         ay[index] = 0;
            //     }
            // }
            // else if ((N % 2 != 0) && (M % 2 == 0))
            // {
            //     if ((i == i0 && j == j0) || (i == i0 && j == j0 - 1))
            //     {
            //         ax[index] = 0;
            //         ay[index] = 0;
            //     }
            // }
            // else
            // {
            //     if ((i == i0 && j == j0))
            //     {
            //         ax[index] = 0;
            //         ay[index] = 0;
            //     }
            // }

        // a=0 at Daphnis
        ax[IX(i0, j0, N)] = 0;
        ax[IX(i0 - 1, j0, N)] = 0;
        ax[IX(i0, j0 - 1, N)] = 0;
        ax[IX(i0 - 1, j0 - 1, N)] = 0;

        ay[IX(i0, j0, N)] = 0;
        ay[IX(i0 - 1, j0, N)] = 0;
        ay[IX(i0, j0 - 1, N)] = 0;
        ay[IX(i0 - 1, j0 - 1, N)] = 0;
        }
    }

    // NxM total gridpoints, N columns, M rows. Meaning Mx(# of cycles) = # of lines
    add_boundary_source(N, M, p, 100);
    write_out(N, M, p);
    //write_out_u(N, M, vx);
    //write_out_v(N, M, vy);
    int n = 100000;
    for (int k = 0; k < n; k++)
    {
        vel_step(N, M, vx, vy, vx0, vy0, ax, ay, S_vx, S_vy, p, T_p, visc, dx, dt);
        dens_step(N, M, p, p0, vx, vy, D, dx, dt);
        if ((k + 1) % 100 == 0)
        {
            write_out(N, M, p);
            //write_out_u(N, M, vx);
            //write_out_v(N, M, vy);
        }
        printf("%i/%i\r", k + 1, n);
    }

    free(p0);
    free(p);
    free(S);
    free(vx);
    free(vy);
    free(vx0);
    free(vy0);
    free(ax);
    free(ay);
    free(T_p);
}
