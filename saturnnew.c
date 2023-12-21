#include <stdio.h>
#include <stdlib.h>
#include <math.h>                           // https://www.tutorialspoint.com/c_standard_library/math_h.htm
#include <time.h>

/* Guide on C commenting etiquette */
// https://users.cs.utah.edu/~germain/PPS/Topics/commenting.html
// https://www.doc.ic.ac.uk/lab/cplus/cstyle.html
// https://improvingsoftware.com/2011/06/27/5-best-practices-for-commenting-your-code/

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif
#define G 6.67e-20                          // big G with km instead of m
#define saturn_equatorial_radius 60268      // km
#define saturn_polar_radius 60268           // km 54364
#define saturn_mass 5.683e26                // kg

#define mimas_radius 198                    // km
#define mimas_mass 3.75e19                  // kg (should be e19)
#define mimas_semi_major_axis 185539        // km
#define mimas_periapsis 181902              // km
#define mimas_apoapsis 189176               // km
#define mimas_ecc 0.0196

#define titan_radius 2575                   // km
#define titan_mass 1.35e23                  // kg
// titan_period = 16 * 24 * 60 * 60         // seconds
#define dist_titan_saturn 1221870           // km

#define huygens_gap_radius 117680           // km
#define huygens_gap_width 350               // km

#define epimetheus_mass 5.266e17            // km
#define janus_mass 1.8975e18                // km

/* Declare functions */
double randfrom(double min, double max);
double *add(double *output, double *arr1, double *arr2);
double *subtract(double *output, double *arr1, double *arr2);
double *scalar_mult(double *output, double *arr, double scalar);
double norm(double *r);
void append(double *twoD_arr, double *arr, int n_orbits);

void moon_positions_to_text_file(double *positions, int twoD_arr_length);
void moon1_positions_to_text_file(double *positions, int twoD_arr_length);
void moon2_positions_to_text_file(double *positions, int twoD_arr_length);
void particle_text_file_append(double *arr, int twoD_arr_length);
void arr1d_text_file_append(double *arr, int arr_length);
void two_arr1d_text_file_append(double *arr, int arr_length);

void verlet_vectorize(double *mimas_positions, double min_radius, double max_radius, int n_orbits, int n_particles, double timestep);
void verlet_vectorize_subcycles(double min_radius, double max_radius, int n_orbits, int n_particles, double timestep, int n_subcycles);
void verlet_vectorize_subcycles_new(double min_radius, double max_radius, int n_orbits, int n_particles, double timestep, int n_subcycles);

void new_mimas_pos(double *new_r, int n_subcycles, int t, int k);
void new_mimas_pos_kepler(double *new_r, double mimas_period, double timestep, int k);
void new_mimas_pos_kepler_parametric(double *new_r, double mimas_period, double timestep, int k);

void add_2darray(double *output, double *arr1, double *arr2, int arr_length);
void subtract_2darray(double *output, double *arr1, double *arr2, int arr_length);
void subtract_2darray_1darray(double *output, double *arr2d, double *arr1d, int arr_length);
void scalar_mult_2darray(double *output, double *arr, double scalar, int arr_length);
void vector_mult_2darray(double *output, double *arr, double *scalars, int arr_length);
void dot_prod_2darray(double *output, double *arr1, double *arr2, int arr_length);
void array_to_2darray(double *output, double *arr1d, int arr_length);
void norm_2darray(double *output, double *arr, int arr_length);
void pow_neg1_array(double *output, double *arr, int arr_length);
void pow_neg3_array(double *output, double *arr, int arr_length);
void acceleration_particle_2darray(double *output, double *r_mimas, double *r_particles, double *arr2d, double *arr1d, int arr_length);
void append_packet_2darray(double *packet, double *arr2d, int step, int n_particles);
void append_packet_1darray(double *packet, double *arr1d, int step, int n_particles);
void append_packet_2_1darrays(double *packet, double *arr1d_1, double *arr1d_2, int step, int n_particles);
void scalar_mult_1darray(double *output, double *arr, double scalar, int arr_length);
void replace(double *arr1, double *arr2, int arr_length);
void new_janus_epimetheus_pos(double *new_r_J, double *new_r_E, double *old_r_J, double *old_r_E, double *old_v_J, double *old_v_E, double *a1, double *a2, double *b1, double *b2, double *b3, double *b4, double timestep, int k);

void acceleration_saturn_2darray(double *output, double *r_particles, double *arr1d, int arr_length);
void acceleration_mimas_2darray(double *output, double *r_mimas, double *r_particles, double *arr1d, int arr_length);
void acceleration_janus_2darray(double *output, double *r_janus, double *r_particles, double *arr1d, int arr_length);
void acceleration_epimetheus_2darray(double *output, double *r_epimetheus, double *r_particles, double *arr1d, int arr_length);
void acceleration_janus_from_epimetheus_1darray(double *output, double *r_janus, double *r_epimetheus);
void acceleration_epimetheus_from_janus_1darray(double *output, double *r_epimetheus, double *r_janus);
void acceleration_saturn_1darray(double *output, double *r_pos);
void zero_2d_array(double *output, int arr_length);

int main()
{
    srand(time(NULL));
    double timestep = 100;
    int n_subcycles = 10000;
    int n_orbits = 10000;
    int n_particles = 10000;
    double min_radius = 115500; // 95400 - 1000;
    double max_radius = 118000; // 95400 + 1000;

    clock_t begin = clock();
    //double a = fmin((double) 11/3, M_PI);
    //printf("a=%0.15f\n", a);
    
    verlet_vectorize_subcycles_new(min_radius, max_radius, n_orbits, n_particles, timestep, n_subcycles);

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Elapsed time: %0.15f seconds\n", time_spent);

    return 0;
}

double randfrom(double min, double max)
{
    //srand(time(NULL));
    double range = (max - min);
    double div = RAND_MAX / range;
    //srand(time(NULL));
    return min + (rand() / div);
}

double *add(double *output, double *arr1, double *arr2)
{
    // Read about what "#pragma GCC ivdep" does here: https://codeforces.com/blog/entry/96344
    #pragma GCC ivdep
    for (int i = 0; i < 3; i++)
    {
        output[i] = arr1[i] + arr2[i];
    }
    return output;
}

double *subtract(double *output, double *arr1, double *arr2)
{
    for (int i = 0; i < 3; i++)
    {
        output[i] = arr1[i] - arr2[i];
    }
    return output;
}

double *scalar_mult(double *output, double *arr, double scalar)
{
    for (int i = 0; i < 3; i++)
    {
        output[i] = arr[i] * scalar;
    }
    return output;
}

double norm(double *r)
{
    double mag_r = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    return mag_r;
}

void append(double *twoD_arr, double *arr, int n_orbits)
{
    #pragma GCC ivdep
    for (int i = 0; i < 3; i++)
    {
        *(twoD_arr + 3 * n_orbits + i) = *(arr + i);
    }
}

void moon_positions_to_text_file(double *positions, int twoD_arr_length)
{
    FILE *f = fopen("E:/Warwick/saturn/mimas_x_y_dt=100_sub=10000.txt", "a");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    for (int i = 0; i < twoD_arr_length; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            fprintf(f, "%0.15f ", *(2 * i + positions + j));
        }
        fprintf(f, "\n");
    }
    fclose(f);
    //printf("DONE TEXT FILE!\n");
}

void moon1_positions_to_text_file(double *positions, int twoD_arr_length)
{
    FILE *f = fopen("E:/Warwick/saturn/janus_x_y_z_dt=100_sub=10000_2.txt", "a");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    for (int i = 0; i < twoD_arr_length; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            fprintf(f, "%0.15f ", *(3 * i + positions + j));
        }
        fprintf(f, "\n");
    }
    fclose(f);
    //printf("DONE TEXT FILE!\n");
}

void moon2_positions_to_text_file(double *positions, int twoD_arr_length)
{
    FILE *f = fopen("E:/Warwick/saturn/epimetheus_x_y_z_dt=100_sub=10000_2.txt", "a");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    for (int i = 0; i < twoD_arr_length; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            fprintf(f, "%0.15f ", *(3 * i + positions + j));
        }
        fprintf(f, "\n");
    }
    fclose(f);
    //printf("DONE TEXT FILE!\n");
}

void particle_text_file_append(double *arr, int twoD_arr_length)
{
    FILE *f = fopen("E:/Warwick/saturn/positions_particles_x_y.txt", "a");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }
    for (int j = 0; j < twoD_arr_length; j++)
    {
        for (int i = 0; i < 3; i++)
        {
            fprintf(f, "%0.15f ", *(3 * j + arr + i));
        }
        fprintf(f, "\n");
        // gap between each particle
        // if ((j+1) % n == 0)
        // {
        //     fprintf(f,"\n");
        // }
    }
    // gap between each packet of particles.
    //fprintf(f,"\n\n\n");
    fclose(f);
}

void arr1d_text_file_append(double *arr, int arr_length)
{
    FILE *f = fopen("E:/Warwick/saturn/positions_particles_huygens_mimas_massx1_115500-118000km_dt=100_n_orb=10000_sub=10000.txt", "a");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }
    for (int j = 0; j < arr_length; j++)
    {
        fprintf(f, "%0.15f ", *(j + arr));
        fprintf(f, "\n");
    }
    fclose(f);
}

void two_arr1d_text_file_append(double *arr, int arr_length)
{
    // "E:/Warwick/saturn/positions_particles_janus_epimetheus_massx1_95400km_dt=100_n_orb=10000_sub=10000.txt"
    FILE *f = fopen("E:/Warwick/saturn/positions_particles_huygens_mimas_massx1_115500-118000km_dt=100_n_orb=10000_sub=10000_ellipse_parametric.txt", "a");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }
    for (int j = 0; j < arr_length; j++)
    {
        for (int i = 0; i < 2; i++)
        {
            fprintf(f, "%0.15f ", *(2 * j + arr + i));
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

void verlet_vectorize(double *mimas_positions, double min_radius, double max_radius, int n_orbits, int n_particles, double timestep)
{
    double *old_r, *new_r;
    double *old_r1, *new_r1;
    double *old_a1, *new_a1;
    double *old_v1, *new_v1;
    double *a1, *a2;

    new_r1 = malloc(n_particles * 3 * sizeof(double));
    new_v1 = malloc(n_particles * 3 * sizeof(double));
    old_a1 = malloc(n_particles * 3 * sizeof(double));
    new_a1 = malloc(n_particles * 3 * sizeof(double));
    a1 = malloc(n_particles * 3 * sizeof(double));
    a2 = malloc(n_particles * 3 * sizeof(double));
    // for the acceleration function
    //a3 = malloc(n_particles * 3 * sizeof(double));
    //a4 = malloc(n_particles * sizeof(double));

    // total number of iterations is (int) (n_orbits * t / timestep)
    // we just need last orbit, so (int) (t/timestep)
    // so packet size should be a fraction of that, maybe 1/10 or 1/100.

    double mimas_period = 2 * M_PI * sqrt(pow(mimas_semi_major_axis, 3) / (G * saturn_mass));
    //int n_positions = (int)mimas_period / timestep;
    int packet_size = 100; // <== need to change this depending on what int(t/timestep) is. Right now it is 815, so 100 is good.
    double *packet = malloc(n_particles * packet_size * 3 * sizeof(double));

    // initialise initial position of n particles.
    old_r1 = malloc(n_particles * 3 * sizeof(double));
    old_v1 = malloc(n_particles * 3 * sizeof(double));
    for (int i = 0; i < n_particles; i++)
    {
        double fraction = ((double)i + 1) / (n_particles + 1);
        double theta = randfrom(0, 2 * M_PI);
        //double epsilon = randfrom(-1, 1);
        //double r = randfrom(min_radius,max_radius);
        double r = min_radius + (max_radius - min_radius) * fraction;
        old_r1[3*i] = r * cos(theta);
        old_r1[3*i + 1] = r * sin(theta);
        old_r1[3*i + 2] = 0;
        fraction = r / mimas_semi_major_axis;
        double v = 2 * M_PI * mimas_semi_major_axis / (mimas_period * sqrt(fraction));
        old_v1[3*i] = -1 * v * sin(theta);
        old_v1[3*i + 1] = v * cos(theta);
        old_v1[3*i + 2] = 0;
        //printf("r = %f  theta = %f  x = %f  y = %f  z = %f\n", r, theta, *(old_r1 + 3 * i), *(old_r1 + 3 * i + 1), *(old_r1 + 3 * i + 2));
    }
    printf("Done intialising.\n");
    // loop in time
    // remember to just save the last orbit, j > n = (n_orbits - 1)*mimas_period/timestep
    old_r = mimas_positions; //we just need first array so 1D array

    int t_end = n_orbits * mimas_period / timestep; // # of positions for n_orbits
    int upper_limit = (n_orbits - 1) * mimas_period / timestep; // equal to zero if n_orbits = 1
    int step = 0;
    for (int j = 1; j < t_end; j++)
    {
        //printf("%i/%i\n",j+1,t_end);
        int p = 0 * 3;
        //printf("%f %f %f\n",*(old_r1+p), *(old_r1+p+1), *(old_r1+p+2)); // first array
        new_r = mimas_positions + 3 * j;
        acceleration_particle_2darray(old_a1, old_r, old_r1, a1, a2, n_particles); // old_a1
        //printf("i= %i/%i r=%f v=%f a=%f\n", j+1, t_end, norm(old_r1 + p), norm(old_v1 + p), norm(old_a1 + p));

        scalar_mult_2darray(a1, old_v1, timestep, n_particles); // a1
        //printf("i= %i/%i r=%f v=%f a=%f\n", j+1, t_end, norm(old_r1 + p), norm(old_v1 + p), norm(old_a1 + p));
        add_2darray(a1, old_r1, a1, n_particles); // reuse a1
        //printf("i= %i/%i r=%f v=%f a=%f\n", j+1, t_end, norm(old_r1 + p), norm(old_v1 + p), norm(old_a1 + p));
        scalar_mult_2darray(a2, old_a1, 0.5 * timestep * timestep, n_particles); // a3
        add_2darray(new_r1, a1, a2, n_particles); // new_r1
        //new_r1 = add(add(old_r1,scalar_mult(old_v1,timestep)),scalar_mult(old_a1,0.5*pow(timestep,2)));
        acceleration_particle_2darray(new_a1, new_r, new_r1, a1, a2, n_particles); // new_a1

        add_2darray(a1, old_a1, new_a1, n_particles); // reuse a1
        scalar_mult_2darray(a1, a1, 0.5 * timestep, n_particles); // reuse a1
        add_2darray(new_v1, old_v1, a1, n_particles); // new_v1
        //new_v1 = add(old_v1,scalar_mult(add(old_a1,new_a1),0.5*timestep));

        if (j >= upper_limit)
        {
            printf("i= %i/%i r=%f v=%f a=%f\n", j+1, t_end, norm(old_r1 + p), norm(old_v1 + p), norm(old_a1 + p));
            step %= packet_size; // step=0,1,2,...,packet_size-1,0,1,2,...
            append_packet_2darray(packet, new_r1, step, n_particles); // length of packet currently will be: step+1
            if (step == packet_size - 1 || j == t_end - 1)
            {
                //append txt file with packet, step will then cycle back to zero on next loop and will reuse packet.
                particle_text_file_append(packet, (step + 1) * n_particles); // n_particles doesn't really matter, wont make sense in this case, and it's not used in text append function anyway.
                // (step+1) above ^ is crucial (don't use packet_size), to get the remainder to be appended properly.
            }
            step += 1; // 1,2,3,...,packet_size
        }
        replace(old_r1, new_r1, n_particles);
        replace(old_v1, new_v1, n_particles);
        // replace function is crucial to avoid memory leaks, can't yse old_x1 = new_x1 since its malloc'd array pointers.
        old_r = new_r;
    }
    free(old_r1);
    free(old_v1);
    free(old_a1);
    free(new_a1);
    free(a1);
    free(a2);
    //free(a3);
    //free(a4);
    free(packet);
}

void verlet_vectorize_subcycles(double min_radius, double max_radius, int n_orbits, int n_particles, double timestep, int n_subcycles)
{
    //double *old_r, *new_r;
    //double old_r[3], new_r[3];
    double *old_r1, *old_v1;
    double *old_a1, *new_a1;
    double *new_r1, *new_v1;
    double *a1, *a2;//, *a3, *a4;
    double *d_r1, *d_v1;

    d_r1 = malloc(n_particles * 3 * sizeof(double));
    d_v1 = malloc(n_particles * 3 * sizeof(double));
    new_r1 = malloc(n_particles * 3 * sizeof(double));
    new_v1 = malloc(n_particles * 3 * sizeof(double));
    old_a1 = malloc(n_particles * 3 * sizeof(double));
    new_a1 = malloc(n_particles * 3 * sizeof(double));
    a1 = malloc(n_particles * 3 * sizeof(double));
    a2 = malloc(n_particles * 3 * sizeof(double));
    // for the acceleration function
    //a3 = malloc(n_particles * 3 * sizeof(double));
    //a4 = malloc(n_particles * 3 * sizeof(double));

    //total number of iterations is (int) (n_orbits * t / timestep)
    //we just need last orbit, so (int) (t/timestep)
    //so packet size should be a fraction of that, maybe 1/10 or 1/100.

    double mimas_period = 2 * M_PI * sqrt(pow(mimas_semi_major_axis, 3) / (G * saturn_mass));
    //int n_positions = (int)mimas_period / timestep;
    int packet_size = 100; // <== need to change this depending on what int(t/timestep) is. Right now it is 815, so 100 is good.
    double *packet = malloc(n_particles * packet_size * 3 * sizeof(double)); // 1 if 1d_array, 3 if 2d_array

    // initialise initial position of n particles.
    old_r1 = malloc(n_particles * 3 * sizeof(double));
    old_v1 = malloc(n_particles * 3 * sizeof(double));
    for (int i = 0; i < n_particles; i++)
    {
        double fraction = ((double)i + 1) / (n_particles + 1);
        double theta = randfrom(0, 2 * M_PI);
        //double epsilon = randfrom(-1, 1);
        //double r = randfrom(min_radius,max_radius);
        double r = min_radius + (max_radius - min_radius) * fraction;
        old_r1[3*i] = r * cos(theta);
        old_r1[3*i + 1] = r * sin(theta);
        old_r1[3*i + 2] = 0;
        fraction = r / mimas_semi_major_axis;
        double v = 2 * M_PI * mimas_semi_major_axis / (mimas_period * sqrt(fraction));
        old_v1[3*i] = -1 * v * sin(theta);
        old_v1[3*i + 1] = v * cos(theta);
        old_v1[3*i + 2] = 0;
        //printf("r = %f  theta = %f  x = %f  y = %f  z = %f\n", r, theta, *(old_r1 + 3 * i), *(old_r1 + 3 * i + 1), *(old_r1 + 3 * i + 2));
    }
    printf("Done intialising.\n");
    // loop in time
    // remember to just save the last orbit, j > n = (n_orbits - 1)*mimas_period/timestep

    //double old_r[3] = {mimas_semi_major_axis, 0, 0};  // circular
    double old_r[3] = {mimas_periapsis, 0, 0};    // elliptic
    double new_r[3];

    int t_end = n_orbits * mimas_period / timestep; // # of positions for n_orbits
    // int t = mimas_period / timestep; // # of positions for 1 orbit
    //int upper_limit = (n_orbits - 1) * mimas_period / (timestep); // equal to zero if n_orbits = 1. divide by: X = # of subcycles
    int step = 0;
    /* the entire code above is the same as without subcycles, just initialising so far */
    
    for (int j = 0; j < t_end/n_subcycles; j++)
    {
        step %= packet_size; // step=0,1,2,...,packet_size-1,0,1,2,...
        
        /* radial vel */
        norm_2darray(a1, old_r1, n_particles);
        pow_neg1_array(a1, a1, n_particles);
        dot_prod_2darray(a2, old_r1, old_v1, n_particles);
        vector_mult_2darray(a2, a2, a1, n_particles);

        /* radial pos */
        norm_2darray(a1, old_r1, n_particles);

        /* length of packet currently will be: step+1 */
        //append_packet_1darray(packet, a1, step, n_particles); // r
        //append_packet_2darray(packet, old_r1, step, n_particles); // x y z
        append_packet_2_1darrays(packet, a1, a2, step, n_particles); // r v
        
        int j0 = t_end/n_subcycles;
        printf("%i/%i ",j+1,j0);
        if (step == packet_size - 1 || j == j0 - 1)
        {
            //append txt file with packet, step will then cycle back to zero on next loop and will reuse packet.
            
            /* (step+1)*n_particles is crucial (don't use packet_size), to get the remainder to be appended properly. */
            //arr1d_text_file_append(packet, (step+1)*n_particles); // r
            //particle_text_file_append(packet, (step+1)*n_particles); // x y z
            two_arr1d_text_file_append(packet, (step + 1) * n_particles); // r v
            printf("WRITE OUT");
        }
        else
        {
            printf("APPENDED %i particles", n_particles);
        }
        step += 1; // 1,2,3,...,packet_size
        printf("\n");
        
        
        zero_2d_array(d_r1, n_particles);
        zero_2d_array(d_v1, n_particles);
        for (int k = 0; k < n_subcycles; k++)
        {
            //new_r = mimas_positions + 3 + 3 * n_subcycles * (j - 1) + 3 * k;
            //new_mimas_pos(new_r, n_subcycles, t, k + j*n_subcycles + 1);
            new_mimas_pos_kepler(new_r, mimas_period, timestep, k + j*n_subcycles + 1);

            // SATURN ONLY
            add_2darray(a1, old_r1, d_r1, n_particles);
            acceleration_saturn_2darray(old_a1, a1, a2, n_particles); // old_a1
            
            scalar_mult_2darray(a1, old_v1, timestep, n_particles); // a1
            add_2darray(a1, old_r1, a1, n_particles); // reuse a1
            scalar_mult_2darray(a2, old_a1, 0.5 * timestep * timestep, n_particles); // a2
            add_2darray(new_r1, a1, a2, n_particles); // new_r1 (reusing old_r1)
            
            add_2darray(a1, new_r1, d_r1, n_particles);
            acceleration_saturn_2darray(new_a1, a1, a2, n_particles); // new_a1

            add_2darray(a1, old_a1, new_a1, n_particles); // reuse a1
            scalar_mult_2darray(a1, a1, 0.5 * timestep, n_particles); // reuse a1
            add_2darray(new_v1, old_v1, a1, n_particles); // new_v1 (reusing old_v1)

            // PERTURBATION FROM MOON ONLY
            add_2darray(a1, old_r1, d_r1, n_particles);
            acceleration_mimas_2darray(old_a1, old_r, a1, a2, n_particles); // reuse old_a1, a4
            
            scalar_mult_2darray(a1, old_a1, 0.5 * timestep * timestep, n_particles);
            add_2darray(d_r1,d_r1,a1, n_particles);
            
            add_2darray(a1, new_r1, d_r1, n_particles);
            acceleration_mimas_2darray(new_a1, new_r, a1, a2, n_particles);
            
            add_2darray(a1, old_a1, new_a1, n_particles);
            scalar_mult_2darray(a1, a1, 0.5 * timestep, n_particles);
            add_2darray(d_v1,d_v1,a1,n_particles);

            //old_r = new_r;
            replace(old_r, new_r, 1);
            replace(old_r1, new_r1, n_particles);
            replace(old_v1, new_v1, n_particles);
        }
        add_2darray(old_r1,new_r1,d_r1, n_particles);
        add_2darray(old_v1,new_v1,d_v1, n_particles);
    }
    free(old_r1);
    free(old_v1);
    free(old_a1);
    free(new_a1);
    free(a1);
    free(a2);
    //free(a3);
    //free(a4);
    free(packet);
}

void verlet_vectorize_subcycles_new(double min_radius, double max_radius, int n_orbits, int n_particles, double timestep, int n_subcycles)
{
    double *old_r, *new_r;
    double *old_r1, *old_v1;
    double *old_a1, *new_a1;
    double *new_r1, *new_v1;
    double *a1, *a2;//, *a3, *a4;
    double *b1, *b2, *b3, *b4;
    double *d_r1, *d_v1;
    // double *old_r_J, *new_r_J;
    // double *old_r_E, *new_r_E;
    // double *old_v_J, *old_v_E;

    d_r1 = malloc(n_particles * 3 * sizeof(double));
    d_v1 = malloc(n_particles * 3 * sizeof(double));
    new_r1 = malloc(n_particles * 3 * sizeof(double));
    new_v1 = malloc(n_particles * 3 * sizeof(double));
    old_a1 = malloc(n_particles * 3 * sizeof(double));
    new_a1 = malloc(n_particles * 3 * sizeof(double));
    a1 = malloc(n_particles * 3 * sizeof(double));
    a2 = malloc(n_particles * 3 * sizeof(double));
    // for Janus-Epimetheus function
    // b1 = malloc(3 * sizeof(double));
    // b2 = malloc(3 * sizeof(double));
    // b3 = malloc(3 * sizeof(double));
    // b4 = malloc(3 * sizeof(double));

    //total number of iterations is (int) (n_orbits * t / timestep)
    //we just need last orbit, so (int) (t/timestep)
    //so packet size should be a fraction of that, maybe 1/10 or 1/100.

    double mimas_period = 2 * M_PI * sqrt(pow(mimas_semi_major_axis, 3) / (G * saturn_mass)); // correct mimas_period
    //int n_positions = (int)mimas_period / timestep;
    int packet_size = 100; // need to change this depending on what n_cycles is. Right now n_cycles=815, so 100 is good.
    double *particles_packet = malloc(n_particles * packet_size * 3 * sizeof(double)); // 1 if 1d_array, 3 if 2d_array
    // double *janus_pos_packet = malloc(packet_size * 3); // x y z
    // double *epimetheus_pos_packet = malloc(packet_size * 3); // x y z
    double *mimas_pos_packet = malloc(packet_size * 3); // x y z

    /* Initialise initial position of n particles. */
    old_r1 = malloc(n_particles * 3 * sizeof(double));
    old_v1 = malloc(n_particles * 3 * sizeof(double));
    double r, v, fraction, theta;
    for (int i = 0; i < n_particles; i++)
    {
        fraction = ((double)i + 1) / (n_particles + 1);
        theta = randfrom(0, 2 * M_PI);
        //r = min_radius + (max_radius - min_radius) * fraction;
        r = randfrom(min_radius, max_radius);
        old_r1[3 * i] = r * cos(theta);
        old_r1[3 * i + 1] = r * sin(theta);
        old_r1[3 * i + 2] = 0;
        fraction = r / mimas_semi_major_axis;
        v = 2 * M_PI * mimas_semi_major_axis / (mimas_period * sqrt(fraction));
        old_v1[3 * i] = -1 * v * sin(theta);
        old_v1[3 * i + 1] = v * cos(theta);
        old_v1[3 * i + 2] = 0;
        //printf("r = %f  theta = %f  x = %f  y = %f  z = %f\n", r, theta, *(old_r1 + 3 * i), *(old_r1 + 3 * i + 1), *(old_r1 + 3 * i + 2));
    }
    printf("Done intialising.\n");
    // loop in time
    // remember to just save the last orbit, j > n = (n_orbits - 1)*mimas_period/timestep

    /* Mimas */
    //double old_r[3] = {mimas_semi_major_axis, 0, 0};  // circular
    // old_r[3] = {mimas_periapsis, 0, 0};    // elliptic
    // double new_r[3];

    old_r[0] = mimas_periapsis;
    old_r[1] = 0;
    old_r[2] = 0; // equivalent to *(old_r + 2) https://stackoverflow.com/a/11625225/7875204

    /* Janus */
    // double old_r_J[3] = {151460, 0, 0};
    v = sqrt(G * saturn_mass / 151460);
    // double old_v_J[3] = {0, v, 0};
    // double new_r_J[3];
    
    // old_r_J[0] = 151460;
    // old_r_J[1] = 0;
    // old_r_J[2] = 0;

    // old_v_J[0] = v;
    // old_v_J[1] = 0;
    // old_v_J[2]  = 0;

    /* Epimetheus */
    // double old_r_E[3] = {-151410, 0, 0};
    v = sqrt(G * saturn_mass / 151410);
    // double old_v_E[3] = {0, -v, 0};
    // double new_r_E[3];
    
    // old_r_E[0] = -151410;
    // old_r_E[1] = 0;
    // old_r_E[2] = 0;

    // old_v_E[0] = v;
    // old_v_E[1] = 0;
    // old_v_E[2]  = 0;


    int t_end = n_orbits * mimas_period / (timestep); // # of positions for n_orbits
    // POSITIONS FOR JANUS AND EPIMETHEUS
    int n_cycles = t_end/n_subcycles; // = 815 positions

    // int t = mimas_period / timestep; // # of positions for 1 orbit
    //int upper_limit = (n_orbits - 1) * mimas_period / (timestep); // equal to zero if n_orbits = 1. divide by: X = # of subcycles
    int step = 0;
    /* the entire code above is the same as without subcycles, just initialising so far */
    
    for (int j = 0; j < n_cycles; j++)
    {
        step = j % packet_size;                    // step=0,1,2,...,packet_size-1,0,1,2,...
        
        /* radial vel */
        norm_2darray(a1, old_r1, n_particles);
        pow_neg1_array(a1, a1, n_particles);
        dot_prod_2darray(a2, old_r1, old_v1, n_particles);
        vector_mult_2darray(a2, a2, a1, n_particles);

        /* radial pos */
        norm_2darray(a1, old_r1, n_particles); // is it better to use less variables or less function calls? (memory vs compute)

        /* length of packet currently will be: step+1 */
        //append_packet_1darray(particles_packet, a1, step, n_particles); // r
        //append_packet_2darray(particles_packet, old_r1, step, n_particles); // x y z
        append_packet_2_1darrays(particles_packet, a1, a2, step, n_particles); // r v

        /* append last Janus and Epimetheus position at the end of each cycle to their packets */
        // append_packet_2darray(janus_pos_packet, old_r_J, step, 1);
        // append_packet_2darray(epimetheus_pos_packet, old_r_E, step, 1);
        append_packet_2darray(mimas_pos_packet, old_r, step, 1);
        
        /* WRITE OUT to text files at the end of each packet (100 cycles) or at the very end */
        printf("%i/%i ", j+1, n_cycles);
        if (step == packet_size - 1 || j == n_cycles - 1)
        {
            //append txt file with packet, step will then cycle back to zero on next loop and will reuse packet.
            
            /* (step+1)*n is crucial (don't use packet_size), to get the remainder to be appended properly. */
            //arr1d_text_file_append(particles_packet, (step + 1) * n); // r
            //particle_text_file_append(particles_packet, (step+1)*n); // x y z
            two_arr1d_text_file_append(particles_packet, (step + 1) * n_particles); // r v

            /* Write out Janus and Epimetheus positions */
            // moon1_positions_to_text_file(janus_pos_packet, (step + 1) * 1);
            // moon2_positions_to_text_file(epimetheus_pos_packet, (step + 1) * 1);

            /* Write out Mimas positions */
            moon_positions_to_text_file(mimas_pos_packet, (step + 1) * 1);

            printf("WRITE OUT");
        }
        else
        {
            printf("APPENDED %i particles", n_particles);
        }
        //step += 1; // 1,2,3,...,packet_size
        printf("\n");
        
        /* Start of NEW CYCLE */
        zero_2d_array(d_r1, n_particles); // Reset d_r1
        zero_2d_array(d_v1, n_particles); // Reset d_v1
        for (int k = 0; k < n_subcycles; k++)
        {
            /* SATURN ONLY */
            //add_2darray(a1, old_r1, d_r1, n_particles);
            acceleration_saturn_2darray(old_a1, old_r1, a2, n_particles); // old_a1
            
            scalar_mult_2darray(a1, old_v1, timestep, n_particles); // a1
            add_2darray(a1, old_r1, a1, n_particles); // reuse a1
            scalar_mult_2darray(a2, old_a1, 0.5 * timestep * timestep, n_particles); // a2
            add_2darray(new_r1, a1, a2, n_particles); // new_r1
            
            //add_2darray(a1, new_r1, d_r1, n_particles);
            acceleration_saturn_2darray(new_a1, new_r1, a2, n_particles); // new_a1

            add_2darray(a1, old_a1, new_a1, n_particles); // reuse a1
            scalar_mult_2darray(a1, a1, 0.5 * timestep, n_particles); // reuse a1
            add_2darray(new_v1, old_v1, a1, n_particles); // new_v1

            /* PERTURBATION FROM MOONS ONLY */

            // MIMAS
            // new_mimas_pos(new_r, n_subcycles, t, k + j*n_subcycles + 1); // USE KEPLER (ELLIPSE)
            // new_mimas_pos_kepler(new_r, mimas_period, timestep, k + j*n_subcycles + 1); // USE PARAMETRIC
            new_mimas_pos_kepler_parametric(new_r, mimas_period, timestep, k + j*n_subcycles + 1);

            acceleration_mimas_2darray(old_a1, old_r, old_r1, a2, n_particles); // reuse old_a1, a4
            
            add_2darray(d_r1, d_r1, old_a1, n_particles);
            
            acceleration_mimas_2darray(new_a1, new_r, new_r1, a2, n_particles);
            
            add_2darray(d_v1, old_a1, new_a1, n_particles);

            /* JANUS-EPIMETHEUS */
            // new_janus_epimetheus_pos(new_r_J, new_r_E, old_r_J, old_r_E, old_v_J, old_v_E, a1, a2, b1, b2, b3, b4, timestep, k + j*n_subcycles + 1); // get new pos, and store new v in old.
            // // JANUS
            // acceleration_janus_2darray(old_a1, old_r_J, old_r1, a2, n_particles);
            // add_2darray(d_r1, d_r1, old_a1, n_particles);
            // acceleration_janus_2darray(new_a1, new_r_J, new_r1, a2, n_particles);
            // add_2darray(d_v1, d_v1, old_a1, n_particles);
            // add_2darray(d_v1, d_v1, new_a1, n_particles);

            // // EPIMETHEUS
            // acceleration_epimetheus_2darray(old_a1, old_r_E, old_r1, a2, n_particles);
            // add_2darray(d_r1, d_r1, old_a1, n_particles);
            // acceleration_epimetheus_2darray(new_a1, new_r_E, new_r1, a2, n_particles);
            // add_2darray(d_v1, d_v1, old_a1, n_particles);
            // add_2darray(d_v1, d_v1, new_a1, n_particles);

            //scalar_mult_2darray(d_r1, d_r1, 0.5 * timestep * timestep, n);
            //scalar_mult_2darray(a1, a1, 0.5 * timestep, n);
            //add_2darray(d_v1,d_v1,a1,n);


            // replace function is crucial to avoid memory leaks, can't use old_x1 = new_x1 since its malloc'd array pointers.
            // replace(old_r, new_r, 1);       // Mimas
            // replace(old_r_J, new_r_J, 1);   // Janus
            // replace(old_r_E, new_r_E, 1);   // Epimetheus
            // replace(old_r1, new_r1, n_particles);
            // replace(old_v1, new_v1, n_particles);

            old_r = new_r;
            // old_r_J = new_r_J;
            // old_r_E = new_r_E;
            old_r1 = new_r1;
            old_v1 = new_v1;
        }
        scalar_mult_2darray(d_r1, d_r1, 0.5 * timestep * timestep, n_particles);
        scalar_mult_2darray(d_v1, d_v1, 0.5 * timestep, n_particles);
        //printf("dr=%0.15f %0.15f, dv=%0.15f %0.15f\n",d_r1[0],d_r1[1],d_v1[0],d_v1[1]);
        add_2darray(old_r1,new_r1,d_r1, n_particles);
        add_2darray(old_v1,new_v1,d_v1, n_particles);
    }
    free(old_r1);
    free(old_v1);
    free(old_a1);
    free(new_a1);
    free(a1);
    free(a2);
    // free(b1);
    // free(b2);
    // free(b3);
    // free(b4);
    free(particles_packet);
    // free(janus_pos_packet);
    // free(epimetheus_pos_packet);
    free(mimas_pos_packet);
}

void new_mimas_pos(double *new_r, int n_subcycles, int t, int k)
{
    double fraction = (double) k / t;
    new_r[0] = mimas_semi_major_axis * cos(2*M_PI*(fraction));
    new_r[1] = mimas_semi_major_axis * sin(2*M_PI*(fraction));
    new_r[2] = 0;
}

void new_mimas_pos_kepler(double *new_r, double mimas_period, double timestep, int k)
{
    double n = 2 * M_PI / mimas_period; // mean motion
    double M = n * k * timestep;        // mean anomaly

    // approx eccentric anomaly using Newton method
    double E = M; // initial guess
    for (int i = 0; i < 20; i++)
    {
        E -= (E - mimas_ecc * sin(E) - M) / (1 - mimas_ecc * cos(E));
    }

    new_r[0] = mimas_semi_major_axis * (cos(E) - mimas_ecc);
    new_r[1] = mimas_semi_major_axis * sqrt(1 - mimas_ecc*mimas_ecc) * sin(E);
    new_r[2] = 0;
}

void new_mimas_pos_kepler_parametric(double *new_r, double mimas_period, double timestep, int k)
{
    double a = mimas_semi_major_axis;
    double b = a * sqrt(1 - mimas_ecc*mimas_ecc);
    double c = mimas_periapsis - a;
    double fraction = k*timestep/mimas_period;
    new_r[0] = a * cos(2*M_PI*fraction) + c;
    new_r[1] = b * sin(2*M_PI*fraction);
    new_r[2] = 0;
}

void new_janus_epimetheus_pos(double *new_r_J, double *new_r_E, double *old_r_J, double *old_r_E, double *old_v_J, double *old_v_E, double *a1, double *a2, double *b1, double *b2, double *b3, double *b4, double timestep, int k)
{
    /* new_r_J */
    acceleration_saturn_1darray(b1, old_r_J); // old_a_J_sat
    acceleration_janus_from_epimetheus_1darray(b2, old_r_J, old_r_E); // old_a_J_E
    add(a2, b1, b2);

    scalar_mult(a1, old_v_J, timestep); // a1
    add(a1, old_r_J, a1); // reuse a1
    scalar_mult(a2, a2, 0.5 * timestep * timestep); // old_a_J * 0.5 * t^2
    add(new_r_J, a1, a2); // new_r1 (reusing old_r1)

    /* new_r_E */
    acceleration_saturn_1darray(b3, old_r_E); // old_a_E_sat
    acceleration_epimetheus_from_janus_1darray(b4, old_r_E, old_r_J); // old_a_E_J
    add(a2, b3, b4);

    scalar_mult(a1, old_v_E, timestep); // a1
    add(a1, old_r_E, a1); // reuse a1
    scalar_mult(a2, a2, 0.5 * timestep * timestep); // old_a_E * 0.5 * t^2
    add(new_r_E, a1, a2); // new_r1 (reusing old_r1)

    /* new_v_J */
    acceleration_saturn_1darray(a2, new_r_J); // new_a_J_sat
    add(a1, a2, b1); // new_a_J_sat + old_a_J_sat
    acceleration_janus_from_epimetheus_1darray(a2, new_r_J, new_r_E); // new_a_J_E
    add(a1, a1, a2);
    add(a1, a1, b2);
    scalar_mult(a1, a1, 0.5 * timestep);
    add(old_v_J, old_v_J, a1);

    /* new_v_E */
    acceleration_saturn_1darray(a2, new_r_E); // new_a_E_sat
    add(a1, a2, b3); // new_a_J_sat + old_a_J_sat
    acceleration_epimetheus_from_janus_1darray(a2, new_r_E, new_r_J); // new_a_J_E
    add(a1, a1, a2);
    add(a1, a1, b4);
    scalar_mult(a1, a1, 0.5 * timestep);
    add(old_v_E, old_v_E, a1);
}

void replace(double *arr1, double *arr2, int arr_length)
{
    #pragma GCC ivdep
    for (int i = 0; i < 3 * arr_length; i++)
    {
        //*(arr1 + i) = *(arr2 + i);
        arr1[i] = arr2[i];
    }
}

void append_packet_2darray(double *packet, double *arr2d, int step, int n_particles)
{
    // arr = [[x1,y1,z1],[x2,y2,z2],...,[xn,yn,zn]] for n different particles
    // need to append arr to packet, but the particles need to be sorted so they are in the right area.
    for (int j = 0; j < n_particles; j++)
    {
        #pragma GCC ivdep
        for (int i = 0; i < 3; i++)
        {
            //*(packet + step * n_particles + 3 * j + i) = *(twoDarr + 3 * j + i);
            packet[step*3*n_particles + 3*j + i] = arr2d[3*j + i];
        }
    }
}

void append_packet_1darray(double *packet, double *arr1d, int step, int n_particles)
{
    #pragma GCC ivdep
    for (int j = 0; j < n_particles; j++)
    {
        packet[step*n_particles + j] = arr1d[j];
    }
}

void append_packet_2_1darrays(double *packet, double *arr1d_1, double *arr1d_2, int step, int n_particles)
{
    #pragma GCC ivdep
    for (int j = 0; j < n_particles; j++)
    {
        packet[step*2*n_particles + 2*j] = arr1d_1[j];
        packet[step*2*n_particles + 2*j + 1] = arr1d_2[j];
    }
}

void add_2darray(double *output, double *arr1, double *arr2, int arr_length)
{
    #pragma GCC ivdep
    for (int j = 0; j < arr_length; j++)
    {
        //ans[j][3];
        // for (int i = 0; i < 3; i++)
        // {
        //     ans[3*j+i] = *(arr1 + 3*j+i) + *(arr2 + 3*j+i);
        //     //ans[3*j+i] = arr1[3*j+i] + arr2[3*j+i];
        // }
        output[3 * j] = arr1[3 * j] + arr2[3 * j];
        output[3 * j + 1] = arr1[3 * j + 1] + arr2[3 * j + 1];
        output[3 * j + 2] = arr1[3 * j + 2] + arr2[3 * j + 2];
    }
}

void subtract_2darray(double *output, double *arr1, double *arr2, int arr_length)
{
    for (int j = 0; j < arr_length; j++)
    {
        #pragma GCC ivdep
        for (int i = 0; i < 3; i++)
        {
            //output[3 * j + i] = *(arr1 + 3 * j + i) - *(arr2 + 3 * j + i);
            output[3*j+i] = arr1[3*j+i] - arr2[3*j+i];
        }
    }
}

void subtract_2darray_1darray(double *output, double *arr2d, double *arr1d, int arr_length)
{
    for (int j = 0; j < arr_length; j++)
    {
        #pragma GCC ivdep
        for (int i = 0; i < 3; i++)
        {
            //output[3 * j + i] = *(arr1 + 3 * j + i) - *(arr2 + 3 * j + i);
            output[3*j+i] = arr2d[3*j+i] - arr1d[i];
        }
    }
}

void scalar_mult_2darray(double *output, double *arr, double scalar, int arr_length)
{
    // scalar_mult_2darray([[1,2,3],[3,0,0]],4,2) --> [[4,8,12],[12,0,0]]
    for (int j = 0; j < arr_length; j++)
    { 
        #pragma GCC ivdep
        for (int i = 0; i < 3; i++)
        {
            //output[3 * j + i] = *(arr + 3 * j + i) * scalar;
            output[3*j+i] = arr[3*j+i] * scalar;
        }
    }
}

void dot_prod_2darray(double *output, double *arr1, double *arr2, int arr_length)
{
    #pragma GCC ivdep
    for (int i = 0; i < arr_length; i++)
    {
        output[i] = arr1[3*i] * arr2[3*i] + arr1[3*i + 1] * arr2[3*i + 1] + arr1[3*i + 2] * arr2[3*i + 2];
    }
}

void vector_mult_2darray(double *output, double *arr, double *scalars, int arr_length)
{
    // vector_mult_2darray([[1,2,3],[3,0,0]], [2,3], 2) --> [[2,4,6],[9,0,0]]
    for (int j = 0; j < arr_length; j++)
    {
        #pragma GCC ivdep
        for (int i = 0; i < 3; i++)
        {
            //output[3 * j + i] = *(arr + 3 * j + i) * *(scalars + j);
            output[3*j+i] = arr[3*j+i] * scalars[j];
        }
    }
}

void array_to_2darray(double *output, double *arr1d, int arr_length)
{
    // array_to_2darray([1,2,3],n) --> [[1,2,3],[1,2,3],...]
    for (int j = 0; j < arr_length; j++)
    {
        #pragma GCC ivdep
        for (int i = 0; i < 3; i++)
        {
            //output[3 * j + i] = *(arr1d + i);
            output[3*j+i] = arr1d[i];
        }
    }
}

void norm_2darray(double *output, double *arr, int arr_length)
{
    // norm_2darray([[3,4,5],[1,1,1],[3,0,0],[0,2,0]],4) --> [sqrt(3^2+4^2+5^2),sqrt(3),3,2]
    for (int i = 0; i < arr_length; i++)
    {
        output[i] = norm(arr + 3 * i);
    }
}

void pow_neg1_array(double *output, double *arr, int arr_length)
{
    // pow_array([1,2,3,4],power=3,arr_length=4) --> [1,8,27,64]
    // pow_neg3_array([1,2,3,4],4) --> [1,1/8,1/27,1/64]
    #pragma GCC ivdep
    for (int i = 0; i < arr_length; i++)
    {
        //output[i] = 1 / ((*(arr + i)) * (*(arr + i)) * (*(arr + i)));
        output[i] = 1/(arr[i]);
    }
}

void pow_neg3_array(double *output, double *arr, int arr_length)
{
    // pow_array([1,2,3,4],power=3,arr_length=4) --> [1,8,27,64]
    // pow_neg3_array([1,2,3,4],4) --> [1,1/8,1/27,1/64]
    #pragma GCC ivdep
    for (int i = 0; i < arr_length; i++)
    {
        //output[i] = 1 / ((*(arr + i)) * (*(arr + i)) * (*(arr + i)));
        output[i] = 1/(arr[i]*arr[i]*arr[i]);
    }
}

void scalar_mult_1darray(double *output, double *arr, double scalar, int arr_length)
{
    // scalar_mult_2darray([[1,2,3],[3,0,0]],4,2) --> [[4,8,12],[12,0,0]]
    #pragma GCC ivdep
    for (int i = 0; i < arr_length; i++)
    {
        //output[i] = *(arr + i) * scalar;
        output[i] = arr[i] * scalar;
    }
}

void acceleration_particle_2darray(double *output, double *r_mimas, double *r_particles, double *arr2d, double *arr1d, int arr_length)
{
    //this function takes 3 array, a 3*n array and an arr_length, and returns a 3*n array; acceleration vectors for all n particle arrays.
    // arr_length = n
    // r_mimas = [x,y,z]
    // r_particle = [[x0,y0,z0],...,...] initial positions for all n particles.
    // r_mimas_2d_array = [[x,y,z],[x,y,z],...] same vector n times.

    // r_mimas is just a size 3 array, r_particle is a 3*n array.
    // need to first scalar mult an array with all ones with r_mimas.
    // then use subtract_2darray(r_particle,r_mimas)

    //double *r2 = malloc(3*arr_length*sizeof(double));
    //double *a1 = malloc(arr_length*sizeof(double));
    
    // r2 = arr2d, a1 = arr1d <-- important
    
    // need to get norms for all arrays (size 3) in these 3*n arrays.
    norm_2darray(arr1d, r_particles, arr_length);
    double scalar1 = -1 * saturn_mass * G;
    pow_neg3_array(arr1d, arr1d, arr_length);
    scalar_mult_1darray(arr1d, arr1d, scalar1, arr_length); // a1 and a2 are NOT 2d arrays, use 1darray method
    vector_mult_2darray(output, r_particles, arr1d, arr_length);

    subtract_2darray_1darray(arr2d, r_particles, r_mimas, arr_length);
    norm_2darray(arr1d, arr2d, arr_length); // reusing a1, now that it wasn't needed anymore
    double scalar2 = -1 * mimas_mass * G; // how could I forget this...
    pow_neg3_array(arr1d, arr1d, arr_length);
    scalar_mult_1darray(arr1d, arr1d, scalar2, arr_length);
    vector_mult_2darray(arr2d, arr2d, arr1d, arr_length); // finished with a1 again
    
    add_2darray(output, output, arr2d, arr_length);

    //free(r2);
    //free(a1); don't need to free anymore
}

void acceleration_saturn_2darray(double *output, double *r_particles, double *arr1d, int arr_length)
{
    //this function takes 3 array, a 3*n array and an arr_length, and returns a 3*n array; acceleration vectors for all n particle arrays.
    // arr_length = n
    // r_mimas = [x,y,z]
    // r_particle = [[x0,y0,z0],...,...] initial positions for all n particles.
    // r_mimas_2d_array = [[x,y,z],[x,y,z],...] same vector n times.

    // r_mimas is just a size 3 array, r_particle is a 3*n array.
    // need to first scalar mult an array with all ones with r_mimas.
    // then use subtract_2darray(r_particle,r_mimas)

    //double *r2 = malloc(3*arr_length*sizeof(double));
    //double *a1 = malloc(arr_length*sizeof(double));
    
    // r2 = arr2d, a1 = arr1d <-- important
    
    // need to get norms for all arrays (size 3) in these 3*n arrays.
    norm_2darray(arr1d, r_particles, arr_length);
    double scalar1 = -1 * saturn_mass * G;
    pow_neg3_array(arr1d, arr1d, arr_length);
    scalar_mult_1darray(arr1d, arr1d, scalar1, arr_length); // a1 and a2 are NOT 2d arrays, use 1darray method
    vector_mult_2darray(output, r_particles, arr1d, arr_length);

    //subtract_2darray_1darray(arr2d, r_particles, r_mimas, arr_length);
    //norm_2darray(arr1d, arr2d, arr_length); // reusing a1, now that it wasn't needed anymore
    //double scalar2 = -1 * mimas_mass * G; // how could I forget this...
    //pow_neg3_array(arr1d, arr1d, arr_length);
    //scalar_mult_1darray(arr1d, arr1d, scalar2, arr_length);
    //vector_mult_2darray(arr2d, arr2d, arr1d, arr_length); // finished with a1 again
    
    //add_2darray(output, output, arr2d, arr_length);

    //free(r2);
    //free(a1); don't need to free anymore
}

void acceleration_mimas_2darray(double *output, double *r_mimas, double *r_particles, double *arr1d, int arr_length)
{
    //this function takes 3 array, a 3*n array and an arr_length, and returns a 3*n array; acceleration vectors for all n particle arrays.
    // arr_length = n
    // r_mimas = [x,y,z]
    // r_particle = [[x0,y0,z0],...,...] initial positions for all n particles.
    // r_mimas_2d_array = [[x,y,z],[x,y,z],...] same vector n times.

    // r_mimas is just a size 3 array, r_particle is a 3*n array.
    // need to first scalar mult an array with all ones with r_mimas.
    // then use subtract_2darray(r_particle,r_mimas)

    //double *r2 = malloc(3*arr_length*sizeof(double));
    //double *a1 = malloc(arr_length*sizeof(double));
    
    // r2 = arr2d, a1 = arr1d <-- important
    
    // need to get norms for all arrays (size 3) in these 3*n arrays.
    //norm_2darray(arr1d, r_particles, arr_length);
    //double scalar1 = -1 * saturn_mass * G;
    //pow_neg3_array(arr1d, arr1d, arr_length);
    //scalar_mult_1darray(arr1d, arr1d, scalar1, arr_length); // a1 and a2 are NOT 2d arrays, use 1darray method
    //vector_mult_2darray(output, r_particles, arr1d, arr_length);

    subtract_2darray_1darray(output, r_particles, r_mimas, arr_length);
    norm_2darray(arr1d, output, arr_length); // reusing a1, now that it wasn't needed anymore
    double scalar2 = -1 * mimas_mass * G; // how could I forget this...
    pow_neg3_array(arr1d, arr1d, arr_length);
    scalar_mult_1darray(arr1d, arr1d, scalar2, arr_length);
    vector_mult_2darray(output, output, arr1d, arr_length); // finished with a1 again
    
    //add_2darray(output, output, arr2d, arr_length);

    //free(r2);
    //free(a1); don't need to free anymore
}

void acceleration_janus_2darray(double *output, double *r_janus, double *r_particles, double *arr1d, int arr_length)
{
    subtract_2darray_1darray(output, r_particles, r_janus, arr_length);
    norm_2darray(arr1d, output, arr_length); // reusing a1, now that it wasn't needed anymore
    double scalar2 = -1 * janus_mass * G; // how could I forget this...
    pow_neg3_array(arr1d, arr1d, arr_length);
    scalar_mult_1darray(arr1d, arr1d, scalar2, arr_length);
    vector_mult_2darray(output, output, arr1d, arr_length); // finished with a1 again
}

void acceleration_epimetheus_2darray(double *output, double *r_epimetheus, double *r_particles, double *arr1d, int arr_length)
{
    subtract_2darray_1darray(output, r_particles, r_epimetheus, arr_length);
    norm_2darray(arr1d, output, arr_length); // reusing a1, now that it wasn't needed anymore
    double scalar2 = -1 * epimetheus_mass * G; // how could I forget this...
    pow_neg3_array(arr1d, arr1d, arr_length);
    scalar_mult_1darray(arr1d, arr1d, scalar2, arr_length);
    vector_mult_2darray(output, output, arr1d, arr_length); // finished with a1 again
}

void acceleration_janus_from_epimetheus_1darray(double *output, double *r_janus, double *r_epimetheus)
{
    // arr_length = 1
    subtract(output, r_janus, r_epimetheus);
    double a1 = norm(output); // reusing a1, now that it wasn't needed anymore
    double scalar2 = -1 * epimetheus_mass * G; // how could I forget this...
    a1 = scalar2/(a1*a1*a1);
    scalar_mult(output, output, a1); // finished with a1 again
}

void acceleration_epimetheus_from_janus_1darray(double *output, double *r_epimetheus, double *r_janus)
{// arr_length = 1
    subtract(output, r_epimetheus, r_janus);
    double a1 = norm(output); // reusing a1, now that it wasn't needed anymore
    double scalar2 = -1 * janus_mass * G; // how could I forget this...
    a1 = scalar2/(a1*a1*a1);
    scalar_mult(output, output, a1); // finished with a1 again
}

void acceleration_saturn_1darray(double *output, double *r_pos)
{
    // arr_length = 1
    double a1 = norm(r_pos); // reusing a1, now that it wasn't needed anymore
    double scalar2 = -1 * saturn_mass * G; // how could I forget this...
    a1 = scalar2/(a1*a1*a1);
    scalar_mult(output, r_pos, a1); // finished with a1 again
}

void zero_2d_array(double *output, int arr_length)
{
    for (int i = 0; i < 3*arr_length; i++)
    {
        output[i] = 0;
    }
}