#include <vector>
#include <cmath>
#include "forces.h"
#include <cstdio>

namespace {
    double theta = falcON::Default::theta;
    int kernel = falcON::Default::kernel;
    bool individual_eps = false;
    double eps = 0.0;  // global softening length
    double code_time = 0.0;
    double dt = 0.0;
    double t0 = 0.0;

    // falcON snapshot class holding particle attributes
    falcON::bodies bodies(falcON::fieldset(  // which particle properties are allocated:
        falcON::fieldset::gravity |          // pos,vel,mass,potential,acceleration,flags;
        falcON::fieldset::e)                 // [if needed], also individual softening lengths
    );

    void compute_force()
    {
        // construct the octtree from all particles
        falcON::forces forces(&bodies,
            /*global softening length*/ eps,
            /*tree opening angle*/ theta,
            /*type of softening kernel*/ static_cast<falcON::kern_type>(kernel),
            /*whether to use individual softening lengths*/ individual_eps);
        forces.grow();

        // compute the force for active particles only
        forces.approximate_gravity(/*use all particles*/ true);
    }

    void leapfrog_step(double dt)
    {
        // 1. update velocity for the first half-step
        for(falcON::body b=bodies.begin_all_bodies(); b; ++b) {
            b.vel() += (0.5*dt) * b.acc();
        }
        // 2. update position for a full step
        for(falcON::body b=bodies.begin_all_bodies(); b; ++b) {
            b.pos() += dt * b.vel();
        }
        // 3. recompute accelerations at the end of the step
        compute_force();
        // 4. update velocity for the second half-step
        for(falcON::body b=bodies.begin_all_bodies(); b; ++b) {
            b.vel() += (0.5*dt) * b.acc();
        }
    }
}  // namespace

int new_particle (int *index, double mass, double x, double y, double z, double vx, double vy, double vz, double radius)
{
    *index = bodies.N_bodies(falcON::bodytype::std);
    falcON::body b = bodies.new_body(falcON::bodytype::std, 1024);
    b.mass() = mass;
    b.pos()[0] = x;
    b.pos()[1] = y;
    b.pos()[2] = z;
    b.vel()[0] = vx;
    b.vel()[1] = vy;
    b.vel()[2] = vz;
    b.eps() = radius;
    return 0;
}

int delete_particle (int index_of_the_particle)
{
    return -1;  // not implemented
}

int get_number_of_particles (int *value)
{
    *value = bodies.N_bodies(falcON::bodytype::std);
    return 0;
}

int get_index_of_first_particle (int *index)
{
    if(bodies.N_bodies(falcON::bodytype::std) > 0) {
        *index = 0;
        return 0;
    }
    return 1;
}

int get_index_of_next_particle (int index, int *next_index)
{
    if(index+1 < bodies.N_bodies(falcON::bodytype::std)) {
        *next_index = index+1;
        return 0;
    }
    return 1;
}

int evolve_model (double tend)
{
    if(dt==0 || (!individual_eps && eps==0)) {   // the two required params are not set
        printf("Timestep and/or softening length are not initialized\n");
        return -1;
    }
    while(code_time + t0 < tend) {
        double step = fmin(tend - code_time, dt);
        code_time += step;
        leapfrog_step(step);
    }
    return 0;
}

int get_time (double *current_time)
{
    *current_time = code_time + t0;
    return 0;
}

int get_mass (int index, double *mass)
{
    if(index >= bodies.N_bodies(falcON::bodytype::std)) { return -1; }
    falcON::body b = bodies.begin_all_bodies(); b+= index;
    *mass = b.mass();
    return 0;
}

int set_mass (int index, double  mass)
{
    if(index >= bodies.N_bodies(falcON::bodytype::std)) { return -1; }
    falcON::body b = bodies.begin_all_bodies(); b+= index;
    b.mass() = mass;
    return 0;
}

int get_radius (int index, double *radius)
{
    if(index >= bodies.N_bodies(falcON::bodytype::std)) { return -1; }
    falcON::body b = bodies.begin_all_bodies(); b+= index;
    *radius = b.eps();
    return 0;
}

int set_radius (int index, double  radius)
{
    if(index >= bodies.N_bodies(falcON::bodytype::std)) { return -1; }
    falcON::body b = bodies.begin_all_bodies(); b+= index;
    b.eps() = radius;
    return 0;
}


int get_position (int index, double *x, double *y, double *z)
{
    if(index >= bodies.N_bodies(falcON::bodytype::std)) { return -1; }
    falcON::body b = bodies.begin_all_bodies(); b+= index;
    *x = b.pos()[0];
    *y = b.pos()[1];
    *z = b.pos()[2];
    return 0;
}

int set_position (int index, double  x, double  y, double  z)
{
    if(index >= bodies.N_bodies(falcON::bodytype::std)) { return -1; }
    falcON::body b = bodies.begin_all_bodies(); b+= index;
    b.pos()[0] = x;
    b.pos()[1] = y;
    b.pos()[2] = z;
    return 0;
}

int get_velocity (int index, double *vx, double *vy, double *vz)
{
    if(index >= bodies.N_bodies(falcON::bodytype::std)) { return -1; }
    falcON::body b = bodies.begin_all_bodies(); b+= index;
    *vx = b.vel()[0];
    *vy = b.vel()[1];
    *vz = b.vel()[2];
    return 0;
}

int set_velocity (int index, double  vx, double  vy, double  vz)
{
    if(index >= bodies.N_bodies(falcON::bodytype::std)) { return -1; }
    falcON::body b = bodies.begin_all_bodies(); b+= index;
    b.vel()[0] = vx;
    b.vel()[1] = vy;
    b.vel()[2] = vz;
    return 0;
}

int get_acceleration (int index, double *ax, double *ay, double *az)
{
    if(index >= bodies.N_bodies(falcON::bodytype::std)) { return -1; }
    falcON::body b = bodies.begin_all_bodies(); b+= index;
    *ax = b.acc()[0];
    *ay = b.acc()[1];
    *az = b.acc()[2];
    return 0;
}

int set_acceleration (int index, double  ax, double  ay, double  az)
{
    if(index >= bodies.N_bodies(falcON::bodytype::std)) { return -1; }
    falcON::body b = bodies.begin_all_bodies(); b+= index;
    b.acc()[0] = ax;
    b.acc()[1] = ay;
    b.acc()[2] = az;
    return 0;
}

int get_potential (int index, double *potential)
{
    if(index >= bodies.N_bodies(falcON::bodytype::std)) { return -1; }
    falcON::body b = bodies.begin_all_bodies(); b+= index;
    *potential = b.pot();
    return 0;
}

int get_state (int index, double *mass, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *radius)
{
    if(index >= bodies.N_bodies(falcON::bodytype::std)) { return -1; }
    falcON::body b = bodies.begin_all_bodies(); b+= index;
    *x  = b.pos()[0];
    *y  = b.pos()[1];
    *z  = b.pos()[2];
    *vx = b.vel()[0];
    *vy = b.vel()[1];
    *vz = b.vel()[2];
    *mass = b.mass();
    *radius = b.eps();
    return 0;
}

int set_state (int index, double  mass, double  x, double  y, double  z, double  vx, double  vy, double  vz, double  radius)
{
    if(index >= bodies.N_bodies(falcON::bodytype::std)) { return -1; }
    falcON::body b = bodies.begin_all_bodies(); b+= index;
    b.pos()[0] = x;
    b.pos()[1] = y;
    b.pos()[2] = z;
    b.vel()[0] = vx;
    b.vel()[1] = vy;
    b.vel()[2] = vz;
    b.mass() = mass;
    b.eps() = radius;
    return 0;
}

int initialize_code () {
    //printf("falcon::initialize_code\n");
    return 0;
}

int cleanup_code () {
    //printf("falcon::cleanup_code\n");
    return 0;
}

int set_individual_epsilon (bool individual_epsilon)
{
    individual_eps = individual_epsilon;
    return 0;
}
int get_eps2 (double *eps2)
{
    *eps2 = eps*eps;
    return 0;
}

int set_eps2 (double  eps2)
{
    eps = sqrt(eps2);
    return 0;
}

int get_epsilon (double *epsilon)
{
    *epsilon = eps;
    return 0;
}

int set_epsilon (double  epsilon)
{
    eps = epsilon;
    return 0;
}

int get_time_step (double *timestep)
{
    *timestep = dt;
    return 0;
}

int set_time_step (double  timestep)
{
    dt = timestep;
    return 0;
}

int get_begin_time (double *begin_time)
{
    *begin_time = t0;
    return 0;
}

int set_begin_time (double  begin_time)
{
    t0 = begin_time;
    return 0;
}

int commit_parameters ()
{
    //printf("falcon::commit_parameters\n");
    return 0;
}

int recommit_parameters () {
    //printf("falcon::recommit_parameters\n");
    return 0;
}

int commit_particles () {
    //printf("falcon::commit_particles\n");
    compute_force();
    return 0;
}

int recommit_particles () {
    //printf("falcon::recommit_particles\n");
    compute_force();
    return 0;
}

int synchronize_model ()
{
    //printf("falcon::synchronize_model\n");
    return 0;
}

int get_total_mass (double *total_mass)
{
    *total_mass = bodies.TotalMass(falcON::bodytype::std);
    return 0;
}

int get_center_of_mass_position (double *com_x, double *com_y, double *com_z)
{
    *com_x = *com_y = *com_z = 0;
    double invmass = 1./bodies.TotalMass(falcON::bodytype::std);
    for(falcON::body b=bodies.begin_all_bodies(); b; ++b) {
        *com_x += invmass * b.mass() * b.pos()[0];
        *com_y += invmass * b.mass() * b.pos()[1];
        *com_z += invmass * b.mass() * b.pos()[2];
    }
    return 0;
}

int get_center_of_mass_velocity (double *com_vx, double *com_vy, double *com_vz)
{
    *com_vx = *com_vy = *com_vz = 0;
    double invmass = 1./bodies.TotalMass(falcON::bodytype::std);
    for(falcON::body b=bodies.begin_all_bodies(); b; ++b) {
        *com_vx += invmass * b.mass() * b.vel()[0];
        *com_vy += invmass * b.mass() * b.vel()[1];
        *com_vz += invmass * b.mass() * b.vel()[2];
    }
    return 0;
}

int get_total_radius (double *radius) { return -2; }

int get_kinetic_energy (double *kinetic_energy)
{
    *kinetic_energy = 0;
    for(falcON::body b=bodies.begin_all_bodies(); b; ++b)
        *kinetic_energy += 0.5 * b.mass() * (b.vel()[0]*b.vel()[0] + b.vel()[1]*b.vel()[1] + b.vel()[2]*b.vel()[2]);
    return 0;
}

int get_potential_energy (double *potential_energy)
{
    *potential_energy = 0;
    for(falcON::body b=bodies.begin_all_bodies(); b; ++b)
        *potential_energy += 0.5 * b.mass() * b.pot();
    return 0;
}
