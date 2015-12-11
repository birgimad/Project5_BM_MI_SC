#include <iostream>
#include "include/armadillo"
using namespace std;
using namespace arma;

void Derivative(double r[3], double v[3], double (&drdt)[3], double (&dvdt)[3], double G, double mass){
    drdt[0] = v[0];
    drdt[1] = v[1];
    drdt[2] = v[2];

    double distance_squared = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
    double newtonian_force = -G*mass/pow(distance_squared,1.5);
    dvdt[0] = newtonian_force*r[0];
    dvdt[1] = newtonian_force*r[1];
    dvdt[2] = newtonian_force*r[2];
}

void updating_dummies(double dt, double drdt[3], double dvdt[3], double (&r_dummy)[3], double (&v_dummy)[3], double number, double (&kr)[3], double (&kv)[3], double r[3], double v[3])
{
    for (int i = 0; i<3; i++){
        kr[i] = dt*drdt[i];
        kv[i] = dt*dvdt[i];
        r_dummy[i] = r[i] +kr[i]/number;
        v_dummy[i] = v[i] +kv[i]/number;
    }
}


int main()
{
    int n=3;

    double G = 2.96e-4;     //Grav const in the units of AU^3 / ( days^3 * mass_sun )
    double mass_earth = 3e-6;       //mass in the units of mass_sun
    double mass = 1;        //mass in the units of mass_sun
    double t0 = 0.0;
    double t_final = 365.0;     //time in the unit of days
    int number_of_time_step = 100000;
    double time = t0;
    double dt = (t_final - t0)/number_of_time_step;

    double r[3], v[3], drdt[3], dvdt[3], r_dummy[3], v_dummy[3];
    double k1r[3], k2r[3], k3r[3], k4r[3];
    double k1v[3], k2v[3], k3v[3], k4v[3];

    //initializing relative position
    r[0] = 1.0;
    r[1] = 0.0;
    r[2] = 0.0;
    //initializing relative velocity
    v[0] = 0;
    v[1] = 0.017;
    v[2] = 0;

    //creating vectors for print file:
    double r_initial[3], v_initial[3];
    for (int i=0; i<3; i++)
    {
        r_initial[i] = r[i];
        v_initial[i] = v[i];
    }

    cout << "initial position:" << endl;
    for (int i = 0; i<n; i++)
    {
        cout << r[i] << endl;
    }
    cout << "initial velocity:" << endl;
    for (int i = 0; i<n; i++)
    {
        cout << v[i] << endl;
    }

    double distance_sun_initial = pow(r[0]*r[0]+r[1]*r[1]+r[2]*r[2],1/2);
    cout << "Distance to sun: " << distance_sun_initial << " AU" << endl;

    // energy calculated in the units of J/mass_earth (orbiting object)
    double energy_initial = v[0]*v[0]+v[1]*v[1]+v[2]*v[2]*0.5+G*mass*pow(r[0]*r[0]+r[1]*r[1]+r[2]*r[2],-0.5);
    cout << "Energy_initial = " << energy_initial << endl;

    //RK4:
    while(time<=t_final){
    Derivative(r,v,drdt,dvdt,G,mass);
    updating_dummies(dt,drdt,dvdt,r_dummy,v_dummy,2,k1r,k1v,r,v);
    Derivative(r_dummy,v_dummy,drdt,dvdt,G,mass);
    updating_dummies(dt,drdt,dvdt,r_dummy,v_dummy,2,k2r,k2v,r,v);
    Derivative(r_dummy,v_dummy,drdt,dvdt,G,mass);
    updating_dummies(dt,drdt,dvdt,r_dummy,v_dummy,1,k3r,k3v,r,v);
    Derivative(r_dummy,v_dummy,drdt,dvdt,G,mass);
    for (int i = 0; i<n; i++){
        k4r[i] = dt*drdt[i];
        k4v[i] = dt*dvdt[i];
    }
    for (int i=0; i<n;i++){
        r[i] = r[i] +(1.0/6.0)*(k1r[i]+2*k2r[i]+2*k3r[i]+k4r[i]);
        v[i] = v[i] +(1.0/6.0)*(k1v[i]+2*k2v[i]+2*k3v[i]+k4v[i]);
    }
    time += dt;
    }

    cout << "final position:" << endl;
    for (int i = 0; i<n; i++)
    {
        cout << r[i] << endl;
    }
    cout << "final velocity:" << endl;
    for (int i = 0; i<n; i++)
    {
        cout << v[i] << endl;
    }

    double distance_sun_final = pow(r[0]*r[0]+r[1]*r[1]+r[2]*r[2],1/2);
    cout << "Distance to sun: " << distance_sun_final << " AU" << endl;

    double energy_final = v[0]*v[0]+v[1]*v[1]+v[2]*v[2]*0.5+G*mass*pow(r[0]*r[0]+r[1]*r[1]+r[2]*r[2],-0.5);
    cout << "Energy_final = " << energy_final << endl;

/*
    ofstream myfile ("RungeKutta4_2body3D.txt");
        if (myfile.is_open())
        {
            myfile << "Runge-Kutta Method, 2 body, 3D" << endl;
            myfile << "Time: " << t_final << " days" << endl;
            myfile << "Number of time steps: " << number_of_time_step << endl;
            myfile << "Time step: " << dt << " days" << endl;

            myfile << "initial position:" << endl;
            for (int i = 0; i<n; i++)
            {
                myfile << r_initial[i] << endl;
            }
            myfile << "initial velocity:" << endl;
            for (int i = 0; i<n; i++)
            {
                myfile << v_initial[i] << endl;
            }
            myfile << "Initial distance to sun: " << distance_sun_initial << " AU" << endl;
            myfile << "Initial energy:  " << energy_initial << endl;

            myfile << "final position:" << endl;
            for (int i = 0; i<n; i++)
            {
                myfile << r[i] << endl;
            }
            myfile << "initial velocity:" << endl;
            for (int i = 0; i<n; i++)
            {
                myfile << v[i] << endl;
            }
            myfile << "Final distance to sun: " << distance_sun_final << " AU" << endl;
            myfile << "Final energy:  " << energy_final << endl;
        }
*/
    return 0;
}



