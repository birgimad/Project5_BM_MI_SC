#include <iostream>
#include "include/armadillo"
using namespace std;
using namespace arma;

void Derivative(double (&r)[3], double (&v)[3], double (&drdt)[3], double (&dvdt)[3], double G, double mass){
    drdt[0] = v[0];
    drdt[1] = v[1];
    drdt[2] = v[2];

    double distance_squared = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
    double newtonian_force = -G*mass/pow(distance_squared,1.5);
    dvdt[0] = newtonian_force*r[0];
    dvdt[1] = newtonian_force*r[1];
    dvdt[2] = newtonian_force*r[2];
}

int main()
{
    int n=3;

    double G = 2.96e-4;    //Grav const in the units of AU^3 / ( days^3 * mass_sun )
    //double mass_of_earth = 3e-6;  //mass in the units of mass_sun
    double t0 = 0.0;
    double mass =  1.0; //mass of object in origo of coordinate system in the units of mass_sun
    double t_final = 182.0;     //time in unit of days
    int number_of_time_step = 100000;
    double time = t0;
    double dt = (t_final - t0)/number_of_time_step;
    double r[3], v[3], drdt[3], dvdt[3], v_partly[3];
    //initializing relative position
    r[0] = 1.0;
    r[1] = 0.0;
    r[2] = 0.0;
    //initializing relative velovity
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

    double energy_initial = v[0]*v[0]+v[1]*v[1]+v[2]*v[2]*0.5+G*mass*pow(r[0]*r[0]+r[1]*r[1]+r[2]*r[2],-0.5);
    cout << "Energy_initial = " << energy_initial << endl;

    //start Velocity-Verlet method
    while(time<=t_final){
    Derivative(r,v,drdt,dvdt,G,mass);

    for(int i=0; i<n ; i++){
    r[i] = r[i]+dt*drdt[i]+0.5*dt*dt*dvdt[i];
    v_partly[i] = drdt[i]+0.5*dt*dvdt[i];
    dvdt[i] = v_partly[i];
    }

    Derivative(r,v,drdt,dvdt,G,mass);

    for(int i=0; i<n ; i++){
    v[i] = v_partly[i]+0.5*dt*dvdt[i];
    }

    time += dt;
    }
    //end Velocity-Verlet Method

    cout << "final position:" << endl;
    for (int i = 0; i<n; i++)
    {
        cout << r[i] << endl;
    }
    cout << "final velocity:" << endl;
    for (int i = 0; i<n; i++)
    {
        cout << drdt[i] << endl;
    }

    double distance_sun_final = pow(r[0]*r[0]+r[1]*r[1]+r[2]*r[2],1/2);
    cout << "Distance to sun: " << distance_sun_final << " AU" << endl;

    double energy_final = v[0]*v[0]+v[1]*v[1]+v[2]*v[2]*0.5+G*mass*pow(r[0]*r[0]+r[1]*r[1]+r[2]*r[2],-0.5);
    cout << "Energy_final = " << energy_final<< endl;

/*
        ofstream myfile ("VelocityVerlet2body2D.txt");
            if (myfile.is_open())
            {
                myfile << "Velocity-Verlet Method, 2 body, 3D" << endl;
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

