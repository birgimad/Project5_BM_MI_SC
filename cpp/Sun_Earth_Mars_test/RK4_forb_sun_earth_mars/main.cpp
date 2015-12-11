#include "include/armadillo"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
using namespace std;
using namespace arma;

// generating mass with gaussian distribution arounf 10 with standard deviation 1
void gaussian_mass_generator(vec (&mass), int number_of_particles)
{
  srand(time(NULL));
  for (int i = 0; i < number_of_particles; i++)
  {
  static int iset = 0;
  static double gset;
  double fac, rsq, v1, v2;
    do{
      v1 = 2.*((double) rand() / (RAND_MAX)) -1.0;
      v2 = 2.*((double) rand() / (RAND_MAX)) -1.0;
      rsq = v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.);
    fac = sqrt(-2.*log(rsq)/rsq);
    gset = v1*fac;
    iset = 1;
    mass(i) = v2*fac;
    mass(i) += 10;
  }
} // end function for gaussian deviates

void uniform_pos_generator(mat (&position), int N)
{
double pi=3.14159;
double c = 2*pi;
double R = 20;
vec phi(N);
vec r(N);
vec theta(N);
vec x(N);
vec y(N);
vec v(N);
srand(time(NULL));

for (int i=0;i<N;i++){

        x(i) = ((double) rand() / (RAND_MAX)); //random numbers generated in the interval(0,1)
        y(i) = ((double) rand() / (RAND_MAX));
        v(i) = ((double) rand() / (RAND_MAX));

    }
for (int i=0;i<N;i++){
   phi(i)=c*x(i);
   r(i)=R*pow(y(i),1.0/3.0);
   theta(i)=acos(1.0-2.0*v(i));
   position(i,0)=r(i)*sin(theta(i))*cos(phi(i));
   position(i,1)=r(i)*sin(theta(i))*sin(phi(i));
   position(i,2)= r(i)*cos(theta(i));
    }
}

void Derivative(mat r, mat v, vec m, mat (&drdt), mat (&dvdt), double G, int number_of_particles){
    double acc_x = 0, acc_y = 0, acc_z = 0;
    double distance_squared = 0;
    double force_between_particles = 0;
    for (int i=0; i<number_of_particles; i++)
    {
        force_between_particles = 0;
        acc_x = 0;
        acc_y = 0;
        acc_z = 0;
        for (int j=0; j<number_of_particles; j++)
        {
            if (j!=i)
            {
                distance_squared = 0;
                for (int k=0; k<3; k++)
                {
                    distance_squared += (r(i,k)-r(j,k))*(r(i,k)-r(j,k));
                }
                force_between_particles = m(j)*pow(distance_squared,-1.5);
                acc_x += G*force_between_particles*(r(j,0)-r(i,0));
                acc_y += G*force_between_particles*(r(j,1)-r(i,1));
                acc_z += G*force_between_particles*(r(j,2)-r(i,2));
            }
        }
        dvdt(i,0) = acc_x;
        dvdt(i,1) = acc_y;
        dvdt(i,2) = acc_z;
    }

    for (int i = 0; i<number_of_particles; i++)
    {
        drdt(i,0) = v(i,0);
        drdt(i,1) = v(i,1);
        drdt(i,2) = v(i,2);
    }
}

void making_ks(mat &pos_temp, mat &vel_temp, mat r, mat v, mat drdt, mat dvdt, double time_step, mat &k, double number, int number_of_particles){
    for (int i = 0; i< number_of_particles; i++)
    {
        for (int j=0; j<3; j++)
        {
            k(i,j) = time_step*drdt(i,j);
            k(i,j+3) = time_step*dvdt(i,j);
            pos_temp(i,j) = r(i,j) + k(i,j)/number;
            vel_temp(i,j) = v(i,j) + k(i,j+3)/number;
        }
    }
}

int main()
{
    double G = 2.96e-4; //1.563e-13; //  986.96/3; //2.96e-4;     //Grav const in the units of AU^3 / ( days^3 * mass_sun )
    double t0 = 0.0;
    double t_final = 10*365.0;     //time in the unit of days
    int number_of_time_step = 10*365;
    double time = t0;
    double dt = (t_final - t0)/number_of_time_step;
    int number_of_particles = 3;
    mat r(number_of_particles,3);
    mat v(number_of_particles,3);
    mat r_initial(number_of_particles,3);
    mat v_initial(number_of_particles,3);
    r.zeros();
    v.zeros();
    //position of Sun
    r(0,0) = 1.0;
    r(0,1) = 1.0;
    r(0,2) = 1.0;
    //position of Earth
    r(1,0) = 2.0;
    r(1,1) = 1.0;
    r(1,2) = 1.0;
    //position of Mars
    r(2,0) = -0.5;
    r(2,1) = 1.0;
    r(2,2) = 1.0;

    //velocity of Sun
    v(0,0) = 0.0;
    v(0,1) = 0.0;
    v(0,2) = 0.0;
    //velocity of Earth
    v(1,0) = 0.0;
    v(1,1) = 0.017;
    v(1,2) = 0.0;
    //velocity of Mars
    v(2,0) = 0.0;
    v(2,1) = 0.014;
    v(2,2) = 0.0;

    vec m(number_of_particles);
    m(0) = 1.0;
    m(1) = 3.0e-6;
    m(2) = 3.2e-7;
    r_initial = r;
    v_initial = v;

    //myfile << "Sun_x" << setw(20) << "Sun_y" << setw(20) << "Sun_z" << setw(20) << "Earth_x" << setw(20) << "Earth_y" << setw(20) << "Earth_z" << setw(20) << "Mars_x" << setw(20) << "Mars_y" << setw(20) << "Mars_z" << endl;
    //myfile << r(0,0) << setw(20) << r(0,1) << setw(20) << r(0,2) << setw(20) << r(1,0) << setw(20) << r(1,1) << setw(20) << r(1,2) << setw(20) << r(2,0) << setw(20) << r(2,1) << setw(20) << r(2,2) << endl;


    mat distance_between_particles(number_of_particles,number_of_particles);
    distance_between_particles.zeros();
    mat distance_between_particles_temp(number_of_particles,number_of_particles);
    distance_between_particles_temp.zeros();


    for (int i = 0; i<number_of_particles; i++)
    {
        for (int j = 0; j<number_of_particles; j++)
        {
            for (int k = 0; k<3; k++)
            {
                distance_between_particles_temp(i,j) += (r(i,k)-r(j,k))*(r(i,k)-r(j,k));
            }
            distance_between_particles(i,j) = pow(distance_between_particles_temp(i,j),0.5);
        }

    }
    cout << "distance between earth and sun: " << distance_between_particles(0,1) << " AU" << endl;
    cout << "distance between mars and sun: " << distance_between_particles(0,2) << " AU" << endl;

    vec kin_en(number_of_particles);
    kin_en.zeros();
    vec pot_en(number_of_particles);
    pot_en.zeros();
    vec tot_en(number_of_particles);
    tot_en.zeros();
    for (int i=0; i<number_of_particles; i++)
    {
        for (int k=0; k<3; k++)
        {
            kin_en(i) += v(i,k)*v(i,k);
        }
        kin_en(i) = 0.5*m(i)*kin_en(i);
        for (int j=0; j<number_of_particles; j++)
        {
            if (j != i)
            {
                pot_en(i) += m(j)/(distance_between_particles(i,j));
            }
        }
        pot_en(i) = pot_en(i)*G*m(i);
        tot_en(i) = kin_en(i)+pot_en(i);
    }
    double tot_en_all = 0;
    for (int j=0;j<number_of_particles;j++)
    {
        tot_en_all += tot_en(j);
    }
    cout << "in en: " << tot_en_all << endl;

    mat pos_temp(number_of_particles,3);
    mat vel_temp(number_of_particles,3);
    mat k1(number_of_particles,7);
    mat k2(number_of_particles,7);
    mat k3(number_of_particles,7);
    mat k4(number_of_particles,7);
    pos_temp.zeros();
    vel_temp.zeros();
    k1.zeros();
    k2.zeros();
    k3.zeros();
    k4.zeros();

    mat drdt(number_of_particles,3);
    drdt.zeros();
    mat dvdt(number_of_particles,3);
    dvdt.zeros();

    while(time<=t_final){
    Derivative(r,v,m,drdt,dvdt,G,number_of_particles);
    making_ks(pos_temp,vel_temp,r,v,drdt,dvdt,dt,k1,2.0,number_of_particles);
    Derivative(pos_temp,vel_temp,m,drdt,dvdt,G,number_of_particles);
    making_ks(pos_temp,vel_temp,r,v,drdt,dvdt,dt,k2,2.0,number_of_particles);
    Derivative(pos_temp,vel_temp,m,drdt,dvdt,G,number_of_particles);
    making_ks(pos_temp,vel_temp,r,v,drdt,dvdt,dt,k3,1,number_of_particles);
    Derivative(pos_temp,vel_temp,m,drdt,dvdt,G,number_of_particles);
    for (int j=0; j<number_of_particles; j++)
    {
        for (int i = 0; i<3; i++){
            k4(j,i) = dt*drdt(j,i);
            k4(j,i+3) = dt*dvdt(j,i);
        }
        for (int i=0; i<3; i++)
        {
            r(j,i) += (1.0/6.0)*(k1(j,i)+2*k2(j,i)+2*k3(j,i)+k4(j,i));
            v(j,i) += (1.0/6.0)*(k1(j,i+3)+2*k2(j,i+3)+2*k3(j,i+3)+k4(j,i+3));
        }
    }

    time += dt;
    }

    mat distance_between_particles_final(number_of_particles,number_of_particles);
    distance_between_particles_final.zeros();
    mat distance_between_particles_final_temp(number_of_particles,number_of_particles);
    distance_between_particles_final_temp.zeros();

    for (int i = 0; i<number_of_particles; i++)
    {
        for (int j = 0; j<number_of_particles; j++)
        {
            for (int k = 0; k<3; k++)
            {
                distance_between_particles_final_temp(i,j) += (r(i,k)-r(j,k))*(r(i,k)-r(j,k));
            }
            distance_between_particles_final(i,j) = pow(distance_between_particles_final_temp(i,j),0.5);
        }

    }

    vec kin_en_final(number_of_particles);
    kin_en_final.zeros();
    vec pot_en_final(number_of_particles);
    pot_en_final.zeros();
    vec tot_en_final(number_of_particles);
    tot_en_final.zeros();
    for (int i=0; i<number_of_particles; i++)
    {
        for (int k=0; k<3; k++)
        {
            kin_en_final(i) += v(i,k)*v(i,k);
        }
        kin_en_final(i) = 0.5*m(i)*kin_en_final(i);
        for (int j=0; j<number_of_particles; j++)
        {
            if (j != i)
            {
                pot_en_final(i) += m(j)/distance_between_particles_final(i,j);
            }
        }
        pot_en_final(i) = pot_en_final(i)*G*m(i);
        tot_en_final(i) = kin_en_final(i)+pot_en_final(i);
    }

    double tot_en_all_final = 0.0, kin_en_tot_final = 0.0, pot_en_tot_final = 0.0;
    for (int j=0;j<number_of_particles;j++)
    {
        tot_en_all_final += tot_en_final(j);
        kin_en_tot_final += kin_en_final(j);
        pot_en_tot_final += pot_en_final(j);
    }
    cout << "fin en: " << tot_en_all_final << endl;
    cout << "kin en: " << kin_en_tot_final << endl;
    cout << "pot en: " << pot_en_tot_final << endl;

    return 0;
}



