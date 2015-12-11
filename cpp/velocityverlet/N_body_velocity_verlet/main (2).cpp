#include <armadillo>
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


int main()
{
    int number_of_particles = 100;
    int n = 3;
    double G = 1.563e-13;     //Grav const in the units of lightyear^3 / ( year^3 * mass_sun )
    double t0 = 0.0;
    double t_final = 1e7;     //time in the unit of years
    int number_of_time_step = 1e5;
    double time = t0;
    double dt = (t_final - t0)/number_of_time_step;
    mat r(number_of_particles,3);
    r.zeros();
    mat v(number_of_particles,3);
    v.zeros();
    mat v_partly(number_of_particles,3);
    v_partly.zeros();
    mat r_initial(number_of_particles,3);
    mat v_initial(number_of_particles,3);
    vec m(number_of_particles);
    m.zeros();
    gaussian_mass_generator(m,number_of_particles);
    uniform_pos_generator(r,number_of_particles);
    r_initial = r;
    v_initial = v;
    mat drdt(number_of_particles,3);
    drdt.zeros();
    mat dvdt(number_of_particles,3);
    dvdt.zeros();
    while(time<=t_final){
         Derivative(r,v,m,drdt,dvdt,G,number_of_particles);

        for(int i=0; i<number_of_particles ; i++){
            for(int j=0; j<3 ; j++){
               r(i,j) = r(i,j)+dt*drdt(i,j)+0.5*dt*dt*dvdt(i,j);
               v_partly(i,j) = drdt(i,j)+0.5*dt*dvdt(i,j);
               dvdt(i,j) = v_partly(i,j);
            }
        }
        Derivative(r,v,m,drdt,dvdt,G,number_of_particles);
        for(int i=0; i<number_of_particles ; i++){
            for(int j=0; j<3 ; j++){
        v(i,j) = v_partly(i,j)+0.5*dt*dvdt(i,j);
            }
        }

    time += dt;
    }

    vec distance_to_center_initial(number_of_particles);
    vec distance_to_center_final(number_of_particles);
    for (int i = 0; i<number_of_particles; i++)
    {
        //cout << m(i) << endl;
        distance_to_center_initial(i) = pow(r_initial(i,0)*r_initial(i,0)+r_initial(i,1)*r_initial(i,1)+r_initial(i,2)*r_initial(i,2),0.5);
        distance_to_center_final(i) = pow(r(i,0)*r(i,0)+r(i,1)*r(i,1)+r(i,2)*r(i,2),0.5);
        cout << distance_to_center_initial(i) << setw(15) << distance_to_center_final(i) << endl;
    }

    ofstream myfile ("100particles_time_1e5.txt");
            if (myfile.is_open())
            {
                 myfile << "value of dt: "<< dt;
                 myfile << "initial distance" << setw(15) << "final distance" << endl;
                for (int i=0;i<number_of_particles; i++){
                 myfile << distance_to_center_initial(i) << setw(15) << distance_to_center_final(i) << endl;
                 }
            }

}
