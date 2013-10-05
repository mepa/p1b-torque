////////////////////////////////////////////////////////////////////////
// torque_parallel.cc: Code to calculate the torque on a star cluster 
//   from a stellar disk.
// 
//          G           = 1
//          M_{galaxy}  = 1 = {the mass enclosed within the cluster's fixed, circular orbit}
//          R_{cluster} = 1
////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <strstream>
#include <mpi.h>
#include <cstdlib>
#include <algorithm>

#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//------------------------------
// Constants
const size_t nbody = 1;
const double G = 1.0;
const double two_pi = 2.0 * M_PI;

//------------------------------
// Number of grid points
const size_t num_pts = 32;
char* error_outfile = "Error_0.1_N32.log" ;
const size_t num_steps = num_pts - 1;

//------------------------------
// Cluster parameters:

const double M_cl = pow( 10.0,-2.0 ); // <--- VARY

const double R_cl = 1.0; 
const double R_theta_cl = R_cl * M_PI; 
const double theta_cl = R_theta_cl / R_cl; 
const double z_cl = 0.0; 
const double G_M_cl = G * M_cl;

//------------------------------
// Host parameters
const double M_host = 1.0;
const double G_M_host = G * M_host;
const double V_circ = 1.0;
const double Omega0 = V_circ / R_cl;
const double Omega0_sq = Omega0 * Omega0;
const double two_Omega0 = 2.0 * Omega0;

const double kappa = 1.0 * Omega0; // <--- VARY

const double nu = 2.0 * kappa;
//const double nu = 2.0 * kappa; // Milky Way: nu0 = (70+-4) km s^(-1) kpc^(-1); kappa0 = (37+-3)  km s^(-1) kpc^(-1)
const double nu_sq = nu * nu;
const double A0 = Omega0  - kappa * kappa / (4.0 * Omega0);
const double B0 = A0 - Omega0;
//const double r_jacobi = pow( G_M_cl / (4.0 * Omega0 * A0), 1.0 / 3.0 );
//const double eps = 1.0e-3 * r_jacobi; // softening length
const double eps = 1.0e-5;

//------------------------------
// Grid parameters
const double sigma_R     = 0.18 * V_circ;                     // Milky Way: sigma_R = (38+-2) km/s;  V_circ = (220+-20) km/s
const double sigma_theta = 0.5 * (kappa / Omega0) * sigma_R;  //            sigma_theta = ( kappa * Omega0 / 2.0 ) * sigma_R (BT, p.686)
const double sigma_z     = 0.5 * sigma_R;                     //            sigma_z = (19+-2) km/s

const double R_min         = - 3.0 * sigma_R / Omega0;
const double R_max         =   3.0 * sigma_R / Omega0;
const double z_min         = - 3.0 * sigma_z / nu;
const double z_max         =   3.0 * sigma_z / nu;
const double v_R_min       = - 3.0 * sigma_R;
const double v_R_max       =   3.0 * sigma_R;
const double v_R_theta_min = - 3.0 * R_cl * sigma_theta;
const double v_R_theta_max =   3.0 * R_cl * sigma_theta;  
const double v_z_min       = - 3.0 * sigma_z;
const double v_z_max       =   3.0 * sigma_z;


//const double nu_host = 0.1 * nu;
//const double nu_host_sq = nu_host * nu_host;
//const double nu_disk = 0.9 * nu;
//const double nu_disk_sq = nu_disk * nu_disk;


using namespace std;


//------------------------------
struct Orbit // all information related to a single orbit
{
  double R_init, R_theta_init, z_init, v_R_init, v_R_theta_init, v_z_init;
  double Lx_exchange, Ly_exchange, Lz_exchange, delta_time;
  double E_J;

  //vector<double> R_coords, R_theta_coords, z_coords;
};

inline double sgn(double x) 
{
  if(x < 0.0) return -1.0; else if(x > 0.0) return 1.0; else return 0.0;
}

inline double sq(double x)
{
  return x * x;
}

//------------------------------
//EQUATIONS OF MOTION: w denotes (q,p) = (x,y,z,vx,vy,vz), dwdt denotes (vx,vy,vz,ax,ay,az)
int func( double t, const double w[], double dwdt[], void *params )
{
  const double & R         = w[0]; 
  const double & R_theta   = w[1];
  const double & z         = w[2];
  const double & v_R       = w[3];
  const double & v_R_theta = w[4];
  const double & v_z       = w[5];
  
  const double R_sq = R * R;
  const double R_cub = R_sq * R;

  const double theta = R_theta / R;
  const double R_dot = v_R;
  const double theta_dot = (v_R_theta - R_dot * theta) / R;

  const double cos_theta = cos(theta);
  const double sin_theta = sin(theta);

  const double xij = R * cos_theta - R_cl * cos(theta_cl); 
  const double yij = R * sin_theta - R_cl * sin(theta_cl); 
  const double zij = z - z_cl;
  const double rij = sqrt( xij * xij + yij * yij + zij * zij + eps * eps); // the distance between star and cluster, where eps is the softening length
                                                                     
  const double rij_3 = rij * rij * rij;

  const double x = R * cos_theta;
  const double y = R * sin_theta;

  const double vx = R_dot * cos_theta - R * theta_dot * sin_theta; 
  const double vy = R_dot * sin_theta + R * theta_dot * cos_theta;

  const double rot_x = Omega0_sq * x + two_Omega0 * vy;      // centrifugal and coriolis terms due to rotating reference frame
  const double rot_y = Omega0_sq * y - two_Omega0 * vx; 
  const double rot_z = 0.0;

  const double r_total = sqrt( R_sq + z * z ); // the distance to the star from the center of the galaxy 
  //const double r_total_3 = r_total * r_total * r_total;

//-------------------------------------------------- DEBUG --------------------------------------------------
//   const double minus_G_M_host_over_r_total_3 = - G_M_host / r_total_3;
//   const double gx_host = minus_G_M_host_over_r_total_3 * x; // debug 1
//   const double gy_host = minus_G_M_host_over_r_total_3 * y;
//   const double gz_host = minus_G_M_host_over_r_total_3 * z;

//   const double gx_host = (- 1.0 + 2.0 * (R - R_cl)) * cos_theta; // debug 2
//   const double gy_host = (- 1.0 + 2.0 * (R - R_cl)) * sin_theta;
//   const double gz_host = 0.0;

  //cout << "XXX" << " " << gx_host << " " <<  - 1.0 / (R * R) * cos_theta << " " << (- 1.0 + 2.0 * (R - R_cl)) * cos_theta / (R * R) << endl;
  //cout << "YYY" << " " << gy_host << " " <<  - 1.0 / (R * R) * sin_theta << " " << (- 1.0 + 2.0 * (R - R_cl)) * sin_theta / (R * R) << endl;
  //cout << t << " " << gx_host << " " <<  - 1.0 / (R * R) * cos_theta << " " << (- 1.0 + 2.0 * (R - R_cl)) * cos_theta << " " << gy_host << " " <<  - 1.0 / (R * R) * sin_theta << " " << (- 1.0 + 2.0 * (R - R_cl)) * sin_theta << endl;
  //cout << t << " " << gx_host - (- 1.0 + 2.0 * (R - R_cl)) * cos_theta << " " << gy_host - (- 1.0 + 2.0 * (R - R_cl)) * sin_theta << endl;
//-----------------------------------------------------------------------------------------------------------

//   const double gx_host = - ( Omega0_sq * R_cl + (kappa * kappa - 3.0 * Omega0_sq) * (r_total - R_cl) + nu_host_sq * z ) * x / r_total; 
//   const double gy_host = - ( Omega0_sq * R_cl + (kappa * kappa - 3.0 * Omega0_sq) * (r_total - R_cl) + nu_host_sq * z ) * y / r_total;
//   const double gz_host = - ( Omega0_sq * R_cl + (kappa * kappa - 3.0 * Omega0_sq) * (r_total - R_cl) + nu_host_sq * z ) * z / r_total;

//   const double gz_disk = - nu_disk_sq * z; //- nu_sq * z;

  const double gx_host = - ( Omega0_sq * R_cl + (kappa * kappa - 3.0 * Omega0_sq) * (r_total - R_cl) ) * x / r_total; 
  const double gy_host = - ( Omega0_sq * R_cl + (kappa * kappa - 3.0 * Omega0_sq) * (r_total - R_cl) ) * y / r_total;
  const double gz_host = - ( Omega0_sq * R_cl + (kappa * kappa - 3.0 * Omega0_sq) * (r_total - R_cl) ) * z / r_total;

//   const double gx_host = - kappa * kappa * (r_total - R_cl) * x / r_total; I think this is the correct expression, but need to try it.
//   const double gy_host = - kappa * kappa * (r_total - R_cl) * y / r_total;
//   const double gz_host = - kappa * kappa * (r_total - R_cl) * z / r_total;

  const double gz_disk = - nu_sq * z;

  const double minus_G_M_cl_over_rij_3 = - G_M_cl / rij_3;
  const double gx_cl = minus_G_M_cl_over_rij_3 * xij;                     // cluster's gravitational force if approximated as a point mass
  const double gy_cl = minus_G_M_cl_over_rij_3 * yij;
  const double gz_cl = minus_G_M_cl_over_rij_3 * zij;

  const double force_x = rot_x + gx_host + gx_cl;
  const double force_y = rot_y + gy_host + gy_cl;
  const double force_z = rot_z + gz_host + gz_cl + gz_disk;

  double & d_R        = dwdt[0];
  double & d_R_theta  = dwdt[1];
  double & d_z        = dwdt[2];
  double & dv_R       = dwdt[3];
  double & dv_R_theta = dwdt[4];
  double & dv_z       = dwdt[5];

  d_R = v_R;
  d_R_theta = v_R_theta;
  d_z = v_z;
  dv_R = (x * force_x + y * force_y) / R + sq(y * vx - x * vy) / R_cub;
  dv_R_theta = dv_R * theta + (x * force_y - y * force_x) / R;
  dv_z = force_z;

  return GSL_SUCCESS;
}


//------------------------------
// MAIN CODE
int main( int argc, char* argv[] )
{
  ofstream errorFile( error_outfile );
  cerr.rdbuf( errorFile.rdbuf() );
  
  int ierr = MPI_Init( &argc, &argv ); //<---------------------------------------------------START MPI
  
  int num_procs, rank;
 
  if (ierr != MPI_SUCCESS)
    {
      if (rank == 0) cout << "Error starting MPI program. Terminating." << endl;
      MPI_Abort( MPI_COMM_WORLD, ierr );
    }
  ierr = MPI_Comm_size( MPI_COMM_WORLD, &num_procs );
  ierr = MPI_Comm_rank( MPI_COMM_WORLD, &rank );

//   // Example code for writing separate error file for each processor.
//   ostrstream log_out_name;
//   log_out_name << "log_" << rank << ends;
//   ofstream log_out(log_out_name.str());
//   log_out_name.rdbuf()->freeze(0); // This is necessary if you have a loop and need to deallocate memory.

  cout.setf(ios::scientific, ios::floatfield);
  cout.precision(16);
  
  const double epsrel = 1.0e-10;
  const double epsabs = 1.0e-10;

  //---SETUP---
  const gsl_odeiv_step_type* integrator = gsl_odeiv_step_rkf45; // "Embedded Runge-Kutta-Fehlberg (4,5) method. This method is a good general-purpose integrator."

  gsl_odeiv_step * s = gsl_odeiv_step_alloc( integrator, 6*nbody );
  gsl_odeiv_control * c = gsl_odeiv_control_y_new( epsabs, epsrel ); // Keep local error on each step within an absolute error (1st arg) and relative error (2nd arg) w.r.t. derivatives of solution y_i(t).
  gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc( 6*nbody );

  //   cout << "# INTEGRATOR INFO:" << endl;
  //   cout << "#   Step method: " << gsl_odeiv_step_name( s ) << endl;
  //   cout << "#   Control method: " << gsl_odeiv_control_name( c ) << endl;

  gsl_odeiv_system sys = { func, 0, 6*nbody, 0 }; // Second arg means no Jacobian function, fourth arg means no parameters

  if (rank == 0)
    {
      //cout << "# Number of steps: num_steps = " << num_steps << endl;
      cout << "# Number of points: num_pts = " << num_pts << endl;
      cout << "#   kappa = " << kappa << endl;
      cout << "#   satellite mass = " << M_cl << endl;
    }  

  //---INITIAL CONDITIONS---
  double w[6*nbody];

  double & R         = w[0]; 
  double & R_theta   = w[1];
  double & z         = w[2];
  double & v_R       = w[3];
  double & v_R_theta = w[4];
  double & v_z       = w[5];

  // Generate deviations in initial conditions  
  const double step_size_R = (R_max - R_min) / num_steps;
  const double step_size_v_R_theta = (v_R_theta_max - v_R_theta_min) / num_steps;
  const double step_size_v_R = (v_R_max - v_R_min) / num_steps;
  const double step_size_z = (z_max - z_min) / num_steps;
  const double step_size_v_z = (v_z_max - v_z_min) / num_steps;

  //cout << step_size_R << " " << step_size_z << " " << step_size_v_R << " " << step_size_v_R_theta << " " << step_size_z << endl;

  vector<double> R_vector(num_pts), v_R_theta_vector(num_pts), v_R_vector(num_pts), z_vector(num_pts), v_z_vector(num_pts);

  for(size_t i = 0; i < num_pts; i ++)
    {
      R_vector[i]         = R_min         + double(i) * step_size_R;
      v_R_theta_vector[i] = v_R_theta_min + double(i) * step_size_v_R_theta;
      v_R_vector[i]       = v_R_min       + double(i) * step_size_v_R;
      z_vector[i]         = z_min         + double(i) * step_size_z;
      v_z_vector[i]       = v_z_min       + double(i) * step_size_v_z;

      if (rank == 0) cout << "# " << R_vector[i] << " " << z_vector[i] << " " << v_R_vector[i] << " " << v_R_theta_vector[i] << " " << v_z_vector[i] << endl;
    }


  //---INTEGRATION---

  vector<Orbit> data;
  Orbit results;

  size_t num_rows = num_pts * num_pts * num_pts * num_pts * num_pts;

  Orbit * all_orbits = new Orbit [num_rows / num_procs];
  
  if (num_pts % num_procs != 0)
    {
      cout << "num_pts % num_procs = " << num_pts % num_procs << " != 0. Terminating." << endl;
      abort();
    }

  //if (rank == 0) cerr << "# Begin orbits loop..." << endl; //DEBUG
  //  for (size_t ix = 0; ix < num_pts; ix++)
  for (size_t ix = rank * num_pts / num_procs; ix < (rank + 1) * num_pts / num_procs; ix++)
    {
      const double R_initial = R_cl + R_vector[ix]; //<--Assign R
      //if (rank == 0) cerr << "# " << R_cl << " " << R_vector[ix] << endl; //DEBUG

      const double inv_sqrt_R_initial = 1.0 / sqrt(R_initial);

      for (size_t iz = 0; iz < num_pts; iz++)
	{
	  const double z_initial = z_vector[iz]; //<--Assign z

	  for (size_t ivx = 0; ivx < num_pts; ivx++)
	    {
	      const double v_R_initial = v_R_vector[ivx]; //<--Assign v_R
	      
	      for (size_t ivy = 0; ivy < num_pts; ivy++)
		{
		  const double v_R_theta_initial = inv_sqrt_R_initial - R_initial + v_R_theta_vector[ivy]; //<--Assign v_R_theta
		  //if (rank ==0) cerr << "# " << inv_sqrt_R_initial << " " << R_initial << " " << v_R_theta_vector[ivy] << endl; //DEBUG

		  for (size_t ivz = 0; ivz < num_pts; ivz++)
		    {
		      //if (rank == 0) cout << "#" << flush;

		      const double v_z_initial = v_z_vector[ivz]; //<--Assign v_z
		  
		      // Assign remaining initial conditions  
		      const double R_theta_initial   = 0.0;
	      
		      // SET INITIAL CONDITIONS
		      R         = R_initial;
		      R_theta   = R_theta_initial;
		      z         = z_initial;
		      v_R       = v_R_initial;
		      v_R_theta = v_R_theta_initial;
		      v_z       = v_z_initial;

		      //if (ix > 1 || iz > 0 || ivx > 0 || ivy > 0 || ivz > 0) abort(); // debug
		      
		      if ((R < eps) && (rank == 0)) cerr << "# WARNING: R close to zero! (R, R_theta, z, v_R, v_R_theta, v_z) = (" << R << ", " << R_theta << ", " << z << ", " << v_R << ", " << v_R_theta << ", " << v_z << ")" << endl;

// 		      if ( (R == R_cl) && (z == 0.0) && (v_R == 0.0) && (v_R_theta == 0.0) && (v_z == 0.0) )
// 			{
// 			  //if (rank == 0 ) cout << "# (" << ix << " " << iz << " " << ivx << " " << ivy << " " << ivz << ") " << "Skip integration since R = z = v_R = v_R_theta = v_z = 0.0" << endl;
// 			  continue;
// 			}
		      
		      //if (rank == 0) cerr << "#\n# Initial conditions: " << R << " " << R_theta << " " << z << " " << v_R << " " << v_R_theta << " " << v_z << endl; //DEBUG
		      
		      double t = 0.0; 
		      double h = 1.0e-2; // initial step-size 
		      
		      // Use these to calculate initial angular momentum
		      const double R_i         = results.R_init         = R;
		      const double R_theta_i   = results.R_theta_init   = R_theta;
		      const double z_i         = results.z_init         = z;
		      const double v_R_i       = results.v_R_init       = v_R;
		      const double v_R_theta_i = results.v_R_theta_init = v_R_theta;
		      const double v_z_i       = results.v_z_init       = v_z;
		      
		      const double t_i         = t;
		      const double theta_i     = R_theta_i / R_i;     
		      const double theta_dot_i = (v_R_theta_i - v_R_i * theta_i) / R_i;

		      const double r_total_i   = sqrt( R * R + z * z );
		      const double cos_theta_i = cos(theta_i);
		      const double sin_theta_i = sin(theta_i);
		      const double xij_i = R_i * cos_theta_i - R_cl * cos(theta_cl); 
		      const double yij_i = R_i * sin_theta_i - R_cl * sin(theta_cl); 
		      const double zij_i = z_i - z_cl;
		      const double rij_i = sqrt( xij_i * xij_i + yij_i * yij_i + zij_i * zij_i + eps * eps); // the distance between star and cluster, where eps is the softening length
		      const double p_sq_i = v_R_i * v_R_i + R_i * R_i * (theta_dot_i  + Omega0) * (theta_dot_i + Omega0) + v_z_i * v_z_i;
		      const double pot_host_i = R_cl * Omega0_sq * (r_total_i - R_cl) + 0.5 * (kappa * kappa - 3.0 * Omega0_sq) * (r_total_i - R_cl) * (r_total_i - R_cl);
		      //const double pot_host_i = - G_M_host / r_total_i; // keplerian potential
		      const double potential_i = pot_host_i - G_M_cl / rij_i + 0.5 * nu_sq * z_i * z_i;
		      
		      double theta     = R_theta / R;     
		      double theta_dot = (v_R_theta - v_R * theta) / R;
		      //		      double Lz =  R * R * (theta_dot + Omega0);
		      
		      double theta_previous = theta;
		      double theta_dot_previous = theta_dot;
		      
		      // 		      ostrstream orbit_file_name;
		      // 		      orbit_file_name << "/data1/r900-3/mepa/Orbit_temp/orbit_" << ix << "_" << iz << "_" << ivx << "_" << ivy << "_" << ivz << ends; 
		      // 		      //orbit_file_name << "orbit_" << ix << "_" << iz << "_" << ivx << "_" << ivy << "_" << ivz << ends;
		      // 		      ofstream orbit_out(orbit_file_name.str()); 
	  
		      // 		      orbit_out << t << " " << R << " " << R_theta << " " << z << " " << v_R << " " << v_R_theta << " " << v_z << " " << Lz << endl;
	  
		      //bool met_the_cluster = false;
		      bool cross_zero = false;
		      bool cross_pos_two_pi = false;
		      bool cross_neg_two_pi = false;	

		      //if (rank == 0) cerr << "# Begin integration..." << endl; //DEBUG
		      while( true )
			{
			  const int status = gsl_odeiv_evolve_apply( e, c, s, &sys, &t, t + 2.0 * h, &h, w );
			  
			  theta = R_theta / R;
			  theta_dot = (v_R_theta - v_R * theta) / R;
			  // double Lz =  R * R * (theta_dot + Omega0);
			  
			  const double theta_current = theta;
			  const double theta_dot_current = theta_dot;
			  
			  if (status != GSL_SUCCESS) abort();
			  
			  cross_zero       = theta_previous * theta_current <= 0.0;
			  cross_pos_two_pi = (theta_previous - two_pi) * (theta_current - two_pi) <= 0.0;
			  cross_neg_two_pi = (theta_previous + two_pi) * (theta_current + two_pi) <= 0.0;

			  if ((theta_previous != 0.0) && (cross_zero || cross_pos_two_pi || cross_neg_two_pi)) break; 

// 			  if (t > 1400.0) //USE ONLY FOR KAPPA = 2 * OMEGA0
// 			    {
// 			      if (rank == 0) cerr << "# (" << ix << " " << iz << " " << ivx << " " << ivy << " " << ivz << ") " << "Integration too long!" << endl;
// 			      break;
// 			    }



// 			  const double met = 0.5 * M_PI;
// 			  met_the_cluster = met_the_cluster || ( (met < theta) && (theta < (two_pi - met)) )  ||  ( ( - met > theta) && (theta > - (two_pi - met)) );
// 			  if ( met_the_cluster && ( cross_zero || cross_pos_two_pi || cross_neg_two_pi ) ) break;

                          
// 			  if (rank ==0) //DEBUG
// 			    {
// 			      const double r_total = sqrt( R * R + z * z );
// 			      const double cos_theta = cos(theta);
// 			      const double sin_theta = sin(theta);
// 			      const double xij = R * cos_theta - R_cl * cos(theta_cl); 
// 			      const double yij = R * sin_theta - R_cl * sin(theta_cl); 
// 			      const double zij = z - z_cl;
// 			      const double rij = sqrt( xij * xij + yij * yij + zij * zij + eps * eps); 
// 			      const double p_sq = v_R * v_R + R * R * (theta_dot  + Omega0) * (theta_dot + Omega0) + v_z * v_z;
// 			      const double potential = - G_M_host / r_total - G_M_cl / rij + 0.5 * nu_sq * z * z;
// 			      const double Lz =  R * R * (theta_dot + Omega0);
// 			      const double E_J = 0.5 * p_sq + potential - Omega0 * Lz;
// 			      cout << t << " " << Lz << " " << E_J << endl;
// 			    }		
			  

			  theta_previous = theta_current;
			  theta_dot_previous = theta_dot_current;
			}
		      //if (rank == 0) cerr << "# ...end integration." << endl; //DEBUG

		      //if (t > 1400.0) continue;

		      // Use these to calculate final angular momentum and Jacobi integral
		      const double R_f         = R;
		      const double R_theta_f   = R_theta;
		      const double z_f         = z;
		      const double v_R_f       = v_R;
		      const double v_R_theta_f = v_R_theta;
		      const double v_z_f       = v_z;
		      
		      const double t_f         = t;
		      const double theta_f     = R_theta_f / R_f;     
		      const double theta_dot_f = (v_R_theta_f - v_R_f * theta_f) / R_f;

		      const double r_total_f = sqrt( R * R + z * z );
		      const double cos_theta_f = cos(theta_f);
		      const double sin_theta_f = sin(theta_f);
		      const double xij_f = R_f * cos_theta_f - R_cl * cos(theta_cl); 
		      const double yij_f = R_f * sin_theta_f - R_cl * sin(theta_cl); 
		      const double zij_f = z_f - z_cl;
		      const double rij_f = sqrt( xij_f * xij_f + yij_f * yij_f + zij_f * zij_f + eps * eps); // the distance between star and cluster, where eps is the softening length
		      const double p_sq_f = v_R_f * v_R_f + R_f * R_f * (theta_dot_f  + Omega0) * (theta_dot_f + Omega0) + v_z_f * v_z_f;
		      const double pot_host_f = R_cl * Omega0_sq * (r_total_f - R_cl) + 0.5 * (kappa * kappa - 3.0 * Omega0_sq) * (r_total_f - R_cl) * (r_total_f - R_cl);
		      //const double pot_host_f = - G_M_host / r_total_f; // keplerian potential
		      const double potential_f = pot_host_f - G_M_cl / rij_f + 0.5 * nu_sq * z_f * z_f;

		      // Calculate angular momentum
		      const double Lx_initial =  (R_i * v_z_i - v_R_i * z_i) * sin(theta_i) - R_i * z_i * (theta_dot_i + Omega0) * cos(theta_i);
		      const double Ly_initial = -(R_i * v_z_i - v_R_i * z_i) * cos(theta_i) - R_i * z_i * (theta_dot_i + Omega0) * sin(theta_i);
		      const double Lz_initial =  R_i * R_i * (theta_dot_i + Omega0);
		      
		      const double Lx_final =  (R_f * v_z_f - v_R_f * z_f) * sin(theta_f) - R_f * z_f * (theta_dot_f + Omega0) * cos(theta_f);
		      const double Ly_final = -(R_f * v_z_f - v_R_f * z_f) * cos(theta_f) - R_f * z_f * (theta_dot_f + Omega0) * sin(theta_f);
		      const double Lz_final =  R_f * R_f * (theta_dot_f + Omega0);
		      
		      results.Lx_exchange = Lx_final - Lx_initial;
		      results.Ly_exchange = Ly_final - Ly_initial;
		      results.Lz_exchange = Lz_final - Lz_initial;

		      // Calculate Jacobi integral
		      const double E_J_initial = 0.5 * p_sq_i + potential_i - Omega0 * Lz_initial;
		      const double E_J_final   = 0.5 * p_sq_f + potential_f - Omega0 * Lz_final;
		      if ( (E_J_final - E_J_initial) > (100.0 * epsabs) ) 
			{
			  if (rank==0) 
			    {
			      cerr << "# WARNING: Jacobi integral not conserved? (" << ix << " " << iz << " " << ivx << " " << ivy << " " << ivz << ") " << E_J_final - E_J_initial << endl; 

			      //break;
			    }
			}
		      results.E_J = E_J_final;

		      //if (rank == 0) cerr << E_J_final << " - " << E_J_initial << " = " << E_J_final - E_J_initial << endl; //CHECK CHANGE IN JACOBI INTERAL OVER EACH ORBIT
		      

		      results.delta_time = t_f - t_i;

// 		      if (results.delta_time < 1.0e-5) //DEBUG
// 			{
// 			  if (rank == 0) 
// 			    {
// 			      cerr << "# delta_time = 0? (" << ix << " " << iz << " " << ivx << " " << ivy << " " << ivz << ") " << results.delta_time << endl;
// 			      abort();
// 			    }
// 			}	

// 		      if (results.R_init < 0.1) //DEBUG
// 			{
// 			  if (rank == 0) 
// 			    {
// 			      cerr << "# R_init = 0? (" << ix << " " << iz << " " << ivx << " " << ivy << " " << ivz << ") " << results.R_init << endl;
// 			      abort();
// 			    }
// 			}		    
		      
      
		      data.push_back( results );
		    }
		}
	    }
	}
    }
  //if (rank == 0) cerr << "# ...end orbits loop." << endl; //DEBUG


  if (rank == 0) cerr << "# CHECK: Is " << data.size() << " = " << (num_rows / num_procs) << "?" << endl;


  //if (rank == 0) cerr << "# Begin data write-out loop..." << endl; //DEBUG
  for (int irank = 0; irank < num_procs; irank ++)
    {
      if(rank == irank) copy( data.begin(), data.end(), all_orbits );
      //if ( rank == 0) cerr << data.R_init << endl;

      double * all_orbits_double = reinterpret_cast<double *>(all_orbits);

      int all_orbits_double_count = num_rows * sizeof( Orbit ) / sizeof( double );

      //cerr << "rank " << rank << " doing gather " << endl; //DEBUG
      //cerr << "arguments of MPI_Bcast are: " << all_orbits_double << " " << all_orbits_double_count / num_procs << endl; //DEBUG
  
      int ierr_1;

      ierr_1 = MPI_Bcast( all_orbits_double, all_orbits_double_count / num_procs, MPI_DOUBLE, irank, MPI_COMM_WORLD );  // does unnecessary work

      //cerr << "rank " << rank << " done gather " << ierr_1 << endl; //DEBUG

      if (rank == 0)
	{
	  ostrstream block_out;
	  block_out.setf(ios::scientific, ios::floatfield);
	  block_out.precision(16);

	  for( size_t i = 0; i < num_rows / num_procs; i++)
	    {
	      Orbit orbit;
	      memcpy( &orbit, &all_orbits_double[i * sizeof(Orbit) / sizeof(double)], sizeof(Orbit) );

	      block_out << orbit.delta_time << " " << orbit.R_init << " " << orbit.R_theta_init << " " << orbit.z_init << " " << orbit.v_R_init << " " << orbit.v_R_theta_init << " " << orbit.v_z_init << " " << orbit.Lx_exchange << " " << orbit.Ly_exchange << " " << orbit.Lz_exchange << " " << orbit.E_J << '\n';

	      //cerr << "i = " << i << endl;

// 	      if ( orbit.R_init < 0.1 )
// 		{
// 		  cerr << i << " --> WARNING: orbit.delta_time = " << orbit.delta_time << " and orbit.R_init = " << orbit.R_init << endl;
// 		  //abort();
// 		}
	    }
	  
	  block_out << ends;

	  cout << block_out.str() << flush;

	  block_out.freeze(false);
	}
    }
  //if (rank == 0) cerr << "# ...end data write-out loop." << endl; //DEBUG 

  MPI_Finalize(); //<-------------------------------------------------END MPI
	
  return 0;
}
