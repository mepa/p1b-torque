#include <iostream>
#include <fstream>
//#include <sstream>
#include <cmath>
#include <vector>
#include <map>
#include <string>
#include <strstream>

//------------------------------
// Constants
const size_t nbody = 1;
const double G = 1.0;
const double two_pi = 2.0 * M_PI;

//------------------------------
// Number of grid points
const size_t num_pts = 28; // <-- CHANGE
const size_t num_steps = num_pts - 1;

//------------------------------
// Cluster parameters:

const double M_cl = 1.0e-4; // <-- CHANGE

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

const double kappa = 1.8 * Omega0; // <-- CHANGE

const double nu = 2.0 * kappa; // <-- LEAVE AS 2\kappa
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

const double R_min         = - 3.0 * sigma_R / kappa;
const double R_max         =   3.0 * sigma_R / kappa;
const double z_min         = - 3.0 * sigma_z / nu;
const double z_max         =   3.0 * sigma_z / nu;
const double v_R_min       = - 3.0 * sigma_R;
const double v_R_max       =   3.0 * sigma_R;
const double v_R_theta_min = - 3.0 * R_cl * sigma_theta;
const double v_R_theta_max =   3.0 * R_cl * sigma_theta;  
const double v_z_min       = - 3.0 * sigma_z;
const double v_z_max       =   3.0 * sigma_z;

const double L_R           = R_max - R_min;
const double L_z           = z_max - z_min;
const double L_v_R         = v_R_max - v_R_min;
const double L_v_R_theta   = v_R_theta_max - v_R_theta_min;
const double L_v_z         = v_z_max - v_z_min;

const double grid_volume   = L_R * L_z * L_v_R * L_v_R_theta * L_v_z;
const double grid_density  = num_pts * num_pts * num_pts * num_pts * num_pts / grid_volume;
//const double grid_density = num_pts / grid_volume;

//------------------------------
// Array size
const size_t num_cols = 11; // t, R, Rtheta, z, vR, vRtheta, vRz, Lx, Ly, Lz, EJ
//const size_t num_cols = 10; // t, R, Rtheta, z, vR, vRtheta, vRz, Lx, Ly, Lz
const size_t num_rows = num_pts * num_pts * num_pts * num_pts * num_pts;


using namespace std;


//------------------------------
// Functions
double DF( const double v_R, const double v_theta, const double v_z, const double z )
{
  //const double factor = nu * Sigma / (4 * M_PI * M_PI * sigma_R * sigma_theta * sigma_z * sigma_z);
  const double factor = nu / (4 * M_PI * M_PI * sigma_R * sigma_theta * sigma_z * sigma_z); // Assume Sigma is constant and gets divided out.

  const double DF_v_R     = exp( -0.5 * v_R * v_R / (sigma_R * sigma_R));
  const double DF_v_theta = exp( -0.5 * v_theta * v_theta / (sigma_theta * sigma_theta));
  const double DF_v_z     = exp( -0.5 * v_z * v_z / (sigma_z * sigma_z));
  const double DF_z       = exp( -0.5 * nu_sq * z * z / (sigma_z * sigma_z));

  return factor * DF_v_R * DF_v_theta * DF_v_z * DF_z;
}

istream & read_row(istream & in, vector<double> & row)
{
  string line;
  do
    { 
      getline(in, line);
    }
  while(line.find('#') != string::npos && in);

  if(in)
    {
      istrstream in_line(line.c_str());
      
      for(size_t i = 0; i < num_cols; i ++) in_line >> row[i];
    }

  return in;
}

//------------------------------
// Main code
int main(int argc, char* argv[])
{
  if (argc != 2) cerr << "usage: torque input_file" << endl;

  ifstream f_in(argv[1]);
  if (!f_in)
    {
      cerr << "cannot open input file" <<  argv[1] << endl;
      exit(1);
    }

  cout.setf(ios::scientific, ios::floatfield);
  cout.precision(16);

//   size_t N;
//   double dR, dz, dv_R, dv_R_theta, dv_z;
//   double dx, dvx, dvy, dvz;

//   string s;
//   getline( f_in, s ); //read first line of file into string s

//   istrstream first_line( s.c_str() );
//   first_line >> N >> dx >> dz >> dvx >> dvy >> dvz;

//   cout << "#" << N << endl;
//   cout << "#" << dx << endl;
//   cout << "#" << dz << endl;
//   cout << "#" << dvx << endl;
//   cout << "#" << dvy << endl;
//   cout << "#" << dvz << endl;

//   cout << "#" << dz * dvx * dvy * dvz << " " << dx << endl;

  map<double, double> dLdt_x, dLdt_y, dLdt_z; // maps x coordinate onto the sum dLdt
  double dLdt_x_total = 0.0;
  double dLdt_y_total = 0.0;
  double dLdt_z_total = 0.0;
  //double Sigma_total = 0.0;

//   for (size_t iblock = 0; iblock < num_pts; iblock++)
//     {
      vector<double> row(num_cols);
      vector< vector<double> > all_rows;
      //vector< vector<double> > block_rows;

      for (size_t i = 0; i < (num_rows); i++)
	{
	  read_row(f_in, row);
	  all_rows.push_back(row);
	}

//       while (read_row(f_in, row))
// 	{
// 	  all_rows.push_back(row);
// 	}

      // Integral
      double sum_x = 0.0;
      double sum_y = 0.0;
      double sum_z = 0.0;
      double f = 0.0;
      
      for (size_t i = 0; i < (num_rows); i ++)
	//for (size_t i = 0; i < (num_rows / num_pts); i ++)
	{
	  const double t           = all_rows[i][0];
	  const double R           = all_rows[i][1];
	  const double z           = all_rows[i][3];
	  const double v_R         = all_rows[i][4];
	  const double v_R_theta   = all_rows[i][5];
	  const double v_z         = all_rows[i][6];
	  const double Lx_exchange = all_rows[i][7];
	  const double Ly_exchange = all_rows[i][8];
	  const double Lz_exchange = all_rows[i][9];
	  
	  //cout << t << " " << R << " " << z << " " << v_R << " " << v_R_theta << " " << v_z << " " << Lx_exchange << " " << Ly_exchange << " " << Lz_exchange << endl;
	  
	  const double theta = 0.0;
	  const double v_theta = (v_R_theta - v_R * theta); // / R; (v_theta \equiv R\dot{\theta})
	  const double vy = v_R * sin(theta) + R * v_theta * cos(theta);
	  
	  const double f = DF( v_R, v_theta, v_z, z);
	  
	  if(dLdt_x.find(R) == dLdt_x.end()) dLdt_x[R] = 0.0;
	  if(dLdt_y.find(R) == dLdt_y.end()) dLdt_y[R] = 0.0;
	  if(dLdt_z.find(R) == dLdt_z.end()) dLdt_z[R] = 0.0;
	  
	  // Calculate step size      
	  const double step_size_R         = (R_max - R_min) / num_steps;
	  const double step_size_v_R_theta = (v_R_theta_max - v_R_theta_min) / num_steps;
	  const double step_size_v_R       = (v_R_max - v_R_min) / num_steps;
	  const double step_size_z         = (z_max - z_min) / num_steps;
	  const double step_size_v_z       = (v_z_max - v_z_min) / num_steps;

	  const double dx = step_size_R;
	  const double dz = step_size_z;
	  const double dvx = step_size_v_R;
	  const double dvy = step_size_v_R_theta;
	  const double dvz = step_size_v_z;
      
	  dLdt_x[R] += f * vy * Lx_exchange * dz * dvx * dvy * dvz; // integral over everything but R
	  dLdt_y[R] += f * vy * Ly_exchange * dz * dvx * dvy * dvz;
	  dLdt_z[R] += f * vy * Lz_exchange * dz * dvx * dvy * dvz; 
	  //dLdt_z[R] += f * (Lz_exchange / t) * dz * dvx * dvy * dvz; 
	  //dLdt_z[R] += Lz_exchange / t;

	  dLdt_x_total += f * vy * Lx_exchange * dx * dz * dvx * dvy * dvz;
	  dLdt_y_total += f * vy * Ly_exchange * dx * dz * dvx * dvy * dvz;
	  dLdt_z_total += f * vy * Lz_exchange * dx * dz * dvx * dvy * dvz;

	  //Sigma_total += nu * nu * abs( z ) / (4 * M_PI * G); 
	}
      
      map<double, double>::const_iterator j = dLdt_y.begin();
      map<double, double>::const_iterator k = dLdt_z.begin();
      
      //const double sigma_sq = sigma_R * sigma_R + sigma_theta * sigma_theta + sigma_z * sigma_z;
      //const double type_I_torque = G * G * M_cl * M_cl * Sigma_total / (sigma_R * sigma_R);
      const double type_I_torque = G * G * M_cl * M_cl / (sigma_R * sigma_R); //Assume Sigma = const, Just use sigma_R^2 for dimensionless units
      
      for (map<double, double>::const_iterator i = dLdt_x.begin(); i != dLdt_x.end(); i ++)
	{
	  const double current_x = i->first;
	  
	  const double current_dLdt_x = i->second;
	  const double current_dLdt_y = j->second;
	  const double current_dLdt_z = k->second;
	  
	  cout << current_x << " " << current_dLdt_x / type_I_torque << " " << current_dLdt_y / type_I_torque << " " << current_dLdt_z / type_I_torque << endl; 
	  //<< current_dLdt_z / (num_pts * num_pts * num_pts * num_pts) << endl;
	  
	  j ++; k ++;
	}

      //}

      cout << "# " << grid_density << " " << dLdt_x_total / type_I_torque << " " << dLdt_y_total / type_I_torque << " " << dLdt_z_total / type_I_torque << endl;
      cout << "# " << grid_volume << endl;
      cout << "# " << type_I_torque << endl;
      return 0;
}
