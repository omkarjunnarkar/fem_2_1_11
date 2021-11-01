#include<iostream>
#include<Eigen/Dense>
#include<iomanip>
#include<fstream>

using namespace std;
using namespace Eigen;

MatrixXd assign_mat(int ele_index_cpp, int numele) {
	//ele_index is current element number, not defined here!! different for each ele !!
	MatrixXd a_element(2, numele + 1);
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < numele + 1; j++) { a_element(i, j) = 0; }
	}
	a_element(0, ele_index_cpp) = 1;
	a_element(1, ele_index_cpp +1) = 1;
	return a_element;
}

double absolute_value(double value) {
	if (value < 0) { return -1 * value; }
	else { return value; }
}

double signum_function(double val) {
	if (val < 0) { return -1; }
	else { return 1; }
}

int element_check(int numele) {
	if (numele == 0) {
		//cout << "Minimum number of elements should be 2, considering =2/n";
		numele = 2;
	}
	else
	{
		//cout << "Number of elements= " << numele << endl;
	}
	return numele;
}

MatrixXd B_matrix(double l) {
	MatrixXd B(1, 2);
	B(0, 0) = -1 / l;
	B(1, 0) = 1 / l;
	return B;
}

MatrixXd k_ele(double tangstiff,double area, double l ) {
	//tangstiff is Ct, not defined here!! different for each elem !!
	double k_factor = (tangstiff * area) / l;
	//cout << tangstiff << "," << area << "," << l;
	MatrixXd k_element(2, 2);
	k_element(0, 0) = k_factor;
	k_element(0, 1) = -k_factor;
	k_element(1, 0) = -k_factor;
	k_element(1, 1) = k_factor;
	return k_element;
}

MatrixXd internal_force_ele(double area, double l, double weight, double J, double sigma_element, MatrixXd B) {
	//sigma is not defined here, it is different for each element !!
	return weight * area * J * sigma_element * B.transpose();
}

MatrixXd external_force_global(double force_value, int force_node_index_cpp, int numele) {
	MatrixXd  ext_force_glob = MatrixXd::Zero(numele +1, 1);
	ext_force_glob(force_node_index_cpp, 0) = force_value;
	return ext_force_glob;
}

double strain(MatrixXd B, MatrixXd u) {
	MatrixXd pr = B * u;
	double val = pr(0, 0);
	return val;
}

double stress_e_ep(double E, double strain, double pl_strain) {
	return E * (strain - pl_strain);
}


void main() {
	
	ofstream mydisplacefile("displacement.csv");
	ofstream mystressfile("stress.csv");
	ofstream mystrainfile("strain.csv");
	int numele = 4;
	double f_total = 10000; //Total force (Newton)
	int steps = 1e4; //nummber of time steps taken for loading

	double weight = 2; //Weights for single point quadrature. 
	double vis = 0.5; //Viscosity
	double max_time = 1; //Total time (seconds)
	double limit = 240; //Yield Stress
	double E = 120000; //young's modulus
	double total_len = 150; //total length of body
	double cs_area = 20; //crosssectional area
	double Ct = E; //Tangent stiffness, initially equal to young's modulus
	double del_t = 1 / double(steps);
	double len_of_ele = total_len / numele;
	double J = len_of_ele / 2;

	MatrixXd u_ele = MatrixXd::Zero(2, 1);
	
	double del_u_red = 0.0;
	 
	double strain_ele = 0.0, stress_ele = 0.0;
	double check_stress = 0.0;
	MatrixXd strain_array_global(numele, 1); //Total Strain
	MatrixXd stress_array_global(numele, 1); //global stress all elements

	MatrixXd B = B_matrix(len_of_ele);
	MatrixXd f_int_ele = MatrixXd::Zero(2, 1);
	
	MatrixXd f_ext_global = MatrixXd::Zero(numele +1, 1);

	MatrixXd force_level_array(steps+1, 1);
	double factor = f_total / steps;
	force_level_array(0, 0) = 0;
	for (int i = 1; i < steps+1; i++) {
		force_level_array(i, 0) = force_level_array(i -1, 0) + factor;
	}

	/*PLASTIC REGION QUANTITIIES*/
	double eps_pl_k = 0.0, eps_pl_k_plus = 0.0;
	double h = 1.0 / (100.0*steps);
	double psuedo_stress = 0.0,psuedo_stress_next=0.0;
	double eps_pl_prev = 0.0;
	//MatrixXd eps_p_array= MatrixXd::Zero(numele, 1);


	cout << "All Values,Matrices defined, entering force level loop" << endl;
	//for (int level = 0; level < force_level_array.rows(); level++) {
	//	//cout << force_level_array(level, 0)<<" , ";
	//	f_ext_global = external_force_global(force_level_array(level, 0), 2, numele);
	//	cout << f_ext_global << endl;
	//	cout << "---------" << endl;
	//}
	//cout << force_level_array.rows()<<endl;


	for (int level = 1; level < force_level_array.rows(); level++)
	{
		cout << "Force Value: " << force_level_array(level, 0) << endl;
		MatrixXd G= MatrixXd::Zero(numele +1, 1);
		MatrixXd u_global_backup = MatrixXd::Zero(numele + 1, 1);
		f_ext_global = external_force_global(force_level_array(level, 0), 2, numele); //2=the cpp_node_index on which force is applied
		double s_prev = 0.0;
		MatrixXd eps_p_array = MatrixXd::Zero(numele, 1);
		//cout << "F_External Defined" << endl;
		//cout << f_ext_global << endl;
		
		for (int nr=0;nr<6;nr++) 
		{
			//cout << "Entered the NR ." << endl;
	
			MatrixXd u_global = MatrixXd::Zero(numele + 1, 1);
			MatrixXd K_Global = MatrixXd::Zero(numele + 1, numele + 1);
			MatrixXd f_int_global = MatrixXd::Zero(numele + 1, 1);

			/*if (force_level_array(level, 0) = 0.0) {
				u_ele = MatrixXd::Zero(2, 1);
				cout << "First step, f=0;thus u=0." << endl;
			}*/

			for (int k = 0; k < numele; k++) {
				//cout << "Entered Element loop." << endl;

				MatrixXd A = assign_mat(k, numele);
				//cout << "Assignment matrix defined for element number: " << k << endl;

				u_ele(0, 0) = u_global_backup(k, 0);
				u_ele(1, 0) = u_global_backup(k +1, 0);

				u_global = u_global + (A.transpose() * u_ele);
				//cout << "U_global computed. " << endl;

				strain_ele = strain(B, u_ele);
				eps_pl_k = eps_p_array(k, 0);
				check_stress = E * (strain_ele- eps_pl_k );
				//cout << "Check stress & strain computed. " << endl;

				if (check_stress <= limit) {
					strain_array_global(k, 0) = strain_ele;
					stress_array_global(k, 0) = check_stress;
					stress_ele = check_stress;
					Ct = E;
					//cout << "Elastic Region." << endl;
				}
				
				else {	
					cout << "Entered Plastic Region. " << endl;
					
					for (int c = 0; c < steps; c++) {
						int x = 0;
						s_prev = stress_ele; //backup of previous stress value
						psuedo_stress = E * (strain_ele - eps_pl_k);
						eps_pl_k_plus = eps_pl_k + h * vis * signum_function(psuedo_stress) * (absolute_value(psuedo_stress / limit) - 1);
						while (x < 1000) {
							//cout << "Inside inner while for x = " << x << endl;
							eps_pl_prev = eps_pl_k_plus;
							psuedo_stress_next = E * (strain_ele - eps_pl_k_plus);
							eps_pl_k_plus = eps_pl_k + h * vis * signum_function(psuedo_stress_next) * (absolute_value(psuedo_stress_next / limit) - 1);
							x++;
							if (pow(eps_pl_k_plus - eps_pl_prev, 2) < 1e-5) { break; }

						}
						eps_pl_k = eps_pl_k_plus;
					}
					strain_array_global(k, 0) = strain_ele;
					stress_array_global(k, 0) = E*(strain_ele-eps_pl_k_plus);
					stress_ele = stress_array_global(k, 0);
					Ct= (E / (1 + (E * del_t * vis) / limit));
					eps_p_array(k, 0) = eps_pl_k_plus;
					//eps_pl_k = eps_pl_k_plus;
				}

				f_int_ele = internal_force_ele(cs_area, len_of_ele, weight, J, stress_ele, B);
				f_int_global = f_int_global + (A.transpose() * f_int_ele);
				//cout << "Computed F_internal_global." << endl;
				MatrixXd K_Element = k_ele(Ct, cs_area, len_of_ele);
				K_Global = K_Global + (A.transpose() * K_Element * A);
				//cout << "Computed K_global." << endl;
				//cout << K_Global << endl;
			}
			//cout << "Exited from Element Loop." << endl;
			//cout << f_int_global << endl;
			G = f_int_global - f_ext_global;
			//cout << G << endl;
			MatrixXd G_red(numele-1, 1); //first node fixed thus omit
			MatrixXd u_red(numele-1, 1);//first node fixed thus omit
			MatrixXd K_Global_red(numele-1, numele-1);//first node fixed thus omit
			for (int z = 0; z < numele-1; z++)
				{
					G_red(z, 0) = G(z + 1, 0);
					u_red(z, 0) = u_global(z + 1, 0);
					for (int z2 = 0; z2 < numele-1; z2++) 
					{
						K_Global_red(z, z2) = K_Global(z + 1, z2 + 1);
					}
				}
			
			MatrixXd del_u_red = K_Global_red.inverse() * G_red;
			u_red = u_red - del_u_red;
			for (int z = 0; z < numele-1; z++) 
				{
					u_global(z + 1, 0) = u_red(z, 0);
				}

			u_global_backup = u_global;
			
		}
		mydisplacefile << u_global_backup(2, 0) << endl;
		mystressfile << stress_array_global(1, 0) << endl;
		mystrainfile << strain_array_global(1, 0) << endl;
	} 
	mydisplacefile.close();
	mystressfile.close();
	mystrainfile.close();
}
