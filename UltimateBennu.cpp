#include <iostream>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <string>
#include <iomanip>

using namespace std;


# define res 200
# define PI 3.14159265 
# define dump 50 //dumping interval, no. of time steps

double weight = 0.5; // used in 2nd  order RK
double xlim = 0.99; // x-dimensions of domain
double xstart = 1 - xlim; // x-dimensions of domain

double grav[res + 4]; // store magnitude of gravity
double eta[res + 4]; // angle made by gravity vector

double h[res + 4]; double htemp[res + 4]; // h
double u[res + 4]; double utemp[res + 4]; // u
double v[res + 4]; double vtemp[res + 4]; // v
double p[res + 4]; double ptemp[res + 4]; double ptemphat[res + 4]; // p
double q[res + 4]; double qtemp[res + 4]; double qtemphat[res + 4]; // q
double r[res + 4]; double rtemp[res + 4]; double rtemphat[res + 4]; // r

double eig[res + 4]; // max speed at every grid boundary

double f1, f2, f3, fp, fm; // fluxes
double Qp_p, Qm_p, Qp_q, Qm_q, Qp_r, Qm_r; // plus and minus values of q vector at current location
double vv, phi;
double x; // current downslope coordinates
double tau; // radial coordinate
double buff = 0;//0.0000001;
double alpha = 0; // angular acceleration

// time is money
int timesteps = 0;
double finalt = 10;

// parameters
double delta = 20 * PI / 180; // angle of friction in degrees
double theta = 1.0; // used in min-mod limiter
double Om = 1.5;
double dx = xlim / res;
double dt = dx / 4;
double Kx = 2 * pow(1 / cos(delta), 2) - 1;
double ke_curr = 0;
double ke_prev = 0;
double mass_total = 0;
double mom_tot = 0;
int flag, flagtemp;
double epsilon = 0.005; // shallowness parameter
double maxspeed;
double momincb, xi, zeta, Gn, Gt, rad1, height1;
double t, timep, dom, L, Ldot, hdot, vdot, mass_total_o = 0, integral1 = 0, integral2 = 0, integral1_o = 0, integral2_o = 0, h_x, v_x;
int j, jj, deep; // loop variables
double shmax, maxheight, factor;
double alpha_momin = 0, alpha_angmom = 0; // contributions of respective quantities to spin change
double Dc = 0.005; // filleting distance or the distance between the two cones
double dregolith = 1.0;
double slant_height = 470.0 / 376.0;
double grav_scale = 1.0;
double rho_bar = 1.2 / 0.75;


// ##########################################################
// function declarations

double pfunc(int);
double qfunc(int);
double rfunc(int);

void Gr();

double ax(int);
double maxeigenx(int, int);

double flux(double, double, double, int, int, int);

void def(int, int);
double Hx(int, int);

double derivative(int, int);

double minmod(double, double, double);
void pfun(int);
double S1(int);
double S2(int);
double S3(int);

// ##########################################################
// function main

int main()
{
	t = 0;j = 0;timep = 0;
	xi = 45 * PI / 180;
	zeta = PI / 2 - xi; // semi-apex angle
	Gr();

	// Initial conditions

	double sigsq = 0.01;
	double mumu = 0.2;
	for (j = 2; j <= res + 1; j++)
	{
		x = xstart + xlim * (j - 1) / res - xlim / 2 / res; tau = x * sin(zeta) + buff;
		if (j > int(2 * res / 3))
			h[j] = dregolith;// 0.1 * exp(-0.5 * pow((x - mumu), 2) / sigsq) / sqrt(sigsq * 2 * PI);//
		else
			h[j] = dregolith;// 0.1 * exp(-0.5 * pow((x - mumu), 2) / sigsq) / sqrt(sigsq * 2 * PI);//
		u[j] = 0;//0.5 * (sin(zeta) - tau);
		v[j] = 0;//+0.5*epsilon*h[j]/cos(xi));
		p[j] = h[j] * tau;
		q[j] = h[j] * u[j] * tau;
		r[j] = h[j] * v[j] * tau;
		mass_total = mass_total + 2 * PI * p[j] * dx;
	}
	ofstream myfileh;
	myfileh.open("hout1.txt");
	ofstream myfileu;
	myfileu.open("uout1.txt");
	ofstream myfilev;
	myfilev.open("vout1.txt");
	ofstream myfilet;
	myfilet.open("tout1.txt");
	ofstream myfiles;
	myfiles.open("sout1.txt");
	ofstream myfilem;
	myfilem.open("mout1.txt");
	ofstream myfilek;
	myfilek.open("kout1.txt");

	for (deep = 0; deep < 25; deep++)
	{
		t = 0;
		//grav_scale = 1 - 2 * epsilon * deep;

		if (deep != 0)
		{
			for (j = 2; j <= res + 1; j++)
			{
				x = xstart + xlim * (j - 1) / res - xlim / 2 / res; tau = x * sin(zeta) + buff;
				h[j] = dregolith + h[j];
				u[j] = 0;//0.5 * (sin(zeta) - tau);;
				v[j] = 0;//+0.5*epsilon*h[j]/cos(xi));
				p[j] = h[j] * tau;
				q[j] = h[j] * u[j] * tau;
				r[j] = h[j] * v[j] * tau;
			}
		}

		// updaing moment of inertia of the core

		rad1 = slant_height * cos(xi);
		height1 = slant_height * sin(xi);
		momincb = (rho_bar - 1) * PI * pow(241.8 / 376.0 * cos(PI / 4), 5) / 10 + PI * pow(rad1, 5) / 10;// landslides in both hemi-tops
		slant_height = slant_height - 2 * epsilon * dregolith;

		// time stepping loop	
		t = 0;
		dom = 0;
		while (t < finalt)
		{
			Om = Om + dom;

			// writing output file 2D
			if (timesteps % dump == 0)
			{
				maxheight = 0;
				shmax = 0;
				for (j = 2; j <= res + 1; j++)
				{
					x = xstart + xlim * (j - 1) / res - xlim / 2 / res; tau = x * sin(zeta) + buff;
					if (maxheight < h[j])
					{
						maxheight = h[j];
						shmax = xstart + xlim * (j - 1) / res - xlim / 2 / res;
					}
					myfileh << h[j] << " " << endl;
					myfileu << u[j] << " " << endl;
					myfilev << v[j] << " " << endl;
				}
				myfilem << mass_total << " " << endl;
				myfilek << mass_total << " " << endl;
				myfilet << Om << " " << endl;
				myfiles << timep << " " << endl;
			}

			// construction of conserved variables and their temporary storage
			for (j = 2; j <= res + 1; j++)
			{
				x = xstart + xlim * (j - 1) / res - xlim / 2 / res; tau = x * sin(zeta) + buff;

				utemp[j] = q[j] / p[j];
				vtemp[j] = r[j] / p[j];
				ptemp[j] = max(p[j], 1e-16 * tau);
				htemp[j] = ptemp[j] / tau;
				qtemp[j] = ptemp[j] * utemp[j];
				rtemp[j] = ptemp[j] * vtemp[j];
				/*ptemp[j] = max(p[j], pow(10, -12) * tau);
				qtemp[j] = p[j] / tau > pow(10, -12) ? q[j] : 0;
				rtemp[j] = p[j] / tau > pow(10, -12) ? r[j] : ptemp[j] * v[j];
				htemp[j] = p[j] / tau;
				utemp[j] = p[j] / tau > pow(10, -12) ? u[j] : 0;
				vtemp[j] = v[j];*/
				ptemphat[j] = ptemp[j];
				qtemphat[j] = qtemp[j];
				rtemphat[j] = rtemp[j];
			}

			//ghost cell updates 
			factor = 1.0;

			//ghost cell updates
		//	ptemp[0] = ptemp[3];qtemp[0] = -qtemp[3];rtemp[0] = rtemp[3];
		//	ptemp[1] = ptemp[2];qtemp[1] = -qtemp[2];rtemp[1] = rtemp[2];

		//	ptemp[res + 2] = ptemp[res + 1] ;qtemp[res + 2] = -qtemp[res + 1];rtemp[res + 2] = rtemp[res + 1];
		//	ptemp[res + 3] = ptemp[res];qtemp[res + 3] = -qtemp[res];rtemp[res + 3] = rtemp[res];
			ptemp[0] = ptemp[3];qtemp[0] = -qtemp[3];rtemp[0] = rtemp[3];htemp[0] = htemp[3];utemp[0] = -utemp[3];vtemp[0] = (1 + 2 * theta) * vtemp[2] - 2 * theta * vtemp[3];
			ptemp[1] = ptemp[2];qtemp[1] = -qtemp[2];rtemp[1] = rtemp[2]; htemp[1] = htemp[2];utemp[1] = -utemp[2];vtemp[1] = (1 + theta) * vtemp[2] - theta * vtemp[3];

			ptemp[res + 2] = ptemp[res + 1];qtemp[res + 2] = -qtemp[res + 1];rtemp[res + 2] = rtemp[res + 1];htemp[res + 2] = htemp[res + 1];utemp[res + 2] = -utemp[res + 1];vtemp[res + 2] = (1 + theta) * vtemp[res + 1] - theta * vtemp[res];
			ptemp[res + 3] = ptemp[res];qtemp[res + 3] = -qtemp[res];rtemp[res + 3] = rtemp[res];htemp[res + 3] = htemp[res];utemp[res + 3] = -utemp[res];vtemp[res + 3] = (1 + 2 * theta) * vtemp[res + 1] - 2 * theta * vtemp[res];

			integral1_o = integral2_o = mass_total_o = 0;
			for (j = 2; j <= res + 1; j++)
			{
				x = xstart + xlim * (j - 1) / res - xlim / 2 / res; tau = x * sin(zeta) + buff;
				mass_total_o = mass_total_o + 2 * PI * p[j] * dx;
				integral1_o = integral1_o + 2 * PI * (pow(x * sin(xi), 3) * epsilon * p[j] / tau) * dx;
				integral2_o = integral2_o + 2 * PI * (pow(x * sin(xi), 2) * epsilon * r[j] / tau) * dx;
			}

			// predictor step for interior
			for (j = 2; j <= res + 1; j++)
			{
				x = xstart + xlim * (j - 1) / res - xlim / 2 / res; tau = x * sin(zeta) + buff;
				p[j] = ptemp[j] - ((Hx(j + 1, 1) - Hx(j, 1)) / dx - S1(j)) * dt;
				q[j] = qtemp[j] - ((Hx(j + 1, 2) - Hx(j, 2)) / dx - S2(j)) * dt;
				r[j] = rtemp[j] - ((Hx(j + 1, 3) - Hx(j, 3)) / dx - S3(j)) * dt;

			}




			for (j = 2; j <= res + 1; j++)
			{
				x = xstart + xlim * (j - 1) / res - xlim / 2 / res; tau = x * sin(zeta) + buff;
				utemp[j] = q[j] / p[j];
				vtemp[j] = r[j] / p[j];
				ptemp[j] = max(p[j], 1e-16 * tau);
				htemp[j] = ptemp[j] / tau;
				qtemp[j] = ptemp[j] * utemp[j];
				rtemp[j] = ptemp[j] * vtemp[j];
				//	cout<<j<<"    "<<htemp[j]<<endl;
			}

			//wall boundary conditions at pole
			ptemp[0] = ptemp[3] * factor;qtemp[0] = -qtemp[3];rtemp[0] = rtemp[3];
			ptemp[1] = ptemp[2] * factor;qtemp[1] = -qtemp[2];rtemp[1] = rtemp[2];

			//wall boundary conditions at equator
			ptemp[res + 2] = ptemp[res + 1];qtemp[res + 2] = -qtemp[res + 1];rtemp[res + 2] = rtemp[res + 1];
			ptemp[res + 3] = ptemp[res];    qtemp[res + 3] = -qtemp[res];    rtemp[res + 3] = rtemp[res];

			// corrector step for interior
			for (j = 2; j <= res + 1; j++)
			{
				x = xstart + xlim * (j - 1) / res - xlim / 2 / res; tau = x * sin(zeta) + buff;
				p[j] = ptemphat[j] * weight + (1 - weight) * (ptemp[j] - ((Hx(j + 1, 1) - Hx(j, 1)) / dx - S1(j)) * dt);
				q[j] = qtemphat[j] * weight + (1 - weight) * (qtemp[j] - ((Hx(j + 1, 2) - Hx(j, 2)) / dx - S2(j)) * dt);
				r[j] = rtemphat[j] * weight + (1 - weight) * (rtemp[j] - ((Hx(j + 1, 3) - Hx(j, 3)) / dx - S3(j)) * dt);
				if (q[j] * qtemphat[j] < 0)
					q[j] = q[j] / abs(q[j]) * p[j] * 1e-8;
				if (r[j] * rtemphat[j] < 0)
					r[j] = r[j] / abs(r[j]) * p[j] * 1e-8;
			}

			ke_curr = mass_total = mom_tot = 0;

			// finding change in rotation rate for angular momentum conservation
			integral1 = 0;
			integral2 = 0;
			mass_total = 0;
			for (j = 2; j <= res + 1; j++)
			{
				x = xstart + xlim * (j - 1) / res - xlim / 2 / res; tau = x * sin(zeta) + buff;
				mass_total = mass_total + 2 * PI * p[j] * dx;
				integral1 = integral1 + 2 * PI * (pow(x * sin(xi), 3) * epsilon * p[j] / tau) * dx;
				integral2 = integral2 + 2 * PI * (pow(x * sin(xi), 2) * epsilon * r[j] / tau) * dx; //
			}
			dom = -(Om * (integral1 - integral1_o) + (integral2 - integral2_o)) / (momincb + integral1);
			alpha = dom / dt;	// angular acceleration
			alpha_momin = -(Om * (integral1 - integral1_o)) / (momincb + integral1) / dt;
			alpha_angmom = -(integral2 - integral2_o) / (momincb + integral1) / dt;

			// mass shedding based on basal presure becoming negative from positive i.e. phi becoming positive from negative
			for (j = 2; j <= res; j++)
			{
				x = xstart + xlim * (j - 1) / res - xlim / 2 / res; tau = x * sin(zeta) + buff;
				def(3, j);
				h[j] = p[j] / tau;
				u[j] = q[j] / p[j];
				v[j] = r[j] / p[j];
				if (phi > 0) // -1e-6
				{
					//cout << "mass shed";
					p[j] = pow(10, -8) * tau;
					q[j] = p[j] * u[j];
					r[j] = p[j] * v[j];//0;			
				}
			}
			// mass shedding from last cell, which has an additional longitudinal curvature
			j = res + 1;
			x = xstart + xlim * (j - 1) / res - xlim / 2 / res; tau = x * sin(zeta) + buff;
			def(3, j);
			if ((phi + 2 * pow(q[res + 1] / p[res + 1], 2) * sin(zeta) / Dc) > 0)
			{
				p[j] = pow(10, -8) * tau;
				q[j] = p[j] * u[j];
				r[j] = p[j] * v[j];
			}

			// reconstruction of main variables
			for (j = 2; j <= res + 1; j++)
			{
				x = xstart + xlim * (j - 1) / res - xlim / 2 / res; tau = x * sin(zeta) + buff;

				h[j] = p[j] / tau;
				u[j] = q[j] / p[j];
				v[j] = r[j] / p[j];

				ke_curr = ke_curr + (pow(u[j], 2) + pow(v[j], 2)) / 2; // estimate of kinetic energy
				mom_tot = mom_tot + abs(h[j] * u[j] * tau) + abs(h[j] * (abs(v[j]) > 1e-5 ? v[j] : 0) * tau);
			}

			maxspeed = 0.0000000001;
			// to evaluate dt from maximum speeds (CFL condition)
			for (j = 1; j <= res + 1; j++)
			{   //pfun(j);
				//x = xstart + xlim * j / res + xlim / 2 / res; tau = x * sin(zeta) + buff;
				eig[j] = ax(j);
				if (maxspeed < abs(eig[j]))
					maxspeed = abs(eig[j]);
			}
			dt = min(dx / 4, dx / 4 / maxspeed);

			// time updation
			t = t + dt;
			timep = timep + dt;
			timesteps = timesteps + 1;
			cout << std::setprecision(10) << timep << "   " << timesteps << "      " << alpha << "    " << Om << "     " << dt << endl;

		} // timesteps loop
	}
	myfileh.close();
	myfileu.close();
	myfilev.close();
	myfilet.close();
	myfiles.close();
	myfilem.close();
	myfilek.close();

	return 0;
}

// use all values for p, q and r from previous step to conserve mass

double pfunc(int jp) // polynomial approximation for p evaluated at these coordinates
{
	return ptemp[jp];
}
double qfunc(int jp) // polynomial approximation for q evaluated at these coordinates
{
	return qtemp[jp];
}
double rfunc(int jp) // polynomial approximation for r evaluated at these coordinates
{
	return rtemp[jp];
}

double ax(int jj)
{    	//	cout<<maxeigenx(jj, 1)<<"   "<< maxeigenx(jj, 2)<<"   "<<"max2"<<endl;

	return (std::max(maxeigenx(jj, 1), maxeigenx(jj, 2))); // 1,2 for plus,minus
}

// calculate eigenvalues
double maxeigenx(int jj, int pm)
// pm is 1,2 for plus and minus of q vector
{
	double e1, e2, e3; //eigenvalues
	double root;
	double xx = xstart + xlim * (jj - 1) / res;
	pfun(jj + 1);
	def(pm, jj + 1);
	if (pm == 1)
	{
		double tau = (xx)* sin(zeta) + buff;

		root = epsilon * Kx * tau * (-Gn * pow(Qp_p, 5) * pow(tau, 4) - pow(Qp_p, 3) * pow(Qp_r, 2) * pow(tau, 3) * cos(zeta) - 2 * pow(Qp_p, 4) * Qp_r * pow(tau, 4) * Om * cos(zeta));
		if (root < 0)
		{//cout<<"non-hperbolic";
			root = 0;
		}
		else
		{
			root = sqrt(root) / (pow(Qp_p, 2) * pow(tau, 3));
		}
		e1 = abs(Qp_q / Qp_p);
		e2 = abs(e1 - root);
		e3 = abs(e1 + root);

		return (std::max(e1, std::max(e2, e3)));
	}
	else
	{
		double tau = (xx)* sin(zeta) + buff;
		root = epsilon * Kx * tau * (-Gn * pow(Qm_p, 5) * pow(tau, 4) - pow(Qm_p, 3) * pow(Qm_r, 2) * pow(tau, 3) * cos(zeta) - 2 * pow(Qm_p, 4) * Qm_r * pow(tau, 4) * Om * cos(zeta));
		//  cout<<phi<<"  "<<jj<<"  "<<root<<endl;
		if (root < 0)
		{//cout<<"non-hperbolic";
			root = 0;
		}
		else
		{
			root = sqrt(root) / (pow(Qm_p, 2) * pow(tau, 3));
		}

		e1 = abs(Qm_q / Qm_p);
		e2 = abs(e1 - root);
		e3 = abs(e1 + root);
		return (std::max(e1, std::max(e2, e3)));
	}
}

double flux(double pp, double qq, double rr, int jj, int no, int pm)
// no is 1,2,3 for f1,f2,f3 
{
	def(pm, jj);
	f1 = qq;
	f2 = qq * qq / pp - epsilon * Kx * pow((pp / tau), 2) * phi * tau / 2.0;
	f3 = qq * rr / pp;

	if (no == 1)
		return f1;
	else if (no == 2)
		return f2;
	else
		return f3;
}

void def(int pm, int jj)
//pm is 1,2,3 for p,m,current
{

	//cout<<jj<<endl;
	double xx = xstart + xlim * (jj - 2) / res;
	double tau = (xx)* sin(zeta) + buff;
	if (pm == 1)
	{
		Gn = Om * Om * cos(zeta) * tau - grav[jj] * cos(xi - eta[jj]);
		vv = Qp_r / Qp_p;
		phi = Gn + 2 * Om * vv * cos(zeta) + cos(zeta) * vv * vv / tau;
		//		cout<<phi<<"   "<<jj<<"    "<<pm<<endl;
	}
	else if (pm == 2)
	{
		Gn = Om * Om * cos(zeta) * tau - grav[jj - 1] * cos(xi - eta[jj - 1]);
		vv = Qm_r / Qm_p;
		phi = Gn + 2 * Om * vv * cos(zeta) + cos(zeta) * vv * vv / tau;
		//		cout<<phi<<"   "<<jj<<"    "<<pm<<endl;
	}
	else
	{
		double tau = (xx + xlim / 2 / res) * sin(zeta) + buff;
		Gn = Om * Om * cos(zeta) * tau - grav[jj] * cos(xi - eta[jj]);
		vv = r[jj] / p[jj];
		phi = Gn + 2 * Om * vv * cos(zeta) + cos(zeta) * vv * vv / tau;
	}
}

void pfun(int jj)
{
	Qp_p = pfunc(jj) - dx * derivative(jj, 1) / 2;
	Qm_p = pfunc(jj - 1) + dx * derivative(jj - 1, 1) / 2;
	Qp_q = qfunc(jj) - dx * derivative(jj, 2) / 2;
	Qm_q = qfunc(jj - 1) + dx * derivative(jj - 1, 2) / 2;
	Qp_r = rfunc(jj) - dx * derivative(jj, 3) / 2;
	Qm_r = rfunc(jj - 1) + dx * derivative(jj - 1, 3) / 2;
}

double Hx(int jj, int no)
{

	pfun(jj);
	fp = flux(Qp_p, Qp_q, Qp_r, jj, no, 1); // first 1 is for f and second 1 is for plus
	fm = flux(Qm_p, Qm_q, Qm_r, jj, no, 2);

	if (no == 1)
	{  //cout<<Qp_p<<"    "<<Qm_p<<"     "<<"   "<<jj<<endl;
		return ((fp + fm) / 2 - ax(jj - 1) * (Qp_p - Qm_p) / 2);
	}
	else if (no == 2)
	{
		return ((fp + fm) / 2 - ax(jj - 1) * (Qp_q - Qm_q) / 2);
	}
	else
	{
		return ((fp + fm) / 2 - ax(jj - 1) * (Qp_r - Qm_r) / 2);
	}
}

double derivative(int jp, int no)
{
	if (jp == res + 3)
		cout << "Daya kuch to gadbad hain" << endl;
	if (no == 1)
	{
		return minmod(theta * (ptemp[jp] - ptemp[jp - 1]) / dx, (ptemp[jp + 1] - ptemp[jp - 1]) / 2 / dx, theta * (ptemp[jp + 1] - ptemp[jp]) / dx);
	}
	else if (no == 2)
	{
		return minmod(theta * (qtemp[jp] - qtemp[jp - 1]) / dx, (qtemp[jp + 1] - qtemp[jp - 1]) / 2 / dx, theta * (qtemp[jp + 1] - qtemp[jp]) / dx);
	}
	else
	{
		return minmod(theta * (rtemp[jp] - rtemp[jp - 1]) / dx, (rtemp[jp + 1] - rtemp[jp - 1]) / 2 / dx, theta * (rtemp[jp + 1] - rtemp[jp]) / dx);
	}

}

double minmod(double a, double b, double c) // calculate minmod
{
	if (a < 0 && b < 0 && c < 0)
		return std::max(a, std::max(b, c));
	else if (a > 0 && b > 0 && c > 0)
		return std::min(a, std::min(b, c));
	else
		return 0;
}

double S1(int jp)
{
	return 0;
}

double S2(int jp)
{
	double ht = htemp[jp];
	double ut = utemp[jp];
	double vt = vtemp[jp];
	def(3, jp);
	int sgn = (ut > 0) ? 1 : ((ut < 0) ? -1 : 0);
	Gt = Om * Om * sin(zeta) * tau + grav[jp] * sin(xi - eta[jp]);
	return sin(zeta) * pow(vt, 2) * ht + ht * tau * phi * tan(delta) * sgn + tau * ht * (Gt + 2 * vt * sin(zeta) * Om);
}

double S3(int jp)
{
	double ht = htemp[jp];
	double ut = utemp[jp];
	double vt = vtemp[jp];
	int sgn = (vt > 0) ? 1 : ((vt < 0) ? -1 : 0);
	return -sin(zeta) * ut * vt * ht + ht * tau * phi * tan(delta) * sgn - tau * ht * (tau * alpha + 2 * ut * sin(zeta) * Om);
}

void Gr()
{
	/*C:\Users\Deepayan Banik\Desktop\So that it looks clean*/
	ifstream etafile("eta45.txt");
	ifstream gravfile("grav45.txt");
	string line;
	grav[0] = eta[0] = grav[res + 3] = eta[res + 3] = 0;

	if (etafile.is_open() || gravfile.is_open())
	{
		int count = 2;
		while (getline(etafile, line))
		{
			eta[count] = stof(line);
			count++;
		}
		count = 2;
		while (getline(gravfile, line))
		{
			grav[count] = grav_scale * stof(line);
			count++;
		}
		grav[1] = grav[2];grav[res + 2] = grav[res + 1];
		eta[1] = eta[2];eta[res + 2] = eta[res + 1];
		etafile.close();
		gravfile.close();
	}
	else cout << "Unable to open file";

}
