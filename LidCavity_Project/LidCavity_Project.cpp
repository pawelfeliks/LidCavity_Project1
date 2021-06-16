#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#define grid 50

using namespace std;

int main()
{   // u, v - velocities, p - pressure
	// uc, vc, pc - for gridpoints
	// un, vn, pn - updated values
	// e - for checking divergence
	// dx, dy - step in x and y direction
	// dt - timestep
	// Re - Reynolds number
	double u[grid][grid + 1], un[grid][grid + 1], uc[grid][grid];
	double v[grid + 1][grid], vn[grid + 1][grid], vc[grid][grid];
	double p[grid + 1][grid + 1], pn[grid + 1][grid + 1], pc[grid][grid];
	double e[grid + 1][grid + 1];
	int i, j, step;
	double dx, dy, dt, tau, delta, error, Re;
	step = 1;
	dx = 1.0 / (grid - 1);
	dy = 1.0 / (grid - 1);
	dt = 0.001;
	delta = 4.5;
	error = 1.0;
	Re = 100.0;

	// Initializing u - velocity in x direction
	for (i = 0; i <= (grid - 1); i++)
	{
		for (j = 0; j <= (grid); j++)
		{
			u[i][j] = 0.0;
			u[i][grid] = 1.0;   // u velocity on the top of domain
			u[i][grid - 1] = 1.0; // u velocity under top of domain
		}
	}

	// Initializing v - velocity in y direction
	for (i = 0; i <= (grid); i++)
	{
		for (j = 0; j <= (grid - 1); j++)
		{
			v[i][j] = 0.0;
		}
	}

	// Initializing p
	for (i = 0; i <= (grid); i++)
	{
		for (j = 0; j <= (grid); j++)
		{
			p[i][j] = 1.0;
		}
	}

	// Iterative loop depending on error
	while (error > 0.001)
	{
		// Solving x-momentum equation for interior points
		for (i = 1; i <= (grid - 2); i++)
		{
			for (j = 1; j <= (grid - 1); j++)
			{
				un[i][j] = u[i][j] - dt * ((u[i + 1][j] * u[i + 1][j] - u[i - 1][j] * u[i - 1][j]) / dx / 2.0
					+ ((u[i][j] + u[i][j + 1]) * (v[i][j] + v[i + 1][j]) - (u[i][j] + u[i][j - 1]) * (v[i + 1][j - 1] + v[i][j - 1])) / 4.0 / dy)
					- (p[i + 1][j] - p[i][j]) * dt / dx
					+ dt * 1.0 / Re * ((u[i + 1][j] - 2.0 * u[i][j] + u[i - 1][j]) / dx / dx + (u[i][j + 1] - 2.0 * u[i][j] + u[i][j - 1]) / dy / dy);
			}
		}

		// Boundary conditions for u velocity
		for (j = 1; j <= (grid - 1); j++)
		{
			un[0][j] = 0.0;                     // Left wall of the domain
			un[grid - 1][j] = 0.0;                // Right wall of the domain
		}

		for (i = 0; i <= (grid - 1); i++)
		{
			un[i][0] = -un[i][1];               // Bottom wall of the domain
			un[i][grid] = 2 - un[i][grid - 1];    // Top wall of the domain
		}


		// Solving y-momentum equation for interior points
		for (i = 1; i <= (grid - 1); i++)
		{
			for (j = 1; j <= (grid - 2); j++)
			{
				vn[i][j] = v[i][j] - dt * (0.25 * ((u[i][j] + u[i][j + 1]) * (v[i][j] + v[i + 1][j]) - (u[i - 1][j] + u[i - 1][j + 1]) * (v[i][j] + v[i - 1][j])) / dx
					+ (v[i][j + 1] * v[i][j + 1] - v[i][j - 1] * v[i][j - 1]) / 2.0 / dy)
					- dt / dy * (p[i][j + 1] - p[i][j])
					+ dt * 1.0 / Re * ((v[i + 1][j] - 2.0 * v[i][j] + v[i - 1][j]) / dx / dx + (v[i][j + 1] - 2.0 * v[i][j] + v[i][j - 1]) / dy / dy);
			}
		}

		// Boundary conditions for v velocity
		for (j = 1; j <= (grid - 2); j++)
		{
			vn[0][j] = -vn[1][j];           // Left wall of the domain
			vn[grid][j] = -vn[grid - 1][j];   // Right wall of the domain
		}

		for (i = 0; i <= (grid); i++)
		{
			vn[i][0] = 0.0;                 // Bottom wall of the domain
			vn[i][grid - 1] = 0.0;            // Top wall of the domain
		}

		// Solving continuity equation
		for (i = 1; i <= (grid - 1); i++)
		{
			for (j = 1; j <= (grid - 1); j++)
			{
				pn[i][j] = p[i][j] - dt * delta * ((un[i][j] - un[i - 1][j]) / dx + (vn[i][j] - vn[i][j - 1]) / dy);
			}
		}


		// Boundary conditions for pressure
		for (j = 0; j <= (grid); j++)
		{
			pn[0][j] = pn[1][j];            // Left wall of the domain
			pn[grid][j] = pn[grid - 1][j];    // Right wall of the domain
		}

		for (i = 1; i <= (grid - 1); i++)
		{
			pn[i][0] = pn[i][1];            // Bottom wall of the domain
			pn[i][grid] = pn[i][grid - 1];    // Top wall of the domain
		}


		// Checking error
		error = 0.0;

		for (i = 1; i <= (grid - 1); i++)
		{
			for (j = 1; j <= (grid - 1); j++)
			{
				e[i][j] = ((un[i][j] - un[i - 1][j]) / dx + (vn[i][j] - vn[i][j - 1]) / dy);
				error = error + abs(e[i][j]);
			}
		}

		if (step % 1000 == 1)
		{
			printf("For the step %d error is equal to %5.3lf \n", step, error);
		}


		// Overwriting u for updated values
		for (i = 0; i <= (grid - 1); i++)
		{
			for (j = 0; j <= (grid); j++)
			{
				u[i][j] = un[i][j];
			}
		}

		// Overwriting v for updated values
		for (i = 0; i <= (grid); i++)
		{
			for (j = 0; j <= (grid - 1); j++)
			{
				v[i][j] = vn[i][j];
			}
		}

		// Overwriting p for updated values
		for (i = 0; i <= (grid); i++)
		{
			for (j = 0; j <= (grid); j++)
			{
				p[i][j] = pn[i][j];
			}
		}
		step = step + 1;
	}

	// Final results with use of linear interpolation
	for (i = 0; i <= (grid - 1); i++)
	{
		for (j = 0; j <= (grid - 1); j++)
		{
			uc[i][j] = 0.5 * (u[i][j] + u[i][j + 1]);
			vc[i][j] = 0.5 * (v[i][j] + v[i + 1][j]);
			pc[i][j] = 0.25 * (p[i][j] + p[i + 1][j] + p[i][j + 1] + p[i + 1][j + 1]);
		}
	}

	// OUTPUT DATA
	ofstream outfile;
	outfile.open("Results_lid_cav.txt");

	// verify the file opnened
	if (outfile.is_open())
	{
		outfile << "Size of mesh: " << grid << "x" << grid << endl;
		outfile << "    X           Y           U           V           P   \n";

		for (j = 0; j < (grid); j++)
		{
			for (i = 0; i < (grid); i++)
			{
				double x, y;
				x = i * dx;
				y = j * dy;
				setprecision(6);

				outfile << fixed << x << "  " << setw(10) << y << "  " << setw(10) << uc[i][j] << "  " << setw(10) << vc[i][j] << "  " << setw(10) << pc[i][j] << endl;
			}
		}
		outfile << "    X           Y           U           V           P   \n";

		outfile.close();
	}

	return 0;
}