
#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<cmath>
using namespace std;

int main(void){
	FILE *fp;
	FILE *fp1;
	FILE *fp2;

	char NAME[50];

	double sumvx;
	double sumvy;
	double sumv2x;
	double sumv2y;
	double sumv2;

	int i,j,I,K;

	int npart=7; // 7 points representing the particles are along each direction
	int N = (npart-1)*(npart-1); // N is the number of real particles

	/* explanation: array of 36*36,
	Since over a length space the first and the last points are the same, out of 7 only 6 will be real particles
	hence the array is of size 36 by 36.
	*/
	double delt;
	double timesteps;
/* Define the vectors of position and velocities */
	double x[N]; // x array
	double y[N]; // y array
	double xm[N]; // approximation of the previous x particle position
	double ym[N]; // approximation of the previous y particle position
	double vx[N]; // x velocity array
	double vy[N]; // y velocity array

	double temp = 0.1; // temperature in scaled units
	double fs; // scaling constant

	double en; // current potential energy
	double fx[N],fy[N]; // array for holding the force
	double xr, yr; // x and y component of r
	double r2; // r^2, square of r
	double box; // box length
	double boxby2; // half the boxlength
	double rc2; // rcut square
	double ff, r2i, r6i; //
	double ecut; // energy cut units

	double a; // Try to keep this close to 1 or more than 1
	double t; // time variable
	double xx, yy;
	double etot; // total energy
	double tmp; //

	/*create a folder named output and*/
	/*delete the previous output folder: print statement is just to show*/
	printf("rm -rf output/*\n");
	system("rm -rf output/*");

	/*create an output folder*/
	system("mkdir output");
	/* ====================================================================   */
	box = 5.5; // total box length: What is the relation with box length and number of particles?
	a = box/(npart-2);  // how much space is available to one molecule
	/*
	Explanation: Suppose 3 molecules occupy a length of 2. Then how much length is available
	for each molecule? This is given by 2/3=0.66. This is in fact the length of freedom for the molecule.
	The reason of taking npart-2 (actual particles are npart-1) is an attempt to get how much length a particle does fully enjoy.
	*/
	box = box + a;

	/* here the box length is increased 'a' step/length just to consider
	enough space for the unconsidered particle in the previous step. Now, this 'box'
	is the full box length */

	boxby2 = box*0.5; // half of the updated full box length
	rc2 = 2.5; /* cut off radius and is manually decided */

	/*
	Now, the cut-off radius cannot be greater than half the system length.
	This is because,cut-off radius is used to minimize the load of calculations. This is
	the disance after which the potential of one molecule does not affect the other. Hence its
	impact on those are not considered. Since the system is symmetric about the centre,
	if a molecule is at the centre of the box, it cannot hve an impacting potential greater than the boxlength/2.
	For the molecules situated at the rigth edge of the box, the impact will be on the molecules of left edge (being periodic system)
	*/

	if(rc2 > boxby2){
		printf("Cut-off radius greater than half-box-length; quitting\n");
		exit(0);
	}

	rc2 = rc2*rc2; // find rcut square.
	ecut = 1./(rc2*rc2*rc2); //
	ecut = 4.*ecut*(ecut-1); // value of the LJ potential at the rcut.
	delt = 1.e-3; // time step
	timesteps = (int) (1./delt); // no of time steps required to reach 1 unit in steps of delt.
	timesteps = timesteps + 10002; // + two more time steps: this is base less anyway and is done just to make the total timesteps to be equal to 1000.

	/** Initialise the MD program **/
	/* Initialize the positions to start the program*/
	/* i and j are the particle indices along x and y direction respectively*/
	for(i=1; i < npart; ++i){
		for(j=1; j < npart; ++j){
		  I = j + (i-1)*(npart-1); // This is the index of each molecules in a square lattice.
			if(i==1 && j==1){
				x[1] = a/2.; // x-position of the first particle
				y[1] = a/2.; // y-position of the first particle
			}
			/* Positions of all the subsequent particles */
			else{
				x[I] = x[1] + (i-1)*a + (double) (rand()/RAND_MAX)*0.0005;
				y[I] = y[1] + (j-1)*a + (double) (rand()/RAND_MAX)*0.0005;
			}
		}
	}
	/* Initialize the velocities to start the program*/
	/* Start with a little random velocity */
	for(i=1; i < npart; ++i){
		for(j=1; j < npart; ++j){
			I = j + (i-1)*(npart-1);
			vx[I] = (double) rand()/RAND_MAX - 0.5;
			vy[I] = (double) rand()/RAND_MAX - 0.5;
		}
	}
	/*Since, initial kinetic energy must be zero, and energy is conserved
	we initiate with zero sum of velocities along x and y direction
	sumvx and sumvy are the component velocity of the centre of mass */
	sumvx = 0.;
	sumvy = 0.;
	/*square sum of velocities*/
	sumv2x = 0.;
	sumv2y = 0.;
	sumv2 = 0.0;

	/*Now, add the initiated velocities to the sum of velocities for all the particles*/
	for(i=1; i < N+1; ++i){
		sumvx = sumvx + vx[i];
		sumvy = sumvy + vy[i];
		sumv2x = sumv2x + vx[i]*vx[i];
		sumv2y = sumv2y + vy[i]*vy[i];
	}
	/* Average sum of velocities */
	sumvx = sumvx/N;
	sumvy = sumvy/N;
	/* Average sum of square velocities */
	sumv2x = sumv2x/N;
	sumv2y = sumv2y/N;
	/* Sum of square velocities */
	sumv2 = sumv2x+ sumv2y;
	/* velocity scaling factor (kT/(0.5mv^2))*/
	fs = sqrt(2.*temp/sumv2);

	/* Subtract the individual velocities from the average sum of
	velocities and scale it by the scaling factor*/

	for(i=1; i < npart; ++i){
		for(j=1; j < npart; ++j){
			I = j + (i-1)*(npart-1); // particle index
			vx[I] = (vx[I] - sumvx)*fs;
			vy[I] = (vy[I] - sumvy)*fs;
			/* Find the previous positions */
			xm[I] = x[I] - vx[I]*delt;
			ym[I] = y[I] - vy[I]*delt;
		}
	}
	/* Start the time loop */
	t = 0.0;

	/* Print data */
	fp = fopen("output/output.dat","w");
	fprintf(fp,"# temp %le\n",temp);
	fprintf(fp,"# N %d\n",npart-1);
	fprintf(fp,"# box %le\n",box);
	fprintf(fp,"# a %le\n",a);
	fprintf(fp,"# delt %le\n",delt);
	fprintf(fp,"# rcut %le\n",sqrt(rc2));
	fprintf(fp,"#time KE PE Tot.E temp\n");

	fp2 = fopen("output/plot.plt","w");
	fprintf(fp2,"set si rat 1.0\n");
	fprintf(fp2,"set view 0,90\n");

	/* Time loop starts here */
	for(K = 0; K < timesteps; ++K){

		/** Calculate the forces **/
		for(i=1; i < npart; ++i){
			for(j=1; j < npart; ++j){
			I = j + (i-1)*(npart-1);
			fx[I] = 0.0;
			fy[I] = 0.0;
			}
		}
		/* start with zero potential energy */
		en = 0.0;

		/*
		loop over all the particles:
		i should start from 1 to N-1 and j from 2 to N
		just to distinguish between the particles
		*/
		for(i=1; i < N; ++i){
			for(j=i+1; j < N+1; ++j){
				/* calculate the difference between the molecules one by one*/
				xr = x[i] - x[j];
				yr = y[i] - y[j];
				/* If the difference happens to be out of the box,
				create an image particle in the box maintaining the periodicity */

				if(xr > boxby2) xr = xr-box;
				else if(xr < -boxby2) xr = xr + box;

				if(yr > boxby2) yr = yr-box;
				else if(yr < -boxby2) yr = yr + box;
				/*square the difference*/
				r2 = xr*xr + yr*yr;

				/* If the difference is more than a certain number and less than
				the cutoff radius then only calulate the forces*/

				if(r2 > 1.e-12 && r2 < rc2){
					r2i = 1./r2;
					r6i = r2i*r2i*r2i;
					ff = 48.*r2i*r6i*(r6i-0.5); // potential energy

					/* By Newton's thrid law if the first particle is
					moved to right the seond in the interaction would move left*/

					fx[i] = fx[i] + ff*xr;
					fy[i] = fy[i] + ff*yr;
					fx[j] = fx[j] - ff*xr;
					fy[j] = fy[j] - ff*yr;
					en = en + 4.*r6i*(r6i-1.) - ecut;
				}
			}
		}

	/** Integrate the equations of motion: push the particles to a new position **/

		/*Once again initiate the centre of mass velocities to be zero*/
			sumvx = 0.;
			sumvy = 0.;
			sumv2x = 0.;
			sumv2y = 0.;

			for(i=1; i < N+1; ++i){
				/*Calculate the new positions from the verlet algorithm */
				xx = 2*x[i] - xm[i] +delt*delt*fx[i];
				yy = 2*y[i] - ym[i] +delt*delt*fy[i];
				if(xx > box){
				xx = xx - box;
				}
				else if (xx < 0){
					xx = xx + box;
				}
				if(yy > box){
					yy = yy - box;
				}
				else if (yy < 0){
					yy = yy + box;
				}
				xr = xx-xm[i];
				yr = yy-ym[i];
				if(xr > boxby2) xr = xr-box;
				else if(xr < -boxby2) xr = xr + box;
				if(yr > boxby2) yr = yr-box;
				else if(yr < -boxby2) yr = yr + box;

				/* velocity verlet */
				vx[i] = xr/(2.*delt);
				vy[i] = yr/(2.*delt);

				/*Update the centre of mass velocities*/
				sumvx = sumvx + vx[i];
				sumvy = sumvy + vy[i];
				sumv2x = sumv2x + vx[i]*vx[i];
				sumv2y = sumv2y + vy[i]*vy[i];

				/*Update positions in previous time*/
				xm[i] = x[i];
				ym[i] = y[i];

				/*Update positions in current time*/
				x[i] = xx;
				y[i] = yy;
			}

			sumv2 = sumv2x + sumv2y;

			/*Instantaneous temperature*/
			temp = sumv2/(2.*N);

			/*Total energy per particle*/
			etot = (en + 0.5*sumv2)/N;

			/*Print the diagnostics*/
			fprintf(fp,"%lf %lf %lf %lf %lf\n",t,0.5*sumv2/N,en/N,etot,temp);

			if(K%100 == 0){
				sprintf(NAME,"time%d.dat",K);
				//fprintf(fp2,"plot \"%s\"\n",NAME);
				fprintf(fp2,"plot \"%s\" pt 7 ps 3\n",NAME);
				fprintf(fp2,"pause 1\n");
				sprintf(NAME,"output/time%d.dat",K);
				fp1 = fopen(NAME,"w");
				for(i=1; i < npart; ++i){
					for(j=1; j < npart; ++j){
					I = j + (i-1)*(npart-1);
					fprintf(fp1,"%le %le 1\n",xm[I],ym[I]);
					}
				fprintf(fp1,"\n");
				}
				fclose(fp1);

				printf("Time=%g\n",t);
			}

	/*Forward the time*/
	t = t + delt;
}
fclose(fp);
fclose(fp2);
return 0;
}
