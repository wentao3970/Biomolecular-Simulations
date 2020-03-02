/*********************************************************
  nummod.c
  -------------------
copyright : (C) 2006 by Ryan Brenke and Philip Yang Shen
email : rbrenke@bu.edu yangshen@bu.edu
 *********************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>

#include _MOL_INCLUDE_

#define PI 3.14159265

void print_short_args (char* app_name);
void print_help_args (char* app_name);
void print_args (char* app_name);
void my_rotate_atomgrp (struct atomgrp* prot, struct atomgrp* rotprot, struct rmatrix* rmatrix, struct tvector* center_of_rotation);
void perturb_atomgrp (struct atomgrp* ag, struct atomgrp* moved_ag, double translate_rangei, double rotate_range, struct tvector* center_of_rotation);
float calEnergy(struct atomgrp* A, struct atomgrp* B, struct prms* prms);
float* calGradients(struct atomgrp* agA, struct atomgrp* agB, struct prms* prms);
#define MAXSLEN 200
char* ATOM_PRM_FILE;

int main (int argc, char* argv[])
{
	char* app_name = argv[0];
	if (argc < 2)
	{
		print_short_args (app_name);
		print_help_args (app_name);
		exit (EXIT_FAILURE);
	}

	char* rec_ifile; // input file
	char* lig_ifile; // input file
	//char* ofile = (char*) mymalloc (MAXSLEN * sizeof (char));
	ATOM_PRM_FILE = atom_prm_file (ATOM_PRM);


	size_t slen; // string length

	int c;
	while (1)
	{

		c = getopt (argc, argv, "hp:");
		if (c == -1)
			break;
		switch (c)
		{
			case 'h':
				print_args (app_name);
				return 0;
			case 'p':
				slen = strlen (optarg);
				if (slen > MAXSLEN)
				{
					fprintf (stderr, "atom parameter file name %s is too long\n", optarg);
					exit (EXIT_FAILURE);
				}
				ATOM_PRM_FILE = optarg;
				break;
			default:
				break;
		}
	}

	if (optind+1 < argc)
	{
		rec_ifile = argv[optind];
		optind++;
		lig_ifile = argv[optind];
		optind++;
		while (optind < argc)
		{
			printf ("ignored argument: %s\n", argv[optind]);
			optind++;
		}
	}
	else
	{
		print_short_args (app_name);
		print_help_args (app_name);
		exit (EXIT_FAILURE);
	}


	struct prms* prms = read_prms (ATOM_PRM_FILE, _MOL_VERSION_);
	//Read receptor and ligand structure
	struct atomgrp* agA = read_file_atomgrp (rec_ifile, prms);
	struct atomgrp* agB = read_file_atomgrp (lig_ifile, prms);


	struct tvector* temp_tv = (struct tvector*) mymalloc (sizeof (struct tvector));
	//temp_tv->X = 10;
	//temp_tv->Y = 0;
	//temp_tv->Z = 0;
	
	// Create two protein A and B
 	struct atomgrp* agA_moved = copy_atomgrp(agA);
	struct atomgrp* agB_copy = copy_atomgrp(agB);
	
	// Caculate the CoM of A and B
	struct tvector* cmA = center_of_mass(agA_moved);
	struct tvector* cmB = center_of_mass(agB_copy);
	
	// Assign the temp_tv coordinate
	temp_tv -> X = 0.5* (cmA->X - cmB->X);
	temp_tv -> Y = 0.5* (cmA->Y - cmB->Y);
	temp_tv -> Z = 0.5* (cmA->Z - cmB->Z);
	
	// Translate agA
	translate_atomgrp (agA_moved, agA_moved,temp_tv); 
	
	char* current_ofile = (char*) mymalloc (100 * sizeof (char)); 
	sprintf (current_ofile, "test%d.ms", 1);
	fprint_file_atomgrp(current_ofile, agA_moved, prms);
	struct tvector* cm=center_of_mass(agA_moved);
	printf("Center of mass %.3f %.3f %.3f\n",cm->X,cm->Y,cm->Z);
	srand(time(NULL));
	double r=(double)(rand())/((RAND_MAX+1.0));
	printf("Random number %.4f\n",r);
	// Evaluate E
	clock_t start, end;
	double elapsed;
	float E=0;
	//E = calEnergy(agA, agB, prms);

	int i;
	start = clock();
	for (i = 0; i < 10; i++)
	{
		E = calEnergy (agA, agB, prms);
	}
	end = clock();
	elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf ("time: %.3f\n", elapsed);
	printf("Energy %.3f\n", E); 
	printf("Gradients array is: ", calGradients(agA, agB, prms));

    sizeofarray = sizeof(calGradients(agA, agB, prms)) / sizeof(calGradients(agA, agB, prms)[0])
    for (int i=0; i 

	return EXIT_SUCCESS;
}

// Q2
float calEnergy(struct atomgrp* agA, struct atomgrp* agB, struct prms* prms ){
	float E = 0;
	float Etemp = 0;
	float x;
	float y;
	float z;
	float dist;
	int i;
	int j;

	for (i = 0; i < agA->natoms; i++){
		for(j = 0; j < agB->natoms; j++){
			x = agA->atoms[i].X - agB->atoms[j].X;
			y = agA->atoms[i].Y - agB->atoms[j].Y;
			z = agA->atoms[i].Z - agB->atoms[j].Z;
			dist = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
			Etemp = (prms->chrgs[agA->atoms[i].atom_typen])*(prms->chrgs[agB->atoms[j].atom_typen])/dist; 
			E = E + Etemp; 
		}
	}
	return E;
}

//Q3

float* calGradients(struct atomgrp* agA, struct atomgrp* agB, struct prms* prms ){
	float E = 0;
	float Etemp = 0;
	float x1,x2;
	float y1,y2;
	float z1,z2;
	float dist;
	int i;
	int j;
	float Fx;
	float Fy;
	float Fz;
	float arr[3*agA->natoms];
	for (i = 0; i < agA->natoms; i++){
		for(j = 0; j < agB->natoms; j++){
			x1 = agA->atoms[i].X; 
			x2 = agB->atoms[j].X;
			y1 = agA->atoms[i].Y;
			y2 = agB->atoms[j].Y;
			z1 = agA->atoms[i].Z;
			z2 = agB->atoms[j].Z;
			dist = sqrt(pow(x1-x2,2)+pow(y1-y2,2)+pow(z1-z2,2));
			Etemp = (prms->chrgs[agA->atoms[i].atom_typen])*(prms->chrgs[agB->atoms[j].atom_typen]);
			Fx = Fx + (-0.5)*3*dist*2*(x1-x2)*Etemp; 
			Fy = Fy + (-0.5)*3*dist*2*(y1-y2)*Etemp;
       		Fz = Fz + (-0.5)*3*dist*2*(z1-z2)*Etemp;
						
		}
		arr[i+0] = Fx;
		arr[i+1] = Fy;
		arr[i+2] = Fz;
	}
	return arr;
}

void print_help_args (char* app_name)
{
	fprintf (stderr, "try '%s -h' for a list of arguments\n", app_name);
}

void print_short_args (char* app_name)
{
	fprintf (stderr, "usage: %s [arguments] RECEPTOR LIGAND\n", app_name);
	fprintf (stderr, "print correlation energy of RECEPTOR and LIGAND\n");
}

void print_args (char* app_name)
{
	print_short_args (app_name);

	printf ("\n");
	printf ("arguments:\n");
	printf ("   %-20s Use <atom.prm> as atom parameter file (default: %s)\n", "-p <atom.prm>", ATOM_PRM_FILE);
}

void my_rotate_atomgrp (struct atomgrp* prot, struct atomgrp* rotprot, struct rmatrix* rmatrix, struct tvector* center_of_rotation)
{
	int atomi;
	for (atomi = 0; atomi < prot->natoms; atomi++)
	{
		float X = prot->atoms[atomi].X;
		float Y = prot->atoms[atomi].Y;
		float Z = prot->atoms[atomi].Z;
		rotprot->atoms[atomi].X = rmatrix->a11*(X - center_of_rotation->X) + rmatrix->a12*(Y - center_of_rotation->Y) + rmatrix->a13*(Z - center_of_rotation->Z) + center_of_rotation->X;
		rotprot->atoms[atomi].Y = rmatrix->a21*(X - center_of_rotation->X) + rmatrix->a22*(Y - center_of_rotation->Y) + rmatrix->a23*(Z - center_of_rotation->Z) + center_of_rotation->Y;
		rotprot->atoms[atomi].Z = rmatrix->a31*(X - center_of_rotation->X) + rmatrix->a32*(Y - center_of_rotation->Y) + rmatrix->a33*(Z - center_of_rotation->Z) + center_of_rotation->Z;
	}
}

void perturb_atomgrp (struct atomgrp* ag, struct atomgrp* moved_ag, double translate_range, double rotate_range, struct tvector* center_of_rotation)
{
	struct tvector* temp_tv = (struct tvector*) mymalloc (sizeof (struct tvector));
	struct rmatrix* temp_rmatrix = (struct rmatrix*) mymalloc (sizeof (struct rmatrix));

	double qw,qx,qy,qz; //quaternions
	double r,theta,phi; //intermediate spherical coordinates for uniform sampling

	if(translate_range<0.0)
	{
		printf ("Input error: translational range should be NONNEGATIVE\n");
		exit (EXIT_FAILURE);
		//printf ("Notice: translational range is forced to be ZERO.\n");
		//translate_range=0;	    
	}
	if(rotate_range<0.0)
	{
		printf ("Input error: rotational range should be NONNEGATIVE\n");
		exit (EXIT_FAILURE);
		//printf ("Notice: rotational range is forced to be ZERO.\n");
		//rotate_range=0;
	}
	else
		if(rotate_range>PI)
		{
			printf ("Input error: maximum rotational range should be PI\n");
			exit (EXIT_FAILURE);
			//printf ("Notice: rotational range is forced to be PI.\n");	    
			//rotate_range=PI;
		}


	//translational perturbation
	/*Uniform sampling in a sphere of radius translate_range*/
	/* intermediate spherical coordinates (r,theta,phi) */

	//random number generator: to modify


	r = translate_range * pow((rand() / ((double)RAND_MAX + 1)),1/3.0);
	phi = acos(1-2*(rand() / ((double)RAND_MAX + 1)));
	theta = 2*PI*(rand() / ((double)RAND_MAX + 1));

	temp_tv->X = r * cos(theta) * sin(phi);
	temp_tv->Y = r * sin(theta) * sin(phi);
	temp_tv->Z = r * cos(phi);


	//rotational perturbation
	//Kuffner paper describes how to generate uniform unit quaternions
	//global uniform sampling: max range PI
	//to modify if need ``local'' orietational perturbation
	//uniform sampling in a sphere of exponential coordinates
	//essential: space of exponential coordinates is similar to the Euclidean space of translations

	r = rotate_range * pow((rand() / ((double)RAND_MAX + 1)),1/3.0);
	phi = acos(1-2*(rand() / ((double)RAND_MAX + 1)));
	theta = 2*PI*(rand() / ((double)RAND_MAX + 1));

	//transform into quaternions
	qw = cos(r/2);
	qx = sin(r/2) * cos(theta) * sin(phi);
	qy = sin(r/2) * sin(theta) * sin(phi);
	qz = sin(r/2) * cos(phi);

	//generate rotation matrix
	temp_rmatrix->a11 = 1 - 2 * (pow(qy,2) + pow(qz,2));
	temp_rmatrix->a12 = 2 * (qx*qy - qz*qw);
	temp_rmatrix->a13 = 2 * (qx*qz + qy*qw);

	temp_rmatrix->a21 = 2 * (qx*qy + qz*qw);
	temp_rmatrix->a22 = 1 - 2 * (pow(qx,2) + pow(qz,2));
	temp_rmatrix->a23 = 2 * (qy*qz - qx*qw);

	temp_rmatrix->a31 = 2 * (qx*qz - qy*qw);
	temp_rmatrix->a32 = 2 * (qy*qz + qx*qw);
	temp_rmatrix->a33 = 1 - 2 * (pow(qx,2) + pow(qy,2));

	my_rotate_atomgrp (ag, moved_ag, temp_rmatrix, center_of_rotation);
	translate_atomgrp (moved_ag, moved_ag, temp_tv);
}
