#include <stdio.h>
#include<cmath>
#include <stdlib.h>


/*

exe arg1 arg2 arg3
arg1->  normalized (to the peak) donor emision spectra file 
arg2->  normalized (to the peak) acceptor absortion spectra file 
arg3->  acceptor maximum extinction coefficient

*/

int main(int argc, char *argv[])
{
	if(argc<4) {
		printf ("Usage: exe arg1 arg2 arg3\n");
		printf ("arg1->  Normalized (to the peak) donor emision spectra file \n");
		printf ("arg2->  Normalized (to the peak) acceptor absortion spectra file \n");
		printf ("arg3->  Acceptor maximum extinction coefficient\n");
		return -1;
	}
	
	FILE  *donor_em, *acceptor_ab;
	
	donor_em = fopen(argv[1],"r");
	if(donor_em == NULL){
		printf("- %s file doesn't exist",argv[1]);
		return -1;
	}
		
	acceptor_ab = fopen(argv[2],"r");
	if(acceptor_ab == NULL){
		printf("- %s file doesn't exist",argv[2]);
		return -1;
	}
	
	
	float e_cof = atof(argv[3]);
	
	double acc_spectra[1000];
	double don_spectra[1000];
	float c;int w;
	float integral=0;
	
	while (fscanf(donor_em,"%d%f",&w,&c)!= EOF){
		don_spectra[w] = c;
		integral += c;
	}
	
	printf("Area: %f\n",integral);
	
	for (int i=0;i<1000;i++){
		don_spectra[i] = don_spectra[i]/integral;
	}
	
	float newintegral = 0;
	
	for (int i=0;i<1000;i++){
		newintegral+=don_spectra[i];
	}
	
	printf("Now area is: %f\n",newintegral);
	
	while (fscanf(acceptor_ab,"%d%f",&w,&c)!= EOF){
		acc_spectra[w] = c;
	}
	
	double e_spectra[1000];
	
	for (int i=300;i<800;i++){
		e_spectra[i] =acc_spectra[i]*e_cof;
	}
	
	double overlap = 0;
	
	for (int i = 300;i<= 800;i++){
		overlap+=(float)pow(i,4)*don_spectra[i]*e_spectra[i];
		//~ printf("overlap: %f %f %f %d\n",overlap,(float)don_spectra[i],(float)e_spectra[i],i);
	}
	
	printf("Overlap: %f\n",(float)overlap);
	
	fclose(donor_em);
	fclose(acceptor_ab);
	return 0;
}

