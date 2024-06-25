// gcc -I/home2/nishizawa/gsl/include -L/home2/nishizawa/gsl/lib -std=c99 membrane_euler.c -o membrane -lgsl -lgslcblas -lm
// icc membrane_euler.c -o membrane -fast -O3 -ip -ipo -unroll
// icx membrane_euler.c -o membrane -fast -O3

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <getopt.h>
#include <omp.h>

#include "gaussrd.c"

#define MAX_LINE_LENGTH       1024
#define MAX_PARAM_NAME_LENGTH 256

#define m 1.0
#define kb 1.0
#define sigma 1.0  // length of reference


int read_param_file(const char* filename, double* params);
int dynamics(char outputFile[], char parameterFile[]);

// Déclaration des variables
	FILE *pf1, *pf2, *pf3, *pf4, *pf5, *pf6;

    int nchains = 1;
    int ntot;
    int n0;
    int nat_mb = 2432;

    double dtscale;
    double fact = 10;

    // independant variables
    double T = 1;

    int nstep;

    int timestepver;
    int nb_int_est;

    int i, j, k ,l, n;
    int xi, yi, zi, xj, yj, zj;
    int l1, l2;
    int count, lcount = 0;
    int ib , jb ;

    // additionnal variables
    double nu, nu2, nu3, nu6, nu12;
    double dlc0, dlc1, dlcmin, dlcmax;
    double d2, d, I;
    double xa, ya, za, xb, yb, zb, ra, rb, Caa, Cbb, Cab;
    double rdgauss, rdbox1, rdbox2, rdphi, rdtheta;
    double dx, dy ,dz;
    double drx, dry, drz;
    double rxa, rya, rza, rxb, ryb, rzb;
    double e0_sc, b_sc, kappa_sc, Delta_sc, var_sc, fact_sc, mb_sc,ch_ind;

    int step1 = 0;  // include or not the initialisation of the system

    double TIME = 0.0;
    long cpu_time;
    double sec;
    double params[10];
    double DT = 0.0;
    double FINISH_TIME = 0.0;
    int OUTPUT_INTERVAL = 0;           //interval step, not time (int)((int)(FINISH_STEP)/(int)(FRAME))

    double rd1, rd2, rd3;
    double fspmax, flnmax, fmbmax, frdmax, fbdmax;
    double fspmoy, flnmoy, fmbmoy, frdmoy, fbdmoy;
    double fspcount=0, flncount=0, fmbcount=0, frdcount=0, fbdcount=0;

int main(int argc, char *argv[]){
//// get param from prameters file in option ///////////////////////////////////////////

    int c;
    int digit_optind = 0;
    // time_t *tp = (time_t *) 0;
    char *parameterFile=NULL;
    char *outputFile=NULL;
    double temp;

    while (1)
    {
        int this_option_optind = optind ? optind : 1;
        int option_index = 0;
        static struct option long_options[] =
        {
            {"file", 1, 0, 'f'},
            {"out", 1, 0, 'o'},
            {0, 0, 0, 0}
        };

        c = getopt_long_only (argc, argv, "",
                              long_options, &option_index);

        if (c == -1)
            break;
        switch (c)
        {
            case 'f':
                //printf("%s\n",optarg);
                parameterFile = (char *)malloc(sizeof(char) * 256);
                strcpy(parameterFile,optarg);
                break;
            case 'o':
                //printf("%s\n",optarg);
                outputFile = (char *)malloc(sizeof(char) * 256);
                strcpy(outputFile,optarg);
		//		outputToFile_init(outputDir);
                break;
            default:
                printf ("?? getopt returned character code 0%o ??\n", c);
        }
    }

    dynamics(outputFile, parameterFile); // commandes de s=run

    // free(parameterFile);
    // free(outputFile);

}


int dynamics(char outputFile[], char parameterFile[]){

    time_t begin = time(NULL); // Initialisation du marqueur de temps
    unsigned long s = getpid();
    init_genrand(s);
    omp_set_num_threads(8);

    double e0_sc, b_sc, kappa_sc, Delta_sc, var_sc, fact_sc, mb_sc;

    // open param file for reading
    if(parameterFile != NULL){
        if (read_param_file(parameterFile, params) != 0) {
            printf("Error: could not read parameter file %s\n", parameterFile);
            return 1;
        }
        dtscale = params[0];
        nstep = (int)params[1];
        e0_sc = params[2];
        b_sc = params[3];
        kappa_sc = params[4];
        Delta_sc = params[5];
        var_sc = params[6];
        fact_sc = params[7];
        mb_sc = params[8];
        ch_ind = params[9];
    }

    int MAX_STEP = (long)(FINISH_TIME/DT);

    printf("\n");
    printf("dtscale = %lf\n", dtscale);
    printf("nstep = %d\n", nstep);
    printf("e0_sc = %lf\n",e0_sc);
    printf("b_sc = %lf\n",b_sc);
    printf("kappa_sc = %lf\n",kappa_sc);
    printf("Delta_sc = %lf\n", Delta_sc);
    printf("Var_sc = %lf\n", var_sc);
    printf("fact_sc = %lf\n", fact_sc);
    printf("mb_sc = %lf\n", mb_sc);
;
    printf("\n");

    pf1 = fopen("conf_in.xyz","r");
    fscanf(pf1,"%d\nchaine\n", &ntot);

    int **links = malloc(ntot * sizeof(int *));
    for(i=0; i<ntot; i++){links[i] = malloc(2 * sizeof(int));}

	double *rx, *ry, *rz;
    rx = (double *)malloc(ntot * sizeof(double)); double *rx0 = rx;
    ry = (double *)malloc(ntot * sizeof(double)); double *ry0 = ry;
    rz = (double *)malloc(ntot * sizeof(double)); double *rz0 = rz;
    for(i=0; i<ntot; i++){*rx=0.0;*ry=0.0;*rz=0.0;}

    double *fx, *fy, *fz;
    fx = (double *)malloc(ntot * sizeof(double)); double *fx0 = fx;
    fy = (double *)malloc(ntot * sizeof(double)); double *fy0 = fy;
    fz = (double *)malloc(ntot * sizeof(double)); double *fz0 = fz;
    for(i=0; i<ntot; i++){*fx=0.0;*fy=0.0;*fz=0.0;}

    int *trinearest;
    trinearest = (int *)malloc(3 * sizeof(int));

    double *d2max;
    d2max = (double *)malloc(3 * sizeof(double));

    int checkflex[ntot-1];
    for(i=0; i<ntot-1; i++){checkflex[i]=0;}
    //printf("tableaux des filaments définis\n");

    for(i=0;i<ntot;i++){fscanf(pf1,"C %lf %lf %lf\n", &*(rx+i), &*(ry+i), &*(rz+i));}
    fclose(pf1);

    // dependant variables

    double e0 = kb*T;                 e0*=e0_sc;
    double e01 = kb*T;                e01 = e01*mb_sc;
    double b = 1000.0*sigma*kb*T;     b*=b_sc;
    double kappa = 1000.0*kb*T;       kappa*=kappa_sc;
    double tau = sqrt( (m*sigma*sigma)/(kb*T) );
    double dt = dtscale*tau;
    double gamma = 0.1*tau;
    double Delta = sigma/2;printf("%lf %lf\n",sigma,Delta);

    double rc2 = (1.122462048309373*sigma)*(1.122462048309373*sigma);
    double rc2m = (1.12246*sigma + Delta)*(1.12246*sigma + Delta);
    double lc0 = 0.5*sigma;
    double lc1 = 0.5*sigma;
    double lcmin = 0.25*sigma;
    double lcmax = 0.75*sigma;

    double var = 2 * kb * T * gamma; var*=var_sc;
    double sqvar = sqrt(var);
    double rdt = sqrt(dt);
    int samechain;

    int lfil = 30;
    int sidecell = 20;
    int nprox = 200;
    double mbsize = 25*sigma;
    double lcell = sqrt(16);
    double midcell = 1.5*mbsize+sqrt(16);

    int step1 = 0;  // include or not the initialisation of the system

    double fact = 10; fact*=fact_sc;
    double fch = ch_ind*fact;
    // to be sure that the force is applied one time only on each particle
    int activeflag[ntot];
    int chiralflag[ntot];

    double vec1x, vec1y, vec1z, vec2x, vec2y, vec2z;
    double vecx, vecy, vecz, norm;
    double normx, normy, normz;
    double inoutx, inouty, inoutz;
    double verifinout;

    pf1 = fopen("membrane.xyz", "r");
    fscanf(pf1, "%d\nmembrane\n", &nat_mb);

	double *rx_mb, *ry_mb, *rz_mb;
    rx_mb = (double *)malloc(nat_mb * sizeof(double));
    ry_mb = (double *)malloc(nat_mb * sizeof(double));
    rz_mb = (double *)malloc(nat_mb * sizeof(double));
    for(i=0; i<nat_mb; i++){rx_mb[i]=0.0;ry_mb[i]=0.0;rz_mb[i]=0.0;}

    double *fx_mb, *fy_mb, *fz_mb;
    fx_mb = (double *)malloc(nat_mb * sizeof(double));
    fy_mb = (double *)malloc(nat_mb * sizeof(double));
    fz_mb = (double *)malloc(nat_mb * sizeof(double));
    for(i=0; i<nat_mb; i++){fx_mb[i]=0.0;fy_mb[i]=0.0;fz_mb[i]=0.0;}

    for(i=0; i<nat_mb; i++){
        fscanf(pf1, "v %lf %lf %lf\n", &*(rx_mb+i), &*(ry_mb+i), &*(rz_mb+i));
    }

for(i=0; i<nat_mb; i++){
    rx_mb[i]*=mbsize;
    ry_mb[i]*=mbsize;
    rz_mb[i]*=mbsize;
}
fclose(pf1);

int ****clist;

    // Allocation de mémoire pour le tableau clist
    clist = malloc(sidecell * sizeof(int***));
    for (i = 0; i < sidecell; i++) {
        clist[i] = malloc(sidecell * sizeof(int**));
        for (j = 0; j < sidecell; j++) {
            clist[i][j] = malloc(sidecell * sizeof(int*));
            for (k = 0; k < sidecell; k++) {
                clist[i][j][k] = malloc(nprox * sizeof(int));
            }
        }
    }

    // Déclaration du tableau clist_mb
    int ****clist_mb;

    // Allocation de mémoire pour le tableau clist
    clist_mb = malloc(sidecell * sizeof(int***));
    for (i = 0; i < sidecell; i++) {
        clist_mb[i] = malloc(sidecell * sizeof(int**));
        for (j = 0; j < sidecell; j++) {
            clist_mb[i][j] = malloc(sidecell * sizeof(int*));
            for (k = 0; k < sidecell; k++) {
                clist_mb[i][j][k] = malloc(nprox * sizeof(int));
            }
        }
    }

printf("start\n");


FILE *fp_outfile;

    // open file for save
    if(outputFile != NULL){
        if((fp_outfile = fopen(outputFile,"w")) == NULL) {
            printf("open error\n");
            return 1;
        }else{
            fprintf(fp_outfile,"%d\nchaine\n", ntot);

        }
    }else{
        fp_outfile=NULL;
        printf("Null\n");

    }

//////////// guiding the filaments on the vesicle ///////////////////////////
printf("%lf %lf %lf %lf %lf\n",b,e0,e01,Delta,kappa);

for(k=0; k<nstep; k++){

    for(j=0; j<ntot; j++){activeflag[j]=0;chiralflag[j]=0;}

    for(xi=0;xi<sidecell;xi++){
        for(yi=0;yi<sidecell;yi++){
            for(zi=0;zi<sidecell;zi++){
                clist[xi][yi][zi][0] = 0;
                clist_mb[xi][yi][zi][0] = 0;
            }
        }
    }

fspmax=0, flnmax=0, fmbmax=0, frdmax=0, fbdmax=0;
fspmoy=0, flnmoy=0, fmbmoy=0, frdmoy=0, fbdmoy=0;
fspcount=0, flncount=0, fmbcount=0, frdcount=0, fbdcount=0;

    for(i=0; i<ntot; i++){

        fx[i] = 0.0;
        fy[i] = 0.0;
        fz[i] = 0.0;

    }

     #pragma omp parallel for private(drx,dry,drz,dlc0,dlcmax,d,d2,i,l,nu,nu2,nu3,nu6,nu12,samechain)

    for(i=0; i<ntot; i++){
        for(l=0; l<i; l++){ // interaction filaments-filaments

                drx = rx[i] - rx[l];
                dry = ry[i] - ry[l];
                drz = rz[i] - rz[l];
                d2 = drx*drx + dry*dry + drz*drz;

                if(d2<rc2 && d2>0.0){
                if(i>l+2 || (int)(i/lfil) != (int)(l/lfil)){;

                    d = sqrt(d2);

                    nu = sigma / (d);
                    nu2 = nu * nu;
                    nu3 = nu2 * nu;
                    nu6 = nu3 * nu3;
                    nu12 = nu6 * nu6;

                    I = e0 * (48.0*nu12 - 24.0*nu6);

                fx[i] += I * (drx/(d2));
                fy[i] += I * (dry/(d2));
                fz[i] += I * (drz/(d2));
                fx[l] -= I * (drx/(d2));
                fy[l] -= I * (dry/(d2));
                fz[l] -= I * (drz/(d2));

                }
                }
        }

        if(i%lfil != 0 && i%lfil != 29){

           drx = 0.5*(rx[i+1] - rx[i-1]);
           dry = 0.5*(ry[i+1] - ry[i-1]);
           drz = 0.5*(rz[i+1] - rz[i-1]);

            }
				// for the extremities
        else{
			if(i%lfil == 29){

				drx = rx[i] - rx[i-1];
				dry = ry[i] - ry[i-1];
				drz = rz[i] - rz[i-1];

            }

			if(i%lfil == 0){

                drx = rx[i+1] - rx[i];
                dry = ry[i+1] - ry[i];
                drz = rz[i+1] - rz[i];

            }}

		d2 = drx*drx + dry*dry + drz*drz;
		d = sqrt(d2);

		fx[i] += (fact * (drx/d) );
		fy[i] += (fact * (dry/d) );
		fz[i] += (fact * (drz/d) );

		activeflag[i]=1;

    }
    #pragma omp parallel for private(xi,yi,zi,n)

    for(n=0;n<ntot;n++){
    xi = (int)((rx[n]+midcell)/lcell);
    yi = (int)((ry[n]+midcell)/lcell);
    zi = (int)((rz[n]+midcell)/lcell);
    if(xi<0 || xi>sidecell || yi<0 || yi>sidecell || zi<0 || zi>sidecell ){
        printf("erreur d'attribution ri : iteration = %d, particle %d \n",k,n);
        printf("erreur d'attribution ri : r = %lf %lf %lf %d %d %d\n",rx[n],ry[n],rz[n],xi,yi,zi);}

    clist[xi][yi][zi][0]++;
    clist[xi][yi][zi][ clist[xi][yi][zi][0] ] = n; // sort the particles into the cells
    }

    #pragma omp parallel for private(xi,yi,zi,n)

    for(n=0;n<nat_mb;n++){
    xi = (int)((rx_mb[n]+midcell)/lcell);
    yi = (int)((ry_mb[n]+midcell)/lcell);
    zi = (int)((rz_mb[n]+midcell)/lcell);

    clist_mb[xi][yi][zi][0]++;
    clist_mb[xi][yi][zi][ clist_mb[xi][yi][zi][0] ] = n; // sort the particles into the cells
    }

    #pragma omp parallel for private(yi,zi,l1,l2,xj,yj,zj,i,j,drx,dry,drz,nu,nu2,nu3,nu6,nu12,d,d2)

    for(xi=1;xi<sidecell-1;xi++){ // crossing all the cells
        for(yi=1;yi<sidecell-1;yi++){
            for(zi=1;zi<sidecell-1;zi++){

                if(clist[xi][yi][zi][0]>0){
                for(l1=1;l1<clist[xi][yi][zi][0]+1;l1++){

                    i = clist[xi][yi][zi][l1];

                    trinearest[0]=0; trinearest[1]=0; trinearest[2]=0;
                    d2max[0]=0.0; d2max[1]=0.0; d2max[2]=0.0;

                    for(xj=xi-1;xj<xi+2;xj++){ // checking the lists around
                        for(yj=yi-1;yj<yi+2;yj++){
                            for(zj=zi-1;zj<zi+2;zj++){

                                for(l2=1;l2<clist_mb[xj][yj][zj][0]+1;l2++){
                                    if(clist[xi][yi][zi][0]>0 && clist_mb[xj][yj][zj][0]>0){

                                    j = clist_mb[xj][yj][zj][l2];

                                    drx = rx[i] - rx_mb[j];
                                    dry = ry[i] - ry_mb[j];
                                    drz = rz[i] - rz_mb[j];
                                    d2 = drx*drx + dry*dry + drz*drz;

                                    if(d2<26 && d2>Delta){

                                        d = sqrt(d2);

                                        nu = sigma / (d-Delta);
                                        nu2 = nu * nu;
                                        nu3 = nu2 * nu;
                                        nu6 = nu3 * nu3;
                                        nu12 = nu6 * nu6;

                                        if(d<rc2m){
                                        fx[i] += (1 * (48.0*nu12 - 24.0*nu6) * (drx/(d*(d-Delta))));
                                        fy[i] += (1 * (48.0*nu12 - 24.0*nu6) * (dry/(d*(d-Delta))));
                                        fz[i] += (1 * (48.0*nu12 - 24.0*nu6) * (drz/(d*(d-Delta))));
                                        }

                                        if(d>=rc2m){
                                        fx[i] += (e01 * (48.0*nu12 - 24.0*nu6) * (drx/(d*(d-Delta))));
                                        fy[i] += (e01 * (48.0*nu12 - 24.0*nu6) * (dry/(d*(d-Delta))));
                                        fz[i] += (e01 * (48.0*nu12 - 24.0*nu6) * (drz/(d*(d-Delta))));
                                        }

                                        if(i%lfil == 0 || i%lfil == 29){

                                        if(d2 >= d2max[0]){
                                            trinearest[2] = trinearest[1]; d2max[2] = d2max[1];
                                            trinearest[1] = trinearest[0]; d2max[1] = d2max[0];
                                            trinearest[0] = j; d2max[0] = d2;
                                        }
                                        if(d2 < d2max[0] && d2 >= d2max[1]){
                                            trinearest[2] = trinearest[1]; d2max[2] = d2max[1];
                                            trinearest[1] = j; d2max[1] = d2;
                                        }
                                        if(d2 < d2max[1] && d2 >= d2max[2]){
                                            trinearest[2] = j; d2max[2] = d2;
                                        }

                                        }
                            }
                        }
                    }

                    }}}

                    vec1x = rx_mb[trinearest[0]] - rx_mb[trinearest[1]];
                    vec1y = ry_mb[trinearest[0]] - ry_mb[trinearest[1]];
                    vec1z = rz_mb[trinearest[0]] - rz_mb[trinearest[1]];

                    vec2x = rx_mb[trinearest[2]] - rx_mb[trinearest[1]];
                    vec2y = ry_mb[trinearest[2]] - ry_mb[trinearest[1]];
                    vec2z = rz_mb[trinearest[2]] - rz_mb[trinearest[1]];

                    normx = vec1y*vec2z - vec1z*vec2y;
                    normy = vec1z*vec2x - vec1x*vec2z;
                    normz = vec1x*vec2y - vec1y*vec2x;

                    inoutx = rx[i] - rx_mb[trinearest[1]];
                    inouty = ry[i] - ry_mb[trinearest[1]];
                    inoutz = rz[i] - rz_mb[trinearest[1]];

                    verifinout = inoutx*normx + inouty*normy + inoutz*normz;

                    if(verifinout<0){

                        normx=-normx;
                        normy=-normy;
                        normz=-normz;

                    }

                    if(i%lfil != 0){

                        drx = rx[i-1] - rx[i];
                        dry = ry[i-1] - ry[i];
                        drz = rz[i-1] - rz[i];

                        vecx = (normy*drz - normz*dry);
                        vecy = (normz*drx - normx*drz);
                        vecz = (normx*dry - normy*drx);

                        norm = sqrt(vecx*vecx + vecy*vecy + vecz*vecz);

                        if(norm != 0.0){

                            vecx /= norm;
                            vecy /= norm;
                            vecz /= norm;

                            fx[i] += fch * vecx;
                            fy[i] += fch * vecy;
                            fz[i] += fch * vecz;

                    }

                    }

                    if(i%lfil != 29){

                        drx = rx[i+1] - rx[i];
                        dry = ry[i+1] - ry[i];
                        drz = rz[i+1] - rz[i];

                        vecx = (normy*drz - normz*dry);
                        vecy = (normz*drx - normx*drz);
                        vecz = (normx*dry - normy*drx);

                        norm = sqrt(vecx*vecx + vecy*vecy + vecz*vecz);

                        if(norm != 0.0){

                            vecx /= norm;
                            vecy /= norm;
                            vecz /= norm;

                            fx[i] += fch * vecx;
                            fy[i] += fch * vecy;
                            fz[i] += fch * vecz;

                    }

                    }



        }// aplying force between the particles i and j
    }
}}}
// end of the listing method

   #pragma omp parallel for private(drx,dry,drz,dlc0,dlcmax,dlc1,dlcmin,d,d2,l,j,i)

    for(i=0; i<ntot-1; i++){

            if(i%lfil != 29){

            drx = rx[i+1] - rx[i];
            dry = ry[i+1] - ry[i];
            drz = rz[i+1] - rz[i];

            d2 = drx*drx + dry*dry + drz*drz;
            d = sqrt(d2);

            if(d2>lc0*lc0  && d2>0.0){

                dlc0 = (lc0 - d)*(lc0 - d);
                dlcmax = (lcmax - d)*(lcmax - d);

                I = (( (dlc0+dlcmax) * -b * exp( 1/(lc0-d) ) ) / (dlc0*dlcmax));

                fx[i+1] += I * (drx/d);
                fy[i+1] += I * (dry/d);
                fz[i+1] += I * (drz/d);
                fx[i] -= I * (drx/d);
                fy[i] -= I * (dry/d);
                fz[i] -= I * (drz/d);

            }

            if(d2<lc1*lc1 && d2>0.0){

                dlc1 = (d - lc1)*(d - lc1);
                dlcmin = (d - lcmin)*(d - lcmin);

                I = (( (dlc1+dlcmin) * b * exp( 1/(d-lc1) ) ) / (dlc1*dlcmin));

                fx[i+1] += I * (drx/d);
                fy[i+1] += I * (dry/d);
                fz[i+1] += I * (drz/d);
                fx[i] -= I * (drx/d);
                fy[i] -= I * (dry/d);
                fz[i] -= I * (drz/d);

            }

    }}

   #pragma omp parallel for private(xa, ya, za, xb, yb, zb, ra, rb, Caa, Cbb, Cab,i)

    for(i=0; i<ntot; i++){

        if(i%lfil > 1){

				xa = rx[i-1] - rx[i-2];
				ya = ry[i-1] - ry[i-2];
				za = rz[i-1] - rz[i-2];
				xb = rx[i] - rx[i-1];
				yb = ry[i] - ry[i-1];
				zb = rz[i] - rz[i-1];
				Caa = xa*xa + ya*ya + za*za;
				Cbb = xb*xb + yb*yb + zb*zb;
				Cab = xa*xb + ya*yb + za*zb;
				fx[i] = fx[i] - kappa * (((Cab/Cbb)*xb - xa) / sqrt(Cbb*Caa));
				fy[i] = fy[i] - kappa * (((Cab/Cbb)*yb - ya) / sqrt(Cbb*Caa));
				fz[i] = fz[i] - kappa * (((Cab/Cbb)*zb - za) / sqrt(Cbb*Caa));
                }

        if(i%lfil > 0 && i%lfil < 29){
				xa = rx[i] - rx[i-1];
				ya = ry[i] - ry[i-1];
				za = rz[i] - rz[i-1];
				xb = rx[i+1] - rx[i];
				yb = ry[i+1] - ry[i];
				zb = rz[i+1] - rz[i];
				Caa = xa*xa + ya*ya + za*za;
				Cbb = xb*xb + yb*yb + zb*zb;
				Cab = xa*xb + ya*yb + za*zb;
				fx[i] = fx[i] + kappa * (((Cab/Cbb)*xb - (Cab/Caa)*xa + xb - xa) / sqrt(Cbb*Caa));
				fy[i] = fy[i] + kappa * (((Cab/Cbb)*yb - (Cab/Caa)*ya + yb - ya) / sqrt(Cbb*Caa));
				fz[i] = fz[i] + kappa * (((Cab/Cbb)*zb - (Cab/Caa)*za + zb - za) / sqrt(Cbb*Caa));
                }

        if(i%lfil < 28){
				xa = rx[i+1] - rx[i];
				ya = ry[i+1] - ry[i];
				za = rz[i+1] - rz[i];
				xb = rx[i+2] - rx[i+1];
				yb = ry[i+2] - ry[i+1];
				zb = rz[i+2] - rz[i+1];
				Caa = xa*xa + ya*ya + za*za;
				Cbb = xb*xb + yb*yb + zb*zb;
				Cab = xa*xb + ya*yb + za*zb;
				fx[i] = fx[i] + kappa * (((Cab/Caa)*xa - xb) / sqrt(Cbb*Caa));
				fy[i] = fy[i] + kappa * (((Cab/Caa)*ya - yb) / sqrt(Cbb*Caa));
				fz[i] = fz[i] + kappa * (((Cab/Caa)*za - zb) / sqrt(Cbb*Caa));

                }

    }

    for(i=0; i<ntot; i++){
        rd1 = gasdev();
        rd2 = gasdev();
        rd3 = gasdev();
        fx[i] += ( sqvar * rd1 )/rdt; if( frdmax < rd1 ){ frdmax = rd1 ;}
        fy[i] += ( sqvar * rd2 )/rdt; if( frdmax < rd2 ){ frdmax = rd2 ;}
        fz[i] += ( sqvar * rd3 )/rdt; if( frdmax < rd3 ){ frdmax = rd3 ;}

    }


    for(i=0; i<ntot; i++){
        rx[i] += ( fx[i] / gamma )*dt;
        ry[i] += ( fy[i] / gamma )*dt;
        rz[i] += ( fz[i] / gamma )*dt;
    }

    if( k % (nstep/1000) == 0 ){
        for(i=0; i<ntot; i++){
            fprintf(fp_outfile,"C %lf %lf %lf\n", rx[i], ry[i], rz[i]);
        }
    }

    fflush(fp_outfile);

    }


    // Affichage du temps d'exécution
	time_t end = time(NULL); // Fin du marqueur de temps
	unsigned long secondes = (unsigned long) difftime(end, begin);
	printf("Temps d'exécution : %ld secondes\n", secondes);

return 0;} // Fin de la simulation

    int read_param_file(const char* filename, double* params) {
    FILE *param_file;
    char line[MAX_LINE_LENGTH];
    char *name, *value;
    int param_idx;

    param_file = fopen(filename, "r");
    if (param_file == NULL) {
        return 1;
    }

    while (fgets(line, MAX_LINE_LENGTH, param_file) != NULL) {
        name = strtok(line, ",");
        value = strtok(NULL, ",");

        if (name == NULL || value == NULL) {
            fclose(param_file);
            return 1;
        }

        if (strcmp(name, "dtscale") == 0) {
            param_idx = 0;
        } else if (strcmp(name, "nstep") == 0) {
            param_idx = 1;
        } else if (strcmp(name, "e0_sc") == 0) {
            param_idx = 2;
        } else if (strcmp(name, "b_sc") == 0) {
            param_idx = 3;
        } else if (strcmp(name, "kappa_sc") == 0) {
            param_idx = 4;
        } else if (strcmp(name, "Delta_sc") == 0) {
            param_idx = 5;
        } else if (strcmp(name, "var_sc") == 0) {
            param_idx = 6;
        } else if (strcmp(name, "fact_sc") == 0) {
            param_idx = 7;
        } else if (strcmp(name, "mb_sc") == 0) {
            param_idx = 8;
        } else if (strcmp(name, "ch_ind") == 0) {
            param_idx = 9;
        } else {
            // Do nothing for unrecognized parameters
            continue;
        }

        params[param_idx] = atof(value);
    }

    fclose(param_file);

    return 0;
}
