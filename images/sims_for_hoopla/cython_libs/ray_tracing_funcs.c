#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <fitsio.h>
//#include <omp.h>
//#include "mycosmology.h"

//--------------------------------------------------------------------
void inverse_cic_omp(double *in_map, double *posy1, double *posy2, double ysc1, double ysc2,double dsi, int nsx, int nsy, int nlx, int nly, double *out_map) {

	int index;
	int i1,j1,i,j;
	double xb1,xb2;
	double ww1,ww2,ww3,ww4,wx,wy;

#pragma omp parallel num_threads(1)	\
	shared(in_map,nlx,nly,ysc1,ysc2,dsi,nsx,nsy,out_map) \
	private(index,i,j,i1,j1,xb1,xb2,wx,wy,ww1,ww2,ww3,ww4)
	{
	double *out_map_sp;
	out_map_sp = (double *)calloc(nlx*nly,sizeof(double));
	#pragma omp for schedule(dynamic,16)

	for(i=0;i<nlx;i++) for(j=0;j<nly;j++) {

		index = i*nlx+j;

		xb1 = (posy1[index]-ysc1)/dsi+(double)nsx/2.0-0.5;
		xb2 = (posy2[index]-ysc2)/dsi+(double)nsy/2.0-0.5;
		printf("%d %lf %lf\n",index,posy1[i],posy2[i]);

		i1 = (int)xb1;
		j1 = (int)xb2;

		wx = 1.-(xb1-(double)(i1));
		wy = 1.-(xb2-(double)(j1));

		ww1 = wx*wy;
		ww2 = wx*(1.0-wy);
		ww3 = (1.0-wx)*wy;
		ww4 = (1.0-wx)*(1.0-wy);

		if (i1<0||i1>nsx-2||j1<0||j1>nsy-2) continue;

		out_map_sp[index] = ww1*in_map[i1*nsx+j1]
					      + ww2*in_map[i1*nsx+j1+1]
					      + ww3*in_map[(i1+1)*nsx+j1]
					      + ww4*in_map[(i1+1)*nsx+j1+1];
	}

	#pragma omp critical
	{
		for(i=0;i<nlx;i++) for(j=0;j<nly;j++) {
			out_map[i*nly+j] += out_map_sp[i*nly+j];
		}
	}
	free(out_map_sp);
	}
}
//--------------------------------------------------------------------
void inverse_cic(double *in_map, double *posy1, double *posy2, double ysc1, double ysc2,double dsi, int nsx, int nsy, int nlx, int nly, double *out_map) {

	int i1,j1,i,j;
	int index;
	double xb1,xb2;
	double ww1,ww2,ww3,ww4,wx,wy;

	for(i=0;i<nlx;i++) for(j=0;j<nly;j++) {

		index = i*nly+j;

		xb1 = (posy1[index]-ysc1)/dsi+(double)nsx/2.0-0.5;
		xb2 = (posy2[index]-ysc2)/dsi+(double)nsy/2.0-0.5;

		i1 = (int)xb1;
		j1 = (int)xb2;

		wx = 1.-(xb1-(double)(i1));
		wy = 1.-(xb2-(double)(j1));

		ww1 = wx*wy;
		ww2 = wx*(1.0-wy);
		ww3 = (1.0-wx)*wy;
		ww4 = (1.0-wx)*(1.0-wy);

		if ((i1<0)||(i1>nsx-2)||(j1<0)||(j1>nsy-2)) continue;

		out_map[index] = ww1*in_map[i1*nsy+j1]
					   + ww2*in_map[i1*nsy+j1+1]
					   + ww3*in_map[(i1+1)*nsy+j1]
					   + ww4*in_map[(i1+1)*nsy+j1+1];
	}
}
//----------------------------------------------------------------------------------
