/*
 * util.c
 *
 *  Created on: Feb 20, 2014
 *      Author: justin
 */
void heapsort(double* R, int* Ni, int* Nj, int n ){
	int l = n/2+1;
	int r = n;
	int i,j;
	double RR = 0;
	int Nii, Njj;
	while(1){
		if (l > 1) {
			l = l-1;
			RR = R[l];
			Nii = Ni[l];
			Njj = Nj[l];
		} else {
			RR = R[r];
			Nii = Ni[r];
			Njj = Nj[r];
			R[r] = R[1];
			Ni[r] = Ni[1];
			Nj[r] = Nj[1];
			r -= 1;
			if (r == 1) {
				R[1] = RR;
				Ni[1] = Nii;
				Nj[1] = Njj;
				break;
			}
		}
		j = l;
		while(1){
			i = j;
			j *= 2;
			if (j>r){
				R[i] = RR;
				Ni[i] = Nii;
				Nj[i] = Njj;
				break;
			}
			if (j<r){
				if(R[j]<R[j+1]){
					j=j+1;
				}
			}
			if(RR>=R[j]){
				R[i]=RR;
				Ni[i] = Nii;
				Nj[i] = Njj;
				break;
			}
			R[i]=R[j];
			Ni[i] = Ni[j];
			Nj[i] = Nj[j];
		}
	}
}



