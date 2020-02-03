#include "asd-dist.h"
#include "pbar.h"

unsigned int *make_thread_partition(int &num_threads, int ncols) {
	if (num_threads > ncols) num_threads = ncols;
	unsigned int *NUM_PER_THREAD = new unsigned int[num_threads];
	unsigned int div = ncols / num_threads;

	for (int i = 0; i < num_threads; i++)
	{
		NUM_PER_THREAD[i] = 0;
		NUM_PER_THREAD[i] += div;
	}

	for (int i = 0; i < ncols % num_threads; i++)
	{
		NUM_PER_THREAD[i]++;
	}

	return NUM_PER_THREAD;
}

void calc_pw_as_dist(void *order)
{
	work_order_t *p = (work_order_t *)order;
	short **data = p->stru_data->data;
	int size = p->stru_data->nind;
	int numThreads = p->threads;

	Bar *pbar = p->bar;
	int step = size / pbar->totalTicks;
	if (step == 0) step = 1;

	short A, B;

	double ps;
	double *row = NULL;
	int *num_loci = NULL;
	int *ibs0 = NULL;
	int *ibs1 = NULL;
	int *ibs2 = NULL;
	for (int j = 0; j < size; j++)
	{
		if (j % step == 0) advanceBar(*pbar, double(step));

		row = new double[size - j];
		num_loci = new int[size - j];
		if (p->CALC_ALL_IBS)
		{
			ibs0 = new int[size - j];
			ibs1 = new int[size - j];
			ibs2 = new int[size - j];
		}
		for (int k = j; k < size; k++)
		{
			row[k - j] = 0;
			num_loci[k - j] = 0;
			if (p->CALC_ALL_IBS)
			{
				ibs0[k - j] = 0;
				ibs1[k - j] = 0;
				ibs2[k - j] = 0;
			}
			for (int l = p->first_index; l < p->last_index; l++)
			{
				if (j == k)
				{
					A = data[j][l];
					//B = data[k][l];
					if (A < 0 /*|| B < 0*/)
					{
						num_loci[k - j]--;
					}
				}
				else
				{
					A = data[j][l];
					B = data[k][l];
					if (A < 0 || B < 0)
					{
						num_loci[k - j]--;
					}
					else
					{
						ps = proportion_shared(A, B);
						row[k - j] += ps;
						if (p->CALC_ALL_IBS)
						{
							if (ps == 1) ibs2[k - j]++;
							if (ps == 0.5) ibs1[k - j]++;
							if (ps == 0) ibs0[k - j]++;
						}
					}
				}
			}
		}

		pthread_mutex_lock(&mutex_dist_mat);
		for (int m = j; m < size; m++)  DIST_MAT[j][m] += double(row[m - j]);
		pthread_mutex_unlock(&mutex_dist_mat);

		pthread_mutex_lock(&mutex_loci_mat);
		for (int m = j; m < size; m++)  NUM_LOCI[j][m] += num_loci[m - j];
		pthread_mutex_unlock(&mutex_loci_mat);

		if (p->CALC_ALL_IBS)
		{
			pthread_mutex_lock(&mutex_ibs_0);
			for (int m = j; m < size; m++)  IBS_0_MAT[j][m] += ibs0[m - j];
			pthread_mutex_unlock(&mutex_ibs_0);

			pthread_mutex_lock(&mutex_ibs_1);
			for (int m = j; m < size; m++)  IBS_1_MAT[j][m] += ibs1[m - j];
			pthread_mutex_unlock(&mutex_ibs_1);

			pthread_mutex_lock(&mutex_ibs_2);
			for (int m = j; m < size; m++)  IBS_2_MAT[j][m] += ibs2[m - j];
			pthread_mutex_unlock(&mutex_ibs_2);
		}

		delete [] num_loci;
		delete [] row;
		if (p->CALC_ALL_IBS)
		{
			delete [] ibs0;
			delete [] ibs1;
			delete [] ibs2;
		}
	}

	//delete [] key_list;
	delete p;
	return;

}

void calc_pw_as_dist2(void *order)
{
	work_order_t *p = (work_order_t *)order;
	short **data = p->stru_data->data;
	int size = (p->stru_data->nind);
	int numThreads = p->threads;

	Bar *pbar = p->bar;
	int step = size / pbar->totalTicks;
	if (step == 0) step = 1;

	short A1, A2, B1, B2;

	double ps;
	double *row = NULL;
	int *num_loci = NULL;
	int *ibs0 = NULL;
	int *ibs1 = NULL;
	int *ibs2 = NULL;
	for (int j = 0; j < size; j++)
	{
		if (j % step == 0) advanceBar(*pbar, double(step));

		row = new double[size - j];
		num_loci = new int[size - j];
		if (p->CALC_ALL_IBS)
		{
			ibs0 = new int[size - j];
			ibs1 = new int[size - j];
			ibs2 = new int[size - j];
		}
		for (int k = j; k < size; k++)
		{
			row[k - j] = 0;
			num_loci[k - j] = 0;
			if (p->CALC_ALL_IBS)
			{
				ibs0[k - j] = 0;
				ibs1[k - j] = 0;
				ibs2[k - j] = 0;
			}
			for (int l = p->first_index; l < p->last_index; l++)
			{
				if (j == k)
				{
					A1 = data[2 * j][l];
					A2 = data[2 * j + 1][l];
					if (A1 < 0 || A2 < 0)
					{
						num_loci[k - j]--;
					}
				}
				else
				{
					A1 = data[2 * j][l];
					A2 = data[2 * j + 1][l];
					B1 = data[2 * k][l];
					B2 = data[2 * k + 1][l];
					if (A1 < 0 || B1 < 0 || A2 < 0 || B2 < 0)
					{
						num_loci[k - j]--;
					}
					else
					{
						//ps = proportion_shared(A, B);
						ps = proportion_shared2(A1, A2, B1, B2);
						row[k - j] += ps;
						if (p->CALC_ALL_IBS)
						{
							if (ps == 1) ibs2[k - j]++;
							if (ps == 0.5) ibs1[k - j]++;
							if (ps == 0) ibs0[k - j]++;
						}
					}
				}
			}
		}

		pthread_mutex_lock(&mutex_dist_mat);
		for (int m = j; m < size; m++)  DIST_MAT[j][m] += double(row[m - j]);
		pthread_mutex_unlock(&mutex_dist_mat);

		pthread_mutex_lock(&mutex_loci_mat);
		for (int m = j; m < size; m++)  NUM_LOCI[j][m] += num_loci[m - j];
		pthread_mutex_unlock(&mutex_loci_mat);

		if (p->CALC_ALL_IBS)
		{
			pthread_mutex_lock(&mutex_ibs_0);
			for (int m = j; m < size; m++)  IBS_0_MAT[j][m] += ibs0[m - j];
			pthread_mutex_unlock(&mutex_ibs_0);

			pthread_mutex_lock(&mutex_ibs_1);
			for (int m = j; m < size; m++)  IBS_1_MAT[j][m] += ibs1[m - j];
			pthread_mutex_unlock(&mutex_ibs_1);

			pthread_mutex_lock(&mutex_ibs_2);
			for (int m = j; m < size; m++)  IBS_2_MAT[j][m] += ibs2[m - j];
			pthread_mutex_unlock(&mutex_ibs_2);
		}

		delete [] num_loci;
		delete [] row;
		if (p->CALC_ALL_IBS)
		{
			delete [] ibs0;
			delete [] ibs1;
			delete [] ibs2;
		}
	}

	//delete [] key_list;
	delete p;
	return;

}
void calc_weighted_asd(void *order){
	work_order_t *p = (work_order_t *)order;
	short **data = p->stru_data->data;
	int size = p->stru_data->nind;
	int numThreads = p->threads;

	Bar *pbar = p->bar;
	int step = size / pbar->totalTicks;
	if (step == 0) step = 1;

	map<short,double> *freq;
	freq = new map<short,double>[p->last_index - p->first_index];
	int *n;
	n = new int[p->last_index - p->first_index];
	for (int j = 0; j < 2*size; j++){
		for (int l = p->first_index; l < p->last_index; l++){
			if(j == 0){
				n[l] = 0;
			}

			if(data[j][l] >= 0){
				if(freq[l].count(data[j][l]) > 0){
					freq[l][data[j][l]]++;
				}
				else{
					freq[l][data[j][l]] = 1;
				}
				n[l]++;
			}

			if(j == 2*size-1){
				for(map<short,double>::iterator it = freq[l].begin(); it != freq[l].end(); ++it){
					freq[l][it->first] /= double(n[l]);
					//cerr << it->first << ": " << freq[l][it->first] << " ";
				}
				//cerr << endl;
			}
		}

	}

	delete [] n;

	short A, B, C, D;

	double ps;
	double *row = NULL;
	int *num_loci = NULL;

	for (int j = 0; j < size; j++)
	{
		if (j % step == 0) advanceBar(*pbar, double(step));

		row = new double[size - j];
		num_loci = new int[size - j];
		
		for (int k = j+1; k < size; k++)
		{
			row[k - j] = 0;
			num_loci[k - j] = 0;

			for (int l = p->first_index; l < p->last_index; l++)
			{
				A = data[2 * j][l];
				B = data[2 * j + 1][l];
				C = data[2 * k][l];
				D = data[2 * k + 1][l];
				if (A < 0 || C < 0 || B < 0 || D < 0){
					num_loci[k - j]--;
				}
				else{
					double fA, fB;
					fA = freq[l][A];
					fB = freq[l][B];
					//ps = proportion_shared(A, B);
					ps = proportion_shared_weighted(A, B, C, D, fA, fB);
					row[k - j] += ps;
				}
			}
		}

		pthread_mutex_lock(&mutex_dist_mat);
		for (int m = j; m < size; m++)  DIST_MAT[j][m] += double(row[m - j]);
		pthread_mutex_unlock(&mutex_dist_mat);

		pthread_mutex_lock(&mutex_loci_mat);
		for (int m = j; m < size; m++)  NUM_LOCI[j][m] += num_loci[m - j];
		pthread_mutex_unlock(&mutex_loci_mat);

		delete [] num_loci;
		delete [] row;
	}

	delete [] freq;
	delete p;
	return;
}

double proportion_shared_weighted(short A, short B, short C, short D, double fA, double fB){
	double Iac, Iad, Ibc, Ibd;
	Iac = (A == C) ? 1.0 : 0.0;
	Iad = (A == D) ? 1.0 : 0.0;
	Ibc = (B == C) ? 1.0 : 0.0;
	Ibd = (B == D) ? 1.0 : 0.0;

	return 0.25*((1.0-fA)*(Iac+Iad)+(1.0-fB)*(Ibc+Ibd));
}

void calc_grm(void *order)
{
	work_order_t *p = (work_order_t *)order;
	short **data = p->stru_data->data;
	int size = p->stru_data->nind;
	int numThreads = p->threads;

	Bar *pbar = p->bar;
	int step = size / pbar->totalTicks;
	if (step == 0) step = 1;

	short A, B;

	double ps;
	double *row = NULL;
	int *num_loci = NULL;
	int *ibs0 = NULL;
	int *ibs1 = NULL;
	int *ibs2 = NULL;

	double *mu;
	mu = new double[p->last_index - p->first_index];
	int *n;
	n = new int[p->last_index - p->first_index];
	for (int j = 0; j < size; j++){
		for (int l = p->first_index; l < p->last_index; l++){
			if(j == 0){
				mu[l] = 0;
				n[l] = 0;
			}

			if(data[j][l] >= 0){
				mu[l] += data[j][l];
				n[l]++;
			}

			if(j == size-1){
				mu[l] /= double(n[l]);	
			}
		}

	}

	delete [] n;

	for (int j = 0; j < size; j++)
	{
		if (j % step == 0) advanceBar(*pbar, double(step));

		row = new double[size - j];
		num_loci = new int[size - j];
		
		for (int k = j; k < size; k++)
		{
			row[k - j] = 0;
			num_loci[k - j] = 0;
			for (int l = p->first_index; l < p->last_index; l++)
			{
				/*
				if (j == k)
				{
					A = data[j][l];
					//B = data[k][l];
					if (A < 0)
					{
						num_loci[k - j]--;
					}
				}
				else
				*/
				{
					A = data[j][l];
					B = data[k][l];
					if (A < 0 || B < 0 || mu[l] == 0 || mu[l] == 1)
					{
						num_loci[k - j]--;
					}
					else
					{
						ps = gen_rel(A, B, mu[l]);
						row[k - j] += ps;
						/*
						if (p->CALC_ALL_IBS)
						{
							if (ps == 1) ibs2[k - j]++;
							if (ps == 0.5) ibs1[k - j]++;
							if (ps == 0) ibs0[k - j]++;
						}
						*/
					}
				}
			}
		}

		pthread_mutex_lock(&mutex_dist_mat);
		for (int m = j; m < size; m++)  DIST_MAT[j][m] += double(row[m - j]);
		pthread_mutex_unlock(&mutex_dist_mat);

		pthread_mutex_lock(&mutex_loci_mat);
		for (int m = j; m < size; m++)  NUM_LOCI[j][m] += num_loci[m - j];
		pthread_mutex_unlock(&mutex_loci_mat);
	/*
		if (p->CALC_ALL_IBS)
		{
			pthread_mutex_lock(&mutex_ibs_0);
			for (int m = j; m < size; m++)  IBS_0_MAT[j][m] += ibs0[m - j];
			pthread_mutex_unlock(&mutex_ibs_0);

			pthread_mutex_lock(&mutex_ibs_1);
			for (int m = j; m < size; m++)  IBS_1_MAT[j][m] += ibs1[m - j];
			pthread_mutex_unlock(&mutex_ibs_1);

			pthread_mutex_lock(&mutex_ibs_2);
			for (int m = j; m < size; m++)  IBS_2_MAT[j][m] += ibs2[m - j];
			pthread_mutex_unlock(&mutex_ibs_2);
		}
	*/
		delete [] num_loci;
		delete [] row;
		/*
		if (p->CALC_ALL_IBS)
		{
			delete [] ibs0;
			delete [] ibs1;
			delete [] ibs2;
		}
		*/
	}

	//delete [] key_list;
	delete p;
	delete [] mu;
	return;

}


double proportion_shared2(short A1, short A2, short B1, short B2) {

	if (A1 == B1)
	{
		if (A2 == B2)
		{
			return 1;
		}
		else
		{
			return 0.5;
		}
	}
	else if (A1 == B2)
	{
		if (A2 == B1)
		{
			return 1;
		}
		else
		{
			return 0.5;
		}
	}
	else // A1 is unique
	{
		if (A2 == B2)
		{
			return 0.5;
		}
		else if (A2 == B1)
		{
			return 0.5;
		}
		else
		{
			return 0;
		}
	}
}


double proportion_shared(short A, short B)
{
	if (abs(A - B) == 0) return 1;
	if (abs(A - B) == 1) return 0.5;
	return 0;
}

double gen_rel(short A, short B, double mu)
{
	return ( (double(A) - 2.0*mu) * (double(B) - 2.0*mu) ) / (2.0*mu * (1-mu));
}
