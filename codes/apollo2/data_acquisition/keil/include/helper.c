#include "helper.h"
#include "trained_model.h"
#include "am_util.h"

/**************************************************/
/*                 Preprocessing                  */
/**************************************************/

// Crop data in time
void preprocessing_crop_data(float32_t * p_data_in, float32_t * p_data_out, uint32_t t_start, uint32_t t_end){
	
	uint16_t i, i_helper=0;
	
	for (i=(T_START-1)*DATA_SIZE*DATA_SIZE; i<T_END*DATA_SIZE*DATA_SIZE; i++){
		p_data_out[i_helper] = p_data_in[i];
		i_helper++;
	};
	
};

// Gaussian smoothing
void preprocessing_gaussian_smooth(float32_t * p_data_in, float32_t * p_data_out){
	
	float32_t data_tmp[DATA_SIZE*DATA_SIZE];
	float32_t gaussian_filter[FILTER_SIZE]; // 1D gaussian kernel
	uint16_t t, i;
	
	// Create the 1D kernel, save it in gaussian_filter
	preprocessing_gaussian_filter(gaussian_filter);
	
	for (t=0; t<DELTA_T; t++){
		
		for (i=0; i<DATA_SIZE*DATA_SIZE; i++){
			data_tmp[i] = p_data_in[t*(DATA_SIZE*DATA_SIZE) + i];
		};
				
		preprocessing_gaussian_smooth_single_frame(data_tmp, gaussian_filter, &p_data_out[t*(DATA_SIZE*DATA_SIZE)]);

	};
	
};

void preprocessing_gaussian_smooth_single_frame(float32_t * p_data_in, float32_t * p_gaussian_filter, float32_t * p_data_out){
	
	int16_t i, j; // row and column index
	float32_t data_replicated[(DATA_SIZE + (FILTER_SIZE - 1)) * (DATA_SIZE + (FILTER_SIZE - 1))];
	float32_t data_tmp[(DATA_SIZE + (FILTER_SIZE - 1))];
	float32_t data_dx[(DATA_SIZE + (FILTER_SIZE - 1))*(DATA_SIZE + (FILTER_SIZE - 1) + (FILTER_SIZE - 1))];
	float32_t data_filtered_T[(DATA_SIZE + (FILTER_SIZE - 1) + (FILTER_SIZE - 1))*(DATA_SIZE + (FILTER_SIZE - 1) + (FILTER_SIZE - 1))];
	uint16_t k=0; // index for p_data_out

	
	// Replicate padding
	preprocessing_replicate_padding(p_data_in, data_replicated);
	
	
	// 1D convolution in x direction and save the partially filtered matrix in data_dx
	for (i=0; i<(DATA_SIZE + (FILTER_SIZE - 1)); i++){
		
		// copy the j-th row of replicated matrix into data_tmp
		for (j=0; j<(DATA_SIZE + (FILTER_SIZE - 1)); j++){
			data_tmp[j] = data_replicated[i*(DATA_SIZE + (FILTER_SIZE - 1)) + j];
		};
		
		// compute 1D convolution in x direction
		arm_conv_f32(data_tmp, (DATA_SIZE + (FILTER_SIZE - 1)), p_gaussian_filter, FILTER_SIZE, &data_dx[i*(DATA_SIZE + (FILTER_SIZE - 1) + (FILTER_SIZE - 1))]);

	};
	
	
	// 1D convolution in y direction performed on data_dx
	for (j=0; j<(DATA_SIZE + (FILTER_SIZE - 1) + (FILTER_SIZE - 1)); j++){
		
		// copy the j-th column of data_dx into data_tmp
		for (i=0; i<(DATA_SIZE + (FILTER_SIZE - 1)); i++){
			data_tmp[i] = data_dx[i*(DATA_SIZE + (FILTER_SIZE - 1) + (FILTER_SIZE - 1)) + j];
		};
		
		// compute 1D convolution in y direction, it gives the transpose of the completely filtered matrix
		arm_conv_f32(data_tmp, (DATA_SIZE + (FILTER_SIZE - 1)), p_gaussian_filter, FILTER_SIZE, &data_filtered_T[j*(DATA_SIZE + (FILTER_SIZE - 1) + (FILTER_SIZE - 1))]);
		
	};
	
	
	// output the filtered matrix, cropped into original size (this case 16-by-16)
	for (j=(FILTER_SIZE - 1); j<(DATA_SIZE + (FILTER_SIZE - 1)); j++){
		for (i=(FILTER_SIZE - 1); i<(DATA_SIZE + (FILTER_SIZE - 1)); i++){
			p_data_out[k] = data_filtered_T[i*(DATA_SIZE + (FILTER_SIZE - 1) + (FILTER_SIZE - 1))+j];
			k++;
		};
	};
	
};

void preprocessing_gaussian_filter(float32_t * p_gaussian_filter){
	// Gaussian filter: it computes the 1D gaussian kernel of size 2*ceil(2*SIGMA)+1
	// uint16_t filter_size = (2*ceil(2*SIGMA)+1); // I'm not using Heap
	
	int16_t i, x;
	int16_t x0 = (FILTER_SIZE+1)/2;
	float32_t twoSigmaSquare = 2*SIGMA*SIGMA;
	float32_t sum=0, oneOverSum;
	float32_t kernel[FILTER_SIZE];
	
	// Create the kernel vector
	for (i = -(FILTER_SIZE-1)/2; i <= (FILTER_SIZE-1)/2; i++){
		x = i+x0;
		kernel[x-1] = exp(-((x-x0)*(x-x0))/twoSigmaSquare);
		sum = sum + kernel[x-1];
	};
	
	// Normalization
	oneOverSum = 1.0f/sum;
	arm_scale_f32(kernel, oneOverSum, p_gaussian_filter, FILTER_SIZE);
	
};

void preprocessing_replicate_padding(float32_t * p_data_in, float32_t * p_data_out){
	// Replicate padding
	uint16_t border = (FILTER_SIZE - 1)/2;
	uint16_t size_replicated; // The size of replicated matrix
	uint16_t i, j;
	
	size_replicated = DATA_SIZE + (FILTER_SIZE-1);
	
	for (i=0; i<size_replicated; i++){  // row
		for (j=0; j<size_replicated; j++){  // col
			
			if (i<border){
				if (j<border){ // top left corner
					p_data_out[(i*size_replicated)+j] = p_data_in[0];
				};
				if (j>=size_replicated-border){ // top right corner
					p_data_out[(i*size_replicated)+j] = p_data_in[DATA_SIZE-1];
				};
				if (j>=border && j<size_replicated-border){ // central top rows
					p_data_out[(i*size_replicated)+j] = p_data_in[j-border];
				};
			};
			
			if (i>=size_replicated-border){
				if (j<border){ // bottom left corner
					p_data_out[(i*size_replicated)+j] = p_data_in[(DATA_SIZE-1)*DATA_SIZE];	
				};
				if (j>=size_replicated-border){ // bottom right corner
					p_data_out[(i*size_replicated)+j] = p_data_in[(DATA_SIZE*DATA_SIZE)-1];	
				};
				if (j>=border && j<size_replicated-border){ // central bottom rows
					p_data_out[(i*size_replicated)+j] = p_data_in[((DATA_SIZE-1)*DATA_SIZE)+(j-border)];
				};
			};

			if (i>=border && i<size_replicated-border){
				if (j<border){  // left columns
					p_data_out[(i*size_replicated)+j] = p_data_in[(i-border)*DATA_SIZE];
				};
				if (j>=size_replicated-border){  // right columns
					p_data_out[(i*size_replicated)+j] = p_data_in[(i-border)*DATA_SIZE + (DATA_SIZE-1)];
				};
				if (j>=border && j<size_replicated-border){ // central matrix
					p_data_out[(i*size_replicated)+j] = p_data_in[(i-border)*DATA_SIZE+(j-border)];
				};
			};

		};
	};
	
};


/**************************************************/
/*              ESD predict function              */
/**************************************************/
uint32_t esd_predict(float32_t * p_data_X){
	
	//*****************************************************************************
// declare variables
//*****************************************************************************

		uint32_t predicted_class = 10;
		float32_t max_score_value;
		float32_t acc_score[N_CLASS] = {0,0,0,0,0,0,0};	// Accumulated score: in each cycle the lda_posterior output is added to lda_posterior of the previous cycles
		float32_t acc_score_tmp; // Temporary acc_score
		float32_t mean_score[N_CLASS]; // Averaged posterior distributions
		float32_t multivar_norm;
		float32_t multivar_norm_prior[N_CLASS]; // == lda_posterior == output of lda_predict
		
		/*Buffer for intermediate matrix results*/
		float32_t buf_SubX[SUB_DIM]; // Subspace X for each cycle
		float32_t buf_Mu[N_CLASS * SUB_DIM]; // Buffer for temporary Mu for each cycle
		float32_t buf_Sigma_pinv[SUB_DIM * SUB_DIM]; // Buffer for temporary Sigma_pinv for each cycle

		float32_t buf_Mu_k[SUB_DIM]; // Buffer for Mu(k)
		float32_t buf_SubX_Mu_k[SUB_DIM]; // Buffer for (X-Mu(k))
		float32_t buf_SubX_Mu_k_T[SUB_DIM]; // Buffer for (X-Mu(k))' which is the transpose of (X-Mu(k))
		float32_t buf_First_mat_mult[SUB_DIM]; // Buffer for first part of (X-Mu(k))*Sigma_pinv*(X-Mu(k))'
		float32_t buf_exponent;

		// Matrix instances
		arm_matrix_instance_f32 mat_buf_Sigma_pinv; // Matrix instance for buf_Sigma_pinv
		arm_matrix_instance_f32 mat_buf_SubX_Mu_k; // Matrix instance for buf_SubX_Mu_k
		arm_matrix_instance_f32 mat_buf_SubX_Mu_k_T; // Matrix instance for buf_SubX_Mu_k_T
		arm_matrix_instance_f32 mat_buf_First_mat_mult; // Matrix instance for first part of matrix multiplication
		arm_matrix_instance_f32 mat_buf_exponent;

		// for loops indeces
		uint16_t c; // Index from 1 to TOT_CYCLES
		uint32_t i_sub_idx; // Index of the element in X to be copied into subX in the subspace
		uint32_t i_helper; // i_helper to calculate the index, i to loop through N_CLASS*SUB_DIM to take Mu of one cycle
		uint16_t k, i; // Index from 1 to N_CLASS
		

		// ESD predict
		for (c=0; c<TOT_CYCLES; c++){
			
//			am_util_stdio_printf("Cycle # %i\n", c+1);
			
			// Copy from data_X to buf_subX for each cycle
			for (i=0; i<SUB_DIM; i++){
				i_helper = c*SUB_DIM + i; // Calculate the index in of the element in subspace_idx to be taken for this cycle
				i_sub_idx = subspace_idx[i_helper]; // Take the element in subspace_idx which is the index in X to be copied into subX
        buf_SubX[i] = p_data_X[i_sub_idx]; // Copy subX from X according to subspace_idx
			}

			// Copy from data_Mu to buf_Mu for this cycle
			for (i=0; i<N_CLASS*SUB_DIM; i++){
				i_helper = c*N_CLASS*SUB_DIM + i; // take the elements from c*N_CLASS*SUB_DIM to (c+1)*N_CLASS*SUB_DIM-1 in data_Mu
				buf_Mu[i] = data_Mu[i_helper];
			}

			// Copy from data_Sigma_pinv to buf_Sigma_pinv for this cycle and initialize the matrix
			for (i=0; i<SUB_DIM*SUB_DIM; i++){
				i_helper = c*SUB_DIM*SUB_DIM + i; // take the elements from c*SUB_DIM*SUB_DIM to (c+1)*SUB_DIM*SUB_DIM-1 in data_Sigma_pinv
				buf_Sigma_pinv[i] = data_Sigma_pinv[i_helper];
			}
			arm_mat_init_f32(&mat_buf_Sigma_pinv, SUB_DIM, SUB_DIM, (float32_t *)buf_Sigma_pinv); // Initialize matrix buf_Sigma_pinv with rows, columns, and data array
			
			
			/********************************************************************/
			/* lda_predict to get the lda_posterior distribution for this cycle */
			/********************************************************************/
			
			for (k=0; k<N_CLASS; k++){
			
				// Create Mu(k) copying from buf_Mu to buf_Mu_k for each class
				for (i=0; i<SUB_DIM; i++){
					i_helper = k*SUB_DIM + i;
					buf_Mu_k[i] = buf_Mu[i_helper];
				}

				// SubX - Mu(k)
				arm_sub_f32((float32_t *)buf_SubX, (float32_t *)buf_Mu_k, (float32_t *)buf_SubX_Mu_k, SUB_DIM);

				// Matrix initialization
				arm_mat_init_f32(&mat_buf_SubX_Mu_k, 1, SUB_DIM, (float32_t *)buf_SubX_Mu_k); // Initialize matrix (SubX-Mu(k))
				arm_mat_init_f32(&mat_buf_SubX_Mu_k_T, SUB_DIM, 1, buf_SubX_Mu_k_T); // Initialize transpose matrix of (SubX-Mu(k))
				arm_mat_init_f32(&mat_buf_First_mat_mult, 1, SUB_DIM, buf_First_mat_mult); // Initialize matrix for the first part of matrix multiplication
				arm_mat_init_f32(&mat_buf_exponent, 1, 1, &buf_exponent);
				
				// transpose of (SubX - Mu(k))
				arm_mat_trans_f32(&mat_buf_SubX_Mu_k, &mat_buf_SubX_Mu_k_T);
				
				// Matrix multiplication: (SubX-Mu(k))*buf_Sigma_pinv*(SubX-Mu(k))'
				arm_mat_mult_f32(&mat_buf_SubX_Mu_k, &mat_buf_Sigma_pinv, &mat_buf_First_mat_mult);
				arm_mat_mult_f32(&mat_buf_First_mat_mult, &mat_buf_SubX_Mu_k_T, &mat_buf_exponent);
				
				multivar_norm = ln_k_det[c] - buf_exponent/2.0f;

				multivar_norm_prior[k] = multivar_norm * PRIOR;

			}
			
			// Accumulate the posterior distribution (called here acc_score) cycle by cycle
			// by adding the lda_posterior to lda_posterior of previous cycles
			for (i=0; i<N_CLASS; i++){
				acc_score_tmp = acc_score[i];
				acc_score[i] = acc_score_tmp + multivar_norm_prior[i]; //lda_posterior[i];
			}

		}
		
		// compute the mean score
		arm_scale_f32((float32_t *)acc_score, TOT_CYCLES_1, (float32_t *)mean_score, N_CLASS); //mean_score[k] = acc_score[k] / TOT_CYCLES;

		// find the index of the max element in the mean score vector --> it's the predicted class ([value pred_class] = max(mean_score);)
    arm_max_f32((float32_t *)mean_score, N_CLASS, &max_score_value, &predicted_class);
		
		return predicted_class +=1;
		
	}
