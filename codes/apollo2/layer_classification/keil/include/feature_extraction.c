#include "feature_extraction.h"
#include "helper.h"
#include "am_util.h"
//
// Feature Extraction
//

/**************************************************/
/* Mean, standard deviation (sd), gradient (grad) */
/**************************************************/

/* Features for 1/3 matrix */

// Extract the mean, std, and gradient of 1/3 matrix in time
void feature_extraction_mean3_std3_grad3(float32_t * p_data_in, float32_t * p_feature_out){
	
	// Take into input the 16-by-16 matrix in time in array form
	// Output is the extracted mean, std, grad of 1/3-matrices, concatenated according to time,
	// i.e. mean3(t_0), sd3(t_0), grad3(t_0), mean3(t_1), sd3(t_1), grad3(t_1), mean3(t_2), and so on...
	
	uint16_t t;
	uint32_t i, i_helper_top, i_helper_center, i_helper_center2;
	float32_t data_tmp_top[DATA_SIZE*(DATA_SIZE/3)], data_tmp_center[DATA_SIZE*(DATA_SIZE/3)], data_tmp_bottom[DATA_SIZE*(DATA_SIZE - 2*(DATA_SIZE/3))];
	float32_t mean3_top[DELTA_T], mean3_center[DELTA_T], mean3_bottom[DELTA_T], grad3_tmp_top[DELTA_T], grad3_tmp_bottom[DELTA_T];
	
	for (t=0; t<DELTA_T; t++){
		i_helper_top = t*DATA_SIZE*DATA_SIZE;
		i_helper_center = i_helper_top + DATA_SIZE*(DATA_SIZE/3);
		i_helper_center2 = i_helper_center + DATA_SIZE*(DATA_SIZE/3);
		
		// copy the flattened 16-by-5 matrix of time t into a temporary array (top part)
		for (i=0; i<DATA_SIZE*(DATA_SIZE/3); i++){
			data_tmp_top[i] = p_data_in[i_helper_top+i];
		}
		// copy the flattened 16-by-5 matrix of time t into a temporary array (central part)
		for (i=0; i<DATA_SIZE*(DATA_SIZE/3); i++){
			data_tmp_center[i] = p_data_in[i_helper_center+i];
		}
		// copy the flattened 16-by-6 matrix of time t into a temporary array (bottom part)
		for (i=0; i<DATA_SIZE*(DATA_SIZE - 2*(DATA_SIZE/3)); i++){
			data_tmp_bottom[i] = p_data_in[i_helper_center2+i];
		}
		// call the arm cmsis mean function using the temporary arrays and save the mean in temporary arrays to compute the gradient later
		arm_mean_f32(data_tmp_top, DATA_SIZE*(DATA_SIZE/3), &mean3_top[t]);
		arm_mean_f32(data_tmp_center, DATA_SIZE*(DATA_SIZE/3), &mean3_center[t]);
		arm_mean_f32(data_tmp_bottom, DATA_SIZE*(DATA_SIZE - 2*(DATA_SIZE/3)), &mean3_bottom[t]);

		// call the arm cmsis std function using the temporary arrays
		arm_std_f32(data_tmp_top, DATA_SIZE*(DATA_SIZE/3), &p_feature_out[8*t+3]);
		arm_std_f32(data_tmp_center, DATA_SIZE*(DATA_SIZE/3), &p_feature_out[8*t+4]);
		arm_std_f32(data_tmp_bottom, DATA_SIZE*(DATA_SIZE - 2*(DATA_SIZE/3)), &p_feature_out[8*t+5]);
	}
	// compute the gradient as mean_top - mean_bottom
	arm_sub_f32(mean3_top, mean3_center, grad3_tmp_top, DELTA_T);
	arm_sub_f32(mean3_center, mean3_bottom, grad3_tmp_bottom, DELTA_T);
	
	for (t=0; t<DELTA_T; t++){
		p_feature_out[8*t] = mean3_top[t];
		p_feature_out[8*t+1] = mean3_center[t];
		p_feature_out[8*t+2] = mean3_bottom[t];
		p_feature_out[8*t+6] = grad3_tmp_top[t];
		p_feature_out[8*t+7] = grad3_tmp_bottom[t];
	}
	
};

/**************************************************/
/*     Cross-correlation of subsequent frames     */
/**************************************************/
void feature_extraction_xcorr2_2frames(float32_t * p_data_in_A, float32_t * p_data_in_B, float32_t * p_xcorr2_out){
	
	// Input: matrices A and B (A smaller than B) in array form to be cross-correlated (shifting B) and their size (n_col = n_row)
	// N.B. Here we deal with square matrix, thus n_col = n_row
	// Output is the cross-correlation matrix in array form
	
	// The idea is to perform 2D cross-correlation using 1D cross-correlation of rows
	
	int16_t k=0;
	int16_t i_helper;
	int16_t i, j; // j index for B, i index for A (index of row of matrix form)
	float32_t xcorr_tmp[SIZE_XCORR_ROW];
	float32_t xcorr_add_tmp = 0;
	float32_t vector_A_tmp[W];
	float32_t vector_B_tmp[DATA_SIZE];

	
	// Initialize p_xcorr2_out to 0
	for (i_helper=0; i_helper<SIZE_XCORR_ROW*SIZE_XCORR_ROW; i_helper++){
		p_xcorr2_out[i_helper] = 0;
	};
	
	
	// Top part of cross-correlation matrix. Shifting B upwards, until the index of A and B are the same.
	for (j=DATA_SIZE-1; j >= 0; j--){
		
		if ((DATA_SIZE-1)-j > W-1){
			for (i=0; i<W; i++){

				// Prepare the vectors to use arm_correlate_f32
				for (i_helper=0; i_helper<W; i_helper++){
					vector_A_tmp[i_helper] = p_data_in_A[(i*W)+i_helper];
				}
				for (i_helper=0; i_helper<DATA_SIZE; i_helper++){
					vector_B_tmp[i_helper] = p_data_in_B[((j+i)*DATA_SIZE)+i_helper];
				}
				// 1D cross-correlation
				arm_correlate_f32(vector_A_tmp, W, vector_B_tmp, DATA_SIZE, xcorr_tmp);
				
				// Add the 1D cross-correlation of rows at every iteration of i
				for (i_helper=0; i_helper<SIZE_XCORR_ROW; i_helper++){
					xcorr_add_tmp = p_xcorr2_out[(k*SIZE_XCORR_ROW)+i_helper];
					p_xcorr2_out[(k*SIZE_XCORR_ROW)+i_helper] = xcorr_add_tmp + xcorr_tmp[i_helper];
				};
				
			}
		}
		else {
			
			for (i=0; i<=(DATA_SIZE-1)-j; i++){
				// Prepare the vectors to use arm_correlate_f32
				for (i_helper=0; i_helper<W; i_helper++){
					vector_A_tmp[i_helper] = p_data_in_A[(i*W)+i_helper];
				}
				for (i_helper=0; i_helper<DATA_SIZE; i_helper++){
					vector_B_tmp[i_helper] = p_data_in_B[((j+i)*DATA_SIZE)+i_helper];
				}
				arm_correlate_f32(vector_A_tmp, W, vector_B_tmp, DATA_SIZE, xcorr_tmp);
				
				// Add the 1D cross-correlation of rows at every iteration of i
				for (i_helper=0; i_helper<SIZE_XCORR_ROW; i_helper++){
					xcorr_add_tmp = p_xcorr2_out[(k*SIZE_XCORR_ROW)+i_helper];
					p_xcorr2_out[(k*SIZE_XCORR_ROW)+i_helper] = xcorr_add_tmp + xcorr_tmp[i_helper];
				};
			}
		};
		k++;
	};
	
	// Bottom part of cross-correlation matrix. Shifting B upwards, when index of B is smaller than index of A.
	for (j=W-2; j>=0; j--){
		for (i=W-1; i>=W-1-j; i--){
			
			// Prepare the vectors to use arm_correlate_f32
			for (i_helper=0; i_helper<W; i_helper++){
				vector_A_tmp[i_helper] = p_data_in_A[(i*W)+i_helper];
			}
			for (i_helper=0; i_helper<DATA_SIZE; i_helper++){
				vector_B_tmp[i_helper] = p_data_in_B[((i-(W-1-j))*DATA_SIZE)+i_helper];
			}
			// 1D cross-correlation
			arm_correlate_f32(vector_A_tmp, W, vector_B_tmp, DATA_SIZE, xcorr_tmp);
			
			// Add the 1D cross-correlation of rows at every iteration of i
			for (i_helper=0; i_helper<SIZE_XCORR_ROW; i_helper++){
				xcorr_add_tmp = p_xcorr2_out[(k*SIZE_XCORR_ROW)+i_helper];
				p_xcorr2_out[(k*SIZE_XCORR_ROW)+i_helper] = xcorr_add_tmp + xcorr_tmp[i_helper];
			};
		};
		k++;
	};
	
};

void feature_extraction_xcorr2_2frames_vector_components(float32_t * p_data_in_A, float32_t * p_data_in_B, int16_t * components){
	
	// Extract the x and y component of cross-correlation vector in the matrix of cross-correlation of two frames A and B
	// Input: matrices A and B to be cross-correlated
	// Output: array with the two vector components of the cross-correlation vector
	
	float32_t xcorr2[SIZE_XCORR_ROW*SIZE_XCORR_ROW];
	uint16_t center;
	uint32_t max_index;
	float32_t max_value;
	uint16_t max_index_col, max_index_row;
	
	// Calculate the cross-correlation map of A and B
	feature_extraction_xcorr2_2frames(p_data_in_A, p_data_in_B, xcorr2);
	
	// Center point of the xcorr map
	center = SIZE_XCORR_ROW/2;
	
	// Find the array-index of the max cross-correlation value
	arm_max_f32(xcorr2, SIZE_XCORR_ROW*SIZE_XCORR_ROW, &max_value, &max_index);
	
	// Corresponding index in matrix
	max_index_col = max_index%SIZE_XCORR_ROW;
	max_index_row = ceil(max_index/SIZE_XCORR_ROW);
	
	// horizontal (x) and vertical (y) of the vector connecting the central point to the maximum point in the xcorr map
	components[0] = max_index_col - center; // x component
	components[1] = center - max_index_row; // y component
	
};

void feature_extraction_xcorr2_vector_components(float32_t * p_data_in, float32_t * p_feature_out){
	
	// Extract the x and y component of cross-correlation vector of 16-by-16 matrices in time
	// Input the 16-by-16 matrix in time in array form
	// Output: the extracted cross-correlation vector components x and y in the following order:
	// x_xcorr_t0, x_xcorr_t1, x_xcorr_t2, ..., x_xcorr_t_end, y_xcorr_t0, y_xcorr_t1, y_xcorr_t2, etc.
	
	uint16_t t, i_row, i_col;
	uint16_t crop = (DATA_SIZE-W)/2;
	int16_t components[2];
	float32_t vector_A[W*W];
	float32_t vector_B[DATA_SIZE*DATA_SIZE];
	uint32_t i_helper;
	
	
	for (t=0; t<DELTA_T-1; t++){
		
		// copy matrix(t) at time t to vector_A
		i_helper = 0;
		for (i_row=crop; i_row<DATA_SIZE-crop; i_row++){
			for (i_col=0; i_col<W; i_col++){
				vector_A[i_helper] = p_data_in[(t*DATA_SIZE*DATA_SIZE)+(i_row*DATA_SIZE)+crop+i_col];
				i_helper++;
			};
		};
		// copy matrix(t+1) at time t+1 to vector_B
		for (i_helper=0; i_helper<DATA_SIZE*DATA_SIZE; i_helper++){
			vector_B[i_helper] = p_data_in[((t+1)*DATA_SIZE*DATA_SIZE)+i_helper];
		};
		
		// Extract the x and y components of the cross-correlation map of frame t and frame t+1
		feature_extraction_xcorr2_2frames_vector_components(vector_A, vector_B, components);
		
		// Save the vector components to p_feature_out according to the order explained in the caption of this function (above)
		p_feature_out[t] = (float32_t) components[0];
		p_feature_out[(DELTA_T-1) + t] = (float32_t) components[1]; 
		
	};
	
};
