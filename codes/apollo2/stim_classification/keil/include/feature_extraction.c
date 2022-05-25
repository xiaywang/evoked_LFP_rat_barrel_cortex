#include "feature_extraction.h"
#include "helper.h"
#include "am_util.h"
#include <stdlib.h>

//
// Feature Extraction
//

/**************************************************/
/* Mean, standard deviation (sd), gradient (grad) */
/**************************************************/

/* Features for entire 16-by-16 matrix */

// Extract the mean of 16-by-16 matrix in time
void feature_extraction_mean1(float32_t * p_data_in, float32_t * p_feature_out){
	
	// Take into input the 16-by-16 matrix in time in array form
	// Output is the array of extracted mean of 16-by-16 matrix in time
	
	uint16_t t;
	uint32_t i, i_helper_top;
	float32_t data_tmp[DATA_SIZE*DATA_SIZE];
	
	for (t=0; t<DELTA_T; t++){
		i_helper_top = t*DATA_SIZE*DATA_SIZE;
		// copy the flattened 16-by-16 matrix of time t into a temporary array
		for (i=0; i<DATA_SIZE*DATA_SIZE; i++){
			data_tmp[i] = p_data_in[i_helper_top+i];
		}
		// call the arm cmsis mean function using the temporary array
		arm_mean_f32(data_tmp, DATA_SIZE*DATA_SIZE, &p_feature_out[t]);
	}

};


/* Features for half-matrix */

// Extract the mean of half-matrix in time
void feature_extraction_mean2(float32_t * p_data_in, float32_t * p_feature_out){
	
	// Take into input the 16-by-16 matrix in time in array form
	// Output is the array of extracted mean of 16-by-8 matrices in time
	// (The order is: mean2_top_t0, mean2_bottom_t0, mean2_top_t1, mean2_bottom_t1, mean2_top_t2, etc)
	
	uint16_t t;
	uint32_t i, i_helper_top, i_helper_center;
	float32_t data_tmp_top[DATA_SIZE*(DATA_SIZE/2)], data_tmp_bottom[DATA_SIZE*(DATA_SIZE/2)];
	
	for (t=0; t<DELTA_T; t++){
		i_helper_top = t*DATA_SIZE*DATA_SIZE;
		i_helper_center = i_helper_top + DATA_SIZE*(DATA_SIZE/2);
		// copy the flattened 16-by-8 matrix of time t into a temporary array (top part)
		for (i=0; i<DATA_SIZE*(DATA_SIZE/2); i++){
			data_tmp_top[i] = p_data_in[i_helper_top+i];
		}
		// copy the flattened 16-by-8 matrix of time t into a temporary array (bottom part)
		for (i=0; i<DATA_SIZE*(DATA_SIZE/2); i++){
			data_tmp_bottom[i] = p_data_in[i_helper_center+i];
		}
		// call the arm cmsis mean function using the temporary arrays
		arm_mean_f32(data_tmp_top, DATA_SIZE*(DATA_SIZE/2), &p_feature_out[2*t]);
		arm_mean_f32(data_tmp_bottom, DATA_SIZE*(DATA_SIZE/2), &p_feature_out[2*t+1]);
	}
	
};

/**************************************************/
/*   New features for stimulation classification  */
/**************************************************/

/* Positive Rebound (PR) */

void feature_extraction_positive_rebound(float32_t * p_mean1, float32_t * p_feature_out, uint32_t *pIndex){

	// Extract the positive rebound. It is the maximum value in mV starting from the timepoint 15 after the stimulation onset.
	// Input is the mean LFP signal
	// Output is the PR in mV and index of PR
	

	float32_t mean1_tmp[DELTA_T-PR_CROP];
	uint16_t i_helper;
	uint32_t idx_crop = 0;

	
	for (i_helper=0; i_helper<DELTA_T-PR_CROP; i_helper++){
		mean1_tmp[i_helper] = p_mean1[PR_CROP+i_helper];
	}
	
	arm_max_f32(mean1_tmp, DELTA_T-PR_CROP, p_feature_out, &idx_crop);
	*pIndex = idx_crop + PR_CROP;
	
};


/* Response Peak */

void feature_extraction_response_peak(float32_t * p_mean1, uint32_t *p_i_PR, float32_t * p_feature_out, uint32_t *pIndex){

	// Extract the response peak. The Response Peak is found as the minimum value between timepoint 4 after the stimulation onset and the positive rebound.
	// Input is the mean LFP signal and index of Positive Rebound
	// Output is the RP in mV and index of RP
	

	float32_t *p_mean1_tmp;
	uint16_t i_helper;
	uint32_t idx_crop = 0;
	uint32_t size_mean1_tmp = *p_i_PR-RP_CROP+1;

	p_mean1_tmp = (float32_t*) malloc(size_mean1_tmp * sizeof(float32_t));
	
	for (i_helper=0; i_helper<size_mean1_tmp; i_helper++){
		p_mean1_tmp[i_helper] = p_mean1[RP_CROP+i_helper];
	}
	
	arm_min_f32(p_mean1_tmp, size_mean1_tmp, p_feature_out, &idx_crop);
	free(p_mean1_tmp);
	*pIndex = idx_crop + RP_CROP;
	
};

/* Response Onset */

void feature_extraction_response_onset(float32_t * p_mean1, uint32_t *p_i_RP, uint32_t t_window_tmp, float32_t * p_feature_out, uint32_t *pIndex){

	// Extract the response onset. The Response Onset is found as the first value between the maximum point between the stimulation onset and the response peak and the response peak
  //	which is followed by 4 (defined as constant t_window) consecutive points with decreasing values.
	// Input is the mean LFP signal, index of Response Peak, and time window
	// Output is the RO in mV and index of RO
	

	float32_t *p_mean1_tmp;
	uint16_t i_helper, count=0, check = 0;
	uint32_t start_t = 0;
	float32_t M; // index of maximum point between the stimulation onset and the response peak
	uint32_t size_mean1_tmp = *p_i_RP+1;

	p_mean1_tmp = (float32_t*) malloc(size_mean1_tmp * sizeof(float32_t));
	
	for (i_helper=0; i_helper<size_mean1_tmp; i_helper++){
		p_mean1_tmp[i_helper] = p_mean1[i_helper];
	}
	
	arm_max_f32(p_mean1_tmp, size_mean1_tmp, &M, &start_t);
	
	free(p_mean1_tmp);
	
	for (i_helper = start_t; i_helper < *p_i_RP; i_helper++){
		if (p_mean1[i_helper+1] < p_mean1[i_helper]){
			count++;
			if (count == t_window_tmp){
				check = 1;
				*p_feature_out = p_mean1[i_helper-t_window_tmp+1];
				*pIndex = i_helper-t_window_tmp+1;
				break;
			}
		} else {
			count=0;
		}
		if (check == 0 && t_window_tmp >= 1){
			feature_extraction_response_onset(p_mean1, p_i_RP, t_window_tmp-1, p_feature_out, pIndex);
		} else if(check == 0 && t_window_tmp < 1) {
			*p_feature_out = p_mean1[start_t];
			*pIndex = start_t;
		}
	}
	
};

/* Trapezoidal method */

void feature_extraction_trapz(float32_t * p_mean1, uint16_t BlockSize, float32_t * A){
	
	// Approximation of the integral via trapezoindal method.
	
	uint16_t i;
	
	for (i=0; i<BlockSize-1; i++){
		*A = *A + ((p_mean1[i]+p_mean1[i+1])/2.0);
	}
	
}

/* time-normalized LFP (tLFP) */

void feature_extraction_time_norm_LFP(float32_t * p_mean1, uint32_t *p_i_RO, float32_t * p_feature_out){

	// Extract the time-normalized LFP = AUC/RD with AUC = Area Under Curve and RD = Response Duration.
	// 
	
	float32_t AUC=0;
	uint16_t t = *p_i_RO + 1;
	uint16_t i_helper;
	
	// Find the first point after RO which has the value higher than RO
	while (t < DELTA_T && p_mean1[t] < p_mean1[*p_i_RO]){
		t++;
	}
	
	// Allocate memory for mean1_tmp to be given to compute AUC of the main response peak
	uint16_t BlockSize = t - *p_i_RO +1;
	float32_t *p_mean1_tmp;
	
	p_mean1_tmp = (float32_t*) malloc(BlockSize * sizeof(float32_t));
	
	for (i_helper=0; i_helper<BlockSize; i_helper++){
		p_mean1_tmp[i_helper] = p_mean1[*p_i_RO + i_helper];
	}
	// Compute AUC using trapezoindal method
	feature_extraction_trapz(p_mean1_tmp, BlockSize, &AUC);
	
	free(p_mean1_tmp);
	
	// Compute tLFP
	*p_feature_out = AUC/(float32_t)(t-*p_i_RO);
	
};


//
// Feature Extraction combined
//

void feature_extraction_mean2_RPA_PR_RPL_tLFP(float32_t * p_data_in, float32_t * p_feature_out){

	// Extract the mean2, RPA, PR, RPL, tLFP and concatenate them into one single array
	
	float32_t mean1[DELTA_T];
	float32_t RP, RO;
	uint32_t i_RP, i_PR, i_RO;

	feature_extraction_mean2(p_data_in, p_feature_out);
	
	// First: extract the mean LFP signal
	feature_extraction_mean1(p_data_in, mean1);
	
	// Extract the Positive Rebound (PR)
	feature_extraction_positive_rebound(mean1, &p_feature_out[2*DELTA_T+1], &i_PR);
	
	// Extract the Response Peak (RP)
	feature_extraction_response_peak(mean1, &i_PR, &RP, &i_RP);
	
	// Extract the Response Onset (RO)
	feature_extraction_response_onset(mean1, &i_RP, t_window, &RO, &i_RO);
	
	// Extract the Response Peak Amplitude (RPA)
	p_feature_out[2*DELTA_T] = RO - RP;
	
//	// Extract the Response Onset Latency (ROL)
//	ROL = (i_RO+1) * frame_rate;

	// Extract the Response Peak Latency (RPL)
	p_feature_out[2*DELTA_T+2] = (i_RP+1) * frame_rate;
	
	// Extract the time-normalized LFP (tLFP)	
	p_feature_out[2*DELTA_T+3] = 0;
	feature_extraction_time_norm_LFP(mean1, &i_RO, &p_feature_out[2*DELTA_T+3]);
	
};
