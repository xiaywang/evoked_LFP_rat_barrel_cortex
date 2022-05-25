#ifndef FEATURE_EXTRAXTION_H_   /* Include guard */
#define FEATURE_EXTRAXTION_H_

#include "apollo2_device.h"
#include "arm_math.h"


//
// Feature Extraction for each t
//

/**************************************************/
/* Mean, standard deviation (sd), gradient (grad) */
/**************************************************/

/* Features for entire 16-by-16 matrix */

// Extract the mean of 16-by-16 matrix in time
void feature_extraction_mean1(float32_t * p_data_in, float32_t * p_feature_out);

/* Features for half-matrix */

// Extract the mean of half-matrix in time
void feature_extraction_mean2(float32_t * p_data_in, float32_t * p_feature_out);


/**************************************************/
/*   New features for stimulation classification  */
/**************************************************/

/* Positive Rebound (PR) */

void positive_rebound(float32_t * p_mean1, float32_t * p_feature_out, uint32_t *pIndex);

/* Response Peak */

void response_peak(float32_t * p_mean1, uint32_t *p_i_PR, float32_t * p_feature_out, uint32_t *pIndex);

/* Response Onset */

void response_onset(float32_t * p_mean1, uint32_t *p_i_RP, uint32_t t_window_tmp, float32_t * p_feature_out, uint32_t *pIndex);

/* Trapezoidal method */

void trapz(float32_t * p_mean1, uint16_t BlockSize, float32_t * A);

/* time-normalized LFP (tLFP) */

void time_norm_LFP(float32_t * p_mean1, uint32_t *p_i_RO, float32_t * p_feature_out);


//
// Feature Extraction - combined
//

// Extract the mean2, RPA, PR, RPL, tLFP and concatenate them into one single array
void feature_extraction_mean2_RPA_PR_RPL_tLFP(float32_t * p_data_in, float32_t * p_feature_out);

#endif // FEATURE_EXTRAXTION_H_
