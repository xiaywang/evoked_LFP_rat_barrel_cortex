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

// Extract the sd of 16-by-16 matrix in time
void feature_extraction_sd1(float32_t * p_data_in, float32_t * p_feature_out);

/* Features for entire half-matrix */

// Extract the mean of half-matrix in time
void feature_extraction_mean2(float32_t * p_data_in, float32_t * p_feature_out);

// Extract the sd of half-matrix in time
void feature_extraction_sd2(float32_t * p_data_in, float32_t * p_feature_out);

// Extract the grad of 16-by-16 matrix in time
void feature_extraction_grad2(float32_t * p_data_in, float32_t * p_feature_out);

/* Features for entire 1/3 matrix */

// Extract the mean of 1/3 matrix in time
void feature_extraction_mean3(float32_t * p_data_in, float32_t * p_feature_out);

// Extract the sd of 1/3 matrix in time
void feature_extraction_sd3(float32_t * p_data_in, float32_t * p_feature_out);

// Extract the grad of 1/3 matrix in time
void feature_extraction_grad3(float32_t * p_data_in, float32_t * p_feature_out);

/**************************************************/
/*     Cross-correlation of subsequent frames     */
/**************************************************/

// Take in input the matrices A and B (A smaller than B) to be cross-correlated (shifting B), output is the cross-correlation matrix (array form)
void feature_extraction_xcorr2_2frames(float32_t * p_data_in_A, float32_t * p_data_in_B, float32_t * p_xcorr2_out);

// Extract the x and y component of cross-correlation vector in the matrix of cross-correlation of two frames A and B
void feature_extraction_xcorr2_2frames_vector_components(float32_t * p_data_in_A, float32_t * p_data_in_B, int16_t * components);

// Extract the x and y component of cross-correlation vector of 16-by-16 matrices in time
void feature_extraction_xcorr2_vector_components(float32_t * p_data_in, float32_t * p_feature_out);

//
// Feature Extraction - combined
//

// Extract the mean, sd, grad of 1/3-matrices, concatenated according to time,
// i.e. mean3(t_0), sd3(t_0), grad3(t_0), mean3(t_1), sd3(t_1), grad3(t_1), mean3(t_2), and so on...
void feature_extraction_mean3_sd3_grad3(float32_t * p_data_in, float32_t * p_feature_out);
// mean3 sd3 grad3 and xcorr
void feature_extraction_mean3_sd3_grad3_xcorr(float32_t * p_data_in, float32_t * p_feature_out);

#endif // FEATURE_EXTRAXTION_H_
