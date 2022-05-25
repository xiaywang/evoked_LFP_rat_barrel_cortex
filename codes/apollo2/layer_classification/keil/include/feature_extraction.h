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

/* Features for 1/3 matrix */

// Extract the mean of 1/3 matrix in time
void feature_extraction_mean3_std3_grad3(float32_t * p_data_in, float32_t * p_feature_out);


/**************************************************/
/*     Cross-correlation of subsequent frames     */
/**************************************************/

// Take in input the matrices A and B (A smaller than B) to be cross-correlated (shifting B), output is the cross-correlation matrix (array form)
void feature_extraction_xcorr2_2frames(float32_t * p_data_in_A, float32_t * p_data_in_B, float32_t * p_xcorr2_out);

// Extract the x and y component of cross-correlation vector in the matrix of cross-correlation of two frames A and B
void feature_extraction_xcorr2_2frames_vector_components(float32_t * p_data_in_A, float32_t * p_data_in_B, int16_t * components);

// Extract the x and y component of cross-correlation vector of 16-by-16 matrices in time
void feature_extraction_xcorr2_vector_components(float32_t * p_data_in, float32_t * p_feature_out);


#endif // FEATURE_EXTRAXTION_H_
