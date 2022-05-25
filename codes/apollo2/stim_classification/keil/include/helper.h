#ifndef HELPER_H_   /* Include guard */
#define HELPER_H_

#include "apollo2_device.h"
#include "arm_math.h"

//#define FULL_T 69  // time width after cropping in time
//#define T_START 8  // to crop the data in time
//#define T_END 45  // to crop the data in time
#define DELTA_T 38  // time width after cropping in time
#define DATA_SIZE 16  // x width = y width of data matrix (this case 16-by-16 matrix)
#define SIGMA 0.6  // sigma of gaussian kernel for preprocessing
#define FILTER_SIZE 5 // Odd number, the filter kernel size is (uint16_t) (2*ceil(2*SIGMA)+1), I'm using only Stack memory.

#define PR_CROP 15 // to find Positive Rebound, start from timepoint 15
#define RP_CROP 4 // to find Response Peak, start from timepoint 4
#define t_window 4 // time window to find the response onset
#define frame_rate 1.602564102 // sampling rate between consecutive frames


/**************************************************/
/*                 Preprocessing                  */
/**************************************************/

// Crop data in time
//void preprocessing_crop_data(float32_t * p_data_in, float32_t * p_data_out, uint32_t t_start, uint32_t t_end);

// Gaussian smoothing
void preprocessing_gaussian_smooth(float32_t * p_data_in, float32_t * p_data_out); //, float32_t sigma
// Gaussian filter: taken sigma in input, it computes the gaussian kernel of size 2*ceil(2*sigma)+1
void preprocessing_gaussian_filter(float32_t * p_gaussian_filter);
// Replicate padding
void preprocessing_replicate_padding(float32_t * p_data_in, float32_t * p_data_out);
// Gaussian smoothing for single frame (single matrix)
void preprocessing_gaussian_smooth_single_frame(float32_t * p_data_in, float32_t * p_gaussian_filter, float32_t * p_data_out);


/**************************************************/
/*              ESD predict function              */
/**************************************************/
uint32_t esd_predict(float32_t * p_data_X);

#endif // HELPER_H_
