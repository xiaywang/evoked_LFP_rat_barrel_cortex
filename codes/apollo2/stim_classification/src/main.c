//*****************************************************************************
//
//! @file main.c
//!
//! @brief Layer classification
//!
//! This code performs the layer classification with ensemble subspace
//! discriminant (ESD) prediction after preprocessing and extracting the
//! features.
//
//*****************************************************************************

//*****************************************************************************
//
// Copyright (c) 2017, Ambiq Micro
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
// 
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived from this
// software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// This is part of revision 1.2.11 of the AmbiqSuite Development Package.
//
//*****************************************************************************

#include "apollo2_device.h"
#include "am_mcu_apollo.h"
#include "am_bsp.h"
#include "am_util.h"
#include "arm_math.h"
#include "config.h"
#include "data.h"
#include "helper.h"
#include "feature_extraction.h"


//*****************************************************************************
//
// Define a circular buffer to hold the ADC samples
//
//*****************************************************************************

#define N_SAMPLES 9728
uint32_t SampleBuffer;
uint32_t samples_counter = 0;
bool Bool_classifier = false;

//*****************************************************************************
//
// Interrupt handler for the ADC.
//
//*****************************************************************************
void am_adc_isr(void)
{
    uint32_t ui32Status;

    //
    // Read the interrupt status.
    //
    ui32Status = am_hal_adc_int_status_get(true);

    //
    // Clear the ADC interrupt.
    //
    am_hal_adc_int_clear(ui32Status);

	
    //
    // If we got a conversion completion interrupt (which should be our only
    // ADC interrupt), go ahead and read the data.
    //
    if (ui32Status & AM_HAL_ADC_INT_CNVCMP)
    {
			
			// Read the value from the FIFO into the circular buffer.
			SampleBuffer = AM_HAL_ADC_FIFO_SAMPLE(am_hal_adc_fifo_pop());
			
			// Increase the counter of number of samples
			samples_counter++;
//			am_util_stdio_printf("counter %i\n", samples_counter);
			
			// Check if the number of converted samples is equal to N_SAMPLES
			if (samples_counter >= N_SAMPLES){
				
//				// Disable interrupts
//        am_hal_interrupt_master_disable();
				
				// Reset the samples counter
				samples_counter = 0;
				
				// Set Bool_classifier to true
				Bool_classifier = true;
				
				// Disable the ADC.
				am_hal_adc_disable();
				
			}
    }
		

}


//*****************************************************************************
//
// Main
//
//*****************************************************************************
int main(void){
//*****************************************************************************
	
    //
    // Set the clock frequency.
    //
    am_hal_clkgen_sysclk_select(AM_HAL_CLKGEN_SYSCLK_MAX);
		am_hal_mcuctrl_bucks_enable();

    //
    // Set the default cache configuration
    //
    am_hal_cachectrl_enable(&am_hal_cachectrl_defaults);
	
		//
    // Allow the XTAL to turn off.
    //
    am_hal_clkgen_osc_stop(AM_HAL_CLKGEN_OSC_XT);

    //
    // Turn off the voltage comparator
    //
    am_hal_vcomp_disable();

    //
    // Configure the board for low power operation.
    //
    am_bsp_low_power_init();

    //
    // Start the ITM interface.
    //
    itm_start();
		
		//
    // Start the CTIMER A3 for timer-based ADC measurements.
    //
    init_timerA3_for_ADC();

    //
    // Enable interrupts.
    //
    am_hal_interrupt_enable(AM_HAL_INTERRUPT_ADC);
    am_hal_interrupt_master_enable();

    //
    // Set a pin to act as our ADC input
    //
    am_hal_gpio_pin_config(16, AM_HAL_PIN_16_ADCSE0);

    //
    // Configure the ADC
    //
    adc_config();

    //
    // Trigger the ADC sampling for the first time manually.
    //
    am_hal_adc_trigger();
		
		//
		// Config gpio pin to output to measure execution time
		//
		am_hal_gpio_pin_config(12, AM_HAL_GPIO_OUTPUT);
		

    /********************************************************************/
						/* ========== Stimulation classification ========== */
		/********************************************************************/
    am_util_stdio_terminal_clear();
    am_util_stdio_printf("--- Stimulation classification ---\n\n");
		
		float32_t stimulations[3] = { // Classes == stimulation amplitudes
			0.5, 1.0, 1.4
		};		

		/**************************************************/
		/*               Declare variables                */
		/**************************************************/

//		float32_t data_cropped[DATA_SIZE*DATA_SIZE*DELTA_T];
		float32_t data_smoothed[DATA_SIZE*DATA_SIZE*DELTA_T];
		float32_t features[(2*DELTA_T) + 4];
		uint32_t predicted_class = 0;

		
		//
    // We are done printing. Disable debug printf messages on ITM.
    //
    am_bsp_debug_printf_disable();
		
		while (1){
			
			//
			// Put the core to sleep.
			//
			sleep();
			
			if (Bool_classifier == true){
			
					/********************************************************************/
												 /* ========== CLASSIFIER ========== */
					/********************************************************************/

					/**************************************************/
					/*                 Preprocessing                  */
					/**************************************************/
					
					// Crop data in time
//					preprocessing_crop_data(data, data_cropped, T_START, T_END);		

// Set the output state high for one GPIO.
					am_hal_gpio_out_bit_set(12);				
						
					// Gaussian smoothing
					preprocessing_gaussian_smooth(data, data_smoothed);
						

					/**************************************************/
					/*               Feature Extraction               */
					/**************************************************/

					

					// Mean2 RPA PR RPL tLFP
					feature_extraction_mean2_RPA_PR_RPL_tLFP(data_smoothed, features);

					

					/**************************************************/
					/*                  ESD predict                   */
					/**************************************************/
					
					// Set the output state high for one GPIO.
					am_hal_gpio_out_bit_set(12);
					
					predicted_class = esd_predict(features);
					
					
					// Sets the output state to low for one GPIO.
					am_hal_gpio_out_bit_clear(12);
					
//					// Delay
//					am_util_delay_ms(500);
//					
//				}
					
					// print the predicted class
//					am_util_stdio_printf("\nThe predicted class is %f V\n", stimulations[predicted_class]);
					
					
					/********************************************************************/
											 /* ========== END CLASSIFIER ========== */
					/********************************************************************/
					
					
							// Set Bool_classifier to false
							Bool_classifier = false;
							
							// Delay
							am_util_delay_ms(500);
							
							//
							// reConfigure the ADC
							//
							adc_config();

							//
							// reTrigger the ADC sampling for the first time manually.
							//
							am_hal_adc_trigger();
					

				}
	}
					
					
					
//    //
//    // We are done printing. Disable debug printf messages on ITM.
//    //
//    am_bsp_debug_printf_disable();

//    //
//    // Loop forever while sleeping.
//    //
//    while (1)
//    {
//        //
//        // Go to Deep Sleep.
//        //
//        am_hal_sysctrl_sleep(AM_HAL_SYSCTRL_SLEEP_DEEP);
//    }
	
}
