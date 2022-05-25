//*****************************************************************************
//
//! @file adc_lpmode0.c
//!
//! @brief Example that takes samples with the ADC at high-speed.
//!
//! This example shows the CTIMER-A3 triggering repeated samples of an external
//! input at 1.2Msps in LPMODE0.  The example uses the CTIMER-A3 to trigger
//! ADC sampling.  Each data point is 128 sample average and is read from the
//! ADC FIFO into an SRAM circular buffer.
//!
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

#include "am_mcu_apollo.h"
#include "am_bsp.h"
#include "am_util.h"
#include "adc_sampling.h"


//*****************************************************************************
//
// Define a circular buffer to hold the ADC samples
//
//*****************************************************************************

#define N_SAMPLES 6400
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
				
				// Disable interrupts
        am_hal_interrupt_master_disable();
				
				// Reset the samples counter
				samples_counter = 0;
				
				// Set Bool_classifier to true
				Bool_classifier = true;
				
//				// Disable the ADC.
//				am_hal_adc_disable();
				
			}
    }
		

}


//*****************************************************************************
//
// Main function.
//
//*****************************************************************************
int
main(void)
{
    //
    // Set the system clock to maximum frequency, and set the default low-power
    // settings for this board.
    //
    am_hal_clkgen_sysclk_select(AM_HAL_CLKGEN_SYSCLK_MAX);
		am_hal_mcuctrl_bucks_enable();
    am_hal_vcomp_disable();

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
    // Print the banner.
    //
    am_util_stdio_terminal_clear();


    //
    // We are done printing. Disable debug printf messages on ITM.
    //
    am_bsp_debug_printf_disable();

    //
    // Loop forever.
    //
    while(1)
    {
			
//			//
//			// Disable interrupts
//			//
//			am_hal_interrupt_master_disable();

			//
			// Put the core to sleep.
			//
			sleep();
		
//			//
//			// Enable interrupts.
//			//
//			am_hal_interrupt_master_enable();

			
			
			if (Bool_classifier == true){

				am_util_stdio_printf("the last value is %f\n", SampleBuffer*3.3/16384.0);
				

				
				// Set Bool_classifier to false
				Bool_classifier = false;
				
				
				
				am_util_delay_ms(1000);
				
//				//
//				// Configure the ADC
//				//
//				adc_config();

//				//
//				// Trigger the ADC sampling for the first time manually.
//				//
//				am_hal_adc_trigger();
				
				// Enable interrupts.
        am_hal_interrupt_master_enable();
			}
			
			

    }
}
