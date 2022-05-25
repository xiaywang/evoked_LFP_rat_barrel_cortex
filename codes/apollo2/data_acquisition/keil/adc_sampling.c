#include "am_mcu_apollo.h"
#include "am_bsp.h"
#include "am_util.h"
#include "adc_sampling.h"


//*****************************************************************************
//
// Start up the ITM interface.
//
//*****************************************************************************
void itm_start(void)
{
    //
    // Initialize the printf interface for ITM/SWO output.
    //
    am_util_stdio_printf_init((am_util_stdio_print_char_t) am_bsp_itm_string_print);

    //
    // Initialize the SWO GPIO pin
    //
    am_bsp_pin_enable(ITM_SWO);

    //
    // Enable the ITM.
    //
    am_hal_itm_enable();

    //
    // Enable debug printf messages using ITM on SWO pin
    //
    am_bsp_debug_printf_enable();
}


//*****************************************************************************
//
// Set up the core for sleeping, and then go to sleep.
//
//*****************************************************************************
void sleep(void)
{
    //
    // Disable things that can't run in sleep mode.
    //
    am_bsp_debug_printf_disable();

    //
    // Go to Deep Sleep.
    //
    am_hal_sysctrl_sleep(AM_HAL_SYSCTRL_SLEEP_DEEP);

    //
    // Re-enable peripherals for run mode.
    //
    am_bsp_debug_printf_enable();
}

//*****************************************************************************
//
// Configure the ADC.
//
//*****************************************************************************
void adc_config(void)
{
    am_hal_adc_config_t sADCConfig;

    //
    // Enable the ADC power domain.
    //
    am_hal_pwrctrl_periph_enable(AM_HAL_PWRCTRL_ADC);

    //
    // Set up the ADC configuration parameters. These settings are reasonable
    // for accurate measurements at a low sample rate.
    //
    sADCConfig.ui32Clock = AM_HAL_ADC_CLOCK_HFRC;
    sADCConfig.ui32TriggerConfig = AM_HAL_ADC_TRIGGER_SOFT;
    sADCConfig.ui32Reference = AM_HAL_ADC_REF_INT_2P0;
    sADCConfig.ui32ClockMode = AM_HAL_ADC_CK_LOW_POWER;
    sADCConfig.ui32PowerMode = AM_HAL_ADC_LPMODE_1;
    sADCConfig.ui32Repeat = AM_HAL_ADC_REPEAT;
    am_hal_adc_config(&sADCConfig);

    //
    // For this example, the samples will be coming in slowly. This means we
    // can afford to wake up for every conversion.
    //
    am_hal_adc_int_enable(AM_HAL_ADC_INT_CNVCMP);

    //
    // Set up an ADC slot
    //
    am_hal_adc_slot_config(0, AM_HAL_ADC_SLOT_AVG_1 |
                              AM_HAL_ADC_SLOT_14BIT |
                              AM_HAL_ADC_SLOT_CHSEL_SE0 |
                              AM_HAL_ADC_SLOT_ENABLE);
    //
    // Enable the ADC.
    //
    am_hal_adc_enable();
}

//*****************************************************************************
//
// Initialize the ADC repetitive sample timer A3.
//
//*****************************************************************************
void init_timerA3_for_ADC(void)
{
    //
    // Start a timer to trigger the ADC periodically (1 second).
    //
    am_hal_ctimer_config_single(3, AM_HAL_CTIMER_TIMERA,
                                   AM_HAL_CTIMER_HFRC_12MHZ |
                                   AM_HAL_CTIMER_FN_REPEAT |
                                   AM_HAL_CTIMER_INT_ENABLE |
                                   AM_HAL_CTIMER_PIN_ENABLE);

    am_hal_ctimer_int_enable(AM_HAL_CTIMER_INT_TIMERA3);

    am_hal_ctimer_period_set(3, AM_HAL_CTIMER_TIMERA, 74, 0);

    //
    // Enable the timer A3 to trigger the ADC directly
    //
    am_hal_ctimer_adc_trigger_enable();

    //
    // Start the timer.
    //
    am_hal_ctimer_start(3, AM_HAL_CTIMER_TIMERA);
}


