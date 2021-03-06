#ifndef ADC_SAMPLING_H_   /* Include guard */
#define ADC_SAMPLING_H_



//*****************************************************************************
//
// Start up the ITM interface.
//
//*****************************************************************************
void itm_start(void);

//*****************************************************************************
//
// Set up the core for sleeping, and then go to sleep.
//
//*****************************************************************************
void sleep(void);

//*****************************************************************************
//
// Configure the ADC.
//
//*****************************************************************************
void adc_config(void);

//*****************************************************************************
//
// Initialize the ADC repetitive sample timer A3.
//
//*****************************************************************************
void init_timerA3_for_ADC(void);



#endif // HELPER_H_
