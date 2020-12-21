//+------------------------------------------------------------------+
//|                                        ehlers_roofing_filter.mq5 |
//|         Copyright © 2013 by John F. Ehlers. All rights reserved. |
//|                             Converted to MQL5 by thetestspecimen |
//|                                  https://www.thetestspecimen.com |
//|                               https://github.com/thetestspecimen |
//+------------------------------------------------------------------+

#property copyright "2013, John F. Ehlers"
#property link "https://github.com/thetestspecimen"
#property version "1.00"

//--------------------------------------------------------------------

#property indicator_separate_window
#property indicator_buffers 3
#property indicator_plots 2

//Plot Line
#property indicator_label1 "Roofing Filter"
#property indicator_type1 DRAW_LINE
#property indicator_color1 clrGoldenrod
#property indicator_style1 STYLE_SOLID
#property indicator_width1 2

//Plot Signal Line
#property indicator_label2  "Roofing Filter Signal"
#property indicator_type2  DRAW_LINE
#property indicator_color2  clrDodgerBlue
#property indicator_style2  STYLE_DASH
#property indicator_width2  2

// set enumerated values for prices
enum enPrices {
    pr_close, // Close
    pr_open, // Open
    pr_high, // High
    pr_low, // Low
    pr_median, // Median
    pr_typical, // Typical
    pr_weighted, // Weighted
    pr_average, // Average (high+low+open+close)/4
    pr_medianbody, // Average median body (open+close)/2
    pr_trendbiased, // Trend biased price
    pr_ha_close, // Heiken ashi close
    pr_ha_open, // Heiken ashi open
    pr_ha_high, // Heiken ashi high
    pr_ha_low, // Heiken ashi low
    pr_ha_median, // Heiken ashi median
    pr_ha_typical, // Heiken ashi typical
    pr_ha_weighted, // Heiken ashi weighted
    pr_ha_average, // Heiken ashi average
    pr_ha_medianbody, // Heiken ashi median body
    pr_ha_trendbiased // Heiken ashi trend biased price
};

// inputs
input enPrices Price = pr_close; // Price Type
input int hpPeriod = 80;        // High Pass Period
input int lpPeriod = 40;        // Low Pass Period

// buffers
double result[]; 	// result buffer
double signal[];    // signal line buffer
double highPass[]; 	// highpass buffer

//////////////////////////////////////////////////////////////////

int OnInit() {
    SetIndexBuffer(0, result, INDICATOR_DATA);
    SetIndexBuffer(1, signal, INDICATOR_DATA);
    SetIndexBuffer(2, highPass, INDICATOR_CALCULATIONS);
    return (0);
}

//////////////////////////////////////////////////////////////////

double work[];

int OnCalculate(
    const int       rates_total, 
    const int       prev_calculated,
    const datetime& time[],
    const double&   open[],
    const double&   high[],
    const double&   low[],
    const double&   close[],
    const long&     tick_volume[],
    const long&     volume[],
    const int&      spread[]) {

    if (Bars(_Symbol, _Period) < rates_total) {
        return (-1);
    }

    // resize the working array to fit the added data on each new bar
    if (ArrayRange(work, 0) != rates_total) {
        ArrayResize(work, rates_total);
    }

    /* most of the time prev_calculated == rates_total == _bars in which case this loop only calculates
    one bar of data. However, additional history or a refresh of data will cause prev_calculated
    to reset to 0, and hence cause a loop through all previous bars */
    for (int i = (int) MathMax(prev_calculated - 1, 0); i < rates_total; i++) {

        // get the price of the current bar based on the Price variable input
        work[i] = getPrice(Price, open, close, high, low, i, rates_total);

        // highpass filter
        highPass[i] = processHighPass(hpPeriod, work, highPass, i);

        // lowpass filter (supersmoother)
        result[i] = processLowPass(lpPeriod, highPass, result, i);

        if(i>3){
         signal[i] = result[i-2];   
        }
        
    }

    // return value of prev_calculated for next call
    return (rates_total);
}


// high pass filter
double processHighPass(double highPassPeriod, double& inputArr[], double& outputArr[], int i) {

    double pi = 3.14159265358979323846264338327950288;
    double alpha1 = (MathCos(MathSqrt(2.0) * pi / highPassPeriod) + MathSin(MathSqrt(2.0) * pi / (double)highPassPeriod) - 1.0) / MathCos(MathSqrt(2.0) * pi / (double)highPassPeriod);
    if (i < 3) {
        return ((1.0 - alpha1 / 2.0) * (1.0 - alpha1 / 2.0) * (inputArr[i]));
    } else {
        return ((1.0 - alpha1 / 2.0) * (1.0 - alpha1 / 2.0) * (inputArr[i] - 2.0 * inputArr[i - 1] + inputArr[i - 2]) + 2.0 * (1.0 - alpha1) * outputArr[i - 1] - (1.0 - alpha1) * (1.0 - alpha1) * outputArr[i - 2]);
    }
}

// low pass filter (supersmoother)
double processLowPass(int lowPassPeriod, double& inputArr[], double& outputArr[], int i) {

    double pi = 3.14159265358979323846264338327950288;
    double a1 = MathExp(-1.0 * MathSqrt(2.0) * pi / (double)lowPassPeriod);
    double b1 = 2.0 * a1 * MathCos(MathSqrt(2.0) * pi / (double)lowPassPeriod);
    double c2 = b1;
    double c3 = -a1 * a1;
    double c1 = 1.0 - c2 - c3;
    if (i < 3) {
        return (c1 * inputArr[i]);
    } else {
        return (c1 * (inputArr[i] + inputArr[i - 1]) / 2.0 + c2 * outputArr[i - 1] + c3 * outputArr[i - 2]);
    }
}

// this method gets allows the selection of which price to use 
// it also includes Heiken Ashi price outputs
double workHa[][4];
double getPrice(int priceType, const double& open[], const double& close[], const double& high[], const double& low[], int i, int bars) {

    if (priceType >= pr_ha_close) {

        if (ArrayRange(workHa, 0) != bars) {
            ArrayResize(workHa, bars);
        }

        double haOpen;
        if (i > 0) {
            haOpen = (workHa[i - 1][2] + workHa[i - 1][3]) / 2.0;
        } else {
            haOpen = (open[i] + close[i]) / 2;
        }

        double haClose = (open[i] + high[i] + low[i] + close[i]) / 4.0;
        double haHigh = MathMax(high[i], MathMax(haOpen, haClose));
        double haLow = MathMin(low[i], MathMin(haOpen, haClose));

        workHa[i][0] = haLow;
        workHa[i][1] = haHigh;
        workHa[i][2] = haOpen;
        workHa[i][3] = haClose;

        switch (priceType) {
        case pr_ha_close:       return (haClose);
        case pr_ha_open:        return (haOpen);
        case pr_ha_high:        return (haHigh);
        case pr_ha_low:         return (haLow);
        case pr_ha_median:      return ((haHigh + haLow) / 2.0);
        case pr_ha_medianbody:  return ((haOpen + haClose) / 2.0);
        case pr_ha_typical:     return ((haHigh + haLow + haClose) / 3.0);
        case pr_ha_weighted:    return ((haHigh + haLow + haClose + haClose) / 4.0);
        case pr_ha_average:     return ((haHigh + haLow + haClose + haOpen) / 4.0);
        case pr_ha_trendbiased:
            if (haClose > haOpen) {
                return ((haHigh + haClose) / 2.0);
            } else {
                return ((haLow + haClose) / 2.0);
            }
        }

    } else {

        switch (priceType) {
        case pr_close:      return (close[i]);
        case pr_open:       return (open[i]);
        case pr_high:       return (high[i]);
        case pr_low:        return (low[i]);
        case pr_median:     return ((high[i] + low[i]) / 2.0);
        case pr_medianbody: return ((open[i] + close[i]) / 2.0);
        case pr_typical:    return ((high[i] + low[i] + close[i]) / 3.0);
        case pr_weighted:   return ((high[i] + low[i] + close[i] + close[i]) / 4.0);
        case pr_average:    return ((high[i] + low[i] + close[i] + open[i]) / 4.0);
        case pr_trendbiased:
            if (close[i] > open[i]) {
                return ((high[i] + close[i]) / 2.0);
            } else {
                return ((low[i] + close[i]) / 2.0);
            }
        }
    }
    return (0);
}
