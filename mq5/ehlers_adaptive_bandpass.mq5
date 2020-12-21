//+------------------------------------------------------------------+
//|                                     ehlers_adaptive_bandpass.mq5 |
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
#property indicator_buffers 5
#property indicator_plots 2
#property indicator_maximum 1.2
#property indicator_minimum -1.2
#property indicator_level1 -0.707
#property indicator_level2 0.707

//Plot Line
#property indicator_label1 "Adaptive BandPass"
#property indicator_type1 DRAW_LINE
#property indicator_color1 clrGoldenrod
#property indicator_style1 STYLE_SOLID
#property indicator_width1 2

//Plot Signal Line
#property indicator_label2  "Bandpass Signal"
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
input double bandwidth = 0.3; 
input int avgLength = 3;        // Averaging Length
input int hpPeriod = 48;        // High Pass Period
input int lpPeriod = 10;        // Low Pass Period

// buffers
double result[]; 	// result buffer
double signal[];    // signal line buffer
double filt[]; 		// filter buffer
double highPass[]; 	// highpass buffer
double bandPass[]; 	// bandPass buffer

//////////////////////////////////////////////////////////////////

int OnInit() {
    SetIndexBuffer(0, result, INDICATOR_DATA);
    SetIndexBuffer(1, signal, INDICATOR_DATA);
    SetIndexBuffer(2, filt, INDICATOR_CALCULATIONS);
    SetIndexBuffer(3, highPass, INDICATOR_CALCULATIONS);
    SetIndexBuffer(4, bandPass, INDICATOR_CALCULATIONS);
    return (0);
}

//////////////////////////////////////////////////////////////////

double work[];
double corr[];
double sqSum[];
double r[][2];
double maxPwr = 0.0;
double pwr[];
double peak = 0.0;

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

    // if the default array size if not chosen then adjust the array sizes
    if(ArrayRange(corr, 0) != (hpPeriod + 1)){
    	ArrayResize(corr, (hpPeriod + 1));
    	ArrayResize(sqSum, (hpPeriod + 1));
    	ArrayResize(r, (hpPeriod + 1));
    	ArrayResize(pwr, (hpPeriod + 1));
    }

    /* Most of the time prev_calculated == rates_total == _bars in which case this loop only calculates
    one bar of data. However, additional history or a refresh of data will cause prev_calculated
    to reset to 0, and hence cause a loop through all previous bars */

    for (int i = (int) MathMax(prev_calculated - 1, 0); i < rates_total; i++) {

        // get the price of the current bar based on the Price variable input
        work[i] = getPrice(Price, open, close, high, low, i, rates_total);

        // highpass filter
        highPass[i] = processHighPass(hpPeriod, work, highPass, i);

        // lowpass filter (supersmoother)
        filt[i] = processLowPass(lpPeriod, highPass, filt, i);

        // pearson correlation on lagged values
        processPearsonCorr(hpPeriod, avgLength, filt, i, corr);

        // discrete fourier transform
        processDFT(lpPeriod, hpPeriod, corr, i, sqSum);

        // EMA to smooth power measurement at each period
        processEMA(lpPeriod, hpPeriod, sqSum, r);

        // fast attack slow decay Auto Gain Control (AGC)
        processAGC(lpPeriod, hpPeriod, maxPwr, r, pwr, i);

        // dominant cycle using the CG of the spectrum
        double dominantCycle = computeDominantCycle(lpPeriod, hpPeriod, pwr);

        // stochastic (Result)
        computeBandPass(lpPeriod, hpPeriod, filt, (int)dominantCycle, bandPass, result, i, bandwidth, peak, signal);
    }

    // return value of prev_calculated for next call
    return (rates_total);
}

void computeBandPass(int lowPassPeriod, int highPassPeriod, double& inputArr[], int domCycle, double& bandpass[], double& finalResult[], int i, double bandWidth, double& peakVal, double& signalLine[]){

    if (i > (highPassPeriod + 1)){

        double pi = 3.14159265358979323846264338327950288;
        double beta1 = MathCos((2*pi) / (0.9 * domCycle));
        double gamma1 = 1 / MathCos((2*pi) * bandWidth / (0.9*domCycle));
        double alpha2 = gamma1 - MathSqrt(pow(gamma1,2) - 1);

        bandpass[i] = 0.5*(1 - alpha2)*(inputArr[i] - inputArr[i-2]) + beta1*(1 + alpha2)*bandpass[i-1] - alpha2*bandpass[i-2];

        peakVal = 0.991*peakVal;
        if(MathAbs(bandpass[i]) > peakVal){
            peakVal = MathAbs(bandpass[i]);
        }
        if(peakVal != 0){
            finalResult[i] = bandpass[i] / peakVal;
        } else {
            finalResult[i] = 0;
        }

        signalLine[i] = 0.9 * finalResult[i-1];


    } else {

        finalResult[i] = 0;

    }
}

// pass input data through a high pass filter

double processHighPass(double highPassPeriod, double& inputArr[], double& outputArr[], int i) {

    double pi = 3.14159265358979323846264338327950288;
    double alpha1 = (MathCos(MathSqrt(2.0) * pi / highPassPeriod) + MathSin(MathSqrt(2.0) * pi / (double)highPassPeriod) - 1.0) / MathCos(MathSqrt(2.0) * pi / (double)highPassPeriod);
    if (i < 3) {
        return ((1.0 - alpha1 / 2.0) * (1.0 - alpha1 / 2.0) * (inputArr[i]));
    } else {
        return ((1.0 - alpha1 / 2.0) * (1.0 - alpha1 / 2.0) * (inputArr[i] - 2.0 * inputArr[i - 1] + inputArr[i - 2]) + 2.0 * (1.0 - alpha1) * outputArr[i - 1] - (1.0 - alpha1) * (1.0 - alpha1) * outputArr[i - 2]);
    }
}

// pass input data through a low pass filter (supersmoother)

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

void processPearsonCorr(int barsLag, int averageLength, double& inputArr[], int i, double& outputArr[]) {

    int avLength = averageLength;
    if(averageLength == 0){
        avLength = barsLag;
    }

    if (i > (barsLag + 1 + avLength)) {

        for (int lag = 0; lag <= barsLag; lag++) {

            double m = (double)averageLength;
            if (averageLength == 0) {
                m = (double)lag;
            }

            double Sx = 0;
            double Sy = 0;
            double Sxx = 0;
            double Syy = 0;
            double Sxy = 0;

            for (int count = 0; count < m; count++) {
                double X = inputArr[i - count];
                double Y = inputArr[i - lag - count];
                Sx = Sx + X;
                Sy = Sy + Y;
                Sxx = Sxx + X * X;
                Sxy = Sxy + X * Y;
                Syy = Syy + Y * Y;
            }

            if (((m * Sxx - Sx * Sx) * (m * Syy - Sy * Sy)) > 0) {
                outputArr[lag] = (m * Sxy - Sx * Sy) / MathSqrt((m * Sxx - Sx * Sx) * (m * Syy - Sy * Sy));
            }

        }
    }
}

// Discrete Fourier Transform (DFT)

void processDFT(int lowPassPeriod, int highPassPeriod, double& inputArr[], int i, double& outputArr[]) {

    double pi = 3.14159265358979323846264338327950288;
    double cosPart[];
    double sinPart[];
    ArrayResize(cosPart, highPassPeriod + 1);
    ArrayResize(sinPart, highPassPeriod + 1);

    for (int period = lowPassPeriod; period <= highPassPeriod; period++) {

        cosPart[period] = 0.0;
        sinPart[period] = 0.0;

        // find cosine and sine correlated components
        for (int n = 3; n <= highPassPeriod; n++) {
            cosPart[period] = cosPart[period] + inputArr[n] * MathCos(2 * pi * (double) n / (double) period);
            sinPart[period] = sinPart[period] + inputArr[n] * MathSin(2 * pi * (double) n / (double) period);
        }

        outputArr[period] = pow(cosPart[period], 2.0) + pow(sinPart[period], 2.0);
    }
}

void processEMA(int lowPassPeriod, int highPassPeriod, double& inputArr[], double& outputArr[][2]) {

    for (int period = lowPassPeriod; period <= highPassPeriod; period++) {
        outputArr[period, 1] = outputArr[period, 0];
        outputArr[period, 0] = 0.2 * pow(inputArr[period], 2.0) + 0.8 * outputArr[period, 1];
    }

}

void processAGC(int lowPassPeriod, int highPassPeriod, double& maxPower, double& inputArr[][2], double& outputArr[], int i) {

    maxPower = 0.991 * maxPower;
    for (int period = lowPassPeriod; period <= highPassPeriod; period++) {
        if (inputArr[period, 0] > maxPower) {
            maxPower = inputArr[period, 0];
        }
    }
    if (maxPower != 0) {
        for (int period1 = 3; period1 <= highPassPeriod; period1++) {
            outputArr[period1] = inputArr[period1, 0] / maxPower;
        }
    } else {
        for (int period1 = 3; period1 <= highPassPeriod; period1++) {
            outputArr[period1] = 0;
        }
    }
}

//Compute the dominant cycle using the CG of the spectrum

double computeDominantCycle(int lowPassPeriod, int highPassPeriod, double& inputArr[]){

    double Spx = 0;
    double Sp = 0;
    for (int period = lowPassPeriod; period <= highPassPeriod; period++) {
        if (inputArr[period] >= 0.5) {
            Spx = Spx + period * inputArr[period];
            Sp = Sp + inputArr[period];
        }
    }

    double domCycle = 0;
    if (Sp != 0) {
        domCycle = Spx / Sp;
    }
    if (domCycle < lowPassPeriod) {
        domCycle = (double)lowPassPeriod;
    }

    return domCycle;
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
