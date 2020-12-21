# What is in this repository?

This repository contains (MetaQuotes Language 5) MQL5 code (used in MetaTrader 5) for some of the indicators featured in "Cycle Analytics for Traders - Advanced Technical Trading Concepts" by John F. Ehlers.

# Why are the files in this repository necessary?

The aforementioned book contains full coded versions of the indicators discussed within the book. However, the book contains code written in EasyLanguage (which 
is typically used with TradeStation). MQL5 on the other hand is based on the object oriented language C++. They are very different.

This repository is just a basic translation of the EasyLanguage code featured in the book into MQL5 code required for MetaTrader 5.

# Can I use these indicators?

The translation to MQL5 is free to use, however, the original EasyLanguage code is copyrighted by John F. Ehlers. 

Therefore, as far as I am aware it would be required to own a copy of the original book before usage of the files in this repository:

Cycle Analytics for Traders - Advanced Technical Trading Concepts by John F. Ehlers\
Copyright © 2013 by John F. Ehlers. All rights reserved.

[Wiley Publishing](https://www.wiley.com/en-us/Cycle+Analytics+for+Traders%3A+Advanced+Technical+Trading+Concepts%2C+%2B+Downloadable+Software-p-9781118728604)\
[Amazon UK](https://www.amazon.co.uk/Cycle-Analytics-Traders-Technical-Downloadable/dp/1118728513/)\
[Amazon US](https://www.amazon.com/Cycle-Analytics-Traders-Downloadable-Software/dp/1118728513)

Furthermore, it is recommended to purhase the book, as without it you will have no insight into the applicability of each of the indicators.

# Indicators

## What is included?

There are two directories included in this repository. The first (mq5) includes the raw MQL5 code for each indicator. The second
(ex5) contains the compiled version of each indicator.

In terms of indicators the following are included:

- Supersmoother
- Roofing Filter
- Even Better Sinewave
- Decycler Oscillator
- Autocorrelation Reversals
- Adaptive Bandpass
- Adaptive Bandpass Cube
- Adaptive CCI
- Adaptive RSI
- Adaptive RSI Fischer
- Adaptive Stochastic
- Adaptive Stochastic Inverse Fischer

## How to use them?

If you just want to use the indicators directly, then all you need is the .ex5 files (i.e. the compiled versions). Which can be
loaded directly into MetaTrader 5 as with any other indicator file.

If you want to make some tweaks to the files, you need the .mq5 code files, and then you will need to re-compile them once your
adjustments have been made.

## Can I use these files with MetaTrader 4?

No.

MetaQuotes (the company who make the MetaTrader program) made some major changes to the coding language between MetaTrader 4 
and MetaTrader 5. They are completely different, and not compatible.

# Disclaimer

I will not be held responsible for any financial losses you may incur as a result of the usage of the files in this repository, 
regardless of whether there are errors in either the translated code, or the original code.

If you choose to download and/or use the code in this repository, you do so at your own risk. If unsure, do not download or use it.

# License

As stated above, the copyright to the original EasyLanguage code is held by John F. Ehlers.\
Copyright © 2013 by John F. Ehlers. All rights reserved.

The translation from EasyLanguage to MQL5 is licensed under the [MIT License](LICENSE.md)