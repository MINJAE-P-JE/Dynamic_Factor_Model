# Canadian, U.S. and Chinese Nowcasting

This project is based on the the [New York Fed Staff Nowcasting Report](https://www.newyorkfed.org/research/policy/nowcast) released every Friday.
The macroeconomic indicators generated by this nowcasting model were used to build commodity investment strategies.

This code implements the nowcasting framework described in "[Macroeconomic Nowcasting and Forecasting with Big Data](https://www.newyorkfed.org/research/staff_reports/sr830.html)" by Brandyn Bok, Daniele Caratelli, Domenico Giannone, Argia M. Sbordone, and Andrea Tambalotti, Staff Reports 830, Federal Reserve Bank of New York (prepared for Volume 10 of the Annual Review of Economics).

## File and folder description
- data : example US data downloaded from FRED
- functions : functions for loading data, estimating model, and updating predictions
- example_DFM.m : example script to estimate a dynamic factor model (DFM) for a panel of monthly and quarterly series
- example_Nowcast.m : example script to produce a nowcast or forecast for a target variable, e.g., real GDP growth
- ResDFM.mat : example DFM estimation output
- Spec_US_example.xls : example model specification for the US

## Reference
Banbura, M., D. Giannone, M. Modugno, and L. Reichlin. 2013. “[Nowcasting and the Real-Time Data Flow.](https://ideas.repec.org/h/eee/ecofch/2-195.html)” In G. Elliott and A. Timmermann, eds., Handbook of Economic Forecasting, Vol. 2. Amsterdam: Elsevier-North Holland.  

Bok, B., D. Caratelli, D. Giannone, A. Sbordone, and A. Tambalotti. 2017. “[Macroeconomic Nowcasting and Forecasting with Big Data.](https://www.newyorkfed.org/research/staff_reports/sr830)” Federal Reserve Bank of New York Staff Reports, no. 830, November.  

Giannone, D., L. Reichlin, and D. Small. 2008. “[Nowcasting: The Real-Time Informational Content of Macroeconomic Data.](https://ideas.repec.org/a/eee/moneco/v55y2008i4p665-676.html)” Journal of Monetary Economics 55, no.4 (May): 665-76.