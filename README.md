This example was given by Montoya (2017) for the study of a group work in computer science.  
Female participants (N = 51) read two syllabi for different computer science classes.
One of the syllabi reported the class would have group projects throughout, and the other syllabi stated that individual projects would be scheduled throughout. 
Measured within-subject variables are “interest in the two classes” (interest) (Y_1,Y_2), “perceptions that the class has a communal environment” (perception) (M_1,M_2).
Each variable is the average of two or three variables measured by 7-point Likert scores
Measured between-subject moderator is “preference for group work” (preference) (W).
This example corresponds to a two-condition within-participant mediation model with a between-participant moderator (W: preference). 
This file is a R code for analyzing the two-condition within-participant mediation design using function "lavaan( )".
1. define SEM models for the natural and rotated condition approaches
2. estimate the path coefficients and indirect effects
3. test the coefficients and indirect effects using the bootstrap confidence intervals
4. probe the interaction effects by plotting the conditional indirect effects against the moderator values.
