# BEPA
Brute-force Exploratory Path Analysis in R

This R package helps perform exploratory path analysis in R by generating an index of all possible model structures.

As the name suggests, this package is not efficient or elegant. As the number of variables increases, the number of possible models increases exponentially. This code has only be tested on up to 6 variables, and it already takes a substantial amount of time to generate the model combinations for this many variables. 

Certain kinds of relationships can also be excluded. For example, you can specify all possible models structures except for those in which x predicts y.


