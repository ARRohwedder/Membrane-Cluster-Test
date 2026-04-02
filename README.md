Membrane Cluster Test plugin technical description.
A temporary copy of the stack is generated to perform necessary changes for the analytics.
For outline detection all images of the stack are summed while values are multiplied by the factor 1.2 and stored in a new image. Auto threshold using the “Huang” method is applied to the image, followed be rank filters Median and Bright_outliers to remove noise. A 2 pixel frame is drawn around the image outline to prevent detection of the entire image. The resulting binary is analysed using the particle analyzer and the largest identified “particle” translated into a poly-ROI as outline.
Coordinates of the ROI outline and centre of mass are stored in integer arrays/variables.
Calculation of vector from outline points to centre:

For each outline position the angle Θ is calculated and stored into an array to provide base for position calculation of points below outline:

Grey values of 4 points along the vector h towards the centre of mass are measured in the 2 compared images and then replaced by background colour in the work copy of the original image:

The averages of the grey values per slice are stored in arrays. Based on the average normalized intensity data a co-localization array is generated if Rcol is within confidence interval:

Other positions are filled with 0 Positive positions are tested for their neighbours and regarded as cluster if when the neighbour is positive as well. 

Start of the clusters are drawn in output image (cluster). 
