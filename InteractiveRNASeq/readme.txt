FH-Differential application.

This folder contains all of the R scripts and data needed to run the FH-Differential application created by Jack Bruton during Jan-March 2023. In order to run the application successfully
from a local machine. open the app.R script located in this folder using R studio and Run the script as a Shiny app. You will need to have installed the Shiny package in order to run this.
Alternativley, you can also access the application through the web page: https://jackbruton.shinyapps.io/fh_differential/

Please note that as of now the deployed app does not work well as it exists on a free plan that limits the memory usage and you may find that frequent disconnects occur. As of now it is 
better to run the application locally to take advantage of your own computers RAM using the steps outlined above.

For Devs and anyone else who wants to work on the application:
The Shiny modules that build each tab are located in the R directory,
normalised counts and log2foldchange data for heatmap generation are stored in the data directory.

