# DNF-API
This is a web API for the DNF algorithm. Namely, you can use this API to get BHK lab's DNF network in JSON format, as well as upload new drug data to integrate into the default DNF network.
# What server does this use?
This project uses OpenCPU's free server. Refer to their documentation on how OpenCPU works.
# How can this API be improved?
1. OpenCPU has issues dealing with any packages that import rJava, so in the future we may have to switch from OpenCPU to a different server, possibly R Shiny or Microsoft Operationalization.
2. The R code for integrating the different DNF data is largely hard-coded for demonstration purpose in the DNF-Explorer webapp. For scalability, it needs to process the different data types (sensitivity, structure, and perturbation) in their raw forms.
# How to develop on this project?
## Requirements
1. R Studio
2. R OpenCPU package
3. Roxygen2
3. Postman (not absolutely necessary but great for debugging)
## Create an API endpoint
All endpoints must be declared in DNF.R in the R folder using the proper documentation conventions (as described in roxygen). The project must be installed/built for a change to take effect.
## Debugging
1. In R Studio load the necessary libraries (DNFAPI and OpenCPU)
2. In R Studio start the server with `ocpu_start_server(port=x)`
3. Use Postman to send requests to server x to test the endpoint (some examples are included in `postman` folder)
## Deploying (assuming you use OpenCPU as server)
You have to do three things to deploy the project. Firstly in R Studio, call `roxygen2::roxygenise()`, which will automatically modify certain miscellaneous files that are necessary for building the project on OpenCPU. Secondly, build the project using R Studio build/install (this is mostly a check to see if the API can be deployed). Lastly, commit the changes and push it onto the GitHub master branch of this project. The push will trigger OpenCPU to deploy the API onto their cloud server. You can check whether the API was deployed if you go to the GitHub repository -> Settings -> Webhook and click on the Webhook URL. It will say that the deployment failed, but you will also receive an email reporting the status of the deployment. If the email says that the deployment build was successful, then you can forget about the fail message on GitHub.
