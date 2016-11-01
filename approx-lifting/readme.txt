Approximately lift the MLN using parametric or non=parametric clustering methods. Please refer to the following papers.
Deepak Venugopal and Vibhav Gogate, "Evidence-based Clustering for Scalable Inference in MLNs", ECML 2014
Deepak Venugopal, Somdeb Sarkhel and Kyle Cherry, "Non-parametric Domain Approximation for Scalable Gibbs sampling in MLNs", UAI, 2016

To run the pre-processor, execute lift_approx.py
For parametric, use lift_approx.py -p mln evidence query compression-ratio
For non-parametric, use lift_approx.py -np mln evidence query upper-bound[real-value between 0 and 1] lower-bound[real-value between 0 and 1, and greater than upper bound]

The output is a new MLN file cmln.txt and cevid.txt that represents the compressed MLN and evidence.