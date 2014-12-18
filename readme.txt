Module: hierarchical motion segments generation
Author: Subarna tripathi

Please cite our paper if you use the code.

Subarna Tripathi, Youngbae Hwang, Serge Belongie, and Truong  Nguyen; "Improving Streaming Video Segmentation with Early and Mid-Level Visual Processing", WACV 2014

Disclaimer: This alpha-version code runs only in windows 32-bit MATLAB due to dependencies on some third party binaries (provided in the package). 

Details:
This uses stability based segmentation on optical flow fields between consecutive frames in a video. Then it hierarchically merges regions based on geometric distances (directed divergence to be specific) in affine motion space. MRF-based smoothing is performed on the final output. Backward warping based segmentation of the previous frame serves as the initialization of the of the forward warping based segmentation of the next video frame, thus enforcing temporal stability. 

this code is little different from what was used for the results generation of the conference paper. I will clean the code in future.