# Adaptive Computer Vision Algorithm for Real Time Moving Vehicle Detection and Instance Segmentation
**Project**: This project describes an adaptive learning approach to detect objects when there is no training data by applying computation geometry and engineering physics. 
## Problem
No training data for object detection</br>
**Assumption**: Stationary camera
## Solution
### Adaptive Learning
- Involves frame subtraction to separate background and foreground
- Applying morphological operation and Laplace of Gaussian Filter to remove noise 
- Implementing Canny Edge detection for corner and boundary 
- Creating the feature vector of blob/object for object classification 
### Polygon Object Instance Segmentation
- Creating a feature vector to store boundary points and centroid of moving object
- Triangulating using boundary points and centroid to implement polygon segmentation on moving object
- Generating graph using boundary points and centroid to store the previous positions of moving object
### Shape Analysis for Object Classification
- Adding more feature using shape analysis
### Object Tracking
Coming Soon
