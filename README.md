# Adaptive Learning for Vehicle Detection and Image Instance Segmentation


**Programmd by**: Juwairiah Zia (javeriazia15@gmail.com) and Enze cui (cuienze2015@gmail.com)<br />
**Data Set**: Stanford Mobile Visual Search Database<br />
**Project**: This project describes an adaptive learning approach to detect object when there is no training data by applying computation geometry and engineering physics. 
## Problem
No training data for object detection
**Project**: Stationary camera
## Solution
### Adaptive Learning
- Involves frame subtraction to separate background and foreground
- Applying morphological operation and laplace of Gaussian Filter to remove noise 
- Implementing Canny Edge detection for corner and boundary 
- Creating the feature vector of blob/object for object classification 
### Polygon Object Instance Segmentation
- Creating a feature vector to store boundary points and centroid of moving object
- Triangulating using boundary points and centroid to implement polygon segmentation on moving object
- Generating graph using boundary points and centroid to store previous position of moving object
