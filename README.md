# CompressedSensing
Recovery of Sparse Matrix in O(N) via Compressed Sensing Techniques

# Background
http://arxiv.org/pdf/1410.4848v1.pdf


![alt tag](https://raw.githubusercontent.com/twhughes/CompressedSensing/master/example.png)

# Brief Explanation:

You have a sparse matrix (top left picture) that you are trying to measure

Instead of sampling the matrix directly (O(N^2) running time), instead sample the same matrix in a basis that is more dense, then transform back into sparse basis (top right picture).  This takes O(N) to sample to arbitrary precision

Code plots number of iterations needed to get to arbitrary accuracy vs. sparsity of matrix (bottom left pictire)

Average of N runs (bottom right)

Red dotted line is amount of iterations needed if you were sampling the sparse matrix directly (sampling all elements)

Green line is the amount of sampling you would need to do if you knew, a priori, which elements were non-zero.



