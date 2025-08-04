import pickle
import numpy as np

class Tools:
    def store_in_pickle(dictionary, tag):
        with open(tag+".pkl", "wb") as file:
            pickle.dump(dictionary,file)
        return
            
    def pickle_to_dict(pkl_file):
        with open(pkl_file + ".pkl", "rb") as file:
            data = pickle.load(file)
        return data
    
    def get_sub(matrix, row_lwb, row_upb, col_lwb, col_upb):
        return matrix[row_lwb:row_upb+1, col_lwb:col_upb+1]
    
    def set_sub(matrix, row_lwb, col_lwb, sub_matrix):
        rows, cols = sub_matrix.shape
        matrix[row_lwb:row_lwb+rows, col_lwb:col_lwb+cols] = sub_matrix

    def make_positive_definite(A: np.ndarray, eps: float = 1e-8) -> np.ndarray:
        """
        Force a symmetric matrix A to be numerically positive definite.

        Parameters
        ----------
        A : np.ndarray
            Input square matrix.
        eps : float
            Small positive number for eigenvalue floor.

        Returns
        -------
        np.ndarray
            Positive definite version of A.
        """
        # Ensure the matrix is symmetric
        A = (A + A.T) / 2

        # Eigendecomposition
        eigenvalues, eigenvectors = np.linalg.eigh(A)

        # Clip the eigenvalues to eps
        eigenvalues_clipped = np.clip(eigenvalues, a_min=eps, a_max=None)

        # Reconstruct the matrix
        A_pd = eigenvectors @ np.diag(eigenvalues_clipped) @ eigenvectors.T

        return A_pd  
    
    def get_near_psd(A):
        C = (A + A.T)/2
        eigval, eigvec = np.linalg.eig(C)
        eigval[eigval < 0] = 0
    
        return eigvec.dot(np.diag(eigval)).dot(eigvec.T)
    