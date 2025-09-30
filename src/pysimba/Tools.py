import pickle
import numpy as np

from . import BASE_DIR


class Tools:

    # Stores Dictionarie in Pickle with given Tag for filename
    def StoreInPickle(dictionary: dict, tag: str):
        with open(BASE_DIR / tag + ".pkl", "wb") as file:
            pickle.dump(dictionary,file)
        return
            
    # Converts Pickle into Dictionaries
    def PickleToDict(pkl_file: str) -> dict:
        pkl_name = pkl_file + ".pkl"
        path = BASE_DIR / pkl_name
        """
        Parameter:
        pkl_file: name of the .pkl file. Doesnt't include the string '.pkl'
        """
        with open(path, "rb") as file:
            data = pickle.load(file)
        return data
    
    # Returns SubMatrix at given indices of a larger Matrix
    def GetSub(matrix, row_lwb: int, row_upb: int, col_lwb: int, col_upb: int) -> list[list[float]]:
        return matrix[row_lwb:row_upb+1, col_lwb:col_upb+1]
    
    # Combines two matrices, by making a larger matrix
    def SetSub(matrix, row_lwb: int, col_lwb: int, sub_matrix):
        rows, cols = sub_matrix.shape
        matrix[row_lwb:row_lwb+rows, col_lwb:col_lwb+cols] = sub_matrix
        return

    # Forces a positive definite Matrice
    def MakePositiveDefinite(A: np.ndarray, eps: float = 1e-8) -> np.ndarray:
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
    
    # Turns every string of a possible Basis/ or Subleading Basis Expansion into the corresponding lambda
    # i.e '0575' -> 0.575 or '10575' -> 0.575
    def StrToLambda(BasisExpansion: str) -> float:
        """
        Examples for Input: '0575' or '10575'
        """

        # Step 1: Remove leading 1 or 2 if present
        if BasisExpansion[0] in ('1', '2'):
            BasisExpansion = BasisExpansion[1:]
        
        # Step 2: If only 1 digit remains -> interpret as X.X
        if len(BasisExpansion) == 1:
            return float(f"{BasisExpansion}.{BasisExpansion}")
        
        # Step 3: Otherwise, take the last up to 3 digits as fraction
        fraction_digits = BasisExpansion[-3:]
        
        # Step 4: Remove leading zeros for correct formatting
        fraction_str = fraction_digits.lstrip('0')
        
        # Step 5: Handle case where all digits were zero
        if not fraction_str:
            fraction_str = '0'
        
        # Step 6: Convert to 0.xxx float
        return float(f"0.{fraction_str}")
    
    def range_2d(x, y, x0 = 0, y0 = 0, xstep = 1, ystep = 1):
        for i in range(x0, x, xstep):
            for j in range(y0, y, ystep):
                yield i, j


    def range_3d(x, y, z, x0 = 0, y0 = 0, z0 = 0, xstep = 1, ystep = 1, zstep = 1):
        for i in range(x0, x, xstep):
            for j in range(y0, y, ystep):
                for k in range(z0, z, zstep):
                    yield i, j, k

    def make_c_string(n: int) -> str:
        # Create string "012...(n-1)"
        subscript = "".join(str(i) for i in range(n))
        # Format as c_{...}
        return f"c_{{{subscript}}}"