import pickle
import numpy as np
from dataclasses import dataclass

import pathlib

# To make the path of the Base directory of the file easily accesible in the whole Code
BASE_DIR = pathlib.Path(__file__).parent.resolve()


# For config File
import yaml


@dataclass
class settings:
    # Open the yaml (config) file and store the information in config (dictionary)
    with open(BASE_DIR / "data/settings.yml", "r") as f:
        config = yaml.safe_load(f)

    # Read out Information of config
    SubLeadCoefficients = config["SubLeadCoefficients"]
    TheoryOrder = config["TheoryOrder"]
    SubLeadTheoryOrder = config["SubLeadTheoryOrder"]
    FitVars = config["FitVars"]
    KeyOrder = config["KeyOrder"]
    BasisExpansion = config["BasisExpansion"]
    SubLeadBasisExpansion = config["SubLeadBasisExpansion"]

    rho2: float = config["Constants"]["rho2"]
    mB: float = config["Constants"]["mB"]
    La2: float = config["Constants"]["La2"]
    N0: float = config["Constants"]["N0"]
    VtbVts: float = config["Constants"]["VtbVts"]
    C2C7: float = config["Constants"]["C2C7"]
    C2C2: float = config["Constants"]["C2C2"]
    C8C7: float = config["Constants"]["C8C7"]
    C8C8: float = config["Constants"]["C8C8"]
    C2C8: float = config["Constants"]["C2C8"]


class Tools:
    # Stores Dictionarie in Pickle with given Tag for filename
    def StoreInPickle(dictionary: dict, tag: str):
        path = tag + ".pkl"
        with open(BASE_DIR / path, "wb") as file:
            pickle.dump(dictionary, file)
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
    def GetSub(
        matrix, row_lwb: int, row_upb: int, col_lwb: int, col_upb: int
    ) -> list[list[float]]:
        return matrix[row_lwb : row_upb + 1, col_lwb : col_upb + 1]

    # Combines two matrices, by making a larger matrix
    def SetSub(matrix, row_lwb: int, col_lwb: int, sub_matrix):
        rows, cols = sub_matrix.shape
        matrix[row_lwb : row_lwb + rows, col_lwb : col_lwb + cols] = sub_matrix
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
        if BasisExpansion[0] in ("1", "2"):
            BasisExpansion = BasisExpansion[1:]

        # Step 2: If only 1 digit remains -> interpret as X.X
        if len(BasisExpansion) == 1:
            return float(f"{BasisExpansion}.{BasisExpansion}")

        # Step 3: Otherwise, take the last up to 3 digits as fraction
        fraction_digits = BasisExpansion[-3:]

        # Step 4: Remove leading zeros for correct formatting
        fraction_str = fraction_digits.lstrip("0")

        # Step 5: Handle case where all digits were zero
        if not fraction_str:
            fraction_str = "0"

        # Step 6: Convert to 0.xxx float
        return float(f"0.{fraction_str}")

    def range_2d(x, y, x0=0, y0=0, xstep=1, ystep=1):
        for i in range(x0, x, xstep):
            for j in range(y0, y, ystep):
                yield i, j

    def range_3d(x, y, z, x0=0, y0=0, z0=0, xstep=1, ystep=1, zstep=1):
        for i in range(x0, x, xstep):
            for j in range(y0, y, ystep):
                for k in range(z0, z, zstep):
                    yield i, j, k

    def make_c_string(n: int) -> str:
        # Create string "012...(n-1)"
        subscript = "".join(str(i) for i in range(n))
        # Format as c_{...}
        return f"c_{{{subscript}}}"

    # def merge_paths(base, rest_path: str):
    #    return base / rest_path

    # Useful functions to add a measurement to the code
    # ATTENTION: These functions will change the data in the dictionary

    def combine_dict_old_with_dict_new(old_path: str, new_path: str, new_key: str):
        old = Tools.PickleToDict(old_path)
        new = Tools.PickleToDict(new_path)

        Tools.StoreInPickle(old, old_path + "_save")
        print(
            "You've stored the old dictionary in a '_save' pickle. You will find it at the location of your old dictionary."
        )

        old[new_key] = new

        Tools.StoreInPickle(old, old_path)
        print(
            "______________________________________________________________________________________"
        )
        print(
            "You just overwrote your old dictionary. The old one and the new one should now be combined. \n \n If you have done this for every nessecary dictionary (Measurements, Leading and Subleading Theory) you should now change the settings.yml file, so that it includes the new measurement. You do that, when you add your 'key' to the 'KeyOrder' \n"
        )
        print(
            "______________________________________________________________________________________"
        )

        test = Tools.PickleToDict(old_path)
        print(
            "The new keys of the dictionary are as follows, please check if your new meaurement is part of it."
        )
        print(
            "______________________________________________________________________________________"
        )
        print(test.keys())
        print(
            "______________________________________________________________________________________"
        )

        return

        # def AddNewMeasurement():
        print("To which dictionary do you want to add a dictionary? Please type")
        print(" '1' - for Measurement")
        print(" '2' - for Leading Theory")
        print(" '3' - for SubLeading Theory")

        choice = input("Enter your choice: ").strip()

        name = input(
            "Enter the name of the pickle file (without .pkl) which you want to add to the chosen dictionary:"
        ).strip()

        new_key = input(
            "Now enter the key of your new measurement (e.g 'belle2'):"
        ).strip()

        if choice == "1":
            Tools.combine_dict_old_with_dict_new(
                settings.config["MeasurementPath"], "data/add/" + name, new_key
            )
        elif choice == "2":
            Tools.combine_dict_old_with_dict_new(
                settings.config["TheoryPath"], "data/add/" + name, new_key
            )
        elif choice == "3":
            Tools.combine_dict_old_with_dict_new(
                settings.config["SubleadingTheoryPath"], "data/add/" + name, new_key
            )
        else:
            print("Not a valid choice!")

        return

    def AddNewMeasurement():
        ##########################################################################
        # This is a test function, remove when real usage starts
        # This one is here to make sure that the dictionaries stay as they are for now
        ##########################################################################
        print("To which dictionary do you want to add a dictionary? Please type")
        print(
            "______________________________________________________________________________________"
        )
        print(" '1' - for Measurement")
        print(" '2' - for Leading Theory")
        print(" '3' - for SubLeading Theory")
        print(
            "______________________________________________________________________________________"
        )

        choice = input("Enter your choice: ").strip()
        print(
            "______________________________________________________________________________________"
        )

        name = input(
            "Enter the name of the pickle file (without .pkl) which you want to add to the chosen dictionary:"
        ).strip()
        print(
            "______________________________________________________________________________________"
        )

        new_key = input(
            "Now enter the key of your new measurement (e.g 'belle2'):"
        ).strip()
        print(
            "______________________________________________________________________________________"
        )

        if choice == "1":
            Tools.combine_dict_old_with_dict_new(
                "data/test_old", "data/" + name, new_key
            )
        else:
            print("Not a valid choice, because this is a test!!!")

        return
