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