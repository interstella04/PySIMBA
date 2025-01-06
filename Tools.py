import pickle


class Tools:
# These two following functions are tools
    def store_in_pickle(dictionary, tag):
        with open(tag+".pkl", "wb") as file:
            pickle.dump(dictionary,file)
        return
            
    def pickle_to_dict(pkl_file):
        with open(pkl_file + ".pkl", "rb") as file:
            data = pickle.load(file)
        return data