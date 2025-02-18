#import uproot
import numpy as np
#from Meas import Meas
import matplotlib.pyplot as plt
from Tools import Tools
#from iminuit import Minuit


class test:
    def __init__(self):
        self.a: int

    def set_a(self, value: int):
        self.a = value

    def get_a(self):
        return self.a
    

t = test()

print(type(t))
#print(t.a)

#print(t.get_a())

t.set_a(42)

print(t.a)

print(t.get_a())