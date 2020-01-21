"""
All variables used in the actual computation should be stored here. Use of dictionary gives flexibility to contents, leaves the door open to changing the iterator between iterations.
"""


"""
Type used to find iteration functions and building initial
model. Extend in file for each algorithm type.
"""
abstract type RowActionMethod end

"""
Type used as supertype of specific model formulations used by 
different algorithms.
"""
abstract type ModelFormulation end

