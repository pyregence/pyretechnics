from pyretechnics.cy_types cimport FuelModel

cdef dict fuel_model_table

cdef FuelModel moisturize(dict fuel_model, list fuel_moisture) noexcept
