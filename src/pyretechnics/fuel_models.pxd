from pyretechnics.cy_types cimport fclaarr, FuelModel

cdef dict fuel_model_table

cdef FuelModel fm_struct(dict fm)

cdef FuelModel moisturize(dict fuel_model, list fuel_moisture) noexcept

cdef FuelModel moisturize_val(FuelModel fuel_model, fclaarr fuel_moisture) noexcept
