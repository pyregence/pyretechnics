# [[file:../../org/pyretechnics.org::fuel-models-pxd][fuel-models-pxd]]
from pyretechnics.cy_types cimport fclaarr, CompactFuelModel, FuelModel

cdef dict[int, CompactFuelModel] compact_fuel_model_table

cdef bint is_burnable_fuel_model_number(int fuel_model_number) noexcept

cdef float compute_exp_A_sigma(float A, float sigma_ij) noexcept

cdef float compute_firemod_size_class(float sigma_i) noexcept

cdef FuelModel expand_compact_fuel_model(int fuel_model_number) noexcept

cdef dict[int, FuelModel] fuel_model_table

cpdef list[int] list_fuel_model_numbers()

cpdef list[FuelModel] list_fuel_models()

cpdef bint fuel_model_exists(int fuel_model_number) noexcept

cpdef FuelModel get_fuel_model(int fuel_model_number) noexcept

cdef FuelModel add_dynamic_fuel_loading(FuelModel fuel_model, fclaarr M_f) noexcept

cdef float compute_gij(fclaarr firemod_size_classes, fclaarr f_ij, float firemod_size_class_ij, bint is_dead) noexcept

cdef FuelModel add_weighting_factors(FuelModel fuel_model) noexcept

cdef FuelModel add_live_moisture_of_extinction(FuelModel fuel_model) noexcept

cpdef FuelModel moisturize(FuelModel fuel_model, fclaarr fuel_moisture) noexcept
# fuel-models-pxd ends here
