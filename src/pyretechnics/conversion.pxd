# [[file:../../org/pyretechnics.org::conversion-pxd][conversion-pxd]]
cimport pyretechnics.cy_types
from pyretechnics.cy_types cimport vec_xy

#==================================================
# Cython functions to cimport into other modules
#==================================================

cdef float rad_to_deg(float radians) noexcept
cdef float deg_to_rad(float degrees) noexcept
cdef float deg_to_ratio(float degrees)
cdef float ratio_to_deg(float ratio)
cdef float F_to_K(float degrees)
cdef float K_to_F(float degrees)
cdef float F_to_C(float degrees)
cdef float C_to_F(float degrees)
cdef float ch_to_m(float ch)
cdef float m_to_ch(float m)
cdef float m_to_ft(float m)
cdef float ft_to_m(float ft)
cdef float mph_to_mps(float mph)
cdef float mps_to_mph(float mps)
cdef float km_hr_to_mps(float km_hr)
cdef float mps_to_km_hr(float mps)
cdef float mph_to_km_hr(float mph)
cdef float km_hr_to_mph(float km_hr) noexcept
cdef float m_min_to_km_hr(float m_min)
cdef float km_hr_to_m_min(float km_hr)
cdef float m_min_to_mph(float m_min)
cdef float mph_to_m_min(float mph)
cdef float mps_to_fpm(float mps)
cdef float fpm_to_mps(float fpm)
cdef float mph_to_fpm(float mph)
cdef float fpm_to_mph(float fpm)
cdef float Btu_ft_s_to_kW_m(float Btu_ft_s)
cdef float kW_m_to_Btu_ft_s(float kW_m)
cdef float Btu_lb_to_kJ_kg(float Btu_lb)
cdef float kJ_kg_to_Btu_lb(float kJ_kg)
cdef float kg_m3_to_lb_ft3(float kg_m3)
cdef float lb_ft3_to_kg_m3(float lb_ft3)
cdef float percent_to_dec(float percent)
cdef float dec_to_percent(float decimal)
cdef float sec_to_min(float seconds)
cdef float min_to_sec(float minutes)
cdef float ms_to_min(float milliseconds)
cdef float min_to_ms(float minutes)
cdef float hour_to_min(float hours)
cdef float min_to_hour(float minutes)
cdef float day_to_min(float days)
cdef float min_to_day(float minutes)
cdef vec_xy cartesian_to_polar(float x, float y)
cdef vec_xy polar_to_cartesian(float r, float theta)
cdef vec_xy cartesian_to_azimuthal(float x, float y) noexcept
cdef vec_xy azimuthal_to_cartesian(float r, float azimuth) noexcept
cdef float opposite_direction(float theta) noexcept
cdef float wind_speed_10m_to_wind_speed_20ft(float wind_speed_10m)
cdef float wind_speed_20ft_to_wind_speed_10m(float wind_speed_20ft)
# conversion-pxd ends here
