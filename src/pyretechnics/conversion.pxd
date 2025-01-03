# [[file:../../org/pyretechnics.org::conversion-pxd][conversion-pxd]]
cimport pyretechnics.cy_types
from pyretechnics.cy_types cimport vec_xy

#==================================================
# Cython functions to cimport into other modules
#==================================================

cpdef float rad_to_deg(float radians)
cpdef float deg_to_rad(float degrees)
cpdef float deg_to_ratio(float degrees)
cpdef float ratio_to_deg(float ratio)
cpdef float F_to_K(float degrees)
cpdef float K_to_F(float degrees)
cpdef float F_to_C(float degrees)
cpdef float C_to_F(float degrees)
cpdef float ch_to_m(float ch)
cpdef float m_to_ch(float m)
cpdef float m_to_ft(float m)
cpdef float ft_to_m(float ft)
cpdef float mph_to_mps(float mph)
cpdef float mps_to_mph(float mps)
cpdef float km_hr_to_mps(float km_hr)
cpdef float mps_to_km_hr(float mps)
cpdef float mph_to_km_hr(float mph)
cpdef float km_hr_to_mph(float km_hr)
cpdef float m_min_to_km_hr(float m_min)
cpdef float km_hr_to_m_min(float km_hr)
cpdef float m_min_to_mph(float m_min)
cpdef float mph_to_m_min(float mph)
cpdef float mps_to_fpm(float mps)
cpdef float fpm_to_mps(float fpm)
cpdef float mph_to_fpm(float mph)
cpdef float fpm_to_mph(float fpm)
cpdef float Btu_ft_s_to_kW_m(float Btu_ft_s)
cpdef float kW_m_to_Btu_ft_s(float kW_m)
cpdef float Btu_lb_to_kJ_kg(float Btu_lb)
cpdef float kJ_kg_to_Btu_lb(float kJ_kg)
cpdef float kg_m3_to_lb_ft3(float kg_m3)
cpdef float lb_ft3_to_kg_m3(float lb_ft3)
cpdef float percent_to_dec(float percent)
cpdef float dec_to_percent(float decimal)
cpdef float sec_to_min(float seconds)
cpdef float min_to_sec(float minutes)
cpdef float ms_to_min(float milliseconds)
cpdef float min_to_ms(float minutes)
cpdef float hour_to_min(float hours)
cpdef float min_to_hour(float minutes)
cpdef float day_to_min(float days)
cpdef float min_to_day(float minutes)
cpdef vec_xy cartesian_to_polar(float x, float y)
cpdef vec_xy polar_to_cartesian(float r, float theta)
cpdef vec_xy cartesian_to_azimuthal(float x, float y)
cpdef vec_xy azimuthal_to_cartesian(float r, float azimuth)
cpdef float opposite_direction(float theta)
cpdef float wind_speed_10m_to_wind_speed_20ft(float wind_speed_10m)
cpdef float wind_speed_20ft_to_wind_speed_10m(float wind_speed_20ft)
# conversion-pxd ends here
