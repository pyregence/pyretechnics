# [[file:../../org/pyretechnics.org::conversion-pxd][conversion-pxd]]
from pyretechnics.cy_types cimport vec_xy

cdef float rad_to_deg(float radians) noexcept
cdef float deg_to_rad(float degrees) noexcept
cdef float deg_to_ratio(float degrees) noexcept
cdef float ratio_to_deg(float ratio) noexcept
cdef float F_to_K(float degrees) noexcept
cdef float K_to_F(float degrees) noexcept
cdef float F_to_C(float degrees) noexcept
cdef float C_to_F(float degrees) noexcept
cdef float ch_to_m(float ch) noexcept
cdef float m_to_ch(float m) noexcept
cdef float m_to_ft(float m) noexcept
cdef float ft_to_m(float ft) noexcept
cdef float mph_to_mps(float mph) noexcept
cdef float mps_to_mph(float mps) noexcept
cdef float km_hr_to_mps(float km_hr) noexcept
cdef float mps_to_km_hr(float mps) noexcept
cdef float mph_to_km_hr(float mph) noexcept
cdef float km_hr_to_mph(float km_hr) noexcept
cdef float m_min_to_km_hr(float m_min) noexcept
cdef float km_hr_to_m_min(float km_hr) noexcept
cdef float m_min_to_mph(float m_min) noexcept
cdef float mph_to_m_min(float mph) noexcept
cdef float mps_to_fpm(float mps) noexcept
cdef float fpm_to_mps(float fpm) noexcept
cdef float mph_to_fpm(float mph) noexcept
cdef float fpm_to_mph(float fpm) noexcept
cdef float Btu_ft_s_to_kW_m(float Btu_ft_s) noexcept
cdef float kW_m_to_Btu_ft_s(float kW_m) noexcept
cdef float Btu_lb_to_kJ_kg(float Btu_lb) noexcept
cdef float kJ_kg_to_Btu_lb(float kJ_kg) noexcept
cdef float kg_m3_to_lb_ft3(float kg_m3) noexcept
cdef float lb_ft3_to_kg_m3(float lb_ft3) noexcept
cdef float percent_to_dec(float percent) noexcept
cdef float dec_to_percent(float decimal) noexcept
cdef float sec_to_min(float seconds) noexcept
cdef float min_to_sec(float minutes) noexcept
cdef float ms_to_min(float milliseconds) noexcept
cdef float min_to_ms(float minutes) noexcept
cdef float hour_to_min(float hours) noexcept
cdef float min_to_hour(float minutes) noexcept
cdef float day_to_min(float days) noexcept
cdef float min_to_day(float minutes) noexcept
cdef vec_xy cartesian_to_polar(float x, float y) noexcept
cdef vec_xy polar_to_cartesian(float r, float theta) noexcept
cdef vec_xy cartesian_to_azimuthal(float x, float y) noexcept
cdef vec_xy azimuthal_to_cartesian(float r, float azimuth) noexcept
cdef float opposite_direction(float theta) noexcept
cdef float wind_speed_10m_to_wind_speed_20ft(float wind_speed_10m) noexcept
cdef float wind_speed_20ft_to_wind_speed_10m(float wind_speed_20ft) noexcept
# conversion-pxd ends here
