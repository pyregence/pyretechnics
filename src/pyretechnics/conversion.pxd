# [[file:../../org/pyretechnics.org::conversion-pxd][conversion-pxd]]
from pyretechnics.cy_types cimport vec_xy

cpdef float rad_to_deg(float radians) noexcept
cpdef float deg_to_rad(float degrees) noexcept
cpdef float deg_to_ratio(float degrees) noexcept
cpdef float ratio_to_deg(float ratio) noexcept
cpdef float F_to_K(float degrees) noexcept
cpdef float K_to_F(float degrees) noexcept
cpdef float F_to_C(float degrees) noexcept
cpdef float C_to_F(float degrees) noexcept
cpdef float ch_to_m(float ch) noexcept
cpdef float m_to_ch(float m) noexcept
cpdef float m_to_ft(float m) noexcept
cpdef float ft_to_m(float ft) noexcept
cpdef float mph_to_mps(float mph) noexcept
cpdef float mps_to_mph(float mps) noexcept
cpdef float km_hr_to_mps(float km_hr) noexcept
cpdef float mps_to_km_hr(float mps) noexcept
cpdef float mph_to_km_hr(float mph) noexcept
cpdef float km_hr_to_mph(float km_hr) noexcept
cpdef float m_min_to_km_hr(float m_min) noexcept
cpdef float km_hr_to_m_min(float km_hr) noexcept
cpdef float m_min_to_mph(float m_min) noexcept
cpdef float mph_to_m_min(float mph) noexcept
cpdef float mps_to_fpm(float mps) noexcept
cpdef float fpm_to_mps(float fpm) noexcept
cpdef float mph_to_fpm(float mph) noexcept
cpdef float fpm_to_mph(float fpm) noexcept
cpdef float Btu_ft_s_to_kW_m(float Btu_ft_s) noexcept
cpdef float kW_m_to_Btu_ft_s(float kW_m) noexcept
cpdef float Btu_lb_to_kJ_kg(float Btu_lb) noexcept
cpdef float kJ_kg_to_Btu_lb(float kJ_kg) noexcept
cpdef float kg_m3_to_lb_ft3(float kg_m3) noexcept
cpdef float lb_ft3_to_kg_m3(float lb_ft3) noexcept
cpdef float percent_to_dec(float percent) noexcept
cpdef float dec_to_percent(float decimal) noexcept
cpdef float sec_to_min(float seconds) noexcept
cpdef float min_to_sec(float minutes) noexcept
cpdef float ms_to_min(float milliseconds) noexcept
cpdef float min_to_ms(float minutes) noexcept
cpdef float hour_to_min(float hours) noexcept
cpdef float min_to_hour(float minutes) noexcept
cpdef float day_to_min(float days) noexcept
cpdef float min_to_day(float minutes) noexcept
cpdef vec_xy cartesian_to_polar(float x, float y) noexcept
cpdef vec_xy polar_to_cartesian(float r, float theta) noexcept
cpdef vec_xy cartesian_to_azimuthal(float x, float y) noexcept
cpdef vec_xy azimuthal_to_cartesian(float r, float azimuth) noexcept
cpdef float opposite_direction(float theta) noexcept
cpdef float wind_speed_10m_to_wind_speed_20ft(float wind_speed_10m) noexcept
cpdef float wind_speed_20ft_to_wind_speed_10m(float wind_speed_20ft) noexcept
# conversion-pxd ends here
