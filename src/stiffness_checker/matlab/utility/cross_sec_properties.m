function [Jx, Iy, Iz] = cross_sec_properties(m_p)

Jx = 0.5 * pi * m_p.r^4;
Iy = pi * m_p.r^4 / 4;
Iz = Iy;

end