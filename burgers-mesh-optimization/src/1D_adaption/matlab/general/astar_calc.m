
gmma = 1.4;
gm1 = gmma -1;
R = 286.9;
p0 = 300000.0;
t0 = 300.0;
pback = 297330.0;
m_e = sqrt( 2/(gm1) * ( (p0/pback)^(gm1/gmma) - 1))

t_e = t0/(1 + (gm1/2)*m_e^2)

mdot = (pback/(R*t_e))* 1 * (m_e*sqrt(gmma*R*t_e))

m_t_old = 0.5;

for  k=1:500,
    
    m_t_new = mdot *sqrt(R*t0)/(p0*0.2*sqrt(gmma)) * (1 + (gm1/2)*m_t_old^2)^((gmma +1)/(2*gm1));
    
    err = abs(m_t_new-m_t_old)
    
    m_t_old = m_t_new;
end

m_t_new

astar = 0.2*m_t_new*( (2/(gmma+1) )*(1+ (gm1/2)*m_t_new^2) )^(-(gmma+1)/(2*gm1))

