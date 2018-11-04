function [K_loc] = local_stiffness_matrix(L, A, Jx, Iy, Iz, dim, m_p)
DA = m_p.E * A;
mu = m_p.mu;

switch dim
    case 2
        DS = m_p.G * A;
        DB = m_p.E * Iy;
        DT = DB / (L^2*DS + 12*DB);
        
        K_loc = zeros(6,6);
        K_loc(1,4) = -DA;
        K_loc(2,3) = 6*L*DT*DS;
        K_loc(2,5) = -12*DT*DS;
        K_loc(2,6) = 6*L*DT*DS;
        K_loc(3,5) = -6*L*DT*DS;
        K_loc(3,6) = DT*(2*L^2-12*DB);
        K_loc(5,6) = -6*L*DT*DS;
        
        K_loc = K_loc + K_loc';
        K_loc(1,1) = DA;
        K_loc(2,2) = 12*DT*DS;
        K_loc(3,3) = DT*(4*L^2*DS + 12*DB);
        K_loc(4,4) = DA;
        K_loc(5,5) = 12*DT*DS;
        K_loc(6,6) = DT*(4*L^2 + 12*DB);
        
        K_loc = K_loc / L;
        
    case 3
        K_loc = zeros(12,12);
        K_block = zeros(6,6);
        
        % block_00 and block_11
        K_block(2,6) = 6*Iz / L^2;
        K_block(3,5) = - 6*Iy / L^2;
        K_block = K_block +  K_block';
        K_loc(1:6,1:6) = K_block;
        K_loc(7:12,7:12) = -K_block;
        
        d_v = zeros(6);
        d_v(1) = A/L;
        d_v(2) = 12*Iz / L^3;
        d_v(3) = 12*Iy / L^3;
        d_v(4) = Jx / (2*(1+mu)*L);
        d_v(5) = 4*Iy / L;
        d_v(6) = 4*Iz / L;
        K_loc(1:6,1:6) = K_loc(1:6, 1:6) + diag(d_v);
        K_loc(7:12,7:12) = K_loc(7:12,7:12) + diag(d_v);
        
        % block_01 and block_10
        K_block = zeros(6,6);
        
        K_block(2,6) = 6*Iz / L^2;
        K_block(3,5) = - 6*Iy / L^2;
        K_block = K_block - K_block';
        K_loc(1:6,7:12) = K_block;
        K_loc(7:12,1:6) = -K_block;
        
        d_v = zeros(6);
        d_v(1) = -A/L;
        d_v(2) = -12*Iz / L^3;
        d_v(3) = -12*Iy / L^3;
        d_v(4) = -Jx / (2*(1+mu)*L);
        d_v(5) = 2*Iy / L;
        d_v(6) = 2*Iz / L;
        K_loc(1:6,7:12) = K_loc(1:6,7:12) + diag(d_v);
        K_loc(7:12,1:6) = K_loc(7:12,1:6) + diag(d_v);
        
        K_loc = K_loc * m_p.E;    
end

end