#ifndef __IG_H_MESH_ 
#define __IG_H_MESH_

#include "util.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>
#include <cassert>

namespace kanba
{
    struct Cell
    {
        std::int32_t m_blank;
        double* m_ptr_vn;
        double* m_ptr_ve;
        double* m_ptr_vs;
        double* m_ptr_vw;
        double* m_ptr_pc;
        Cell* m_ptr_cn;
        Cell* m_ptr_ce;
        Cell* m_ptr_cs;
        Cell* m_ptr_cw;
        Cell() : 
                 m_blank(1),
                 m_ptr_vn(nullptr),
                 m_ptr_ve(nullptr),
                 m_ptr_vs(nullptr),
                 m_ptr_vw(nullptr),
                 m_ptr_pc(nullptr),
                 m_ptr_cn(nullptr),
                 m_ptr_ce(nullptr),
                 m_ptr_cs(nullptr),
                 m_ptr_cw(nullptr) {}
    };

    // Square Cylinder must have 4 Cells of offset for all directions!
    struct SquareCylinder
    {
        Pair<std::size_t> m_origin;
        Pair<std::size_t> m_origin_ext;
        Pair<std::size_t> m_length;
        SquareCylinder(const Pair<std::size_t> origin, const Pair<std::size_t> length) : 
            m_origin(origin),
            m_origin_ext(origin.x + 2, origin.y + 2),
            m_length(length) {}
    };

    class Mesh
    {
        private:
        // Size : MeshSize (Exclude outer boundary condition area)
        Pair<std::size_t> m_size;
        // Size : MeshSize (Include outer boudary condition area)
        Pair<std::size_t> m_size_ext;
        // Size : ElemSize (Max vector element size)
        Pair<std::size_t> m_size_ext_pls;
        // Length : MeshLength (Exclude outer boundary condition area)
        Pair<double> m_length;

        double m_delta_t;

        // X direction flow velocity on staggered mesh edge Matrix (u)
        // DIM : (m_size_ext.x + 1, m_size_ext.y)
        std::vector<std::vector<double>> m_vec_u;
        // Y direction flow velocity on staggered mesh edge Matrix (v)
        // DIM : (m_size_ext.x, m_size_ext.y + 1)
        std::vector<std::vector<double>> m_vec_v;
        // Pressure on center of staggered mesh Matrix(p)
        // DIM : (m_size_ext.x, m_size_ext.y)
        std::vector<std::vector<double>> m_vec_p;

        // Cell of references Matrix (c)
        // DIM : (m_size_ext.x, m_size_ext.y)
        std::vector<std::vector<Cell>> m_vec_c;

        // Square Cylinder registry Vector (s)
        std::vector<SquareCylinder> m_vec_s;
 
        public:
        // Constructor (Setup sizes, initializers)
        Mesh(const Pair<std::size_t> size, const Pair<double> length, double delta_t) : 
            m_size(size), 
            m_size_ext(size.x + 2, size.y + 2),
            m_size_ext_pls(size.x + 3, size.y + 3),
            m_length(length),
            m_delta_t(delta_t)
        {
            // Initialize vectors
            m_vec_u = std::vector<std::vector<double>>(m_size_ext_pls.x, std::vector<double>(m_size_ext.y, 0.0));
            m_vec_v = std::vector<std::vector<double>>(m_size_ext.x, std::vector<double>(m_size_ext_pls.y, 0.0));
            m_vec_p = std::vector<std::vector<double>>(m_size_ext.x, std::vector<double>(m_size_ext.y, 0.0));
            m_vec_c = std::vector<std::vector<Cell>>(m_size_ext.x, std::vector<Cell>(m_size_ext.y, Cell()));
            m_vec_s = std::vector<SquareCylinder>();

            // FOR ALL CELLS
            for(std::size_t xi = 0; xi < m_size_ext.x; xi++)
            {
                for(std::size_t yi = 0; yi < m_size_ext.y; yi++)
                {
                    // Setup Cells references
                    // IF IT IS NOT WEST EDGE
                    if(xi > 0)
                    {
                        m_vec_c[xi][yi].m_ptr_cw = &m_vec_c[xi - 1][yi];
                    }
                    // IF IT IS NOT EAST EDGE
                    if(xi < m_size_ext.x - 1)
                    {
                        m_vec_c[xi][yi].m_ptr_ce = &m_vec_c[xi + 1][yi];
                    }
                    // IF IT IS NOT SOUTH EDGE
                    if(yi > 0)
                    {
                        m_vec_c[xi][yi].m_ptr_cs = &m_vec_c[xi][yi - 1];
                    }
                    // IF IT IS NOT NORTH EDGE
                    if(yi < m_size_ext.y - 1)
                    {
                        m_vec_c[xi][yi].m_ptr_cn = &m_vec_c[xi][yi + 1];
                    }
                    // Setup flow velocities, pressure references
                    m_vec_c[xi][yi].m_ptr_pc = &m_vec_p[xi][yi];
                    m_vec_c[xi][yi].m_ptr_vw = &m_vec_u[xi][yi];
                    m_vec_c[xi][yi].m_ptr_ve = &m_vec_u[xi + 1][yi];
                    m_vec_c[xi][yi].m_ptr_vs = &m_vec_v[xi][yi];
                    m_vec_c[xi][yi].m_ptr_vn = &m_vec_v[xi][yi + 1];
                }
            }
        }
        void add_square_cylinder(SquareCylinder object)
        {
            // Register object to the vector
            m_vec_s.push_back(object);
            // Setup cell information
            for(std::size_t xi = object.m_origin.x + 2; xi < object.m_origin.x + 2 + object.m_length.x; xi++)
            {
                for(std::size_t yi = object.m_origin.y + 2; yi < object.m_origin.y + 2 + object.m_length.y; yi++)
                {
                    m_vec_c[xi][yi].m_blank = 0; 
                }
            }
        }
        void set_uniform_flow(const Pair<double> velocity)
        {
            // X direction velocity
            for(std::size_t xi = 0; xi < m_size_ext_pls.x; xi++)
            {
                for(std::size_t yi = 0; yi < m_size_ext.y; yi++)
                {
                    m_vec_u[xi][yi] = velocity.x;
                }
            }
            // Y direction velocity
            for(std::size_t xi = 0; xi < m_size_ext.x; xi++)
            {
                for(std::size_t yi = 0; yi < m_size_ext_pls.y; yi++)
                {
                    m_vec_v[xi][yi] = velocity.y;
                }
            }
        }
        void process_outer_boudary_condition()
        {
            // Inflow Outflow Bounday (West and East outer BC)
            for(std::size_t yi = 0; yi < m_size_ext.y; yi++)
            {
                // X direction Velocity
                m_vec_u[0][yi] = 1.0;
                m_vec_u[1][yi] = 1.0;
                m_vec_u[m_size_ext_pls.x - 1][yi] = m_vec_u[m_size_ext_pls.x - 3][yi];
                m_vec_u[m_size_ext_pls.x - 2][yi] = m_vec_u[m_size_ext_pls.x - 3][yi];
                // Center Pressure
                m_vec_p[0][yi] = m_vec_p[1][yi];
                m_vec_p[m_size_ext.x - 1][yi] = m_vec_p[m_size_ext.y - 2][yi];
            }
            for(std::size_t yi = 0; yi < m_size_ext_pls.y; yi++)
            {
                // Y direction Velocity
                m_vec_v[0][yi] = 0.0;
                m_vec_v[1][yi] = 0.0;
                m_vec_v[m_size_ext.x - 1][yi] = m_vec_v[m_size_ext.x - 3][yi];
                m_vec_v[m_size_ext.y - 2][yi] = m_vec_v[m_size_ext.x - 3][yi];
            }
            // Slip Boundary (North and South outer BC)
            for(std::size_t xi = 2; xi < m_size_ext_pls.x - 2; xi++)
            {
                // X direction velocity (South)
                m_vec_u[xi][0] = m_vec_u[xi][2];
                m_vec_u[xi][1] = m_vec_u[xi][2];
                // X direction velocity (North)
                m_vec_u[xi][m_size_ext.y - 1] = m_vec_u[xi][m_size_ext.y - 3];
                m_vec_u[xi][m_size_ext.y - 2] = m_vec_u[xi][m_size_ext.y - 3];
            }
            for(std::size_t xi = 2; xi < m_size_ext.x - 2; xi++)
            {
                // Y direction velocity (South)
                m_vec_v[xi][0] = -m_vec_v[xi][3];
                m_vec_v[xi][1] = -m_vec_v[xi][2];
                // Y direction velocity (North)
                m_vec_v[xi][m_size_ext_pls.y - 1] = -m_vec_v[xi][m_size_ext_pls.y - 4]; 
                m_vec_v[xi][m_size_ext_pls.y - 2] = -m_vec_v[xi][m_size_ext_pls.y - 3];
            }
            for(std::size_t xi = 1; xi < m_size_ext_pls.x - 2; xi++)
            {
                // Pressure (South)
                m_vec_p[xi][0] = m_vec_p[xi][1];
                // Pressure (North)
                m_vec_p[xi][m_size_ext.y - 1] = m_vec_p[xi][m_size_ext.y - 2];
            }
        }
        void process_wall_boudary_condition()
        {
            // FOR ALL CELLS IN COMPUTING AREA
            for(std::size_t xi = 2; xi < m_size_ext.x - 2; xi++)
            {
                for(std::size_t yi = 2; yi < m_size_ext.y - 2; yi++)
                {
                    Cell& cell = m_vec_c[xi][yi];

                    //Clear cell's velocities
                    *(cell.m_ptr_vn) = 0.0;
                    *(cell.m_ptr_ve) = 0.0;
                    *(cell.m_ptr_vs) = 0.0;
                    *(cell.m_ptr_vw) = 0.0;

                    // IF CELL IS IN OBJECT
                    if(cell.m_blank == 0)
                    {
                        // IF SURROUDED BY SQUARE CYLINDER CELL
                        if(cell.m_ptr_cn->m_blank == 0 && cell.m_ptr_ce->m_blank == 0 && cell.m_ptr_cs->m_blank == 0 && cell.m_ptr_cw->m_blank == 0)
                        {
                            // Set pressure to 0.0
                            *(cell.m_ptr_pc) = 0.0;
                        }
                        // IF SOUTH CELL IS NOT IN OBJECT (SOUTH EDGE OF SQUARE CYLINDER)
                        if(cell.m_ptr_cs->m_blank == 1)
                        {
                            // Pressure
                            *(cell.m_ptr_pc) = *(cell.m_ptr_cs->m_ptr_pc);
                            // Insert adjacent flow velocities
                            if(cell.m_ptr_cw->m_blank != 1)
                            {
                                *(cell.m_ptr_vw) = -*(cell.m_ptr_cs->m_ptr_vw);
                                *(cell.m_ptr_cn->m_ptr_vw) = -*(cell.m_ptr_cs->m_ptr_cs->m_ptr_vw);
                            }
                            *(cell.m_ptr_vn) = -*(cell.m_ptr_cs->m_ptr_vs);
                            *(cell.m_ptr_cn->m_ptr_vn) = -*(cell.m_ptr_cs->m_ptr_cs->m_ptr_vs);
                            // Wall Surface flow velocity
                            *(cell.m_ptr_vs) = 0.0;
                        }
                        // IF NORTH CELL IS NOT IN OBJECT (NORTH EDGE OF SQUARE CYLINDER)
                        if(cell.m_ptr_cn->m_blank == 1)
                        {
                            // Pressure
                            *(cell.m_ptr_pc) = *(cell.m_ptr_cn->m_ptr_pc);
                            // Insert adjacent flow velocities
                            if(cell.m_ptr_cw->m_blank != 1)
                            {
                                *(cell.m_ptr_vw) = -*(cell.m_ptr_cn->m_ptr_vw);
                                *(cell.m_ptr_cs->m_ptr_vw) = -*(cell.m_ptr_cn->m_ptr_cn->m_ptr_vw);
                            }
                            *(cell.m_ptr_vs) = -*(cell.m_ptr_cn->m_ptr_vn);
                            *(cell.m_ptr_cs->m_ptr_vs) = -*(cell.m_ptr_cn->m_ptr_cn->m_ptr_vn);
                            // Wall Surface flow velocity
                            *(cell.m_ptr_vn) = 0.0;
                        }
                        // IF WEST CELL IS NOT IN OBJECT (WEST EDGE OF SQUARE CYLINDER)
                        if(cell.m_ptr_cw->m_blank == 1 && cell.m_ptr_cs->m_blank != 1 && cell.m_ptr_cn->m_blank != 1)
                        {
                            // Pressure
                            *(cell.m_ptr_pc) = *(cell.m_ptr_cw->m_ptr_pc);
                            // Insert adjacent flow velocities
                            *(cell.m_ptr_vs) = -*(cell.m_ptr_cw->m_ptr_vs);
                            *(cell.m_ptr_ce->m_ptr_vs) = -*(cell.m_ptr_cw->m_ptr_cw->m_ptr_vs);
                            *(cell.m_ptr_ve) = -*(cell.m_ptr_cw->m_ptr_vw);
                            *(cell.m_ptr_ce->m_ptr_ve) = -*(cell.m_ptr_cw->m_ptr_cw->m_ptr_vw);
                            // Wall Surface flow velocity
                            *(cell.m_ptr_vw) = 0.0;
                        }
                        // IF EAST CELL IS NOT IN OBJECT (EAST EDGE OF SQUARE CYLINDER)
                        if(cell.m_ptr_ce->m_blank == 1 && cell.m_ptr_cs->m_blank != 1 && cell.m_ptr_cn->m_blank != 1)
                        {
                            // Pressure
                            *(cell.m_ptr_pc) = *(cell.m_ptr_ce->m_ptr_pc);
                            // Insert adjacent flow velocities
                            *(cell.m_ptr_vs) = -*(cell.m_ptr_ce->m_ptr_vs);
                            *(cell.m_ptr_cw->m_ptr_vs) = -*(cell.m_ptr_ce->m_ptr_ce->m_ptr_vs);
                            *(cell.m_ptr_vw) = -*(cell.m_ptr_ce->m_ptr_ve);
                            *(cell.m_ptr_cw->m_ptr_vw) = -*(cell.m_ptr_ce->m_ptr_ce->m_ptr_ve);
                            // Wall Surface flow velocity
                            *(cell.m_ptr_ve) = 0.0;
                        }
                        // OVERWRITE SOUTH WEST CORNER PRESSURE
                        if(cell.m_ptr_cw->m_blank == 1 && cell.m_ptr_cs->m_blank == 1)
                        {
                            *(cell.m_ptr_pc) = *(cell.m_ptr_cs->m_ptr_cw->m_ptr_pc);
                        }
                        // OVERWRITE NORTH WEST CORNER PRESSURE
                        if(cell.m_ptr_cw->m_blank == 1 && cell.m_ptr_cn->m_blank == 1)
                        {
                            *(cell.m_ptr_pc) = *(cell.m_ptr_cn->m_ptr_cw->m_ptr_pc);
                        }
                        // OVERWRITE SOUTH EAST CORNER PRESSURE
                        if(cell.m_ptr_ce->m_blank == 1 && cell.m_ptr_cs->m_blank == 1)
                        {
                            *(cell.m_ptr_pc) = *(cell.m_ptr_cs->m_ptr_ce->m_ptr_pc);
                        }
                        // OVERWRITE NORTH EAST CORNER PRESSURE
                        if(cell.m_ptr_ce->m_blank == 1 && cell.m_ptr_cn->m_blank == 1)
                        {
                            *(cell.m_ptr_pc) = *(cell.m_ptr_cn->m_ptr_ce->m_ptr_pc);
                        }
                    }
                }
            }
        }
        void process_HSMAC(std::size_t max_iter, double relax_coef)
        {
            double dx = m_length.x / static_cast<double>(m_size.x);
            double dy = m_length.y / static_cast<double>(m_size.y);
            for(std::size_t i = 0; i < max_iter; i++)
            {
                process_outer_boudary_condition();
                process_wall_boudary_condition();

                double max_divergence = 0.0;
                double coef_0 = (-relax_coef) / (2.0 * m_delta_t * (1.0 / (dx * dx) + 1.0 / (dy * dy)));
                for(std::size_t xi = 1; xi < m_size_ext.x - 1; xi++)
                {
                    for(std::size_t yi = 1; yi < m_size_ext.y - 1; yi++)
                    {
                        Cell& cell = m_vec_c[xi][yi];
                        // Calculation
                        double cell_divergence = static_cast<double>(cell.m_blank) * (*cell.m_ptr_ve - *cell.m_ptr_vw) / dx + (*cell.m_ptr_vn - *cell.m_ptr_vs) / dy;
                        double delta_p = coef_0 * cell_divergence;
                        // Relaxation
                        *(cell.m_ptr_pc) += delta_p;
                        *(cell.m_ptr_vw) += (m_delta_t / dx) * delta_p * static_cast<double>(cell.m_ptr_cw->m_blank);
                        *(cell.m_ptr_ve) -= (m_delta_t / dx) * delta_p * static_cast<double>(cell.m_ptr_ce->m_blank);
                        *(cell.m_ptr_vs) += (m_delta_t / dy) * delta_p * static_cast<double>(cell.m_ptr_cs->m_blank);
                        *(cell.m_ptr_vn) -= (m_delta_t / dy) * delta_p * static_cast<double>(cell.m_ptr_cn->m_blank);

                        // Update
                        max_divergence = std::max(max_divergence, cell_divergence);
                    } 
                }
                //if(max_divergence < EPS){}
            }
        }
        void process_forward_rk4()
        {
            Mesh k_1 = Mesh(m_size, m_length, m_delta_t);
            Mesh k_2 = Mesh(m_size, m_length, m_delta_t);
            Mesh k_3 = Mesh(m_size, m_length, m_delta_t);
            Mesh k_4 = Mesh(m_size, m_length, m_delta_t);

            get_rhs(this, &k_1);
            k_1.store_added_velocities(this, &k_1);
            k_1.store_multiplied_velocity(&k_1, m_delta_t / 2.0);
            get_rhs(&k_1, &k_2);
            k_2.store_added_velocities(this, &k_2);
            k_2.store_multiplied_velocity(&k_2, m_delta_t / 2.0);
            get_rhs(&k_2, &k_3);
            k_3.store_added_velocities(this, &k_3);
            k_3.store_multiplied_velocity(&k_3, m_delta_t);
            get_rhs(&k_3, &k_4);
        }
        void get_rhs(Mesh* base, Mesh* store)
        {
            for(std::size_t xi = 2; xi < m_size_ext.x - 2; xi++)
            {
                for(std::size_t yi = 2; yi < m_size_ext.y - 2; yi++)
                {
                    Cell& cell_b = base->m_vec_c[xi][yi];
                    Cell& cell_s = store->m_vec_c[xi][yi];
                }
            }
        }
        void store_added_velocities(Mesh* a, Mesh* b)
        {
            for(std::size_t xi = 0; xi < m_size_ext_pls.x; xi++)
            {
                for(std::size_t yi = 0; yi < m_size_ext.y; yi++)
                {
                    m_vec_u[xi][yi] = a->m_vec_u[xi][yi] + b->m_vec_u[xi][yi];
                }
            }
            for(std::size_t xi = 0; xi < m_size_ext.y; xi++)
            {
                for(std::size_t yi = 0; yi < m_size_ext_pls.y; yi++)
                {
                    m_vec_v[xi][yi] = a->m_vec_v[xi][yi] + b->m_vec_v[xi][yi];
                }
            }
        }
        void store_multiplied_velocity(Mesh* a, double coef)
        { 
            for(std::size_t xi = 0; xi < m_size_ext_pls.x; xi++)
            {
                for(std::size_t yi = 0; yi < m_size_ext.y; yi++)
                {
                    m_vec_u[xi][yi] = coef * a->m_vec_u[xi][yi];
                }
            }
            for(std::size_t xi = 0; xi < m_size_ext.y; xi++)
            {
                for(std::size_t yi = 0; yi < m_size_ext_pls.y; yi++)
                {
                    m_vec_v[xi][yi] = coef * a->m_vec_v[xi][yi];
                }
            }
        }
    };
}

#endif
