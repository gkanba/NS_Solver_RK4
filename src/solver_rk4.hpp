#ifndef __IG_H_SOLVER_
#define __IG_H_SOLVER_

#include "mesh.hpp"

namespace kanba
{

    class SolverRK4
    {
        private:
            Mesh& m_mesh;
            Mesh m_mesh_updated;
        public:
            SolverRK4(const Mesh& mesh): m_mesh(mesh)
            {
            }
    };
}

#endif
