#include "stiffness_checker/Material.h"

namespace conmech
{
namespace material
{
// https://github.com/jpanetta/MeshFEM/blob/master/src/lib/MeshFEM/Materials.cc
void ConstantMaterial::setFromJson(const nlohmann::json &config) {
    std::string type = config["type"];
    if (type == "isotropic_material")        { parseIsotropic<_N>(config, m_E);   }
    else if (type == "isotropic")            { parseIsotropic<_N>(config, m_E);   }
    else if (type == "orthotropic_material") { parseOrthotropic<_N>(config, m_E); }
    else if (type == "orthotropic")          { parseOrthotropic<_N>(config, m_E); }
    else if (type == "symmetric_material")   { parseAnisotropic<_N>(config, m_E); }
    else if (type == "anisotropic")          { parseAnisotropic<_N>(config, m_E); }
    else { throw std::runtime_error("Invalid type."); }
}

void ConstantMaterial::setFromFile(const std::string &materialPath) {
    std::ifstream is(materialPath);
    if (!is.is_open()) {
        throw std::runtime_error("Couldn't open material " + materialPath);
    }
    nlohmann::json config;
    is >> config;
    setFromJson(config);
}

nlohmann::json ConstantMaterial::getJson() const {
    nlohmann::json config;

    std::array<std::array<Real, flatLen(N)>, flatLen(N)> mat;
    for (size_t i = 0; i < flatLen(N); ++i) {
        for (size_t j = 0; j < flatLen(N); ++j) {
            mat[i][j] = m_E.D(i, j);
        }
    }

    config["type"] = "anisotropic";
    config["material_matrix"] = mat;

    return config;
}

} // end ns material
} // end ns conmech