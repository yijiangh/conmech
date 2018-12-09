# this file is modified from conmech:
# https://github.com/libconmech/libconmech/blob/master/cmake/LibconmechDownloadExternal.cmake
################################################################################
include(DownloadProject)

# Shortcut function
function(conmech_download_project name)
    download_project(
            PROJ         ${name}
            SOURCE_DIR   ${CONMECH_EXTERNAL}/${name}
            DOWNLOAD_DIR ${CONMECH_EXTERNAL}/.cache/${name}
            QUIET
            ${CONMECH_EXTRA_OPTIONS}
            ${ARGN}
    )
endfunction()

################################################################################

### Eigen
#function(conmech_download_eigen)
#    conmech_download_project(eigen
#            URL           http://bitbucket.org/eigen/eigen/get/3.2.10.tar.gz
#            URL_MD5       8ad10ac703a78143a4062c9bda9d8fd3
#            )
#endfunction()
#
### pybind11
#function(conmech_download_pybind11)
#    conmech_download_project(pybind11
#            GIT_REPOSITORY https://github.com/pybind/pybind11.git
#            GIT_TAG        2d0507db43cd5a117f7843e053b17dffca114107
#            )
#endfunction()

## GLFW
function(conmech_download_glfw)
    conmech_download_project(glfw
            GIT_REPOSITORY https://github.com/glfw/glfw.git
            GIT_TAG        58cc4b2c5c2c9a245e09451437dd6f5af4d60c84
            )
endfunction()

## ImGui
function(conmech_download_imgui)
    conmech_download_project(imgui
            GIT_REPOSITORY https://github.com/ocornut/imgui.git
            GIT_TAG        bc6ac8b2aee0614debd940e45bc9cd0d9b355c86
            )
    conmech_download_project(libigl-imgui
            GIT_REPOSITORY https://github.com/libigl/libigl-imgui.git
            GIT_TAG        a37e6e59e72fb07bd787dc7e90f72b9e1928dae7
            )
endfunction()