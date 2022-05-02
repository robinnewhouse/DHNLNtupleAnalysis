################################################################################
# Package: DHNLNtupleAnalysis
################################################################################

# Declare the package name:
atlas_subdir( DHNLNtupleAnalysis )

# External dependencies:
find_package( ROOT COMPONENTS Core MathCore Hist Graf Gpad File Tree Graf2d RIO )


atlas_add_library( DHNLNtupleAnalysis
  DHNLNtupleAnalysis/*.h Root/*.cxx
  PUBLIC_HEADERS DHNLNtupleAnalysis
  INCLUDE_DIRS ${ROOT_INCLUDE_DIRS} 
  PRIVATE_INCLUDE_DIRS Root 
  LINK_LIBRARIES ${ROOT_LIBRARIES} pthread NtupleAnalysisUtils
  # ${extra_libs}
  )


# Build all executables in util/ 
file(GLOB_RECURSE files "util/*.cxx")

foreach(_exeFile ${files})
    get_filename_component(_theExec ${_exeFile} NAME_WE)
    get_filename_component(_theLoc ${_exeFile} DIRECTORY)
    # we specify a folder for programs we do not want to compile. Useful during release transition...
    if(${_theLoc} MATCHES "DoNotBuild" OR ${_theLoc} MATCHES "ToUpdate")
    else()
      # build an executable program matching the cxx file name, 
      # linked against NTAU
      atlas_add_executable( ${_theExec}
      ${_exeFile}
      INCLUDE_DIRS ${ROOT_INCLUDE_DIRS} 
      LINK_LIBRARIES ${ROOT_LIBRARIES} DHNLNtupleAnalysis NtupleAnalysisUtils )
    endif()
endforeach()

atlas_install_python_modules( python/*.py python/generated_files_info )