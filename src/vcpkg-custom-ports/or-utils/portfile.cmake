vcpkg_from_gitlab(
  GITLAB_URL
  https://gitlab.kuleuven.be
  OUT_SOURCE_PATH
  SOURCE_PATH
  REF
  v1.0.3
  SHA512
  d872ceaad2d8a0adffc716885a69e90e0c8d851ade53aadcc869a14502ee3fa63ea9c4f41158af7734a4993b464443e0150d164c99eb7e9501aead71e4240f37
  REPO
  u0056096/or-utils
  HEAD_REF
  develop
)

if(VCPKG_TARGET_IS_WINDOWS)
  set(ENV{GUROBI_HOME} "C:/gurobi912/win64")
endif()
vcpkg_cmake_configure(SOURCE_PATH ${SOURCE_PATH})

vcpkg_cmake_install()
file(
  INSTALL "${SOURCE_PATH}/LICENSE.md"
  DESTINATION "${CURRENT_PACKAGES_DIR}/share/${PORT}"
  RENAME copyright
)
vcpkg_cmake_config_fixup()
file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include")
vcpkg_copy_pdbs()
