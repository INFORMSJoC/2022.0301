vcpkg_from_gitlab(
  GITLAB_URL
  https://gitlab.kuleuven.be
  OUT_SOURCE_PATH
  SOURCE_PATH
  REF
  v0.1.3
  SHA512
  fe96982c47765b8cfe1883dd49823b86ef4f284823dc70f1a1ffc2ba4afbbd91903297ee9eb51070816e641f9cffdc26c82ce7451e7080ff381a2e44c2453a0b
  REPO
  u0056096/branch-and-bound
)

vcpkg_cmake_configure(SOURCE_PATH ${SOURCE_PATH})

vcpkg_cmake_install()
vcpkg_cmake_config_fixup(PACKAGE_NAME "branch-and-bound")
file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include")
vcpkg_copy_pdbs()
